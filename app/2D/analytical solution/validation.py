import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.integrate as integrate
from scipy.integrate import dblquad


import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import scrubbers_canisters as sg
import model_params as model

gas_canister = sg.GasCan2D(x_loc=0.5, y_loc=0.5, radius=0.1, concentration=1.0)

def run_simulation(gas_canister = gas_canister):
    parameters   = model.parameters()

    def build_matrix(parameters: model.parameters, scrubbers: list, x, y):
        """
        Builds the large block matric necessary for 2d calculations.

        Args:
            parameters        Uses parameters defined in the model to compute the large block Matrix
                Nx_points     number of x spatial points
                Ny_points     number of y spatial points
                dt, dx, dy    temporal and spatial step size
                Cx, Cy        C[_] = dt / (d[_]**2)

        Outputs:
            A           Large block matrix described on page __  
        """
        b = 1 + 2*parameters.diffusion_coeff*parameters.dt*(1/parameters.dx**2 + 1/parameters.dy**2)
        B = np.diag(b)

        for i in range(parameters.Nx_points - 1):
            B[i+1, i] = -parameters.Cx
            B[i, i+1] = -parameters.Cx

        A = sparse.bmat([[B if i == j else -parameters.Cy * np.eye(parameters.Nx_points) if abs(i-j)==1
                        else None for i in range(parameters.Nx_points)]
                        for j in range(parameters.Ny_points)], format='csr')
        
        return A

    x = np.linspace(0, parameters.len_x, parameters.Nx_points)
    y = np.linspace(0, parameters.len_y, parameters.Ny_points)

    A = build_matrix(parameters=parameters, x=x, y=y)
    U = np.zeros((parameters.Nt_points, parameters.Nx_points * parameters.Ny_points))

    u_init = np.flipud(gas_canister.get_initial_concentration(x, y))
    u_init = u_init.reshape((parameters.Nx_points * parameters.Ny_points),)
  
    u = u_init.copy()
    U[0, :] = u_init

    for n in range(1, parameters.Nt_points):
        u = sp.linalg.minres(A, u)[0]
        U[n, :] = u
        print(f"On iteration {n}")

    print('COMPLETED ALL ITERATIONS')

    return U


parameters = model.parameters()

N, M = 20, 20

lamb_nm = np.zeros((N,M))

for i in range(N):
    for j in range(M):
        D = parameters.diffusion_coeff
        lamb_nm[i][j] = (((i * np.pi)**2) + ((j * np.pi)**2))


# Initialize matrices
A_nm = np.zeros((N, M))

# Calculate A_nm using numerical integration
for n in range(N - 1):
    for m in range(M - 1):
        y_term = (-np.sin(m * np.pi * gas_canister.y_u) + np.sin(m * np.pi * gas_canister.y_l))/(m+1)
        x_term = (-np.sin(n * np.pi * gas_canister.x_u) + np.sin(n * np.pi * gas_canister.x_l))/(m+1)

        A_nm[n][m] = 4 * y_term * x_term


save_A_nm = os.path.join(os.path.dirname(__file__), 'A_nm.npy')
save_lamb_nm = os.path.join(os.path.dirname(__file__), 'lamb_nm.npy')

np.save(save_A_nm, A_nm)
np.save(save_lamb_nm, lamb_nm)


print('COMPLETED BUILDING A_nm AND lamb_nm')