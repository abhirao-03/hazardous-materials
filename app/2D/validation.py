import scrubbers_canisters as sg
import numpy as np
import scipy as sp
import scipy.sparse as sparse
import model_params as model
import scipy.integrate as integrate
from scipy.integrate import dblquad

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

n, m = parameters.Nx_points, parameters.Ny_points

lamb_nm = np.zeros((n,m))

for i in range(n):
    for j in range(m):
        D = parameters.diffusion_coeff
        lamb_nm[i][j] = D * (((i * np.pi/parameters.len_x)**2) + ((j * np.pi/parameters.len_x)**2))

# Parameters (example values for a, b, and D)
a = parameters.len_x  # Length in x-direction
b = parameters.len_y  # Length in y-direction
D = parameters.diffusion_coeff  # Constant in lambda_nm formula

# Define grid size for n and m
N, M = parameters.Nx_points, parameters.Ny_points  # Size of the matrices (e.g., 5x5 matrix for A_nm and lambda_nm)

# Initialize matrices
A_nm = np.zeros((N, M))

# Calculate A_nm using numerical integration
for n in range(N - 1):
    for m in range(M - 1):
        y_term = (-np.sin(m * np.pi * gas_canister.y_u) + np.sin(m * np.pi * gas_canister.y_l))/(m+1)
        x_term = (-np.sin(n * np.pi * gas_canister.x_u) + np.sin(n * np.pi * gas_canister.x_l))/(m+1)

        A_nm[n][m] = 4 * y_term * x_term

def u_exact(x, y, t):
    for i in range(N):
        for i in range(M):
            



print()