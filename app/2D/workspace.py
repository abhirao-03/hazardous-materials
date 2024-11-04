import scrubbers_canistors as sg
import numpy as np
import scipy.sparse as sp
import model_params as model

gas_canistor = sg.GasCan2D(x_loc=0.5, y_loc=0.5, radius=0.1, concentration=1.0)
scrubbers = [sg.Scrubber2D(x_loc=0.5, y_loc=0.7, radius=0.08, efficiency=300.0),
                sg.Scrubber2D(x_loc=0.5, y_loc=0.3, radius=0.08, efficiency=300.0)]


def run_simulation(gas_canistor = gas_canistor, scrubbers = scrubbers):
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
        b = 1 + 2*parameters.diffusion_coeff*parameters.dt*(1/parameters.dx**2 + 1/parameters.dy**2) + parameters.dt*(scrubbers[0].sink(x, y))

        for i in range(1, len(scrubbers)):
            b += parameters.dt*(scrubbers[i].sink(x, y))

        B = np.diag(b)

        for i in range(parameters.Nx_points - 1):
            B[i+1, i] = -parameters.Cx
            B[i, i+1] = -parameters.Cx

        A = sp.bmat([[B if i == j else -parameters.Cy * np.eye(parameters.Nx_points) if abs(i-j)==1
                        else None for i in range(parameters.Nx_points)]
                        for j in range(parameters.Ny_points)], format='csr')
        
        return A

    x = np.linspace(0, parameters.len_x, parameters.Nx_points)
    y = np.linspace(0, parameters.len_y, parameters.Ny_points)

    A = build_matrix(parameters=parameters, scrubbers=scrubbers, x=x, y=y)
    U = np.zeros((parameters.Nt_points, parameters.Nx_points * parameters.Ny_points))

    u_init = gas_canistor.get_initial_concentration(x, y).reshape((parameters.Nx_points * parameters.Ny_points),)
    u = u_init.copy()
    U[0, :] = u_init

    for n in range(1, parameters.Nt_points):
        u = sp.linalg.minres(A, u)[0]
        U[n, :] = u
        print(f"On iteration {n}")

    print('COMPLETED ALL ITERATIONS')

    return U