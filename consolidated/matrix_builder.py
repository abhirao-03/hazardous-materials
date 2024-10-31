import numpy as np
import scipy.sparse as sp
import model_params as model

def build_matrix(parameters: model.parameters):
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

    B = np.zeros((parameters.Nx_points, parameters.Ny_points))
    a = 1 + 2*parameters.dt*(1/parameters.dx**2 + 1/parameters.dy**2)

    for i in range(parameters.Nx_points - 1):
        B[i, i] = a
        B[i+1, i] = -parameters.Cx
        B[i, i+1] = -parameters.Cx

    B[-1, -1] = a

    A = sp.bmat([[B if i == j else -parameters.Cy * np.eye(parameters.Nx_points) if abs(i-j)==1
                    else None for i in range(parameters.Nx_points)]
                    for j in range(parameters.Ny_points)], format='csr')
    
    return A