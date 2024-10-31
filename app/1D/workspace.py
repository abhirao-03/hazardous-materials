import scrubbers_canistors as sg
import model_params as model
import numpy as np

parameters = model.parameters()
canistor = sg.gas_canistor()
scrubber = sg.scrubber(0.2, 0.05, 0.2)

def build_matrix(parameters):
    A = np.zeros((parameters.Nx_points, parameters.Nx_points))

    for i in range(1, parameters.Nx_points - 1):
        A[i, i - 1] = -parameters.C
        A[i, i] = 1 + 2 * parameters.C
        A[i, i + 1] = -parameters.C

    A[0, 0] = 1 + 2 * parameters.C
    A[0, 1] = -2 * parameters.C
    A[parameters.Nx_points - 1, parameters.Nx_points - 2] = -2 * parameters.C
    A[parameters.Nx_points - 1, parameters.Nx_points - 1] = 1 + 2 * parameters.C

    return A

A = build_matrix(parameters=parameters)

x = np.linspace(0, parameters.length, parameters.Nx_points)
u_initial = canistor.get_initial_concentration(x)

U = np.zeros((parameters.Nt_points, parameters.Nx_points))
U[0, :] = u_initial

u = u_initial.copy()

for n in range(1, parameters.Nt_points):
    u = np.linalg.solve(A, u)
    u = scrubber.sink(x, u)
    U[n, :] = u
    print(f'On iteration {n}')

print('COMPLETED ALL ITERATIONS')