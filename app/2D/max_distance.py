import numpy as np
from workspace import *
import model_params as model
parameters   = model.parameters()

x_loc = gas_canister.x_loc
y_loc = gas_canister.y_loc
len_x = parameters.len_x
len_y = parameters.len_y
Nx = parameters.Nx_points
Ny = parameters.Ny_points

x = np.linspace(0, len_x, Nx)
y = np.linspace(0, len_y, Ny)
# Initialize matrix to store maximum concentration at each point
U_max = gas_canister.get_initial_concentration(x, y)

######
U = run_simulation(gas_canister = gas_canister, scrubbers = scrubbers)     

# Unpack the flattened U matrix into a 3D matrix
# ( with dimensions Nt_points, Nx_points, Ny_points )
U_3D = np.zeros((parameters.Nt_points, parameters.Nx_points, parameters.Ny_points))

for i in range(parameters.Nt_points):
    U_3D[i] = np.flipud(U[i].reshape((parameters.Nx_points, parameters.Ny_points)))
######

for time in range(0, parameters.Nt_points, 1000):
    print(f'On time {time}')
    U_current = U_3D[time]

# Compare max concentration to each time step's concentration
    for i in range(Nx):
        for j in range(Ny):
            # Updates maximum concentration at index i,j
            if U_current[i,j] > U_max[i,j]:
                U_max[i,j] =  U_current[i,j]

    # Initializes matrix to record which points would have the detector detect
    binary = np.zeros((Nx, Ny))
    for i in range(Nx):
        for j in range(Ny):
            # If maximum is over 1% then it would detect
            if U_max[i,j] >= 0.01:
                binary[i,j] = 1
        

# Finds the perimeter of the area of detected concentrations
# Intializes array to store maximum distance index from canister in each direction
dire = []

#Finds the actual length of between each x and y points
x_l = len_x/(Nx-1)
y_l = len_y/(Ny-1)

# Loops through every point
# Corners
if binary[0][0] == 1:
    if binary[0][1]*binary[1][0]*binary[1][1] != 0:
        dire.append([0,0])

if binary[0][-1] == 1:
    if binary[0][-2]*binary[1][-2]*binary[0][-1] != 0:
        dire.append([len_x,len_y])

if binary[-1][0] == 1:
    if binary[-2][0]*binary[-2][1]*binary[-1][1] != 0:
        dire.append([0,0])

if binary[-1][-1] == 1:
    if binary[-1][-2]*binary[-2][-2]*binary[-2][-1] != 0:
        dire.append([0,0])
        

# Goes through the top and bottom edges
for j in range(1,Ny-1):
    # Top
    if binary[0][j] == 1:
        if binary[0][j-1]*binary[0][j+1]*binary[1][j-1]*binary[1][j]*binary[1][j+1] != 1:
            dire.append([0,j*y_l])

    # Bottom
    if binary[-1][j] == 1:
        if binary[-1][j-1]*binary[-1][j+1]*binary[-2][j-1]*binary[-2][j]*binary[-2][j+1] != 1:
            dire.append([len_y,j*y_l])

for i in range(1, Nx-1):
    # Left and Right edges
    # Left
    if binary[i][0] == 1:
        if binary[i-1][0]*binary[i+1][0]*binary[i-1][1]*binary[i][1]*binary[i+1][1] != 1:
            dire.append([i*x_l, 0])

    # Bottom
    if binary[i][-1] == 1:
        if binary[i-1][-1]*binary[i+1][-1]*binary[i-1][-2]*binary[i][-2]*binary[i+1][-2] != 1:
            dire.append([j*y_l, len_x])

    # Other points
    for j in range(1, Ny-1):
        if binary[i][j] == 1:
        # Checks if any points surrounding binary[i,j] are not 1 (to find perimeter)
            if binary[i-1][j-1]*binary[i-1][j]*binary[i-1][j+1]*binary[i][j-1]*binary[i][j+1]*binary[i+1][j-1]*binary[i+1][j]*binary[i+1][j+1] != 1:
                # Appends real position 
                dire.append([i*x_l,j*y_l])

# Finds distance from canister
for i in range(len(dire)):
    dire[i] = np.sqrt((x_loc - dire[i][0])**2+(y_loc - j)**2)

# Finds the minimum distance from canister of the perimeter
max_distance = dire[0]

for i in range(len(dire)):
    if dire[i] < max_distance:
        max_distance == dire[i]

print()