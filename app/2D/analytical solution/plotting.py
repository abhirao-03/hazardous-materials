import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import os



load_A_nm = os.path.join(os.path.dirname(__file__), 'A_nm.npy')
load_lam_nm = os.path.join(os.path.dirname(__file__), 'lamb_nm.npy')

A_nm = np.load(load_A_nm)
lamb_nm = np.load(load_lam_nm)

N, M = A_nm.shape

def u(x, y, t):

    sum = 0
    for i in range(N):
        print(f'Iteration {i+1}')
        for j in range(M):
            A = A_nm[i][j]
            lam = lamb_nm[i][j]

            sum += A * np.cos(i*np.pi*x) * np.cos(j*np.pi*y) * np.exp(-lam*t)

    return sum


Nx_points, Ny_points = 1000, 1000
x = np.linspace(0, 1, Nx_points)
y = np.linspace(0, 1, 1000)

X, Y = np.meshgrid(x, y)

u_init = u(X, Y, 0)



min_concentration = 0
max_concentration = 0.2

fig, ax = plt.subplots()
heatmap = ax.imshow(u_init, cmap='magma', vmin=min_concentration, vmax=max_concentration, interpolation='nearest')
ax.set_title("Analytical Solution for Initial Concentration $u(x, y, 0) = f(x, y)$")

plt.colorbar(heatmap, ax=ax)
plt.show()

print()