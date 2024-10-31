import numpy as np
import scipy.sparse as sparse

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.animation import FuncAnimation


class gas_canistor:
    def __init__(self, loc = (0.5, 0.5), radius = 0.05, concentration = 1.0):
        self.radius = radius

        self.x = loc[0]
        self.y = loc[1]

        self.lb_x = self.x - self.radius
        self.ub_x = self.x + self.radius
        self.lb_y = self.y - self.radius
        self.ub_y = self.y + self.radius

        self.concentration = concentration

def f(canistor: gas_canistor, x: np.array, y: np.array):
    concentrations = np.where((x >= canistor.lb_x) & (x <= canistor.ub_x) & (y >= canistor.lb_y) & (y <= canistor.ub_y), canistor.concentration, 0)
    return concentrations

can_1 = gas_canistor()
x = np.linspace(0.0, 1.0, 100)
y = np.linspace(0.0, 1.0, 100)
concentrations = f(can_1, x, x)



L = 1.0
T = 10.0                 # Total time
Nt = 1000                # Number of time steps
Nx = 100
Ny = 100
dx = L / (Nx - 1)       # Spatial step size
dy = L / (Ny - 1)
dt = T / Nt             # Time step size

Cx = dt / dx**2
Cy = dt / dy**2

x = np.linspace(0, L, Nx)
y = np.linspace(0, L, Ny)

X, Y = np.meshgrid(x, y)
Z = f(can_1, X, Y)

u_initial = Z.reshape((Nx*Ny,))

B = np.zeros((Nx, Ny))

a = 1 + 2*dt*(1/dx**2 + 1/dy**2)

for i in range(Nx - 1):
    B[i, i] = a
    B[i+1, i] = -Cx
    B[i, i+1] = -Cx

B[-1, -1] = a

A = sparse.bmat([[B if i == j else -Cy * np.eye(Nx) if abs(i-j)==1
                else None for i in range(Nx)]
                for j in range(Ny)], format='csr')

u = u_initial.copy()

for n in range(Nt - 1):
    u = sparse.linalg.spsolve(A, u)
    print(f"On iteration: {n}")


X, Y = np.meshgrid(x, y)

fig = plt.figure()
ax = plt.axes(projection='3d')

def update(n):
    global u
    u = sparse.linalg.spsolve(A, u)
    ax.clear()
    ax.plot_surface(X, Y, u.reshape((Nx, Ny)), rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_title(f"Time step: {n}")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('Concentration')

anim = FuncAnimation(fig, update, frames=Nt, interval=50)

plt.show()