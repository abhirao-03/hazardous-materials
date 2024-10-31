import numpy as np
import model_params as model
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from workspace import *

U_tracked = U

min_concentration = 0
max_concentration = 0.2

fig, ax = plt.subplots()
heatmap = ax.imshow(U_tracked[0].reshape(parameters.Nx_points, parameters.Ny_points), cmap='magma', vmin=min_concentration, vmax=max_concentration, interpolation='nearest')
heatmap = ax.imshow(U_tracked[0].reshape(parameters.Nx_points, parameters.Ny_points), cmap='magma', vmin=min_concentration, vmax=max_concentration, interpolation='nearest')
ax.set_title("Heatmap over Time")

fixed_point, = ax.plot(30, 50, 'o', color='cyan', markersize=8)
fixed_point, = ax.plot(70, 50, 'o', color='cyan', markersize=8)

def update(frame):
    heatmap.set_array(U_tracked[frame].reshape(parameters.Nx_points, parameters.Ny_points))
    ax.set_title(f"Time Step: {frame}")
    return heatmap,

ani = FuncAnimation(fig, update, frames=parameters.Nt_points, interval=100, blit=False, repeat=True)

plt.colorbar(heatmap, ax=ax)
plt.show()