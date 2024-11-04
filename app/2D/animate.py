import numpy as np
import model_params as model
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from workspace import *

parameters = model.parameters()

U_tracked = run_simulation()

min_concentration = 0
max_concentration = 0.2

fig, ax = plt.subplots()
heatmap = ax.imshow(U_tracked[0].reshape(parameters.Nx_points, parameters.Ny_points), cmap='magma', vmin=min_concentration, vmax=max_concentration, interpolation='nearest')
heatmap = ax.imshow(U_tracked[0].reshape(parameters.Nx_points, parameters.Ny_points), cmap='magma', vmin=min_concentration, vmax=max_concentration, interpolation='nearest')
ax.set_title("Heatmap over Time")

for scrubber in scrubbers:
    fixed_point, = ax.plot(scrubber.x_loc*100, scrubber.y_loc*100, 'o', color='cyan', markersize=8, alpha=0.2)

def update(frame):
    heatmap.set_array(U_tracked[frame].reshape(parameters.Nx_points, parameters.Ny_points))
    ax.set_title(f"Time Step: {frame}")
    return heatmap,

ani = FuncAnimation(fig, update, frames=parameters.Nt_points, interval=100, blit=False, repeat=True)

plt.colorbar(heatmap, ax=ax)
plt.show()