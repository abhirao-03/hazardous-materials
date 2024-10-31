import numpy as np
import matplotlib.pyplot as plt
import model_params as model
from workspace import *
from matplotlib.animation import FuncAnimation

U_tracked = U

# Set up the figure and axis
fig, ax = plt.subplots()
x_values = np.linspace(0, 1, parameters.Nx_points)  # Replace with your actual x values
line, = ax.plot(x_values, U_tracked[0, :], color='r')

# Set up the title and labels
ax.set_title('Concentration Over Time')
ax.set_xlabel('x')
ax.set_ylabel('Concentration')

# Function to update the line for each frame
def update(frame):
    line.set_ydata(U_tracked[frame, :])  # Update the y data of the line
    return line,

# Create the animation
ani = FuncAnimation(fig, update, frames=parameters.Nt_points, blit=False)

# Display the animation
plt.show()