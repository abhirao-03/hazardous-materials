import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize

class gas_canistor:
    def __init__(self, loc = 0.5, radius = 0.05, concentration = 1.0):
        self.radius = radius
        self.loc = loc
        self.lower_bound = self.loc - self.radius
        self.upper_bound = self.loc + self.radius
        self.concentration = concentration

def f(canistor: gas_canistor, x: np.array):
    return np.where((x >= canistor.lower_bound) & (x <= canistor.upper_bound), 1, 0)

class scrubber():
    def __init__(self, location, radius, efficiency):
        self.loc = location
        self.radius = radius
        self.efficiency = efficiency
        self.lower_bound = self.loc - self.radius
        self.upper_bound = self.loc + self.radius

def sink(scrub:scrubber, x, u):
    """
    loc: location of sink
    x: position
    cot: concentration
    """
    
    scrub_loc = np.where((x >= scrub.lower_bound) & (x <= scrub.upper_bound))
    u[scrub_loc] = (1 - scrub.efficiency) * u[scrub_loc]
    
    return u

scrub = scrubber(0.2, 0.05, 0.2)


# Parameters
L = 1.0                 # Length of interval
T = 10.0                 # Total time
Nx = 100                # Number of spatial points
Nt = 10000                # Number of time steps
dx = L / (Nx - 1)       # Spatial step size
dt = T / Nt             # Time step size
C = dt / dx**2

x = np.linspace(0, L, Nx)

canistor = gas_canistor()
u_initial = f(canistor, x)

u = np.copy(u_initial)

U_tracked = np.zeros((Nt, Nx))
U_tracked[0, :] = u

A = np.zeros((Nx, Nx))

for i in range(1, Nx - 1):
    A[i, i - 1] = -C
    A[i, i] = 1 + 2 * C
    A[i, i + 1] = -C

A[0, 0] = 1 + 2 * C
A[0, 1] = -2 * C   
A[Nx - 1, Nx - 2] = -2 * C 
A[Nx - 1, Nx - 1] = 1 + 2 * C


for n in range(Nt - 1):
    u = np.linalg.solve(A, u)
    u = sink(scrub, x, u)
    U_tracked[n+1, :] = u


U_tracked_x = np.copy(U_tracked)
U_tracked_y = np.copy(U_tracked)


# Initialize the figure and axis
fig, ax = plt.subplots()

# Create an initial heatmap
x_curr = U_tracked_x[0,:, np.newaxis]
y_curr = U_tracked_y[0,:, np.newaxis]
A = x_curr @ y_curr.T

norm = Normalize(vmin=0, vmax=0.7)

heatmap = ax.imshow(A, cmap='hot', interpolation='nearest', norm=norm)

cbar = fig.colorbar(heatmap, ax=ax)
cbar.set_label('Concentration of Gas')

ax.set_title('Gas Containment Failure')


# Function to update the heatmap
def update(data):
    heatmap.set_data(data)
    return heatmap,

# Function to generate the next matrix (modify as per your update logic)
def generate_data():
    for i in range(Nt-1):
        x_curr = U_tracked_x[i+1,:, np.newaxis]
        y_curr = U_tracked_y[i+1,:, np.newaxis]
        new_A = x_curr @ y_curr.T
        yield new_A


# Create the animation
ani = animation.FuncAnimation(fig, update, generate_data, blit=True, interval=100, repeat=True)

plt.show()