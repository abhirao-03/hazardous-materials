# Hazardous Materials Industrial Mathematics
A Python-based simulation tool for modeling gas diffusion in 2D space with multiple gas sources, scrubbers, and optimal detector placement calculation.

## Overview

This project simulates the diffusion of gas from multiple sources in a 2D space, accounting for gas scrubbers, and calculates optimal positions for gas detectors. It uses numerical methods to solve the diffusion equation and implements sophisticated algorithms for detector placement optimization.


## Project Structure

```txt
app/
├─ 1D/
│  ├─ animate.py
│  ├─ model_params.py
│  ├─ workspace.py
│  ├─ scrubbers_canistors.py
├─ 2D/
│  ├─ animate.py
│  ├─ scrubbers_canistors.py
│  ├─ workspace.py
│  ├─ model_params.py
│  ├─ max_distance_multiple.py
│  ├─ max_distance.py
exploration/
visualisation/
```

For both 1D and 2D the following files are the same:

- `model_params.py`: Model parameters and configuration
- `scrubbers_canisters.py`: Classes for gas canisters and scrubbers
- `workspace.py`: Core simulation functionality
- `animate.py`: Visualization and animation utilities

For 2D the following files implement our detector optimization algorithms:
- `max_distance.py`
- `max_distance_multiple.py`

## Usage

### Animating Gas Failure

To animate gas containment failure simply run `animate.py`. This animation has a gas canistor with radius 0.23 metres, centered in a 4 metres by 4 metres room.

You can change the parameters of the room and the gas canistor by changing `model_params.py` and `workspace.py` respectively.

### Gas Canistors and Scrubbers

You can specify gas canistor(s), and scrubbers using the `GasCan2d` and `Scrubber2D` object class we created. Sample usage below:

```python

# Define gas canisters
gas_canisters = [
    GasCan2D(x_loc=0.9, y_loc=0.9, radius=0.05, concentration=1.0),
    GasCan2D(x_loc=0.1, y_loc=0.9, radius=0.05, concentration=1.0)
]

# Define scrubbers
scrubbers = [
    Scrubber2D(x_loc=0.5, y_loc=0.7, radius=0.08, efficiency=500.0),
    Scrubber2D(x_loc=0.5, y_loc=0.3, radius=0.08, efficiency=500.0)
]
```

### Basic Simulation

```python
import model_params as model
from workspace import run_simulation

# Create model parameters
parameters = model.parameters()

# Run simulation
U_tracked = run_simulation()
```

### Detector Optimization

```python
from max_distance_multiple import max_distance
from scrubbers_canisters import GasCan2D, Scrubber2D

# Define gas canisters
gas_canisters = [
    GasCan2D(x_loc=0.9, y_loc=0.9, radius=0.05, concentration=1.0),
    GasCan2D(x_loc=0.1, y_loc=0.9, radius=0.05, concentration=1.0)
]

# Define scrubbers
scrubbers = [
    Scrubber2D(x_loc=0.5, y_loc=0.7, radius=0.08, efficiency=500.0),
    Scrubber2D(x_loc=0.5, y_loc=0.3, radius=0.08, efficiency=500.0)
]

# Run optimization
U_3D = max_distance(parameters, gas_canisters, scrubbers, run_simulation)
```

## Model Parameters

The simulation can be configured through the `parameters` class:

- `len_x`, `len_y`: Domain dimensions
- `Nx_points`, `Ny_points`: Spatial discretization points
- `t_end`: Simulation end time
- `Nt_points`: Number of time steps
- `diffusion_coeff`: Gas diffusion coefficient

## Visualization

The project visualizes the following:
- Real-time concentration heatmaps
- Gas spread boundaries
- Canister and scrubber locations
- Optimal detector positions
- Animated diffusion process

An overview of these can be seen in the "visualisation" notebook.
