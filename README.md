# Hazardous Materials Industrial Mathematics
A Python-based simulation tool for modeling gas diffusion in 2D space with multiple gas sources, scrubbers, and optimal detector placement calculation.

## Overview

This project simulates the diffusion of gas from multiple sources in a 2D space, accounting for gas scrubbers, and calculates optimal positions for gas detectors. It uses numerical methods to solve the diffusion equation and implements sophisticated algorithms for detector placement optimization.


## Project Structure


- `model_params.py`: Model parameters and configuration
- `scrubbers_canisters.py`: Classes for gas canisters and scrubbers
- `workspace.py`: Core simulation functionality
- `animate.py`: Visualization and animation utilities
- `max_distance_multiple.py`: Detector optimization algorithms

## Usage

### Gas Canistors and Scrubbers

You can specify gas canistor(s), and scrubbers using the `GasCan2d` and `Scrubber2D` object class we created. Sample usage below:

```python

# Define gas canisters
gas_canisters = [
    GasCan2D(x_loc=0.9, y_loc=0.9, radius=0.05, concentration=1.0),
    GasCan2D(x_loc=0.1, y_loc=0.9, radius=0.05, concentration=1.0)
]

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

The project includes several visualization features:
- Real-time concentration heatmaps
- Gas spread boundaries
- Canister and scrubber locations
- Optimal detector positions
- Animated diffusion process