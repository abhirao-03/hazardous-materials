# Hazardous Materials Industrial Mathematics
A Python-based simulation tool for modeling gas diffusion in 2D space with multiple gas sources, scrubbers, and optimal detector placement calculation.

## Overview

This project simulates the diffusion of gas from multiple sources in a 2D space, accounting for gas scrubbers, and calculates optimal positions for gas detectors. It uses numerical methods to solve the diffusion equation and implements sophisticated algorithms for detector placement optimization.

## Features

- 2D gas diffusion simulation
- Multiple gas source (canister) support
- Gas scrubber simulation with configurable efficiency
- Optimal detector placement calculation
- Real-time visualization with animated heatmaps
- Comprehensive analysis of detector coverage zones

## Project Structure

- `model_params.py`: Model parameters and configuration
- `scrubbers_canisters.py`: Classes for gas canisters and scrubbers
- `workspace.py`: Core simulation functionality
- `animate.py`: Visualization and animation utilities
- `max_distance_multiple.py`: Detector optimization algorithms

## Usage

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

## Detection Algorithm

The detector placement algorithm:
1. Calculates 1% concentration threshold zones for each gas canister
2. Finds overlapping zones for multiple canisters
3. Optimizes detector placement to maximize coverage
4. Handles cases where multiple detectors are required

## Output

The optimization process provides:
- Number of required detectors
- Grid coordinates for each detector
- Physical coordinates (length, width) for installation
- Visualization of detector coverage zones

## Notes

- The simulation assumes identical gas canisters (concentration, dimensions)
- Scrubber efficiency can be adjusted to model different absorption rates
- The coordinate system origin is at a room corner
- Grid indexing is modified to match Python's default system
