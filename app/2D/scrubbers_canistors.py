import numpy as np


class GasCan2D:
    def __init__(self, x_loc, y_loc, radius, concentration):
        self.concentration = concentration        # Initial concentration in the canister
        self.x_loc = x_loc                        # x-coordinate of the canister's center (in meters)
        self.y_loc = y_loc                        # y-coordinate of the canister's center (in meters)
        self.radius = radius                      # Radius of the canister


    def get_initial_concentration(self, x, y):
        # Create a meshgrid of the x-y plane to apply masks over
        X, Y = np.meshgrid(x, y)

        Y = np.flipud(Y)

        # Create a boolean mask for equation of a circle, i.e the equation of the points where canister sits on the xy-plane
        mask = (X - self.x_loc)**2 + (Y - self.y_loc)**2 <= self.radius**2

        # Create a grid to hold concentration values
        concentration_grid = np.zeros_like(X, dtype=float)

        # Set the concentration in the masked area
        concentration_grid[mask] = self.concentration

        # Returns a grid which is essentially an Nx_points * 
        return concentration_grid


class Scrubber2D():
    def __init__(self, x_loc, y_loc, radius, efficiency):
        self.x_loc = x_loc                          # Central x location of the scrubber
        self.y_loc = y_loc                          # Central y location of the scrubber
        self.radius = radius                        # Radius of the scrubber
        self.efficiency = efficiency                # Efficiency of the scrubber to absorb gas

    def get_affected_indices(self, x, y):
        # Create a meshgrid of the x-y plane to apply masks over
        X, Y = np.meshgrid(x, y)

        Y = np.flipud(Y)

        # Create a boolean mask for equation of a circle, i.e the equation of the points where canister sits on the xy-plane
        mask = (X - self.x_loc)**2 + (Y - self.y_loc)**2 <= self.radius**2
    
        return mask

    def sink(self, x, y):
        exp_term = np.exp(-((x - self.x_loc)**2 + (y - self.y_loc)**2)/(2*self.radius**2))
        return self.efficiency * exp_term