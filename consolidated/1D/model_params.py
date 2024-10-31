import numpy as np

class parameters():
    def __init__(self, length = 1.0, Nx_points = 100, t_end=10.0, Nt_points = 1000):
        self.length = length
        self.Nx_points = Nx_points
        self.t_end = t_end
        self.Nt_points = Nt_points

        self.dx = self.length / self.Nx_points
        self.dt = self.t_end / self.Nt_points

        self.C = self.dt / self.dx**2