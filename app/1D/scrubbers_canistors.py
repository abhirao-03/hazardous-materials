import numpy as np

class gas_canistor:
    def __init__(self, loc = 0.5, radius = 0.05, concentration = 1.0):
        self.radius = radius
        self.loc = loc
        self.lower_bound = self.loc - self.radius
        self.upper_bound = self.loc + self.radius
        self.concentration = concentration
    
    def get_initial_concentration(self, x):
        return np.where((x >= self.lower_bound) & (x <= self.upper_bound), 1, 0)

class scrubber():
    def __init__(self, location, radius, efficiency):
        self.loc = location
        self.radius = radius
        self.efficiency = efficiency
        self.lower_bound = self.loc - self.radius
        self.upper_bound = self.loc + self.radius

    def sink(self, x, u):
        """
        loc: location of sink
        x: position
        cot: concentration
        """
        
        scrub_loc = np.where((x >= self.lower_bound) & (x <= self.upper_bound))
        u[scrub_loc] = (1 - self.efficiency) * u[scrub_loc]
        
        return u