class model_parameters():
    """
    Defines the parameters for our model

    Args:
        len_x          Length of x interval. Doesn't matter since we can non-dimesionalise to get 1.0
        len_y          '''
        Nx_points      Number of x spatial points
        Ny_points      Number of y spatial points
        t_end          End of time horizon
        t_num_steps    Number of t temporal points

    Output:
        model_parameters Object
            Containing all the defined parameters
    """
    def __init__(self, len_x=1.0, len_y=1.0, Nx_points=100, Ny_points=100, t_end=10.0, t_num_steps=1000):
        self.len_x = len_x
        self.len_y = len_y
        self.Nx_points = Nx_points
        self.Ny_points = Ny_points
        self.dx = self.len_x / self.Nx_points
        self.dy = self.len_y / self.Ny_points

        self.t_end = t_end
        self.num_steps = t_num_steps
        self.dt = self.t_end / self.num_steps

        self.Cx = self.dt / (self.dx ** 2)
        self.Cy = self.dt / (self.dt ** 2)