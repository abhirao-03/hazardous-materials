import scrubbers_canistors as sg
from model_params import *
from matrix_builder import build_matrix

parameters   = model_parameters()
gas_canistor = sg.GasCan2D(x_loc=.5, y_loc=.5, radius=0.05, concentration=1.0)
scrubbers_2d = [sg.Scrubber2D(x_loc=0.8, y_loc=0.8, radius=0.02, efficiency=0.40, cap=5.0),
                sg.Scrubber2D(x_loc=0.2, y_loc=0.2, radius=0.02, efficiency=0.60, cap=5.0)
                ]

A = build_matrix(parameters=parameters)

