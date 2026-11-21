import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import rayleigh

# Importing data for the green LED
greenled = np.loadtxt('interferometer/greenled_second_trial.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
gd_dist, g_N, g_err = greenled

for i in gd_dist:
    

# Known values
plastic_square_length = 50.06           #mm
plastic_square_thickness = 7.82         #mm

