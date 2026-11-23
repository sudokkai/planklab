import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import rayleigh

# Importing data for the green LED
greenled = np.loadtxt('interferometer/greenled_second_trial.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
g_dist, g_N, g_err = greenled
delta_x = []
for i in range (len(g_dist)):
    delta_x.append( g_dist[i] - g_dist[i-1] )
print(delta_x)
print(len(delta_x))
print ("yo")
# Known values
plastic_square_length = 50.06           #mm
plastic_square_thickness = 7.82         #mm

# -- Activity 1 -- (green LED)

# Model for thw wavelength of the green LED (N = 2 delta_x / lamda)
def wvl_model (x, lamda):
    return (2 * x) / lamda

# Curve fit parameters
popt1, pcov1 = curve_fit(wvl_model, delta_x, g_N, sigma=g_err, absolute_sigma=True)

fitted_N = popt1