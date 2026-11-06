import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data1 = np.loadtxt('thermal_motion/1.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
data2 = np.loadtxt('thermal_motion/2.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
data3 = np.loadtxt('thermal_motion/3.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
data4 = np.loadtxt('thermal_motion/4.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
data5 = np.loadtxt('thermal_motion/5.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
data6 = np.loadtxt('thermal_motion/6.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
data7 = np.loadtxt('thermal_motion/7.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
data8 = np.loadtxt('thermal_motion/8.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
data9 = np.loadtxt('thermal_motion/9.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
data10 = np.loadtxt('thermal_motion/10.csv', delimiter=',', comments="#", unpack=True, skiprows=1)

#known constants
b_diameter = 1.9 #um (micrometers)
b_err = 0.1
visc = 1.00 #centipoise, at 20oC but decreases 2% with each degree increase in temperature.
visc_err = 0.05 #centipoise
temp = 296.5 #kelvin
temp_err = 0.5 #kelvin

# creating model for <r^2> = slope_D * 6t
def model(t, A, V0):
    return A / (voltage - V0)

# using curve_fit to generate slope and intercept values from our data
# taking wavelength as the x axis and voltage average values as the y, we take the wavelength uncertanties since they're much higher
popt, pcov = curve_fit(model, voltage_avgs, wavelength_values, sigma=wavelength_uncertainties, absolute_sigma=True)
slope_opt, intercept_opt = popt

def step_distance(x, y):
    return np.sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2)

r_values = []
