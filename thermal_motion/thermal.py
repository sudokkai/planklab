import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

x1, y1 = np.loadtxt('thermal_motion/1.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x2, y2 = np.loadtxt('thermal_motion/2.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x3, y3 = np.loadtxt('thermal_motion/3.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x4, y4 = np.loadtxt('thermal_motion/4.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x5, y5 = np.loadtxt('thermal_motion/5.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x6, y6 = np.loadtxt('thermal_motion/6.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x7, y7 = np.loadtxt('thermal_motion/7.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x8, y8 = np.loadtxt('thermal_motion/8.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x9, y9 = np.loadtxt('thermal_motion/9.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x10, y10 = np.loadtxt('thermal_motion/10.csv', delimiter=',', comments="#", unpack=True, skiprows=1)

num_of_image = len(x1[0])
image_per_sec = 2
t_interval = 1/image_per_sec  # each image is taken 0.5 seconds apart

list = {}
def step_dist(x):
    current = 0
    for i in range (len(x)):
        current = x[i+1] - x[i]
        list.append(current)
    return(current)

def averages(x):
    total = 0
    for i in x:
        total += i
    return total


# function that accepts a list and returns the standard deviation of the mean
#def standard_dev_mean (x):
#    sigma = 0
#    average = 0
#    n = len(x)
#    for i in x:
#        average += i
#    average = average / n
#    for i in x:
#        sigma += (i - average)**2
#    sigma = (np.sqrt(sigma / (n-1))) / np.sqrt(n)
#    return sigma


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
