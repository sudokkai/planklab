import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

x1,  y1  = np.loadtxt('thermal_motion/1.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x2,  y2  = np.loadtxt('thermal_motion/2.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x3,  y3  = np.loadtxt('thermal_motion/3.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x4,  y4  = np.loadtxt('thermal_motion/4.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x5,  y5  = np.loadtxt('thermal_motion/5.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x6,  y6  = np.loadtxt('thermal_motion/6.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x7,  y7  = np.loadtxt('thermal_motion/7.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x8,  y8  = np.loadtxt('thermal_motion/8.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x9,  y9  = np.loadtxt('thermal_motion/9.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
x10, y10 = np.loadtxt('thermal_motion/10.csv', delimiter=',', comments="#", unpack=True, skiprows=1)


# facts from the experimantal procedure
scale = 0.1204 * 1e-6               # m/pixel (0.1204 um/px)
fps = 2.0                           # frames per second (2 images/sec)
dt = 1.0 / fps                      # time per frame (0.5 s)
num_frames = 120                    # 120 frames total
num_trials = 10                     # 10 particles/CSVs

# known constants
b_err = 0.1                         # idk this one
visc = 1.00                         # centipoise, at 20oC but decreases 2% with each degree increase in temperature.
visc_err = 0.05                     # centipoise
temp = 296.5                        # kelvin
temp_err = 0.5                      # kelvin
radius = (1.9e-6) / 2.0             # bead radius r = 0.95 um
eta = 0.001                         # viscosity of water in Pa.s (1 cP at 20°C)
gamma = 6 * np.pi * eta * radius    # stokes drag coefficient



r_sq_trials = 0
def r_squared(x, y):
    x_sq = 0
    y_sq = 0
    r_sq = 0
    for i in range(len(x)):
        # Convert pixels to meters
        x_m = x[i] * scale
        y_m = y[i] * scale
        # Compute squared displacement relative to initial position: (x(t) - x0)^2 + (y(t) - y0)^2
        r_sq = (x_m - x_m[0])**2 + (y_m - y_m[0])**2
        r_sq_trials.append(r_sq)
    return r_sq_trials

r_trial1_sq = r_squared(x1,y1)
print(r_trial1_sq)

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



