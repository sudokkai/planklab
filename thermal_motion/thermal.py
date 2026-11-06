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

<<<<<<< HEAD
num_of_image = len(x1[0])
image_per_sec = 2
t_interval = 1/image_per_sec  # each image is taken 0.5 seconds apart
=======

# facts from the experimantal procedure
scale = 0.1204 * 1e-6               # m/pixel (0.1204 um/px)
fps = 2.0                           # frames per second (2 images/sec)
dt = 1.0 / fps                      # time per frame (0.5 s)
num_frames = 120                    # 120 frames total
num_trials = 10                     # 10 particles/CSVs
time_intervals = np.arange(0, 60, 0.5, dtype=None)  # interval of time for each frame

# known constants
b_err = 0.1                         # idk this one
visc = 1.00                         # centipoise, at 20oC but decreases 2% with each degree increase in temperature.
visc_err = 0.05                     # centipoise
temp = 296.5                        # kelvin
temp_err = 0.5                      # kelvin
radius = (1.9e-6) / 2.0             # bead radius r = 0.95 um
eta = 0.001                         # viscosity of water in Pa.s (1 cP at 20°C)
gamma = 6 * np.pi * eta * radius    # stokes drag coefficient

>>>>>>> kai's-branch

# Function to square the relative postion (with the initial position being the first data point) squared of each data point
def r_squared(x, y):
    r_sq_trials = []
    x_sq = 0
    y_sq = 0
    r_sq = 0
    for i in range(len(x)):
        # Convert pixels to meters
        x_m = x * scale
        y_m = y * scale
        # Compute squared displacement relative to initial position: (x(t) - x0)^2 + (y(t) - y0)^2
        r_sq = (x_m[i] - x_m[0])**2 + (y_m[i] - y_m[0])**2
        r_sq_trials.append(r_sq)
    return r_sq_trials

# Calling the function to find all squared positions
r_trial1_sq = r_squared(x1,y1)
r_trial2_sq = r_squared(x2,y2)
r_trial3_sq = r_squared(x3,y3)
r_trial4_sq = r_squared(x4,y4)
r_trial5_sq = r_squared(x5,y5)
r_trial6_sq = r_squared(x6,y6)
r_trial7_sq = r_squared(x7,y7)
r_trial8_sq = r_squared(x8,y8)
r_trial9_sq = r_squared(x9,y9)
r_trial10_sq = r_squared(x10,y10)

# Putting every r^2 into vertical columns in another array
r = np.array([r_trial1_sq, r_trial2_sq, r_trial3_sq, r_trial4_sq, r_trial5_sq, r_trial6_sq, r_trial7_sq, r_trial8_sq, r_trial9_sq, r_trial10_sq])
# Averaging each r^2(t) for each individual t
R_mean = np.mean(r,axis=0)



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

<<<<<<< HEAD
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
for i in range():
=======

plt.scatter(time_intervals, R_mean)
plt.show()
>>>>>>> kai's-branch
