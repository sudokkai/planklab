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

