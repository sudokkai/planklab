import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# importing experimentally collected data from csv 
data = np.loadtxt('plankLab.csv', delimiter=',', comments="#", unpack=True, skiprows=1)

# taking csv data into named variables
voltage_data = data[0]
wavelength_nm = data[1]
err_data = data[2]

# defiining constants
e = 1.602176634e-19       # elementary charge (C)
c = 2.99792458e8          # speed of light (m/s)
true_h = 6.62607015e-34   # real value of plank's constant (m^2 kg / s)

# function to add and average all of the voltages per a single wavelength
def average(x):
    total = 0
    counter = 0
    applicounter = 0
    for i in wavelength_nm:
        if i == x:
            total += voltage_data[counter]
            applicounter += 1
        counter += 1
    average = total / applicounter
    return average

# output a list of all measured voltages that share the same wavelength
def collection(x):
    dummylist = []
    counter = 0
    for i in wavelength_nm:
        if i == x:
            dummylist.append(voltage_data[counter])
        counter += 1
    return dummylist

# function that accepts a list and returns the standard deviation of the mean
def standard_dev_mean (x):
    sigma = 0
    average = 0
    n = len(x)
    for i in x:
        average += i
    average = average / n
    for i in x:
        sigma += (i - average)**2
    sigma = (np.sqrt(sigma / (n-1))) / np.sqrt(n)
    return sigma

# adding all the different wavelength values to an array
wavelength_values = []
for i in wavelength_nm:
    if not i in wavelength_values:
        wavelength_values.append(i)

# appending a list of uncertainties for each of the wavelengths
wavelength_uncertainties = []
for i in wavelength_values:
    wavelength_uncertainties.append(err_data[list(wavelength_nm).index(i)])

# all the average values of voltage per wavelength
voltage_avgs = []
for i in wavelength_values:
    voltage_avgs.append(average(i))

# puts all voltage values for a given wavelength in each slot in the array
grouped_values = []
for i in wavelength_values:
    grouped_values.append(collection(i))

# uncertainty for each average
err_stat = []
for i in grouped_values:
    err_stat.append(standard_dev_mean(i))

# converting into np arrays so that operations like * and ** work
voltage_avgs = np.array(voltage_avgs)
err_stat = np.array(err_stat)

# uncertainty from the multimeter. could do this more accurately i suppose but later
meter_limit = (voltage_avgs * 0.00005) + 0.0002

err_meter = meter_limit         # could add the / np.sqrt(3) for more techincal accuracy but unsure if thats needed

# propogating uncertainty from multimeter and averages
err_per_avg = np.sqrt(err_stat**2 + err_meter**2)

# creating model for voltage = slope * (1 / wavelength) + intercept
def model (wavelength_nm, m, b):
    return m*(1/wavelength_nm) + b

# using curve_fit to generate slope and intercept values from our data
# taking wavelength as the x axis and voltage average values as the y, we take the wavelength uncertanties since they're much higher
popt, pcov = curve_fit(model, voltage_avgs, wavelength_values, sigma=wavelength_uncertainties, absolute_sigma=True)
m_opt, b_opt = popt
slope_uncertainty = np.sqrt(pcov[0, 0])  # first element in pcov = uncertainty in slope

# calculating plank's constant using the fitted data
h = m_opt *(1e-9) * e / c    # accounting for units by converting from (nm) to (m)
h_uncertainty = slope_uncertainty * (1e-9) * e / c

print("Calcualted plank's constant: ", h, "+- ", h_uncertainty, "\nReal plank's constant: ", true_h)
print("Max h = ", h + h_uncertainty, "\nMin h = ", h - h_uncertainty)
# the value is a little different than what is ideal unfortunately. this is likely due to non-ideal leds that have energy loss and turn-on behaviour
# and the uncertainty in the wavelength being very high 

# creating a new array that has 100 points inbetween the start and end of our different wavelengths 
# so that the curve can be smooth
fitted_voltage_avgs = np.linspace(voltage_avgs[0], voltage_avgs[-1], 100)

# plotting
plt.figure(1)
plt.plot(fitted_voltage_avgs, model(fitted_voltage_avgs, m_opt, b_opt), color='pink', label='Graph of stuff')
# scatter plot only using the average values 
plt.scatter(voltage_avgs, wavelength_values, label='Points of stuff')
plt.errorbar(voltage_avgs, wavelength_values, xerr=err_per_avg, yerr=wavelength_uncertainties, fmt='none', color='red', label="Uncertainty in wavelengths")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Voltage (V)")
plt.show()