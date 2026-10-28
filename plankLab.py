import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# importing experimentally collected data from csv 
data = np.loadtxt('plankLab.csv', delimiter=',', comments="#", unpack=True, skiprows=1)

# taking csv data into named variables
voltage_data = data[0]
wavelength_nm = data[1]
err_data = data[2]
data = np.array(data)

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
wavelength_uncertainties = np.array(wavelength_uncertainties)

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

# uncertainty from the multimeter. this could be done more accurately if it was calculated for each individual voltage, rather than 
# the averages i suppose. can change later
meter_limit = (voltage_avgs * 0.0005) + 0.0002

# treating the manufacturer's accuracy limit as a rectangular distribution, we divide by sqrt(3)
err_meter = meter_limit / np.sqrt(3) 

# propogating uncertainty from multimeter and averages
err_per_avg = np.sqrt(err_stat**2 + err_meter**2)

# creating model for wavelength = slope * (voltage - intercept)
def model(voltage, A, V0):
    return A / (voltage - V0)

# using curve_fit to generate slope and intercept values from our data
# taking wavelength as the x axis and voltage average values as the y, we take the wavelength uncertanties since they're much higher
popt, pcov = curve_fit(model, voltage_avgs, wavelength_values, sigma=wavelength_uncertainties, absolute_sigma=True)
slope_opt, intercept_opt = popt
slope_uncertainty = np.sqrt(pcov[0, 0])  # first element in pcov = uncertainty in slope

# propagating uncertainty from the wavelengths and voltages using the partial derivative method
err_wavelength_from_voltage = (np.abs(slope_opt) / (voltage_avgs - intercept_opt)**2) * err_per_avg
err_wavelength_from_voltage = np.array(err_wavelength_from_voltage)
err_effective = np.sqrt(wavelength_uncertainties**2 + err_wavelength_from_voltage**2)

# refiting the curve, now using the propagated uncertainty for sigma
popt_refit, pcov_refit = curve_fit(model, voltage_avgs, wavelength_values, sigma=err_effective, absolute_sigma=True)
slope_opt_refit, intercept_opt_refit = popt_refit
slope_uncertainty_refit = np.sqrt(pcov_refit[0, 0])

# calculating plank's constant using the fitted data
h = slope_opt_refit * (1e-9) * e / c                        # accounting for units by converting from (nm) to (m)
h_uncertainty = slope_uncertainty_refit * (1e-9) * e / c   # this is the actual unertainty of our plank's constant, bc we needed the factors

# calculating the reduced chi squared using the refitted 
chiSquared = 0
N = len(voltage_avgs)
m = 2 
v = N - m
for i in range(N):
    chiSquared += ( (wavelength_values[i] - model(voltage_avgs[i], slope_opt_refit, intercept_opt_refit)) / err_effective[i] )**2
reduced_chiSquared = chiSquared / v
print("Reduced Chi Squared is: ", reduced_chiSquared)

# calculating residuals = real values - predicted values
predicted_value = model(voltage_avgs, * popt_refit)
residuals = wavelength_values - predicted_value

print("Calcualted plank's constant: ", h, "+- ", h_uncertainty, "\nReal plank's constant: ", true_h)
print("Max h = ", h + h_uncertainty, "\nMin h = ", h - h_uncertainty)
# the value is a little different than what is ideal unfortunately. this is likely due to non-ideal leds that have energy loss and 
# turn-on behaviour, and the uncertainty in the wavelength being very high 

# creating a new array that has 100 points inbetween the start and end of our different voltages so that the fitted curve can be smooth
fitted_voltage_avgs = np.linspace(voltage_avgs[0], voltage_avgs[-1], 100)

# for plotting the theoretical h value alongside our calculated one
k_theoretical = true_h * c / e * 1e9  # units: V nm
# note that this uses the same intercept at the fitted function so that the comparison between slopes can be isolated
theoretical_wavelengths = k_theoretical / (fitted_voltage_avgs - intercept_opt_refit)


# plotting the curve of voltage (V) vs wavelength (nm)
plt.figure(1)
plt.plot(fitted_voltage_avgs, model(fitted_voltage_avgs, slope_opt_refit, intercept_opt_refit), color='pink', 
                label='Voltage (V) vs Wavelength (nm) where λ = hc / e (V - V0)')
plt.plot(fitted_voltage_avgs, theoretical_wavelengths, color='gray', linestyle='--',
                label="Theoretical curve (using accepted h)")
plt.scatter(voltage_avgs, wavelength_values, 
                label='Average voltage threshold for each LED to turn on')
plt.errorbar(voltage_avgs, wavelength_values, xerr=err_per_avg, yerr=wavelength_uncertainties, fmt='none', color='red', 
                label="Uncertainty in wavelengths (vertical) and voltage (horizontal)")
plt.xlabel("Voltage (V))")
plt.ylabel("Wavelength (nm)")
plt.title("Plotting Wavelength (nm) of LEDs vs the turn on Voltage (V)")
plt.figtext(0.01, 0.01, "Experimental Planck's Constant: "+str(h), fontsize=12, style="italic", color="gray",)
plt.legend()

# plot of the residuals 
plt.figure(2)
plt.scatter(voltage_avgs, residuals, color='blue', label='Voltages (V)')
plt.errorbar(voltage_avgs, residuals, yerr=wavelength_uncertainties,fmt='none', color='blue', 
                label='Residuals')
plt.axhline(0, color='black', linestyle='--', linewidth=1)   # zero-reference line
plt.title("Plot of the Residuals", fontsize=14) 
plt.figtext(0.01, 0.01, "Reduced Chi Squared: "+str(reduced_chiSquared), fontsize=12, style="italic", color="gray",)
plt.legend()
plt.show()