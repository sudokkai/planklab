import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import rayleigh
# Libraries

# Known values
p_sq_thickness = 7.82         # Measured in (mm) with caliper
p_sq_thickness_err = 0.01     # Uncertainty from the caliper (** CHECK THIS **)

# -- Activity 1 -- (green LED)


# Importing data for the green LED
greenled = np.loadtxt('interferometer/greenled_second_trial.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
g_dist, greenled_N, greenled_err = gre9enled


# Separating the two trials (The mirror distance was reset at the 13th place in the csv since fringes became too hard to see)
trial1_dist = g_dist[1:12]
trial2_dist = g_dist[13:21]

# Displacement (dx) between consecutive readings in absolute value
dx1 = np.abs(np.diff(trial1_dist))
dx2 = np.abs(np.diff(trial2_dist))
delta_x = np.concatenate((dx1, dx2))

# Corresponding fringes counted of each dx
N1 = greenled_N[2:12]
N2 = greenled_N[14:21]
N = np.concatenate((N1, N2))

# and their respective error(s)
err1 = greenled_err[2:12]
err2 = greenled_err[14:21]
N_err = np.concatenate((err1, err2))

# delta_x, N, and N_err are just the lists of values from our data, converted into change in x, change in fringes, and uncertainty respectively

# For plotting purposes only : Making the points be of the cumulative displacement and fringe count so that points don't overlap.
x1_cumulative = np.cumsum(dx1)
x2_cumulative = np.cumsum(dx2)
N1_cumulative = np.cumsum(N1)
N2_cumulative = np.cumsum(N2)

err1_cumulative = np.sqrt(np.cumsum(err1**2))
err2_cumulative = np.sqrt(np.cumsum(err2**2))


# -- Uncertainty in the distance readings: -- 
# Uncertainty in each micrometer postion
sigma_position = 0.5  # Estimating to the nearest 1um right now, but we did measure as if we could distinguish a 0.5um
# displacement = x2 - x1
# independent uncertainties propagated give 
sigma_dx = np.sqrt(2) * sigma_position   # Uncertainty in each displacement
#UNFINISHED


# Wavelength model (N lambda = 2 deltaX) 
# Fitted using a linear model of the form y=mx
def wvl_model(x, m):
    return m*x         

print("Number of: ") # Checking if the lengths match up for curve fitting
print("Displacements: ", len(delta_x))
print("Fringe counts: ", len(N))
print("Uncertainties: ", len(N_err))

# Curve fit for the wavelength
popt1, pcov1 = curve_fit(wvl_model, delta_x, N, sigma=N_err, absolute_sigma=True)

# Getting the optimized values
m_opt = popt1[0]
# Uncertainties = diagonals of pcov1
m_err = np.sqrt( pcov1[0, 0] )

# Calculate wavelength
Lambda_opt = 2/m_opt
Lambda_err = 2*m_err/m_opt**2  # Not sure if I should use this uncertainty propagation (need to verify)

#UNFINISHED
# Calculating residuals = real values - predicted values
predicted_value1 = wvl_model(delta_x, * popt1)
residuals = N - predicted_value1

# Calculating reduced chi squared
chi_squared = np.sum((residuals / N_err)**2)
degrees_freedom = len(N) - len(popt1)
reduced_chi_squared = chi_squared / degrees_freedom


# Printing results
print("Slope: ", m_opt, "+/- ",  m_err, "fringes/um")
print("Wavelength: ", Lambda_opt, "+/- ", Lambda_err, "um")
print("Wavelength: ", Lambda_opt*1000, "+/- ",  Lambda_err*1000, "nm")
print("Residuals:", residuals)
print("Fringe uncertainties: ", N_err)
print("Normalized residuals: ", residuals / N_err)

print("Trial 1:")
print("Displacement:", dx1)
print("Fringes:", N1)

print("\nTrial 2:")
print("Displacement:", dx2)
print("Fringes:", N2)

print("\nFringes per micrometre:")
print(N / delta_x)


# Testing results by fitting each trial independently
popt_trial1, pcov_trial1 = curve_fit(wvl_model, dx1, N1, sigma=err1, absolute_sigma=True)
popt_trial2, pcov_trial2 = curve_fit(wvl_model, dx2, N2, sigma=err2, absolute_sigma=True)
# Calculate wavelengths
lambda_trial1 = 2 / popt_trial1[0]
lambda_trial2 = 2 / popt_trial2[0]
print("Trial 1 wavelength:", lambda_trial1 * 1000, "nm")
print("Trial 2 wavelength:", lambda_trial2 * 1000, "nm")

residuals1 = N1 - wvl_model(dx1, *popt_trial1)
residuals2 = N2 - wvl_model(dx2, *popt_trial2)

chi1 = np.sum((residuals1 / err1)**2)
chi2 = np.sum((residuals2 / err2)**2)

red_chi1 = chi1 / (len(N1) - len(popt_trial1))
red_chi2 = chi2 / (len(N2) - len(popt_trial2))

print("Trial 1 reduced chi-squared:", red_chi1)
print("Trial 2 reduced chi-squared:", red_chi2)

# Ploting activity 1 (green LED)

# Creating fitted time so that the curve can be plotted all the way through the data points
fitted_x = np.linspace(0, max(max(x1_cumulative), max(x2_cumulative)), 200)
plt.figure(figsize=(8, 6))
# Plotting trials 1 and 2's data points seperately 
plt.errorbar(
        x1_cumulative,          # Note that this plot uses the CUMULATIVE change in x. This is not directly what the fit is based off of, but it is here so that the data's alingment with the fit can be assessed
        N1_cumulative,
        yerr=err1_cumulative,
        color='blue',
        fmt='o',
        label='Trial 1' )
# Trial 2 CUMULATIVE x
plt.errorbar(
        x2_cumulative,
        N2_cumulative,
        yerr=err2_cumulative,
        color='purple',
        fmt='o',
        label='Trial 2' )
# Plotting the function
plt.plot(
        fitted_x,
        wvl_model(fitted_x, m_opt), # This uses the fit from the non cumulative calculation
        color='green',
        label='Linear fit' )
plt.xlabel(r'Cumulative Mirror Displacement $\Delta x_{\mathrm{total}}$ ($\mu$m)')
plt.ylabel('Cumulative Number of Fringes Disappeared, $N_{\mathrm{total}}$')
plt.title('Green LED: Total Fringe Count vs. Total Mirror Displacement')
plt.legend()
plt.grid(alpha=0.3)

# Plot of the residuals
plt.figure(figsize=(8, 6))
plt.errorbar(
        delta_x,
        residuals,
        yerr=N_err,
        fmt='o',
        color='blue',
        capsize=3,
        label='Residuals ± uncertainty' )
plt.axhline(0, color='black', linestyle='--', linewidth=1)   # zero-reference line
plt.title("Green LED: Residuals of Linear Fit")
plt.figtext(0.01, 0.01, "Reduced Chi Squared: "+str(reduced_chi_squared), fontsize=12, style="italic", color="gray",)  
plt.xlabel(r'Mirror Displacement Interval $\Delta x$ ($\mu$m)')
plt.ylabel('Residuals (fringes)')
plt.legend()
plt.grid(alpha=0.3)



plt.show()