import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import rayleigh
# Libraries

# Known values
p_sq_thickness = 7.82         # Measured in (mm) with caliper
t = p_sq_thickness * 1000     # Putting thickness into (um) to match wavelength(s)
p_sq_thickness_err = 0.01     # Uncertainty from the caliper (** CHECK THIS **)
intial_temp_of_alum = 24.1    # degrees celcius


# -- Activity 1 -- (green LED)


# Importing data for the green LED
greenlaser = np.loadtxt('interferometer/greenled_second_trial.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
g_dist, greenled_N, greenlaser_err = greenlaser




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
err1 = greenlaser_err[2:12]
err2 = greenlaser_err[14:21]
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


# Ploting activity 1 (green laser wavelength)

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


# -- Activity 2 -- (plastic square index of refraction)

# Importing data
act2 = np.loadtxt('interferometer/plastic_square_2.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
angle_deg, angle_min, N_dis, N_dis_err = act2

# Cutting out first data point and making the N points cumulative
N_dis = np.cumsum(N_dis[1:11])
# Same for cumulartive uncertainty except its the squared product
N_dis_err = np.sqrt(np.cumsum(N_dis_err[1:11]**2))
print(N_dis)
print(len(N_dis))

angle_min = angle_min / 60                 # Converting minutes into decimals of degrees
full_angles = angle_deg + angle_min        # Adding the minutes to the degrees
full_angles = full_angles * np.pi / 180    # Converting to radians
reference_angle = full_angles[0]

# Calculate angles relative to the reference angle
theta = full_angles[1:11] - reference_angle

# Quadratic model: N = a * ((theta**2 - 2 * theta * theta_0)
# where theta is the rotation from the 320 degree starting position, and
# theta_0 is the signed angular offset the starting position to perpendicular incidence
def ref_model(theta, a, theta_0):
    return a * (theta**2 - 2*theta*theta_0)

# N = t/lambda theta^2 (1 - 1/n)
# a = t / lambda (1 - 1/n)

# Optimized curve coefficient and uncertainty
popt2, pcov2 = curve_fit(ref_model, theta, N_dis, sigma=N_dis_err, absolute_sigma=True)
a_opt, theta_0 = popt2
a_err, theta_0_err = np.sqrt(np.diag(pcov2))

# Calculate index of refraction
n_opt = 1 / (1 - (a_opt * Lambda_opt / t))

# And error in n (considering fit uncertainty only)
n_err = (n_opt**2) * (Lambda_opt / t) * a_err
print("Fitted coefficient:", a_opt, "+/-", a_err)
print("Index of refraction:", n_opt, "+/-", n_err)

# Calculating residuals = real values - predicted values
predicted_value2 = ref_model(theta, * popt2)
residuals2 = N_dis - predicted_value2

# Calculating reduced chi squared
chi_squared2 = np.sum((residuals2 / N_dis_err)**2)
degrees_freedom2 = len(N_dis) - len(popt2)
reduced_chi_squared2 = chi_squared2 / degrees_freedom2

print("\n-- ACTIVITY 1: GREEN LASER WAVELENGTH --")

print(f"Fitted slope: {m_opt:.5f} ± {m_err:.5f} fringes/um")
print(f"Wavelength: {Lambda_opt*1000:.2f} ± {Lambda_err*1000:.2f} nm")
print(f"Reduced chi-squared: {reduced_chi_squared:.3f}")

print("\n-- ACTIVITY 2: REFRACTIVE INDEX --")   

print(f"Fitted coefficient (a): {a_opt:.2f} ± {a_err:.2f} rad^-2")
print(f"Fitted angular offset (theta_0): {np.degrees(theta_0):.3f} ± {np.degrees(theta_0_err):.3f} degrees")
print(f"Refractive index (n): {n_opt:.4f} ± {n_err:.4f}")
print(f"Reduced chi-squared: {reduced_chi_squared2:.3f}")


# Ploting activity 2 (index of refraction)

# Creating fitted time so that the curve can be plotted all the way through the data points
fitted_x2 = np.linspace(0, max(theta), 200)
plt.figure(figsize=(8, 6))
plt.errorbar(
        theta,          
        N_dis,
        yerr=N_dis_err,
        color='blue',
        fmt='o',
        label='Trial 1' )
plt.plot(
        fitted_x2,
        ref_model(fitted_x2, a_opt, theta_0),
        color='green',
        label='Quadratic fit' )
plt.title('Plastic Square: Cumulative Fringe Count vs. Angular Displacement')
plt.xlabel(r'Angular Displacement from Initial Position, $\theta$ (rad)')
plt.ylabel(r'Cumulative Number of Fringes Disappeared, $N_{\mathrm{total}}$')
plt.legend()
plt.grid(alpha=0.3)

# Plot of the residuals
plt.figure(figsize=(8, 6))
plt.errorbar(
        theta,
        residuals2,
        yerr=N_dis_err,
        fmt='o',
        color='blue',
        capsize=3,
        label='Residuals ± uncertainty' )
plt.axhline(0, color='black', linestyle='--', linewidth=1)   # zero-reference line
plt.title('Plastic Square: Residuals of Quadratic Fit')
plt.figtext(0.01, 0.01, "Reduced Chi Squared: "+str(reduced_chi_squared2), fontsize=12, style="italic", color="gray",)  
plt.xlabel(r'Angular Displacement from Initial Position, $\theta$ (rad)')
plt.ylabel('Residuals in Fringe Count (fringes)')
plt.legend()
plt.grid(alpha=0.3)

# -- Activity 3 -- (thermal expansion coefficient)

# Known values
rod_length = 8.9         # cm
rod_length_err = 0.05    # cm

l0 = 89000
l000


plt.show()