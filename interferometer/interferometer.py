import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import rayleigh
# Libraries

# -- Activity 1 -- (green LED)


# Importing data for the green LED
greenled = np.loadtxt('interferometer/greenled_second_trial.csv', delimiter=',', comments="#", unpack=True, skiprows=1)
g_dist, greenled_N, greenled_err = greenled


# Separating the two trials (The mirror distance was reset at the 13th place in the csv since fringes became too hard to see)
# Also, appending the starting point of 50um to each trial
trial1_dist = np.concatenate(([50], g_dist[1:12]))
trial2_dist = np.concatenate(([50], g_dist[13:21]))

# Displacement (dx) between consecutive readings
dx1 = np.abs(np.diff(trial1_dist))
dx2 = np.abs(np.diff(trial2_dist))

# Corresponding fringe counts of each dx
N1 = greenled_N[1:12]
N2 = greenled_N[13:21]
# and their respective error(s)
err1 = greenled_err[1:12]
err2 = greenled_err[13:21]

# Cumulative displacement and fringe counts
# Used to graph the results so that it is interpretable
x1 = np.cumsum(dx1) 
x2 = np.cumsum(dx2)

N1_cumulative = np.cumsum(N1)
N2_cumulative = np.cumsum(N2)

# Combining the two trials into one list
delta_x = np.concatenate((x1, x2))
N_cumulative = np.concatenate((N1_cumulative, N2_cumulative))

# Cumulative uncertainties
N_err = np.concatenate((
    np.sqrt(np.cumsum(err1**2)),
    np.sqrt(np.cumsum(err2**2))
))

# Known values
p_sq_thickness = 7.82         # Measured in (mm) with caliper
p_sq_thickness_err = 0.01     # Uncertainty from the caliper (** CHECK THIS **)


# Wavelength model | N lambda = 2 deltaX
# Using linear model, with lambda = 2 / m
def wvl_model(x, m, b):
    return m*x + b

print("Displacement:", len(delta_x))
print("Fringe counts:", len(N_cumulative))
print("Uncertainties:", len(N_err))

# Curve fit for the wavelength
popt1, pcov1 = curve_fit(wvl_model, delta_x, N_cumulative, sigma=N_err, absolute_sigma=True)

# Getting the optimized values
m_opt, b_opt = popt1
m_err, b_err = np.sqrt(np.diag(pcov1))

# Calculate wavelength
Lambda_opt = 2/m_opt
Lambda_err = 2*m_err/m_opt**2

# Printing results
print(f"Slope: {m_opt:.4f} +/- {m_err:.4f} fringes/um")
print(f"Intercept: {b_opt:.2f} +/- {b_err:.2f} fringes")
print(f"Wavelength: {Lambda_opt:.4f} +/- {Lambda_err:.4f} um")
print(f"Wavelength: {Lambda_opt*1000:.1f} +/- {Lambda_err*1000:.1f} nm")

# Ploting activity 1 (green LED)

# Creating fitted time so that the curve can be plotted smoothly
fitted_x = np.linspace(0, max(delta_x), 200)

plt.figure(figsize=(8, 6))
plt.errorbar(
    delta_x,
    N_cumulative,
    yerr=N_err,
    fmt='o',
    color='blue',
    capsize=3,
    label='Experimental data'
)

plt.plot(
    fitted_x,
    wvl_model(fitted_x, *popt1),
    color='pink',
    label='Weighted linear fit'
)

plt.xlabel(r'Cumulative Mirror Displacement $\Delta x$ ($\mu$m)')
plt.ylabel('Cumulative Number of Fringes, N')
plt.title('Green LED: Fringe Count vs. Mirror Displacement')
plt.legend()
plt.grid(alpha=0.3)

plt.show()