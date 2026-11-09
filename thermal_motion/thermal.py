import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Loading data
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

# Facts from the experimantal procedure
scale = 0.1204 * 1e-6               # m/pixel (0.1204 um/px)
fps = 2.0                           # frames per second (2 images/sec)
dt = 1.0 / fps                      # time per frame (0.5 s)
num_frames = 120                    # 120 frames total
num_trials = 10                     # 10 particles; 1 csv for each
time_intervals = np.arange(0, 60, 0.5, dtype=None)                      # interval of time for each frame
fitted_time = np.linspace(time_intervals[0], time_intervals[-1], 100)   # fitted time for graphing the plot (not necessary ?..)

# Known constants
b_err = 0.1                         # idk this one
visc = 1.00                         # centipoise, at 20oC but decreases 2% with each degree increase in temperature.
visc_err = 0.05                     # centipoise
temp = 296.5                        # kelvin
temp_err = 0.5                      # kelvin
radius = (1.9e-6) / 2.0             # bead radius r = 0.95 um
diameter = 1.9e-6                   # meters (1.9 um)
diameter_err = 0.1e-6               # meters (0.1 um uncertainty)
eta = 0.001                         # viscosity of water in Pa.s (1 cP at 20°C)
eta_err = 0.00005                   # Pa.s (0.05 cP uncertainty based on 5% in manual)
gamma = 6 * np.pi * eta * radius    # Stokes drag coefficient

# Function to square the relative postion (with the initial position being the first data point) squared of each data point
def r_squared(x, y):
    r_sq_trials = []
    x_m = x * scale
    y_m = y * scale
    r_sq = 0
    for i in range(len(x)):
        # Convert pixels to meters
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
# Mean squared displacement across trials using previouly defined array
R_sq_mean = np.mean(r,axis=0)

# Calculating the standard error of the mean for each time step
# ddof=1 gives sample standard deviation, so divide by sqrt(num_trials) for standard error after
std_dev = np.std(r, axis=0, ddof=1)
std_dev_mean = std_dev / np.sqrt(num_trials)

# propagating uncertainty (currently just from the SEM)
total_err = std_dev_mean

# Prevent division by zero at t = 0 by replacing zero errors with a small minimum uncertainty (1e-15)
actual_err_fixed = np.where(total_err == 0, 1e-15, total_err)

# 4Dt model for r^2
def Rsquared_model(t, D):
    return 4 * D * t

# fitting the curve for D
poptR, pcovR = curve_fit(Rsquared_model, time_intervals, R_sq_mean, sigma=actual_err_fixed, absolute_sigma=True)
D_opt = poptR[0]
D_unc = np.sqrt(pcovR[0][0])
k = (D_opt * gamma) / temp
k_unc = (D_unc * gamma) / temp              # Uncertainty from just the curve

# Fractional Error Propagation - expand on this more and why it is the way it is (maybe in the report mainly but here too a bit)
frac_err_D = D_unc / D_opt
frac_err_d = diameter_err / diameter
frac_err_eta = eta_err / eta
frac_err_T = temp_err / temp

# Combined fractional uncertainty
frac_err_k = np.sqrt(frac_err_D**2 + frac_err_d**2 + frac_err_eta**2 + frac_err_T**2)

# Total uncertainty in k from all sources
k_unc_total = k * frac_err_k
accepted_k = 1.380649e-23

# Calculating the reduced chi squared for poptR fit
chiSquared = 0
N = len(time_intervals)
m = 2 
v = N - m
for i in range(N):
    chiSquared += ( (R_sq_mean[i] - Rsquared_model(time_intervals[i], D_opt)) / actual_err_fixed[i] )**2
reduced_chiSquared = chiSquared / v
print("Reduced Chi Squared is: ", reduced_chiSquared)

# Calculating residuals = real values - predicted values
predicted_value = Rsquared_model(time_intervals, * poptR)
residuals = R_sq_mean - predicted_value

print("Calcualted k value: ", k, "+- ", k_unc, "\nMax: ", k + k_unc, "\nMin: ", k - k_unc, "\nAccepted k: ", accepted_k)

# Plot of the r^2(t) = 4Dt relation
plt.figure(figsize=(8, 6))
plt.errorbar(time_intervals, R_sq_mean * 1e12, yerr=actual_err_fixed * 1e12, fmt='o', color='blue', 
                label='Averaged Data ± SEM')
plt.plot(fitted_time, Rsquared_model(fitted_time, D_opt) * 1e12, color='pink', 
                label=f'Fit: 4Dt (D={D_opt:.3e})')
plt.xlabel('Time (s)')
plt.ylabel(r'Mean Squared Distance $\langle r^2 \rangle$ ($\mu\text{m}^2$)')
plt.title('Mean-Squared Distance vs. Time')
plt.legend()

# Plot of the residuals
plt.figure(figsize=(8, 6))
plt.scatter(time_intervals, residuals * 1e12, color='blue', label='Time (s)')
plt.errorbar(time_intervals, residuals * 1e12, yerr=actual_err_fixed * 1e12,fmt='none', color='blue', # accounting for conversion from m to um to match other plot
                label=r'Residuals ($\mu\text{m}^2$)')
plt.axhline(0, color='black', linestyle='--', linewidth=1)   # zero-reference line
plt.title("Plot of the Residuals", fontsize=14) 
plt.figtext(0.01, 0.01, "Reduced Chi Squared: "+str(reduced_chiSquared), fontsize=12, style="italic", color="gray",)
plt.xlabel("Time (s)")
plt.ylabel(r'Residuals ($\mu\text{m}^2$)')
plt.legend()
plt.show()