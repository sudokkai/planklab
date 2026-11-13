import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Loading data from csv
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

# Known constants
temp = 296.5                        # Kelvin
temp_err = 0.5                      # Kelvin
radius = (1.9e-6) / 2.0             # Bead radius (0.95 um)
radius_err = (0.1e-6) / 2           # Meters (0.05 um)
diameter = 1.9e-6                   # Meters (1.9 um)
diameter_err = 0.1e-6               # Meters (0.1 um uncertainty)
eta = 0.001                         # Viscosity of water in Pa.s (1 cP at 20°C)
eta_err = 0.00005                   # Pa.s (0.05 cP uncertainty based on 5% in manual)
gamma = 6 * np.pi * eta * radius    # Stokes drag coefficient

# -- Method 1 -- ( <r^2 (t)> = 4Dt)

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

# 4Dt model for r^2
def Rsquared_model(t, D):
    return 4 * D * t

# Sclicing out the first data point, since it will always be (0,0) with 0 error, causing issues with the curve_fit accuracy
t_fit = time_intervals[1:]
r_fit = R_sq_mean[1:]
err_fit = std_dev_mean[1:]

# fitting the curve for D
poptR, pcovR = curve_fit(Rsquared_model, t_fit, r_fit, sigma=err_fit, absolute_sigma=True)
D_opt = poptR[0]
D_unc = np.sqrt(pcovR[0][0])
k = (D_opt * gamma) / temp
k_unc = (D_unc * gamma) / temp              # Uncertainty from just the curve
print("Uncertainty from ONLY SEM: ", k_unc)

# Calculating the reduced chi squared for the poptR fit
chiSquared = 0
N = len(t_fit)
m = 1
v = N - m
for i in range(N):
    chiSquared += ( (r_fit[i] - Rsquared_model(t_fit[i], D_opt)) / err_fit[i] )**2
reduced_chiSquared = chiSquared / v

# Function is k = (6 pi eta radius D) / temp
# Finding the squared partial derivaties * their squared uncertainties for error propagation
dk_dD = np.abs(( (6 * np.pi * eta * radius) / temp )**2 * D_unc**2)
dk_dr = np.abs(( (6 * np.pi * eta * D_opt ) / temp )**2 * radius_err**2)
dk_deta = np.abs(( (6 * np.pi * radius * D_opt) / temp )**2 * eta_err**2)
dk_dT = np.abs(-( ( (6 * np.pi * eta * radius * D_opt) / (temp**2) ) )**2 * temp_err**2)

# Total uncertainty in k from all sources
k_unc_total = np.sqrt(dk_dD + dk_dr + dk_deta + dk_dT)
accepted_k = 1.380649e-23

# Displaying results
print(" --- Method 1: ---\nCalcualted k value: ", k, "+- ", k_unc_total, "\nMax: ", k + k_unc_total, "\nMin: ", k - k_unc_total, "\nAccepted k: ", accepted_k,
      "\nReduced Chi Squared is: ", reduced_chiSquared)

# Calculating residuals = real values - predicted values
predicted_value = Rsquared_model(t_fit, * poptR)
residuals = r_fit - predicted_value

# -- Method 2 -- sum ri^2 = NDdt

# dt is step size (not total time)
# ri is the step distance between one t at some arbitrary time
# sum ri^2 is the sum of the squares of those for each trial at the same time 

# Defining function to find the squared step distance of each point and sum them
def find_D(x, y):
    total = 0
    x_m = x * scale
    y_m = y * scale
    for i in range(len(x) - 1):
        ri_sq =  (x_m[i] - x_m[i+1])**2 + (y_m[i] - y_m[i+1])**2
        total += ri_sq
    return   total / (dt * num_frames)

# Calling the function for each particle
D1 = find_D(x1, y1)
D2 = find_D(x2, y2)
D3 = find_D(x3, y3)
D4 = find_D(x4, y4)
D5 = find_D(x5, y5)
D6 = find_D(x6, y6)
D7 = find_D(x7, y7)
D8 = find_D(x8, y8)
D9 = find_D(x9, y9)
D10 = find_D(x10, y10)

# Calculating the mean
mean_D = (D1 + D2 + D3 + D4 + D5 + D6 + D7 + D8 + D9 + D10) / num_trials
all_Ds = np.array([D1, D2, D3, D4, D5, D6, D7, D8, D9, D10])

# Finding the standard error between all values 
unc_mean_D =  (np.std(all_Ds, axis=0, ddof=1)) / np.sqrt(num_trials)

# Converting from D to k
k_2 = (mean_D * gamma) / temp
k_unc_2 = (unc_mean_D * gamma) / temp              
# Displaying results 
print(" --- Method 2: ---\nk = ", k_2, "+- ", k_unc_2,
      "\nMax: ", k_2 + k_unc_2,
      "\nMin: ", k_2 - k_unc_2)


# -- Plotting --

# Plot of Method 1:
fitted_time = np.linspace(t_fit[0], t_fit[-1], 100)                               # fitted time for graphing the plot (not necessary since its linear?..)
plt.figure(figsize=(8, 6))
plt.errorbar(t_fit, r_fit * 1e12, yerr=err_fit * 1e12, fmt='o', color='blue',     # Scaling the R values by 1e12 so that they match graph units (um)
                label='Averaged Data ± SEM')
plt.plot(fitted_time, Rsquared_model(fitted_time, D_opt) * 1e12, color='pink', 
                label=f'Fit: 4Dt (D={D_opt:.3e})')
plt.xlabel('Time (s)')
plt.ylabel(r'Mean Squared Distance $\langle r^2 \rangle$ ($\mu\text{m}^2$)')
plt.title('Mean-Squared Distance vs. Time')
plt.legend()

# Plot of the residuals
plt.figure(figsize=(8, 6))
plt.scatter(t_fit, residuals * 1e12, color='blue', label='Time (s)')
plt.errorbar(t_fit, residuals * 1e12, yerr=err_fit * 1e12,fmt='none', color='blue',   # Accounting for conversion from m to um to match other plot
                label=r'Residuals ($\mu\text{m}^2$)')
plt.axhline(0, color='black', linestyle='--', linewidth=1)   # zero-reference line
plt.title("Plot of the Residuals", fontsize=14) 
plt.figtext(0.01, 0.01, "Reduced Chi Squared: "+str(reduced_chiSquared), fontsize=12, style="italic", color="gray",)  
plt.xlabel("Time (s)")
plt.ylabel(r'Residuals ($\mu\text{m}^2$)')
plt.legend()

# Plot of Method 3: 




plt.show()

