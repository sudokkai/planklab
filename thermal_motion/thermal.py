import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import rayleigh

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
#num_of_image = len(x1[0])
image_per_sec = 2
t_interval = 1/image_per_sec  # each image is taken 0.5 seconds apart

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
scale = 0.1204 * 1e-6
b_diameter = 1.9 #um (micrometers)
b_err = 0.1
visc = 1.00 #centipoise, at 20oC but decreases 2% with each degree increase in temperature.
visc_err = 0.05 #centipoise
temp = 296.5 #kelvin
temp_err = 0.5 #kelvin

#creating steps
def step_distance(x, y):
    temp=[]
    for i in range(len(x)-1):
        temp.append(np.sqrt((x[i] - x[i+1])**2 + (y[i] - y[i+1])**2))
    return temp


#need to fix so it works for each csv, or its run for each csv
combined_csv = np.vstack(((x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5), (x6, y6), (x7, y7), (x8, y8), (x9, y9), (x10, y10)))
step1 = step_distance(x1, y1)
step2 = step_distance(x2, y2)
step3 = step_distance(x3, y3)
step4 = step_distance(x4, y4)
step5 = step_distance(x5, y5)
step6 = step_distance(x6, y6)
step7 = step_distance(x7, y7)
step8 = step_distance(x8, y8)
step9 = step_distance(x9, y9)
step10 = step_distance(x10, y10)
r_values = np.vstack((step1, step2, step3, step4, step5, step6, step7, step8, step9, step10))
print(r_values)
#print(combined_csv)
#r_values=[]
#for n in range(10):
    #for i in range(119):
    #    r_values.append(step_distance(combined_csv[n][i], combined_csv[n][i+1], combined_csv[n+1][i], combined_csv[n+1][i+1]))
#print(len(r_values))

r_values = np.array(r_values)
r_mean = np.mean(r_values)
#print(r_values)
err = np.std(r_values)
mean_err = err / np.sqrt(10)

#model for distribution, step distances

def stepmodel(r, D):
    return ((r)/(2*D*0.5)) * np.exp((-(r**2))/(4*D*0.5))

#spiderman result
loc, scale = rayleigh.fit(r_values, floc=0)
fitted_pdf = rayleigh.pdf(r_values, loc=loc, scale=scale)

t = np.linspace(0, 119, 119)    
r=np.linspace(0,120,100)

popt1, pcov1 = curve_fit(stepmodel, t, r_values[0],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value1, pcov_val1 = popt1[0], pcov1
popt2, pcov2 = curve_fit(stepmodel, t, r_values[1],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value2, pcov_val2 = popt2[0], pcov2
popt3, pcov3 = curve_fit(stepmodel, t, r_values[2],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value3, pcov_val3 = popt3[0], pcov3
popt4, pcov4 = curve_fit(stepmodel, t, r_values[3],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value4, pcov_val4 = popt4[0], pcov4
popt5, pcov5 = curve_fit(stepmodel, t, r_values[4],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value5, pcov_val5 = popt5[0], pcov5
popt6, pcov6 = curve_fit(stepmodel, t, r_values[5],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value6, pcov_val6 = popt6[0], pcov6
popt7, pcov7 = curve_fit(stepmodel, t, r_values[6],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value7, pcov_val7 = popt7[0], pcov7
popt8, pcov8 = curve_fit(stepmodel, t, r_values[7],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value8, pcov_val8 = popt8[0], pcov8
popt9, pcov9 = curve_fit(stepmodel, t, r_values[8],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value9, pcov_val9 = popt9[0], pcov9
popt10, pcov10 = curve_fit(stepmodel, t, r_values[9],  sigma=mean_err, absolute_sigma=True) #fix sigma
D_value10, pcov_val10 = popt10[0], pcov10

D_val_average = np.mean([D_value1, D_value2, D_value3, D_value4, D_value5, D_value6, D_value7, D_value8, D_value9, D_value10])
r_values_comb = np.concatenate((step1, step2, step3, step4, step5, step6, step7, step8, step9, step10))
plt.hist(r_values_comb, bins=30, density=True, alpha=0.6, color='b', label='Step Distance Histogram')
plt.plot(r, stepmodel(r, D_val_average), 'r-', label='Fitted Rayleigh Distribution1',color='red')
#plt.plot(r, stepmodel(r, D_value2), 'r-', label='Fitted Rayleigh Distribution2',color='blue')
#plt.plot(r, stepmodel(r, D_value3), 'r-', label='Fitted Rayleigh Distribution3',color='green')
#plt.plot(r, stepmodel(r, D_value4), 'r-', label='Fitted Rayleigh Distribution4',color='orange')
#plt.plot(r, stepmodel(r, D_value5), 'r-', label='Fitted Rayleigh Distribution5',color='purple')
#plt.plot(r, stepmodel(r, D_value6), 'r-', label='Fitted Rayleigh Distribution6',color='brown')
#plt.plot(r, stepmodel(r, D_value7), 'r-', label='Fitted Rayleigh Distribution7',color='pink')
#plt.plot(r, stepmodel(r, D_value8), 'r-', label='Fitted Rayleigh Distribution8',color='gray')
#plt.plot(r, stepmodel(r, D_value9), 'r-', label='Fitted Rayleigh Distribution9',color='olive')
#plt.plot(r, stepmodel(r, D_value10), 'r-', label='Fitted Rayleigh Distribution10',color='cyan')
#plt.plot(r_values, fitted_pdf, 'r-', label='Fitted Rayleigh Distribution') LOL
plt.xlabel('Step Distance (fix)')
plt.ylabel('# of Steps')
plt.title('Distribution of Step Distances')
plt.legend()
plt.show()


print("D values: ", D_value1*scale, D_value2*scale, D_value3*scale, D_value4*scale, D_value5*scale, D_value6*scale, D_value7*scale, D_value8*scale, D_value9*scale, D_value10*scale)
#print(combined_csv, r_values, D_value)

