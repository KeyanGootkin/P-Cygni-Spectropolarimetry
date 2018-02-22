import P_Cyg as pc
import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
from astropy.io import fits,ascii
#Grab all the fits files seperated by data type
datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
retfiles = glob(datadir+'hpolret*.fits')
ccdfiles = glob(datadir+'hpolccd*.fits')
#For every file read in the time, average polarization, and error
'''times = []
average_pol = []
average_err = []
for file in retfiles:
    hdu = fits.open(file)
    infohdu = hdu[0].header
    tablehdu = hdu[1].header
    ave_pol,ave_err = tablehdu["Ave-pol"],tablehdu["Ave-erro"]
    obs_time = infohdu["MJD_OBS"]
    times.append(obs_time)
    average_pol.append(ave_pol)
    average_err.append(ave_err)
#Sort all of these boys because SOME TIMES DON"T PLAY NICE
times,average_pol,average_err = np.array(times),np.array(average_pol),np.array(average_err)
times_sort = np.argsort(times)
times = times[times_sort]
average_pol = average_pol[times_sort]
average_err = average_err[times_sort]'''
times,average_pol,average_err = pc.ave_pol_curve(retfiles)
#Plot average polarization vs. time
plt.figure(figsize = [10,5])
plt.errorbar(times,average_pol,yerr=average_err)
plt.show()
#Periodogram
pol_period,pol_power = pc.hybrid_periodogram(times,average_pol,average_err)
plt.figure(figsize = [10,5])
plt.plot(pol_period,pol_power)
plt.xlim(0,1000)
plt.show()
#Find best periods
pol_period,pol_power = np.array(pol_period),np.array(pol_power)
power_sort = np.argsort(pol_power)
pol_power = pol_power[power_sort]
pol_period = pol_period[power_sort]
print(pol_period[len(pol_period)-5:])
print(pol_power[len(pol_period)-5:])
#phase fold at best period
#period = pol_period[len(pol_period)-1]
period=130.431612
phase_time = times%period
plt.scatter(phase_time,average_pol)
