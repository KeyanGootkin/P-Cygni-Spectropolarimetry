import P_Cyg as pc
import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
from astropy.io import fits,ascii
#Grab all the fits files seperated by data type
datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
ccdbfiles = glob(datadir+'hpolccd*b_hw.fits')
ccdrfiles = glob(datadir+'hpolccd*r_hw.fits')

times = []
ave_err = []
ave_pol = []
#for each blue file, run through each red file and see if they are the same date
for bfile in ccdbfiles:
    for rfile in ccdrfiles:
        if rfile[:len(rfile)-9] == bfile[:len(bfile)-9]:
            hdu_r = fits.open(rfile)
            tablehdu_r = hdu_r[1].header
            infohdu_r = hdu_r[0].header
            ave_r = tablehdu_r['ave-pol']
            ave_err_r = tablehdu_r['ave-erro']
            time_r = infohdu_r['mjd-obs']
            
            hdu_b = fits.open(bfile)
            tablehdu_b = hdu_b[1].header
            infohdu_b = hdu_b[0].header
            ave_b = tablehdu_b['ave-pol']
            ave_err_b = tablehdu_b['ave-erro']
            time_b = infohdu_b['mjd-obs']
            
            ave_pol.append((ave_b*ave_r)/2)
            #Pretty sure this is wrong. sqrt of sum of squares?
            ave_err.append((ave_err_b*ave_err_r)/2)
            times.append((time_r*time_b)/2)
plt.figure(figsize = [10,5])
plt.errorbar(times,ave_pol,yerr=ave_err)

times,ave_pol,ave_err = np.array(times),np.array(ave_pol),np.array(ave_err)

#Periodogram
pol_period,pol_power = pc.hybrid_periodogram(times,ave_pol,ave_err)
plt.figure(figsize = [10,5])
plt.plot(pol_period,pol_power)
plt.show()
#Find best periods
pol_period,pol_power = np.array(pol_period),np.array(pol_power)
power_sort = np.argsort(pol_power)
pol_power = pol_power[power_sort]
pol_period = pol_period[power_sort]
print(pol_period[len(pol_period)-5:])
print(pol_power[len(pol_period)-5:])
#phase fold at best period
period = pol_period[len(pol_period)-1]
phase_time = times%period
plt.scatter(phase_time,ave_pol)










