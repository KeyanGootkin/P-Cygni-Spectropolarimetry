import P_Cyg as pc
import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
from astropy.io import fits,ascii
#Grab all the fits files seperated by data type
datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
bfiles = glob(datadir+'hpolccd*b_hw.fits')
rfiles = glob(datadir+'hpolccd*r_hw.fits')

#first do the red redfiles
rtimes,rpol,rerr=pc.ccd_ave_pol_curve(rfiles)
pc.make_figure(rtimes,rpol,rerr,"Average Polarizaiton (red)")
plt.show()
rperiod,rpower=pc.hybrid_periodogram(rtimes,rpol,rerr)
pc.make_periodogram_figure(rperiod,rpower)
plt.show()
best_rperiod= pc.find_best_period(rperiod,rpower)
pc.plot_phase_time(rtimes,rpol,best_rperiod)
plt.show()


#then the blue
btimes,bpol,berr=pc.ccd_ave_pol_curve(bfiles)
pc.make_figure(btimes,bpol,berr,"Average Polarizaiton (blue)")
plt.show()
bperiod,bpower=pc.hybrid_periodogram(btimes,bpol,berr)
pc.make_periodogram_figure(bperiod,bpower)
plt.show()
best_bperiod=pc.find_best_period(bperiod,bpower)
pc.plot_phase_time(btimes,bpol,best_bperiod)
plt.show()
