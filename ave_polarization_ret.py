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

times,pol,pos,ave_err = pc.calculate_average_polarization_ret(retfiles)
pc.make_figure(times,pol,ave_err,"recalculated")
plt.show()
period,power = pc.hybrid_periodogram(times,pol,ave_err)
pc.make_periodogram_figure(period,power)
plt.show()
best_period = pc.find_best_period(period,power)
pc.plot_phase_time(times,pol,best_period)

pc.make_figure(times,pos,ave_err,"position angle")
plt.show()
pos_period,pos_power=pc.hybrid_periodogram(times,pos,ave_err)
pc.make_periodogram_figure(pos_period,pos_power)
plt.show()
best_pos_period = pc.find_best_period(pos_period,pos_power)
pc.plot_phase_time(times,pos,best_pos_period)
plt.show()









