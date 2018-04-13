import P_Cyg as pc
import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
from astropy.io import fits, ascii
from astropy.stats import LombScargle

# Grab all the fits files seperated by data type
datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
bfiles = glob(datadir + 'hpolccd*b_hw.fits')
rfiles = glob(datadir + 'hpolccd*r_hw.fits')

# first do the red redfiles
rtimes, rpol, rerr = pc.ccd_ave_pol_curve(rfiles)
pc.make_figure(rtimes, rpol, rerr, "Average Polarizaiton (red)")
plt.show()
rperiod, rpower, rls = pc.hybrid_periodogram(rtimes, rpol, rerr, give_ls=True)
pc.make_periodogram_figure(rperiod, rpower)
plt.show()
best_rperiod = pc.find_best_period(rperiod, rpower)
pc.plot_phase_time(rtimes, rpol, best_rperiod)
plt.show()


# then the blue
btimes, bpol, berr = pc.ccd_ave_pol_curve(bfiles)
pc.make_figure(btimes, bpol, berr, "Average Polarizaiton (blue)")
plt.show()
bperiod, bpower, bls = pc.hybrid_periodogram(btimes, bpol, berr, give_ls=True)
pc.make_periodogram_figure(bperiod, bpower)
plt.show()
best_bperiod = pc.find_best_period(bperiod, bpower)
pc.plot_phase_time(btimes, bpol, best_bperiod)
plt.show()


rlsperiod, rlspower = rls.autopower()
r_false_alarm = rls.false_alarm_probability(rlspower.max())
blsperiod, blspower = bls.autopower()
b_false_alarm = bls.false_alarm_probability(blspower.max())

print(pc.find_best_period((1 / rlsperiod), rlspower))
print(r_false_alarm * 100)
print(pc.find_best_period((1 / blsperiod), blspower))
print(b_false_alarm * 100)

plt.plot((1 / rlsperiod), rlspower)
plt.show()

plt.plot((1 / blsperiod), blspower)
plt.show()

pc.plot_phase_time(btimes, bpol, pc.find_best_period(
    (1 / blsperiod), blspower))
