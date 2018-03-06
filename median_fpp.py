import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
from astropy.io import fits, ascii
from astropy.table import Table, Column
import P_Cyg as pc
import scipy.stats

datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
figdir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Figures/'
txt_list = glob(datadir + 'hpolret*.txt')
fits_list = glob(datadir + 'hpolret*.fits')
rfitsfiles = glob(datadir + 'hpolccd*r_hw.fits')
rtxtfiles = glob(datadir + 'hpolccd*r_hw.fits.txt')


median_wavelength, median_flux, median_pol, median_pos, median_err = pc.median_flux_pol_pos(
    rtxtfiles, rfitsfiles, 250, radial_velocity=-8.9)

plt.figure(figsize=[15, 20])
plt.subplot(3, 1, 1)
pc.make_figure(median_wavelength, median_flux, 0, title="Median Flux")
plt.xlim(6200, 7000)
# plt.plot([vert_line,vert_line],[min(median_flux),max(median_flux)])
# plt.xlim(vert_line-300,vert_line+300)
plt.subplot(3, 1, 2)
pc.make_figure(median_wavelength, median_pol,
               median_err, title="Median Polarization")
plt.xlim(6200, 7000)
plt.ylim(0.95, 1.2)
# plt.plot([vert_line,vert_line],[0,2])
# plt.xlim(vert_line-300,vert_line+300)
plt.subplot(3, 1, 3)
pc.make_figure(median_wavelength, median_pos, 0, title="Median Position Angle")
plt.xlim(6200, 7000)
plt.ylim(65, 75)
# plt.plot([vert_line,vert_line],[0,100])
# plt.xlim(vert_line-300,vert_line+300)
plt.savefig(figdir + 'Spectral_Lines/Median_Halpha_ccd.eps', overwrite=True)
plt.show()
