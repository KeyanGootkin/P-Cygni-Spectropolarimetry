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
bfitsfiles = glob(datadir + 'hpolccd*b_hw.fits')
btxtfiles = glob(datadir + 'hpolccd*b_hw.fits.txt')

interesting_wavelengths = [4860, 5875, 6560]
w = interesting_wavelengths[2]

median_wavelength, median_flux, median_pol, median_pos, median_err = pc.median_flux_Q_U(
    rtxtfiles, rfitsfiles, 1000, radial_velocity=-8.9)
paerr = pc.position_angle_error(median_pol, median_err)
plt.figure(figsize=[15, 20])
plt.subplot(3, 1, 1)
pc.make_figure(median_wavelength, median_flux, 0, title="Median Flux")
plt.xlim(w - 200, w + 200)
# plt.plot([vert_line,vert_line],[min(median_flux),max(median_flux)])
# plt.xlim(vert_line-300,vert_line+300)
plt.subplot(3, 1, 2)
pc.make_figure(median_wavelength, median_pol,
               median_err, title="Median Q")
plt.xlim(w - 200, w + 200)
plt.ylim(0.105, 0.55)
# plt.plot([vert_line,vert_line],[0,2])
#plt.xlim(vert_line-300,vert_line+300)
plt.subplot(3, 1, 3)
pc.make_figure(median_wavelength, median_pos,
               median_err, title="Median U")
plt.xlim(w - 200, w + 200)
plt.ylim(0.65, 1.405)
# plt.plot([vert_line,vert_line],[0,100])
# plt.xlim(vert_line-300,vert_line+300)
plt.savefig(figdir + 'Spectral_Lines/Halpha_median_Q_U.eps', overwrite=True)
plt.show()
