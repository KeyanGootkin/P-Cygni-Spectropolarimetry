import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
from astropy.io import fits,ascii
from astropy.table import Table,Column
import P_Cyg as pc
import scipy.stats

datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
figdir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Figures/'
txt_list = glob(datadir+'hpolret*.txt')
fits_list = glob(datadir+'hpolret*.fits')






plt.figure(figsize=[15,20])
plt.subplot(3,1,1)
median_wavelength,median_flux,median_pol,median_pos = pc.median_flux_pol_pos(txt_list,fits_list,250,radial_velocity=-8.9)
pc.make_figure(median_wavelength,median_flux,0,title="Median Flux")
#plt.plot([vert_line,vert_line],[min(median_flux),max(median_flux)])
#plt.xlim(vert_line-300,vert_line+300)
plt.subplot(3,1,2)
pc.make_figure(median_wavelength,median_pol,0,title="Median Polarization")
#plt.plot([vert_line,vert_line],[0,2])
#plt.xlim(vert_line-300,vert_line+300)
plt.subplot(3,1,3)
pc.make_figure(median_wavelength,median_pos,0,title="Median Position Angle")
#plt.plot([vert_line,vert_line],[0,100])
#plt.xlim(vert_line-300,vert_line+300)
plt.show()
#plt.savefig(figdir+'MedianFPP.eps',overwrite=True)
