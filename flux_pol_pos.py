import P_Cyg as pc
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as cm
from glob import glob

datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
bfitsfiles = glob(datadir + 'hpolccd*b_hw.fits')
btxtfiles = glob(datadir + 'hpolccd*b_hw.fits.txt')
rfitsfiles = glob(datadir + 'hpolccd*r_hw.fits')
rtxtfiles = glob(datadir + 'hpolccd*r_hw.fits.txt')
rettxtfiles = glob(datadir + 'hpolret*.txt')
retfitsfiles = glob(datadir + 'hpolret*.fits')
allfitsfiles = bfitsfiles + rfitsfiles + retfitsfiles
alltxtfiles = btxtfiles + rtxtfiles + rettxtfiles

vert_line = 4850

pc.vert_stack_ccd_txt_pol_data(rtxtfiles, rfitsfiles, 250, radial_velocity=-8.9,window=[5000,11000])
plt.show()

pc.stack_txt_pol_data(rtxtfiles, rfitsfiles, 250, radial_velocity=-8.9,window=[5000,11000])
plt.show()

pc.stack_txt_pol_data(btxtfiles, bfitsfiles, 250, radial_velocity=-8.9)
plt.show()

pc.stack_ccd_txt_pol_data(alltxtfiles,allfitsfiles,250,radial_velocity=-8.9,window=[3000,11000])
plt.show()

bwave,bflux,bpol,bpos,berr = pc.median_flux_pol_pos(btxtfiles,bfitsfiles,250,radial_velocity=-8.9)
rwave,rflux,rpol,rpos,rerr = pc.median_flux_pol_pos(rtxtfiles,rfitsfiles,250,radial_velocity=-8.9)
plt.figure(figsize=[15,20])
plt.subplot(3,1,1)
plt.plot(bwave,bflux)
plt.plot(rwave,rflux)
plt.plot([6575,6575],[0,10e-11])
plt.xlim(6000,7000)
plt.subplot(3,1,2)
plt.plot(bwave,bpol)
plt.plot(rwave,rpol)
plt.plot([6575,6575],[0,2])
plt.xlim(6000,7000)
plt.subplot(3,1,3)
plt.plot(bwave,bpos)
plt.plot(rwave,rpos)
plt.xlim(6000,7000)
plt.plot([6575,6575],[0,90])
plt.show()
