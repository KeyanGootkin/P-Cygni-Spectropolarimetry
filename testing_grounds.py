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

mtxtfiles = pc.match_rb_txt_files(rtxtfiles, btxtfiles)
mfitsfiles = pc.match_rb_fits_files(rfitsfiles, bfitsfiles)
pc.stack_txt_pol_data(mtxtfiles[0], mfitsfiles[0], 250)
plt.show()
