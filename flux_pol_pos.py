import P_Cyg as pc
import numpy as np
import matplotlib as plt

datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
bfiles = glob(datadir+'hpolccd*b_hw.fits')
rfiles = glob(datadir+'hpolccd*r_hw.fits')


#pc.stack_txt_pol_data(rfiles,250)
#plt.show()
pc.stack_txt_pol_data(bfiles,250)
plt.show()
