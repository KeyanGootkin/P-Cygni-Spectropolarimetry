import P_Cyg as pc
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as cm
from glob import glob

datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
bfiles = glob(datadir + 'hpolccd*b_hw.fits')
rfiles = glob(datadir + 'hpolccd*r_hw.fits')
retfiles = glob(datadir + 'hpolret*.fits')

times, p_list, pos_list, err_list = pc.calculate_average_polarization_ret(
    retfiles)
p = np.mean(p_list)



def find_shape(polarization, i):
    pol = polarization / 100
    h_gamma = ((np.sin(np.deg2rad(i))**2) - pol*(np.sin(np.deg2rad(i))**2))/pol
    shape = (h_gamma - 2) / (2 + (3 * h_gamma))
    return(shape)


x = np.linspace(0, 90, 100)
shapes_list = []
for i in x:
    shape = find_shape(p, i)
    shapes_list.append(shape)

plt.plot(x, shapes_list)
plt.ylim(0,0.4)
plt.show()

print(find_shape(p,90))
