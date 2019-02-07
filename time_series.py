import numpy as np
import matplotlib.pyplot as plt
import P_Cyg as pc
from glob import glob

datadir = 'C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/'
fits = glob(datadir + 'hpolccd*r_hw.fits')
txt = glob(datadir + 'hpolccd*r_hw.fits.txt')

all_wavelength, all_flux, all_pol, all_pos, all_err = pc.get_all_fpp(txt, fits, 1000, radial_velocity=-8.9)

interp_flux = [all_flux[0]]
interp_pol = [all_pol[0]]
interp_pos = [all_pos[0]]
interp_err = [all_err[0]]
for i in range(1, len(all_wavelength)):
    interpolated_flux = np.interp(
        all_wavelength[0], all_wavelength[i], all_flux[i])
    interp_flux.append(interpolated_flux)
    interpolated_pol = np.interp(
        all_wavelength[0], all_wavelength[i], all_pol[i])
    interp_pol.append(interpolated_pol)
    interpolated_pos = np.interp(
        all_wavelength[0], all_wavelength[i], all_pos[i])
    interp_pos.append(interpolated_pos)
    interpolated_err = np.interp(
        all_wavelength[0], all_wavelength[i], all_err[i])
    interp_err.append(interpolated_err)
tspol = []
tspos = []
for i in range(len(all_wavelength[0])):
    pol = []
    pos = []
    for x in range(len(interp_pol)):
        pol.append(interp_pol[x][i])
        pos.append(interp_pos[x][i])
    tspol.append(pol)
    tspos.append(pos)

gind = []
for i in range(len(all_wavelength[0])):
    if all_wavelength[0][i] >= 6400 and all_wavelength[0][i] <= 6700:
        gind.append(i)
x = np.linspace(1,len(tspol[0]),len(tspol[0]))
for i in [111,112,113,114,115,116,117]:
    plt.plot(x,tspol[i])
plt.plot(x,np.median(tspol[gind[0]:gind[len(gind)-1]],axis=0),linewidth=5,color='black')
plt.show()
