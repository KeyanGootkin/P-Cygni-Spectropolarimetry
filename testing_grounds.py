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
'''
hdu = fits.open(datadir+'hpolret_p-cyg_19941204_hw.fits')

tablehdu = hdu[1]

table = tablehdu.data
wavelength = table["wavelength"][0]
Q = np.array(table["Q"][0])
U = np.array(table["U"][0])
fits_w_bin_edges,bin_Q = pc.easy_bin_mean(wavelength,Q,250)
fits_w_bin_edges,bin_U = pc.easy_bin_mean(wavelength,U,250)
pol,pos = pc.polarization(bin_Q,bin_U)
err = table["Error"]




tabletxt = ascii.read(file,names=["Wavelength","Flux","q",'u','err'],data_start=1)

flux = tabletxt["Flux"]
table_wavel = tabletxt["Wavelength"]


plt.plot(table_wavel,flux)
plt.show()

table_w_bin_edges,bin_flux = pc.easy_bin_mean(table_wavel,flux,250)

plt.plot(table_w_bin_edges[1:],bin_flux)
#plt.xlim(4700,5100)
plt.show()


plt.plot(fits_w_bin_edges[1:],pol)
plt.show()
plt.plot(fits_w_bin_edges[1:],pos)
plt.show()
print(min(pos),max(pos))
print(len(wavelength),len(table_wavel))
'''


'''plt.plot(fits_w_bin_edges[1:],bin_Q)
plt.show()

plt.plot(fits_w_bin_edges[1:],bin_U)
plt.show
'''


'''
wavelength,flux,pol,pos = pc.txt_pol_data(file,250)

astropy_table = Table.read(file, format = 'ascii',delimiter='\s',data_start=1,names=['Wavelength','Flux','q','u','err'])

print(astropy_table)






plt.plot(wavelength,flux)
plt.show()
plt.plot(wavelength,pol)
plt.show()
plt.plot(wavelength,pos)
plt.show()
'''

tables_list = glob(datadir+'hpolret*.txt')

pc.stack_txt_pol_data(tables_list,250,radial_velocity=0,window=[6500,7000])

'''
plt.subplot(3,1,2)
for i in range(len(tables_list)):
    files = tables_list[i]
    wavelengthx,fluxx,polx,posx = pc.txt_pol_data(files,250)
    plt.plot(wavelengthx,polx,label=(files[len(files)-20:len(files)-12]))
plt.xlim(6450,6650)
plt.title("% Polarization")


plt.subplot(3,1,3)
for i in range(len(tables_list)):
    files = tables_list[i]
    wavelengthx,fluxx,polx,posx = pc.txt_pol_data(files,250)
    plt.plot(wavelengthx,posx,label=(files[len(files)-20:len(files)-12]))
plt.xlim(6450,6650)
plt.title("Position Angle")
'''
