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
txtretlist = glob(datadir + 'hpolret*.txt')
fitsretlist = glob(datadir + 'hpolret*.fits')
ccdrfitsfiles = glob(datadir + 'hpolccd*r_hw.fits')
ccdbfitsfiles = glob(datadir + 'hpolccd*b_hw.fits')
ccdrtxtfiles = glob(datadir + 'hpolccd*r_hw.fits.txt')
ccdbtxtfiles = glob(datadir + 'hpolccd*b_hw.fits.txt')

ccd_fits_match = pc.match_rb_fits_files(ccdrfitsfiles,ccdbfitsfiles)
ccd_txt_match = pc.match_rb_txt_files(ccdrtxtfiles,ccdbtxtfiles)

ccd_fits = ccd_fits_match[0]
ccd_txt = ccd_txt_match[0]

ccd_r = ccd_txt[0]
ccd_b = ccd_txt[1]

waver,fluxr,polr,posr,errr = pc.txt_pol_data(ccd_r, 250)
waveb,fluxb,polb,posb,errb = pc.txt_pol_data(ccd_b, 250)

plt.plot(waver,polr)
plt.plot(waveb,polb)
plt.show()



wavec = np.array(list(waver)+list(waveb))
fluxc = np.array(list(fluxr)+list(fluxb))
polc = np.array(list(polr)+list(polb))
wavesort = np.argsort(wavec)
wavec = wavec[wavesort]
fluxc = fluxc[wavesort]
polc = polc[wavesort]
print(wavec)

plt.plot(wavec,polc)
plt.show()
print(len(wavec))


print(len(wavec),len(fluxc),len(polc))
wavec,fluxc,polc = list(wavec),list(fluxc),list(polc)
for i in range(500):
    print(i)
    if polc[i] == 0 or fluxc[i] == 0:
        wavec.remove(wavec[i])
        fluxc.remove(fluxc[i])
        polc.remove(polc[i])

print(len(wavec),len(polc))
print(polc)

plt.plot(wavec,polc)
plt.show()
