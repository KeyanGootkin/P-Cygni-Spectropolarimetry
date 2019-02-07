import P_Cyg as pc
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as cm
from glob import glob
from astropy.io import fits, ascii

datadir = 'C:/Users/Keyan/Desktop/Science/Data/P-Cygni/wuppe/'
files= glob(datadir+"*.fits")
all_w = []
all_q = []
all_u = []
all_e = []

for f in files:
    hdu = fits.open(f)
    tablehdu = hdu[1].header
    t = hdu[1].data
    w = t["WAVELENGTH"][0]
    all_w.append(w)
    q = t["Q"][0]
    all_q.append(q)
    u = t["U"][0]
    all_u.append(u)
    e = t["ERROR"][0]
    all_e.append(e)
"""
mq=np.median(all_q,axis=0)
mu=np.median(all_u,axis=0)

mp,mpa = pc.polarization(mq,mu)
"""
for q,u in zip(all_q,all_u):
    p,pa = pc.polarization(q,u)
    plt.plot(w,p)


def NordISP(wavelength):
    wavelength = np.array(wavelength)
    ISP = 1.06*np.exp(-1.15*(np.log(5500/wavelength)**2))
    return(ISP)
isp = NordISP(w)
plt.plot(w,isp)
plt.show()
