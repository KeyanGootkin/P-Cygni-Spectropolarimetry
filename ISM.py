import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from glob import glob
import os
import astropy.coordinates as coord
from astropy.io import fits, ascii
from astropy.table import Table, Column
#import P_Cyg as pc
import scipy.stats

#Only run this part if you need to remake the table
#IMPORTANT remember to take out weird stars if you do remake table :-)
"""file = "C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/ISM_stars/cyg.fits"

hdu = fits.open(file)
datahdu = hdu[1]
data = hdu[1].data

ra_raw = data["RA2000"]
dec_raw = data["DE2000"]
ra = []
for x in ra_raw:
    h = int(x[0:2])
    m = int(x[3:5])
    s = float(x[6:12])
    hms = (h,m,s)
    ra.append(hms)
dec = []
for x in dec_raw:
    d = int(x[:3])
    m = int(x[4:6])
    s = float(x[7:])
    dms = (d,m,s)
    dec.append(dms)
pol = data["Pol"]
err = data['e_Pol']
pa = data["PA"]
paerr = data['e_PA']
ra = coord.Angle(ra,unit='hourangle').deg
dec = coord.Angle(dec, unit = 'deg').deg
spec = list(data["SpType"])
hdnum = data['HD']
dist = data["Dist"]

table = Table([spec,pol,err,pa,paerr,ra,dec,hdnum,dist], names = ['Spectral Class', '% Pol', 'Pol Error', 'PA', 'PA Error', "RA", "DEC","HD","Distance"])

ascii.write(table, "C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/ISM_stars/ISM_stars.dat", delimiter = ':', overwrite = True)
"""
def pol_line(ra, dec, pol, pa, rastretch=1, decstretch=1, color='b'):
    """
    Creates a line centered on the ra and dec of a star which shows the % polarization
    (length of line), and position angle (angle of line) of that star. NOTE: you
    must have the figure open before you call this

    Parameters
    ----------
    ra : Float
        Right Ascension of the star

    dec : Float
        Declination of the star

    pol : Float
        % Polarization of the star

    pa : Float
        Position Angle of the star

    rastretch : Float
        The factor by which the line must be spread in ra in order to preserve the
        correct angle visually on the graph

    decstretch : Float
        The factor by which the line must be stretched in dec in order to preserve
        the correct angle visually on the graph

    color : string
        The matplotlib color you want the line to be
    """
    center = [ra, dec]
    l = pol / 7
    deltadec = l * np.cos(np.deg2rad(pa)) * decstretch
    deltara = l * np.sin(np.deg2rad(pa)) * rastretch
    plt.plot([center[0] - deltara, center[0] + deltara], [center[1] -
             deltadec, center[1] + deltadec], color=color, linewidth=0.5)
    return()



ism = ascii.read("C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/ISM_stars/ISM_stars.dat", delimiter = ':')

pdist = 1700
distcut = 300

ra = ism["RA"]
dec = ism["DEC"]
dist = ism["Distance"]
pra = coord.Angle((20, 15, 56.53) ,unit='hourangle').deg
pdec = coord.Angle((+37, 52, 35.3) , unit="deg").deg

cygism = ism[ism["Distance"] >= (pdist-distcut)]
cygism = cygism[cygism["Distance"] <= (pdist+distcut)]
foreism = ism[ism["Distance"] < (pdist-distcut)]
backism = ism[ism["Distance"] > (pdist+distcut)]

fig = plt.figure()
#ax1 = fig.add_subplot(111,projection='mollweide')
ax = plt.gca()
ax.scatter(cygism["RA"],cygism["DEC"],s=6)

'''ax.scatter(foreism["RA"],foreism["DEC"],s=6,color='g')

ax.scatter(backism["RA"],backism["DEC"],s=6,color='purple')'''

ax.scatter(pra,pdec,color='r',marker='+',label='P Cygni')




xlim = ax.get_xlim()
ylim = ax.get_ylim()
for r,d,pol,pa in zip(cygism["RA"],cygism["DEC"],cygism['% Pol'],cygism['PA']):
    pol_line(r,d,pol,pa,decstretch= (ylim[1]-ylim[0])/(xlim[1]-xlim[0]))

'''for r,d,pol,pa in zip(foreism["RA"],foreism["DEC"],foreism['% Pol'],foreism['PA']):
    pol_line(r,d,pol,pa,decstretch= (ylim[1]-ylim[0])/(xlim[1]-xlim[0]),color='g')

for r,d,pol,pa in zip(backism["RA"],backism["DEC"],backism['% Pol'],backism['PA']):
    pol_line(r,d,pol,pa,decstretch= (ylim[1]-ylim[0])/(xlim[1]-xlim[0]),color='purple')'''

ax.grid()
ax.set_xlabel("RA (Deg)")
ax.set_ylabel("DEC (Deg)")
ax.annotate(xy=[20.22,37.95],s="P Cygni")
ax.set_title("Polarization Map Around P Cygni")
plt.show()
