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
             deltadec, center[1] + deltadec], color=color, linewidth=3)
    return()



ism = ascii.read("C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/ISM_stars/ISM_stars.dat", delimiter = ':')
print(ism)
pdist = 1700
distcut = 1000000

ra = ism["RA"]
dec = ism["DEC"]
dist = ism["Distance"]
pra = coord.Angle((20, 15, 56.53) ,unit='hourangle').deg
pdec = coord.Angle((+37, 52, 35.3) , unit="deg").deg

cygism = ism[ism["Distance"] >= (pdist-distcut)]
cygism = cygism[cygism["Distance"] <= (pdist+distcut)]
foreism = ism[ism["Distance"] < (pdist-distcut)]
backism = ism[ism["Distance"] > (pdist+distcut)]

fig = plt.figure(figsize=(15,15),dpi=300)
#ax1 = fig.add_subplot(111,projection='mollweide')
ax = plt.gca()
ax.scatter(cygism["RA"],cygism["DEC"],s=250,color='w',marker='*')
'''ax.scatter(foreism["RA"],foreism["DEC"],s=6,color='g')

ax.scatter(backism["RA"],backism["DEC"],s=6,color='purple')'''

ax.scatter(pra,pdec,color="#FFEF89",marker='*',label='P Cygni',s=2000)
ax.set_facecolor("#283655")



xlim = ax.get_xlim()
ylim = ax.get_ylim()
for r,d,pol,pa in zip(cygism["RA"],cygism["DEC"],cygism['% Pol'],cygism['PA']):
    pol_line(r,d,pol,pa,decstretch= (ylim[1]-ylim[0])/(xlim[1]-xlim[0]),color='w')
ax.scatter(306.8,35,s=250,color='r',marker='*')
pol_line(306.8,35,3,90,decstretch= (ylim[1]-ylim[0])/(xlim[1]-xlim[0]),color='r')
plt.text(306.15,35.15,"Example Star",color='w',fontsize=30,fontname='Arial')
plt.text(306,34.75,"3% Polarization",color='w',fontsize=30,fontname="Arial")
'''for r,d,pol,pa in zip(foreism["RA"],foreism["DEC"],foreism['% Pol'],foreism['PA']):
    pol_line(r,d,pol,pa,decstretch= (ylim[1]-ylim[0])/(xlim[1]-xlim[0]),color='g')

for r,d,pol,pa in zip(backism["RA"],backism["DEC"],backism['% Pol'],backism['PA']):
    pol_line(r,d,pol,pa,decstretch= (ylim[1]-ylim[0])/(xlim[1]-xlim[0]),color='purple')'''


plt.xticks(fontsize=30, fontname='Arial', color='darkslategrey')
plt.yticks(fontsize=30, fontname='Arial', color='darkslategrey')

ax.set_xlabel(r"Right Ascension ($^\circ$)",size=30,fontname="Arial")
ax.set_ylabel(r"Declination ($^\circ$)",size=30,fontname="Arial")
#ax.annotate(xy=[303.6,38],s="P Cygni")
#ax.set_title("Polarization Map Around P Cygni",size=26,fontname="Garamond")
plt.savefig("C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Figures/Poster/map.png",overwrite=True,transparent=False)
