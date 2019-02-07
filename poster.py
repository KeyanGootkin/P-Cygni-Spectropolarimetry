import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import P_Cyg as pc
datadir = 'C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/'
files = glob(datadir + 'hpolccd_p-cyg_*r_hw.fits.txt')
fits = glob(datadir + 'hpolccd_p-cyg_*r_hw.fits')
font="Arial"
red="#fe8a71ff"
yellow="#ffdf0045"
blue = "#3da4abff"
"""Q U"""

ws,qs,us,es = [],[],[],[]
for f in files:
    w,f,q,u,e = pc.txt_QU_data(f, 1000, radial_velocity=-8.9)
    ws.append(w)
    qs.append(q)
    us.append(u)
    es.append(e)
qq,uu,ee = [],[],[]
for q,u,e in zip(qs,us,es):
    qq.append(np.mean(q))
    uu.append(np.mean(u))
    ee.append(np.mean(e))
f = plt.figure(figsize=(15,15),dpi=300)
plt.scatter(qq,uu,s=250)
plt.scatter(np.mean(qq),np.mean(uu),color=red,marker='*',label="Mean Q vs Mean U",s=1000)
plt.xlabel(r"Stokes Q",size=30,fontname=font)
plt.ylabel("Stokes U",size=30,fontname=font)
plt.xticks(fontsize=30, fontname=font, color='darkslategrey')
plt.yticks(fontsize=30, fontname=font, color='darkslategrey')
plt.savefig("C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Figures/Poster/qquu.png",overwrite=True,transparent=True)

"""Map"""
"""See ISM.py"""

"""Stacked H alpha"""

fig1, (polfig) = plt.subplots(1,dpi=600,figsize=(10,20))
plt.xticks(fontsize=20, fontname=font, color='darkslategrey')
plt.yticks(fontsize=20, fontname=font, color='darkslategrey')

fig2, (pafig) = plt.subplots(1,dpi=600,figsize=(10,20))

all_wavelength = []
all_flux = []
all_pol = []
all_pos = []
for i in range(len(files)):
    txt = files[i]
    wavelength, flux, pol, pos, err = pc.txt_pol_data(txt, 1000)
    wavelength = pc.dedopler(wavelength, -8.9)
    all_wavelength.append(wavelength)
    all_flux.append(flux)
    all_pol.append(pol)
    all_pos.append(pos)
interp_pol = []
interp_pos = []
for i in range(len(all_wavelength)):
    interp_pol.append(np.interp(
        all_wavelength[0], all_wavelength[i], all_pol[i]))
    interp_pos.append(np.interp(
        all_wavelength[0], all_wavelength[i], all_pos[i]))
median_pol = np.median(interp_pol,axis=0)
median_pos = np.median(interp_pos,axis=0)
count = 1
for p,pa in zip(interp_pol,interp_pos):

    polfig.plot(all_wavelength[0],(np.array(p)+count*0.5),color='black',linewidth=3)
    pafig.plot(all_wavelength[0],(np.array(pa)+count*37), color = "black",linewidth=3)
    count+=1
polfig.plot(all_wavelength[0],median_pol,color="r",linewidth=3)
pafig.plot(all_wavelength[0],median_pos,color="r",linewidth=3)
# vertical line plt.plot([6563,6563],[-100,100],linestyle=':',color='black',alpha=0.5)
polfig.set_xlim(6563-100,6563+100)
polfig.set_xlabel(r"Wavelength ($\AA$)",size=30,fontname=font)
pafig.set_xlim(6563-100,6563+100)
pafig.set_xlabel(r"Wavelength ($\AA$)",size=30,fontname=font)
polfig.set_ylim(0,19)
polfig.set_ylabel("% Polarization + Constant",size=30,fontname=font)
pafig.set_ylim(0,1400)
pafig.set_ylabel(r"Position Angle + Constant ($^\circ}$)",size=30,fontname=font)
plt.xticks(fontsize=20, fontname=font, color='darkslategrey')
plt.yticks(fontsize=20, fontname=font, color='darkslategrey')

fig1.savefig("C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Figures/Poster/stackedpol.png",overwrite=True,transparent=True)
fig2.savefig("C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Figures/Poster/stackedpos.png",overwrite=True,transparent=True)

rwave,rflux,rpol,rpos,rerr = pc.median_flux_pol_pos(files,fits,1000,radial_velocity=-8.9)
rpaerr = pc.position_angle_error(rpol,rerr)
fw,fp,fpa,fe,fpae = [],[],[],[],[]
for w,p,pa,e,pae in zip(rwave,rpol,rpos,rerr,rpaerr):
    if w >= 6563-70 and w <= 6563+70:
        fw.append(w)
        fp.append(p)
        fpa.append(pa)
        fe.append(e)
        fpae.append(pae)
f1,polf=plt.subplots(1,figsize=(10,6),dpi=600)
f1.tight_layout(pad=10)
polf.plot(fw,fp,color='r',linewidth=3)
polf.errorbar(fw,fp,yerr=fe,color="r",elinewidth=2)
polf.set_xlabel(r"Wavelength ($\AA$)",size=30,fontname=font)
polf.set_ylabel("% Polarization",size=30,fontname=font)
plt.xticks(fontsize=20, fontname=font, color='darkslategrey')
plt.yticks(fontsize=20, fontname=font, color='darkslategrey')
plt.savefig("C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Figures/Poster/alphapol.png",overwrite=True,transparent=True)
f1.clf()
f2,paf = plt.subplots(1,figsize=(10,6),dpi=600)
f2.tight_layout(pad=10)
paf.plot(fw,fpa,color='r',linewidth=3)

paf.errorbar(fw,fpa,yerr=fpae,color="r",elinewidth=2)
paf.set_xlabel(r"Wavelength ($\AA$)",size=30,fontname=font)
paf.set_ylabel(r"Position Angle ($^\circ}$)",size=30,fontname=font)
plt.xticks(fontsize=20, fontname=font, color='darkslategrey')
plt.yticks(fontsize=20, fontname=font, color='darkslategrey')
plt.savefig("C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Figures/Poster/alphapa.png",overwrite=True,transparent=True)

