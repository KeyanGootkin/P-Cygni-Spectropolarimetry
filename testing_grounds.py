import P_Cyg as pc
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from glob import glob

datadir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Data/'
figdir = 'C:/Users/Keyan/Desktop/Data/P-Cygni/Figures/'

bfitsfiles = glob(datadir + 'hpolccd*b_hw.fits')
btxtfiles = glob(datadir + 'hpolccd*b_hw.fits.txt')
rfitsfiles = glob(datadir + 'hpolccd*r_hw.fits')
rtxtfiles = glob(datadir + 'hpolccd*r_hw.fits.txt')
rettxtfiles = glob(datadir + 'hpolret*.txt')
retfitsfiles = glob(datadir + 'hpolret*.fits')
allfitsfiles = bfitsfiles + rfitsfiles + retfitsfiles
alltxtfiles = btxtfiles + rtxtfiles + rettxtfiles
"""
mtxtfiles = pc.match_rb_txt_files(rtxtfiles, btxtfiles)
mfitsfiles = pc.match_rb_fits_files(rfitsfiles, bfitsfiles)
pc.stack_txt_pol_data(mtxtfiles[0], mfitsfiles[0], 250)
plt.show()
""""""


all_wavelength, all_flux, all_Q, all_U, all_err = pc.get_all_QU(
    rtxtfiles, rfitsfiles, 1000, radial_velocity=-8.9)
count = 0
for w, f, q, u, e in zip(all_wavelength, all_flux, all_Q, all_U, all_err):
    count += 1
    plt.figure(figsize=[15, 20])
    plt.subplot(3, 1, 1)
    pc.make_figure(w,f,0,title='Flux '+str(count))
    plt.xlim(halpha - width, halpha + width)
    plt.subplot(3, 1, 2)
    pc.make_figure(w,q,e,title='Stokes Q '+str(count))
    plt.ylim(-0.25,1)
    plt.xlim(halpha - width, halpha + width)
    plt.subplot(3, 1, 3)
    pc.make_figure(w,u,e,title='Stokes U '+str(count))
    plt.ylim(0.3, 2.1)
    plt.xlim(halpha - width, halpha + width)

    plt.savefig(figdir+'QU/Halpha_QU_'+str(count)+'.eps',overwrite=True)
""""""
count = 0
cmap = cm.get_cmap('magma')
for file in rtxtfiles:
    count += 1
    w,f,q,u,e = pc.txt_QU_data(file,1000,radial_velocity=-8.9)
    good_ind = []
    for i in w:
        if i >= 6540 and i <= 6610:
            good_ind.append(i)
    low_ind = list(w).index(min(good_ind))
    high_ind = list(w).index(max(good_ind))
    plt.figure()
    for i in good_ind:
        plt.plot(q[list(w).index(i)],u[list(w).index(i)],'-o',c=((cmap((w[list(w).index(i)]-w[low_ind])/(w[high_ind]-w[low_ind])))))
    plt.plot(q[low_ind:high_ind],u[low_ind:high_ind])
    plt.title(str(count))
    plt.savefig(figdir+"QU/Halpha_Q_vs_U/Q_vs_U "+str(count)+'.jpg',overwrite=True,dpi=1000)
"""
"""w,f,q,u,e = pc.mean_flux_Q_U(rtxtfiles, rfitsfiles, 50, radial_velocity=-8.9)
plt.plot(q,u,'-o')
plt.savefig(figdir+"QU/Q_vs_U/Mean_Q_vs_U.eps",overwrite=True,dpi=2000)"""
"""
w,f,q,u,e = pc.median_flux_Q_U(rtxtfiles, rfitsfiles, 50, radial_velocity=-8.9)
plt.plot(q,u,'-o')
plt.show()
plt.savefig(figdir+"QU/Q_vs_U/Median_Q_vs_U.eps",overwrite=True,dpi=2000)"""
"""
w,f,p,a,e = pc.fpp_sub_mean(rtxtfiles, rfitsfiles, 1000, radial_velocity=-8.9)
ae = pc.position_angle_error(p,e)
halpha = 6560
width = 100

plt.figure(figsize=[15, 20])
plt.subplot(3, 1, 1)
pc.make_figure(w,f,0,title='Flux ')
plt.xlim(halpha - width, halpha + width)
plt.subplot(3, 1, 2)
pc.make_figure(w,p,e,title='% Polarization')
#plt.ylim(-0.25,1)
plt.xlim(halpha - width, halpha + width)
plt.subplot(3, 1, 3)
pc.make_figure(w,a,ae,title='Position Angle')
#plt.ylim(0.3, 2.1)
plt.xlim(halpha - width, halpha + width)
plt.show()"""
count = 0
for file in rtxtfiles:
    count += 1
    w,f,pol,pos,e = pc.txt_pol_data(file, 1000, radial_velocity = -8.9)
    paerr = pc.position_angle_error(pol,e)
    plt.figure(figsize=[15,20])
    plt.subplot(3, 1, 1)
    pc.make_figure(w, f, 0, title=str(count))
    plt.xlim(6560-200,6560+200)
    plt.subplot(3, 1, 2)
    pc.make_figure(w, pol, e, title = "Polarization")
    plt.xlim(6560-200,6560+200)
    plt.ylim(0.5,2.5)
    plt.subplot(3, 1, 3)
    pc.make_figure(w, pos,paerr, title = "Position Angle")
    plt.xlim(6560-200,6560+200)
    plt.ylim(40,100)
    plt.savefig(figdir+"individual fpp/Halpha_fpp_"+str(count)+".jpg",overwrite=True)
