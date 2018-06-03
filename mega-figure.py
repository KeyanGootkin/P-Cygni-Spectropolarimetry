import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from glob import glob
import os
from astropy.io import fits, ascii
from astropy.table import Table, Column
import P_Cyg as pc
import scipy.stats

datadir = 'C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/'
figdir = 'C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Figures/'
rfitsfiles = glob(datadir + 'hpolccd*r_hw.fits')
rtxtfiles = glob(datadir + 'hpolccd*r_hw.fits.txt')
bfitsfiles = glob(datadir + 'hpolccd*b_hw.fits')
btxtfiles = glob(datadir + 'hpolccd*b_hw.fits.txt')


'''
median_wavelength, median_flux, median_pol, median_pos, median_err = pc.median_flux_pol_pos(
    rtxtfiles, rfitsfiles, 1000, radial_velocity=-8.9)
paerr = pc.position_angle_error(median_pol, median_err)



figshape = (23,21)

plt.figure(figsize=[25,25])
flux = plt.subplot2grid(figshape, (0, 0), colspan=19,rowspan=6)
pc.make_figure(median_wavelength, median_flux, 0, title="Median Flux")
plt.xlim(w - 100, w + 100)

pol = plt.subplot2grid(figshape, (7, 11), colspan=8,rowspan=6)
pc.make_figure(median_wavelength, median_pol,
               median_err, title="Median Polarimetry")
plt.xlim(w - 100, w + 100)
plt.ylim(0.75,1.5)

pos = plt.subplot2grid(figshape, (14, 11), colspan=8,rowspan=6)
pc.make_figure(median_wavelength, median_pos,
               paerr, title="Median Position Angle")
plt.xlim(w - 100, w + 100)
plt.ylim(60,85)

QU = plt.subplot2grid(figshape, (8, 0), rowspan=10,colspan=10)
plt.xlabel("Q",size=17, fontname='Times New Roman')
plt.ylabel("U",size=17, fontname='Times New Roman')

time = plt.subplot2grid(figshape, (21,0), colspan=19, rowspan=3)
time.plot(np.linspace(0,10000),)
plt.show()
#plt.savefig(figdir+"mega-fig_example.eps", overwrite=True, dpi = 1000)
'''
interesting_wavelengths = [4860, 5875, 6560]
gw = interesting_wavelengths[2]
count = 0
cmap = cm.get_cmap('magma')
for txtfile, fitsfile in zip(rtxtfiles, rfitsfiles):
    count += 1
    w, f, pol, pos, e = pc.txt_pol_data(txtfile, 1000, radial_velocity=-8.9)
    paerr = pc.position_angle_error(pol, e)
    cal_f = pc.divide_continuum(w,f,6500,6540,6610,6650)

    figshape = (23, 21)
    fig = plt.figure(figsize=[25, 25])

    plt.subplot2grid(figshape, (0, 0), colspan=19, rowspan=6)
    pc.make_figure(w, cal_f, 0, title="Flux")
    plt.xlim(gw - 100, gw + 100)
    plt.ylim(0,5.5)

    pol_fig = plt.subplot2grid(figshape, (7, 11), colspan=8, rowspan=6)
    pc.make_figure(w, pol, e, title="% Polarization")
    plt.xlim(gw - 100, gw + 100)
    plt.ylim(0.3, 2.5)

    plt.subplot2grid(figshape, (14, 11), colspan=8, rowspan=6)
    pc.make_figure(w, pos, paerr, title="Position Angle")
    plt.xlim(gw - 100, gw + 100)
    plt.ylim(30, 120)

    w, f, q, u, e = pc.txt_QU_data(txtfile, 1000, radial_velocity=-8.9)
    good_ind = []
    for i in w:
        if i >= 6540 and i <= 6610:
            good_ind.append(i)
    low_ind = list(w).index(min(good_ind))
    high_ind = list(w).index(max(good_ind))

    QU = plt.subplot2grid(figshape, (8, 0), rowspan=10, colspan=10)
    for i in good_ind:
        QU.plot(q[list(w).index(i)], u[list(w).index(i)], '-o',
                c=((cmap((w[list(w).index(i)] - w[low_ind]) / (w[high_ind] - w[low_ind])))))
    QU.plot(q[low_ind:high_ind], u[low_ind:high_ind])
    plt.xlabel("Q", size=17, fontname='Times New Roman')
    plt.ylabel("U", size=17, fontname='Times New Roman')
    plt.xlim(-0.5, 1)
    plt.ylim(0.5, 2)

    time = plt.subplot2grid(figshape, (21, 0), colspan=19, rowspan=3)
    t, ave_pol, ave_err = pc.get_ccd_ave_pol(fitsfile)
    time.axvline(x=t, linewidth=4)
    plt.xlim(49500, 53000)

    plt.savefig(figdir + "Mega_Fig/Halpha/Halpha_fig_" +
                str(count) + ".png", overwrite=True)
    plt.close(fig)
