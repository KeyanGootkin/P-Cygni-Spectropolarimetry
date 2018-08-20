import P_Cyg as pc
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from glob import glob

datadir = 'C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/'
figdir = 'C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Figures/'

bfitsfiles = glob(datadir + 'hpolccd*b_hw.fits')
btxtfiles = glob(datadir + 'hpolccd*b_hw.fits.txt')
rfitsfiles = glob(datadir + 'hpolccd*r_hw.fits')
rtxtfiles = glob(datadir + 'hpolccd*r_hw.fits.txt')
rettxtfiles = glob(datadir + 'hpolret*.txt')
retfitsfiles = glob(datadir + 'hpolret*.fits')
allfitsfiles = bfitsfiles + rfitsfiles + retfitsfiles
alltxtfiles = btxtfiles + rtxtfiles + rettxtfiles


def HalphaContinuum(wavelength,data,error=None):
    continuum = []
    c_wavelength = []
    c_error = []
    for i in wavelength:
        if i >= 6000 and i <= 6400:
            continuum.append(data[list(wavelength).index(i)])
            if error != None:
                c_error.append(error[list(wavelength).index(i)])
            c_wavelength.append(i)
        elif i >= 6700 and i <= 7000:
            continuum.append(data[list(wavelength).index(i)])
            if error != None:
                c_error.append(error[list(wavelength).index(i)])
            c_wavelength.append(i)
    continuum = np.interp(wavelength,c_wavelength,continuum)
    if error != None:
        interp_error = np.interp(wavelength,c_wavelength,c_error)
    ha_continuum = []
    ha_error = []
    ha_wavelength = []
    for i in wavelength:
        if i >= 6540 and i <= 6600:
            ha_continuum.append(continuum[list(wavelength).index(i)])
            if error != None:
                ha_error.append(interp_error[list(wavelength).index(i)])
            ha_wavelength.append(i)
    ave_continuum = np.mean(ha_continuum)
    if error != None:
        ave_error = np.mean(ha_error)
    if error == None:
        return ave_continuum
    else:
        return ave_continuum,ave_error

def HalphaLine(wavelength,data,error=None):
    line = []
    l_wavelength = []
    l_error = []
    for i in wavelength:
        if i >= 6540 and i <= 6600:
            line.append(data[list(wavelength).index(i)])
            if error != None:
                l_error.append(error[list(wavelength).index(i)])
    ave_line = np.mean(line)
    if error != None:
        ave_error = np.mean(l_error)
    if error == None:
        return(ave_line)
    else:
        return(ave_line,ave_error)

def meanCONTtoLINE(q,u,error):
    weighted_q_sum_list = []
    weighted_u_sum_list = []
    weight_list = []
    for sq,su,e in zip(q,u,error):
        q_weight = 1/(e**2)
        u_weight = 1/(e**2)
        weight_list.append(q_weight)
        weighted_q_sum_list.append(sq*q_weight)
        weighted_u_sum_list.append(su*u_weight)
    mean_q = sum(weighted_q_sum_list)/(sum(weight_list))
    mean_u = sum(weighted_u_sum_list)/(sum(weight_list))
    return(mean_q,mean_u)


def QUcontinuumTOline(wavelengths,qs,us,errors=None,mean=False):
    qconlist = []
    qlinelist= []
    uconlist = []
    ulinelist=[]

    if not mean:
        for w,q,u in zip(wavelengths,qs,us):
            q_cont = HalphaContinuum(w,q)
            qconlist.append(q_cont)
            q_line = HalphaLine(w,q)
            qlinelist.append(q_line)
            u_cont = HalphaContinuum(w,u)
            uconlist.append(u_cont)
            u_line = HalphaLine(w,u)
            ulinelist.append(u_line)
            plt.plot([q_cont,q_line],[u_cont,u_line],color = 'black')
        plt.scatter(qconlist,uconlist,color='b',label="Continuum")
        plt.scatter(qlinelist,ulinelist,color='orange',label = "Line")
        plt.legend()

    elif mean:
        if errors == None:
            print("No error provided")
        cont_error = []
        line_error = []
        for w,q,u,e in zip(wavelengths,qs,us,errors):
            q_cont,q_cont_err = HalphaContinuum(w,q,error=e)
            cont_error.append(q_cont_err)
            qconlist.append(q_cont)
            q_line,q_line_err = HalphaLine(w,q,error=e)
            line_error.append(q_line_err)
            qlinelist.append(q_line)
            u_cont,u_cont_err = HalphaContinuum(w,u,error=e)
            uconlist.append(u_cont)
            u_line,u_line_err = HalphaLine(w,u,error=e)
            ulinelist.append(u_line)
        q_cont,u_cont = meanCONTtoLINE(qconlist,uconlist,cont_error)
        q_line,u_line = meanCONTtoLINE(qlinelist,ulinelist,line_error)
        plt.plot([q_cont,q_line],[u_cont,u_line],'black')
        plt.scatter(q_cont,u_cont,color='b',label="Continuum")
        plt.scatter(q_line,u_line,color='orange',label = "Line")
        plt.legend()

w,f,q,u,e = pc.get_all_QU(rtxtfiles, rfitsfiles,1000, radial_velocity=-8.9)

QUcontinuumTOline(w,q,u,errors=e,mean=True)
plt.xlim(0.3,0.4)
plt.ylim(0.9,1.1)
plt.show()


"""
rwave,rflux,rpol,rpos,rerr = pc.median_flux_pol_pos(rtxtfiles,rfitsfiles,1000,radial_velocity=-8.9)
plt.plot(rwave,rflux)
plt.xlim(5900,7100)
boundaries = [6000,6400,6540,6600,6700,7000]
for b in boundaries:
    plt.axvline(b)
plt.grid()
plt.show()
"""
"""
mtxtfiles = pc.match_rb_txt_files(rtxtfiles, btxtfiles)
mfitsfiles = pc.match_rb_fits_files(rfitsfiles, bfitsfiles)
pc.stack_txt_pol_data(mtxtfiles[0], mfitsfiles[0], 250)
plt.show()
"""
"""


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
for file in btxtfiles:
    count += 1
    w, f, q, u, e = pc.txt_QU_data(file, 1000, radial_velocity=-8.9)
    good_ind = []
    for i in w:
        if i >= 5860 and i <= 5895:
            good_ind.append(i)
    low_ind = list(w).index(min(good_ind))
    high_ind = list(w).index(max(good_ind))
    fig = plt.figure()
    for i in good_ind:
        plt.plot(q[list(w).index(i)], u[list(w).index(i)], '-o',
                 c=((cmap((w[list(w).index(i)] - w[low_ind]) / (w[high_ind] - w[low_ind])))))
    plt.plot(q[low_ind:high_ind], u[low_ind:high_ind])
    plt.title(str(count))
    plt.savefig(figdir + "QU/He_Q_vs_U/He_Q_vs_U " +
                str(count) + '.png', overwrite=True)
    plt.close(fig)
""""""
w,f,q,u,e = pc.mean_flux_Q_U(rtxtfiles, rfitsfiles, 50, radial_velocity=-8.9)
plt.plot(q,u,'-o')
plt.show()
""""""
w,f,q,u,e = pc.median_flux_Q_U(rtxtfiles, rfitsfiles, 50, radial_velocity=-8.9)
plt.plot(q,u,'-o')
plt.show()
plt.savefig(figdir+"QU/Q_vs_U/Median_Q_vs_U.eps",overwrite=True,dpi=2000)
"""
"""
w,f,p,a,e = pc.fpp_sub_mean(rtxtfiles, rfitsfiles, 1000, radial_velocity=-8.9)
ae = pc.position_angle_error(p,e)
halpha = 6560
width = 100

fig = plt.figure(figsize=[15, 20])
plt.subplot(3, 1, 1)
pc.make_figure(w,f,0,title='Flux ')
plt.xlim(halpha - width, halpha + width)
polfig = plt.subplot(3, 1, 2)
pc.make_figure(w,p,e,title='% Polarization')
#plt.ylim(-0.25,1)
plt.xlim(halpha - width, halpha + width)
pafig = plt.subplot(3, 1, 3)
pc.make_figure(w,a,ae,title='Position Angle')
#plt.ylim(0.3, 2.1)
plt.xlim(halpha - width, halpha + width)
polfig.axvline(6500)
plt.show()
"""
"""count = 0
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
    plt.savefig(figdir+"individual fpp/Halpha_fpp_"+str(count)+".jpg",overwrite=True)"""
"""
count = 0
cmap = cm.get_cmap('magma')
w, f, q, u, e = pc.median_flux_Q_U(
    rtxtfiles, rfitsfiles, 1000, radial_velocity=0)

good_ind = []
for i in w:
    if i >= 6540 and i <= 6610:
        good_ind.append(i)
low_ind = list(w).index(min(good_ind))
high_ind = list(w).index(max(good_ind))
frame = plt.figure(figsize=[20, 10])
plt.subplot2grid((1, 2), (0, 0))
for i in good_ind:
    plt.plot(q[list(w).index(i)], u[list(w).index(i)], '-o',
             c=((cmap((w[list(w).index(i)] - w[low_ind]) / (w[high_ind] - w[low_ind])))))
plt.plot(q[low_ind:high_ind], u[low_ind:high_ind], 'blue')
plt.xlim(0, 0.5)
plt.ylim(0, 1.3)
plt.title("H alpha")
plt.subplot2grid((1, 2), (0, 1))
plt.plot(q, u, "r")
plt.xlim(0, 0.5)
plt.ylim(0, 1.3)
plt.title("ALL")
plt.savefig(figdir+"is halpha polarized.eps",overwrite=True)"""
