import P_Cyg as pc
from glob import glob
import matplotlib.pyplot as plt

txtlist = glob(
    "C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/*r_hw.fits.txt")
fitslist = glob("C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Data/*r_hw.fits")
#all_w,all_f,all_q,all_u,all_e = pc.get_all_QU(txtlist, fitslist, 1000, radial_velocity=-8.9)
testfile = txtlist[4]
m_w, m_f, m_q, m_u, m_e = pc.median_flux_Q_U(
    txtlist, fitslist, 1000, radial_velocity=-8.9)
'''
plt.figure(figsize=[15,15])
pc.make_figure(m_q,m_u,0,title='Median Q vs. U',axis = ['Stokes Q','Stokes U'])
plt.scatter(m_q,m_u,s=3)
plt.savefig('C:/Users/Keyan/Desktop/Science/Data/P-Cygni/Figures/QU/Q_vs_U/Median_Q_vs_U.eps',overwrite=True)
plt.show()
'''
pc.make_figure(m_w, m_f, 0)
plt.show()

all_w, all_f, all_q, all_u, all_e = pc.get_all_QU(
    txtlist, fitslist, 1000, radial_velocity=-8.9)
for w, q, u in zip(all_w, all_q, all_u):
    H_q = pc.find_continuum_w_gap(w, q, 5650, 5680, -2, -1)
    H_u = pc.find_continuum_w_gap(w, u, 5650, 5680, -2, -1)
