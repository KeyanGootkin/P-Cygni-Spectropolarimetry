import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
from astropy.io import fits,ascii
from astropy.stats import LombScargle
from scipy import stats
from astropy.table import Table,Column

def Theta(period, times, magnitudes, errors):
    """
    computes the power of periodicity (Laffler-Kinnman method, also known
    as Phase Dispersion Minimization) of a set of magnitude data (in one
    filter) for a given period.
    
    Parameters
    ---------
    period : float
        the period of time which may or may not match the period of the
        given data set
    
    times : numpy array or list
        the list of times for the data
    
    magnitudes : numpy array or list
        the list of magnitudes for the data
    
    errors : numpy array or list
        the list of errors for the data
    
    
    Returns
    -------
    Theta : float
        power or periodicity of that period for the given data set. In
        the Laffler-Kinnman equation this is represented by the symbol
        upper-case Theta
        
    """
    #Phase folds the times array
    phased_times = (times/period)%1
    #sorts the magnitude and error data by phase-folded time instead of linear time
    phase_mag = magnitudes[np.argsort(phased_times)]
    phase_err = errors[np.argsort(phased_times)]
    top_sum = 0
    bottom_sum = 0
    w_i_sum = 0
    #This runs through the the sums in the Laffler-Kinnman equation
    for i in range(1,len(phase_mag)):
        top = phase_mag[i]-phase_mag[i-1]
        bottom = phase_mag[i]-np.mean(phase_mag)
        w_i = (phase_err[i]**2+phase_err[i-1]**2)**(-1)
        top_sum += (top**2)*w_i
        bottom_sum += bottom**2
        w_i_sum += w_i
    #Computes final Theta
    Theta = top_sum/(bottom_sum*w_i_sum)
    return(Theta)
    
    
def hybrid_periodogram(times,magnitudes,errors):
    """
    computes the power or periodicity (a hybrid method which uses both
    the Lomb-Scargle method and the Laffler Kinnman method) of a set of 
    magnitude data (in one filter) over a series of periods.
    
    Parameters
    ---------    
    times : numpy array or list
        the list of times for the data
    
    magnitudes : numpy array or list
        the list of magnitudes for the data
    
    errors : numpy array or list
        the list of errors for the data
    
    
    Returns
    -------
    period : list
        a list of periods for which periodicity was computed
        
    menorah : list
        a list of periodicities based on both the periodicities returned 
        by the Lomb-Scargle and Laffler-Kinnman methods at each period.
        Called menorah because the symbol used to represent this value is an upper-case
        psi, which I believe looks like a menorah with most of the candle
        holders broken off.
        
    """
    #Lomb-Scargle
    LS_frequency,LS_power = LombScargle(times,magnitudes).autopower()
    #Converts frequency to period
    period = 1/LS_frequency
    #Laffler-Kinman
    theta_list = []
    #for each period that was used in LombScargle, this computes Theta
    for p in period:
        theta_list.append(Theta(p,times,magnitudes,errors))
    theta_array = np.array(theta_list)
    #Finding Menorah (otherwise known as upper-case psi)
    menorah = (2*LS_power)/theta_array
    return period,menorah


def polarization(Q,U):
    Q,U = np.array(Q),np.array(U)
    polarization = np.sqrt(Q**2+U**2)
    position_angle = np.degrees(np.arctan2(U,Q))
    return(polarization,position_angle)
    
def easy_bin_mean(x,data,bins_num):
    """
    I REALLY DIDN'T LIKE THE FORMAT SO I JUST STOLE SCIPY'S STATS THING BUT CHANGED
    THE FORMAT SO I COULD STAND IT. ALL CREDIT TO SCIPY"""
    bin_data, bin_x, binnumber = stats.binned_statistic(x,data,statistic='mean',bins=bins_num)
    return(bin_x,bin_data)



def ret_pol_data(txt_file_name,bin_num):
    
    table = Table.read(txt_file_name, format = 'ascii',delimiter='\s',data_start=1,names=['Wavelength','Flux','q','u','err'])
    wavelength = np.array(table["Wavelength"])
    Q = np.array(table["q"])
    U = np.array(table["u"])
    err = np.array(table["err"])
    flux = table["Flux"]
    wavelength_bin_edges,bin_Q = easy_bin_mean(wavelength,Q,bin_num)
    wavelength_bin_edges,bin_U = easy_bin_mean(wavelength,U,bin_num)
    wavelength_bin_edges,bin_flux = easy_bin_mean(wavelength,flux,bin_num)
    pol,pos = polarization(bin_Q,bin_U)
    return(wavelength_bin_edges[1:],bin_flux,pol,pos)

    
def get_ret_ave_pol(fits_file):
    hdu = fits.open(fits_file)
    infohdu = hdu[0].header
    tablehdu = hdu[1].header
    ave_pol,ave_err = tablehdu["Ave-pol"],tablehdu["Ave-erro"]
    obs_time = infohdu["MJD_OBS"]
    return(obs_time,ave_pol,ave_err)
    
def ret_ave_pol_curve(ret_fits_file_list):
    times_list = []
    average_pol_list = []
    average_err_list = []
    for file in ret_fits_file_list:
        obs_time,ave_pol,ave_err = get_ret_ave_pol(file)
        times_list.append(obs_time)
        average_pol_list.append(ave_pol)
        average_err_list.append(ave_err)
    #Sort all of these boys because SOME TIMES DON"T PLAY NICE
    times_list,average_pol_list,average_err_list = np.array(times_list),np.array(average_pol_list),np.array(average_err_list)
    times_sort = np.argsort(times_list)
    times_list = times_list[times_sort]
    average_pol_list = average_pol_list[times_sort]
    average_err_list = average_err_list[times_sort]
    return(times_list,average_pol_list,average_err_list)
        
        
        
        
        
    
    
    
    
    
    
    


