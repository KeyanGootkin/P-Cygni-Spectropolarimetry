import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
from astropy.io import fits,ascii
from astropy.stats import LombScargle

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
