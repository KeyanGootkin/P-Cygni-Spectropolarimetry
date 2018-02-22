import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
from astropy.io import fits,ascii
from astropy.stats import LombScargle
from scipy import stats
from astropy.table import Table,Column

'''TOOLS'''

def dedopler(raw_wavelength,radial_velocity):
    """
    This function takes an array of wavelengths and returns an array of wavelengths
    which have been corrected to account for doppler shift. This uses a value of
    z which assumes speeds much much lower than c.

    Parameters
    ----------
    raw_wavelength: numpy array
        array of observed wavelengths which have some doppler shift.

    radial_velocity: float
        the radial velocity of the object given in km/s

    Returns
    ----------
    calibrated_wavelength : numpy array
        array of wavelengths which have been corrected to account for doppler
        shift caused by the target's radial velocity
    """
    if radial_velocity == 0:
        return(raw_wavelength)
    else:
        z=radial_velocity/(3*10**5)
        calibrated_wavelength = raw_wavelength*(1+z)
        return(calibrated_wavelength)

def Theta(period, times, magnitudes, errors):
    """
    computes the periodicity (Laffler-Kinnman method, also known
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

def make_figure(x,y,err,title):
    """
    Creates an errorbar figure for a single filter (does not work for
    tables with multiple filters) of Magnitude vs Time (in minutes)

    Parameters
    ---------
    data_table : Table
        the table containing the info for the figure

    title : STRING
        name of the figure
    """
    #Crank up the size of the figure, for beauty purposes.
    plt.figure(figsize=[10,5])
    #Creates actual plot
    plt.errorbar(x,y,yerr=err)
    #Sets the title of the plot
    plt.title(title, size = 23, fontname='Times New Roman')
    #Labels the x and y axis.
    plt.gca().set_xlabel('Time [MJD]', size=17, fontname='Times New Roman')
    plt.ylabel('Polarization', size = 17, fontname='Times New Roman')
    #Makes the axis look nice (Thank you Cayenne)
    plt.xticks(fontsize = 15, fontname='Times New Roman', color = 'darkslategrey')
    plt.yticks(fontsize = 15,fontname='Times New Roman', color = 'darkslategrey')
    #Creates a grid on the plot to make it more readable
    plt.grid()
    return True


def make_periodogram_figure(period,power):
    """
    Creates a periodogram figure, from 0 to 30 minutes in period, given
    an array of periods and power for a single filter.

    Parameters
    ---------
    period : numpy array
        an array of periods for which the periodicity was calculated

    power : numpy array
        an array of computed periodicities
    """
    plt.figure(figsize=[10,5]),
    #creates the actual figure
    plt.errorbar(period,power)
    #sets title
    plt.title('Periodogram', size = 23, fontname='Times New Roman')
    #Labels the x-axis then the y-axis
    plt.gca().set_xlabel('Period [MJD]', size=17, fontname='Times New Roman')
    plt.ylabel('Peroidicity', size = 17, fontname='Times New Roman')
    #Makes the axis look nice (Thank you Cayenne)
    plt.xticks(fontsize = 15, fontname='Times New Roman', color = 'darkslategrey')
    plt.yticks(fontsize = 15,fontname='Times New Roman', color = 'darkslategrey')
    #Sets up grid
    plt.grid()
    return True

def make_scatter(x,y,title):
    #Crank up the size of the figure, for beauty purposes.
    plt.figure(figsize=[10,5])
    #Creates actual plot
    plt.scatter(x,y)
    #Sets the title of the plot
    plt.title(title, size = 23, fontname='Times New Roman')
    #Makes the axis look nice (Thank you Cayenne)
    plt.xticks(fontsize = 15, fontname='Times New Roman', color = 'darkslategrey')
    plt.yticks(fontsize = 15,fontname='Times New Roman', color = 'darkslategrey')
    #Creates a grid on the plot to make it more readable
    plt.grid()
    return True


def plot_phase_time(times,variable,period):
    phase_times = (times/period)%1
    make_scatter(phase_times,variable,'Time-Folded - Period='+str(period))
    plt.show()
    return True


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

'''FLUX, POLARIZATION AND POSITION ANGLE FIGURES'''

def txt_pol_data(txt_file_name,bin_num,radial_velocity=0,window=[3000,8000]):
    #extract all info from txt file
    table = ascii.read(txt_file_name, delimiter='\s',encoding='utf8',
                       data_start=1,names=['Wavelength','Flux','q','u','err'])
    wavelength = np.array(table["Wavelength"])
    Q = np.array(table["q"])
    U = np.array(table["u"])
    err = np.array(table["err"])
    flux = table["Flux"]
    #bin up all of the data
    wavelength_bin_edges,bin_Q = easy_bin_mean(wavelength,Q,bin_num)
    wavelength_bin_edges,bin_U = easy_bin_mean(wavelength,U,bin_num)
    wavelength_bin_edges,bin_flux = easy_bin_mean(wavelength,flux,bin_num)
    #calculate polarization and position angle
    pol,pos = polarization(bin_Q,bin_U)
    #get rid of doppler shift
    wavelength = dedopler(wavelength,radial_velocity)
    return(wavelength_bin_edges[1:],bin_flux,pol,pos)

def stack_txt_pol_data(txt_file_name_list,bin_num,radial_velocity=0,window=[3000,8000]):
    plt.figure(figsize=[15,20])
    for files in txt_file_name_list:
        wavelength,flux,pol,pos = txt_pol_data(files,bin_num)
        wavelength = dedopler(np.array(wavelength),radial_velocity)
        plt.subplot(3,1,1)
        plt.plot(wavelength,flux/np.median(flux),label=(files[len(files)-20:len(files)-12]))
        plt.title("Flux")
        plt.xlim(window[0],window[1])
        plt.subplot(3,1,2)
        plt.plot(wavelength,pol,label=(files[len(files)-20:len(files)-12]))
        plt.title("% Polarization")
        plt.xlim(window[0],window[1])
        plt.subplot(3,1,3)
        plt.plot(wavelength,pos,label=(files[len(files)-20:len(files)-12]))
        plt.title("Position Angle")
        plt.xlim(window[0],window[1])

'''AVERAGE POLARIZATION CURVES'''

def get_ave_pol(fits_file):
    hdu = fits.open(fits_file)
    infohdu = hdu[0].header
    tablehdu = hdu[1].header
    ave_pol,ave_err = tablehdu["Ave-pol"],tablehdu["Ave-erro"]
    obs_time = infohdu["MJD_OBS"]
    return(obs_time,ave_pol,ave_err)

def get_ccd_ave_pol(fits_file):
    hdu = fits.open(fits_file)
    infohdu = hdu[0].header
    tablehdu = hdu[1].header
    ave_pol,ave_err = tablehdu["Ave-pol"],tablehdu["Ave-erro"]
    obs_time = infohdu["MJD-OBS"]
    return(obs_time,ave_pol,ave_err)

def ave_pol_curve(fits_file_list):
    #Create empty lists for later
    times_list = []
    average_pol_list = []
    average_err_list = []
    #For every file given, find time,average polarization, and average error
    for files in fits_file_list:
        obs_time,ave_pol,ave_err = get_ave_pol(files)
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

def ccd_ave_pol_curve(fits_file_list):
    #Create empty lists for later
    times_list = []
    average_pol_list = []
    average_err_list = []
    #For every file given, find time,average polarization, and average error
    for files in fits_file_list:
        obs_time,ave_pol,ave_err = get_ccd_ave_pol(files)
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
