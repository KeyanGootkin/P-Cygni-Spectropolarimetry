import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from glob import glob
import os
from astropy.io import fits, ascii
from astropy.stats import LombScargle
from scipy import stats
from astropy.table import Table, Column

'''TOOLS'''


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


def find_continuum_w_gap(wavelength, data, left_1, left_2, right_1, right_2):
    """
    Takes wavelength dependent data and finds the continuum mean in a range with
    a cut-out for a line.

    Parameters
    ------------
    wavelength : list
        list or array of wavelengths

    data : list
        list or array of the wavelength dependent data (ie flux, % pol or pa)

    left_1 : integer or float
        The value of wavelength for the left-most boundary for the region you
        wish to consider.

    left_2 : integer or float
        The value of wavelength for the left edge of the line (or region you wish
        to ignore).

    right_1 : integer or float
        The value of wavelength for the right edge of the line (or region you wish
        to ignore).

    right_2 : integer or float
        The value of wavelength for the right-most boundary for the region you
        wish to consider.

    Returns
    ------------

    ave_continuum : float
        The average value of the data between left_1, left_2 and right_1, right_2
    """
    continuum = []
    w = wavelength
    for i in w:
        if i >= left_1 and i <= left_2:
            continuum.append(data[list(w).index(i)])
        elif i >= right_1 and i <= right_2:
            continuum.append(data[list(w).index(i)])
    ave_continuum = np.mean(np.array(continuum))
    return(ave_continuum)


def divide_continuum(wavelength, data, left_1, left_2, right_1, right_2):
    """
    Takes wavelength dependent data and divides it by the average continuum value.

    Parameters
    ------------

    wavelength : list
        list or array of wavelengths

    data : list
        list or array of the wavelength dependent data (ie flux, % pol or pa)

    left_1 : integer or float
        The value of wavelength for the left-most boundary for the region you
        wish to consider.

    left_2 : integer or float
        The value of wavelength for the left edge of the line (or region you wish
        to ignore).

    right_1 : integer or float
        The value of wavelength for the right edge of the line (or region you wish
        to ignore).

    right_2 : integer or float
        The value of wavelength for the right-most boundary for the region you
        wish to consider.

    Returns
    ------------
    new_data : numpy array
        The parameter array "data", with each value divided by the average data
        value between left_1, left_2 and right_1, right_2. The data divided by
        the continuum value.
    """
    continuum = find_continuum_w_gap(
        wavelength, data, left_1, left_2, right_1, right_2)
    new_data = np.array(data) / continuum
    return(new_data)


def sub_mean(list):
    """
    Takes a list and subtracts the mean from each indice. NOTE: the mean is done
    at axis=0, this may result in the mean being an array instead of a numerical
    value.

    Parameters
    ------------
    list : list or array
        a list or array of integers or floats values.

    Returns
    ------------
    subbed_list : numpy array
        original list where each value has had the average value of the list
        subtracted from it.

    """
    list = np.array(list)
    subbed_list = list - np.mean(list, axis=0)
    return(subbed_list)


def cut_list(list, low=0, high=1):
    """
    Takes a list of values and returns only the values of the list between the
    low and high (inclusive) given as kwargs.

    Parameters
    ------------
    list : list
        A list of numerical values.

    low : integer or float
        The lowest value you wish to have in the new list.

    high : integer or float
        The highest value you wish to have in the new list.

    Returns
    ------------
    newlist : list
        list, but with only values between low and high (inclusive)
    """
    low = low
    high = high
    newlist = []
    for i in list:
        if i >= low and i <= high:
            newlist.append(i)
    return(newlist)


def bin_errors(errors, wavelength, wavelength_bin_edges):
    """
    When binning data with errors, you can't just take the average error in the
    bin and call it good. Instead you have to take the root of the sum of the
    squares value within the bin. This does that! NOTE: Make sure that the
    wavelength values haven't been binned or anything like that, these Should
    be compleetely raw wavelength values.

    Parameters
    ------------
    errors : list or array
        Your error values for each wavelength

    wavelength : list or array
        Your wavelength values (before using the bin edges)

    wavelenght_bin_edges : list or array
        The wavelength values of the edges of each bin.

    Returns
    ------------
    binned_errors : numpy array
        An array of error values which you can use for the binned QU values!
    """

    binned_errors = []
    # For loop cycles through each bin edge
    for i in range(len(wavelength_bin_edges) - 1):
        # the lower edge of each region in which to calculate error
        low = wavelength_bin_edges[i]
        # if the next edge is not the last edge, then that next edge is the upper limit
        if i + 1 <= len(wavelength_bin_edges) - 1:
            high = wavelength_bin_edges[i + 1]
        # If the next edge is the last one, or greater (somehow), then the last one is the upper limit
        # This protects against trying to make the upper limit beyond the range of the list if the lower limit is set as the last item
        else:
            high = wavelength_bin_edges[len(wavelength_bin_edges) - 1]
        temp_err_storage = []
        # For loop cycles through each individual wavelength measurement
        for x in range(len(wavelength)):
            # If that wavelength is between the low and high (low inclusive), add the corresponding error measurement to the temp err temp_err_storage
            if wavelength[x] >= low and wavelength[x] < high:
                temp_err_storage.append(errors[x])
        # Take everything in temp storage, and add the sqrt of the sum of the squares to the binned bin_errors
        # This is the error measurement for one bin
        binned_errors.append(
            (np.sum(np.array(temp_err_storage)**2) / len(temp_err_storage))**0.5)
    # Return the list of errors for each binned region.
    return(binned_errors)


def error_of_median(all_err):
    """
    Takes a list of lists of errors, and correctly calculates the median error
    at each indice. Uses formula for median error from WR_Extraction which was
    from some statistics website which I'm too lazy to find again. But the
    SHPASWRS repository has the link in there somewhere.

    Parameters
    ------------
    all_err : list or array
        A list or array of all error list. This is a nested list with many lists
        of errors across different observations
    Returns
    ------------
    median_err : array
        The median error over a series of different observations.
    """
    median_err = (1.253 * ((sum(np.array(all_err)**2))**0.5)) / \
        np.sqrt(len(all_err))
    return(median_err)


def position_angle_error(pol, err):
    """
    Finds the error of our polarimetric position angle measurements based on
    % polarization and that values error.

    Parameters
    ------------
    pol : float
        % polarization calculated from the Stokes Q and U parameters

    err : float
        The error of pol

    Returns
    ------------
    paerr : float
        The error of our position angle measurements
    """
    paerr = (err / pol) * (90 / np.pi)
    return(paerr)


def dedopler(raw_wavelength, radial_velocity):
    """
    This function takes an array of wavelengths and returns an array of wavelengths
    which have been corrected to account for doppler shift. This uses a value of
    z which assumes speeds much much lower than c.

    Parameters
    ----------
    raw_wavelength : numpy array
        array of observed wavelengths which have some doppler shift.

    radial_velocity : float
        the radial velocity of the object given in km/s

    Returns
    -------
    calibrated_wavelength : numpy array
        array of wavelengths which have been corrected to account for doppler
        shift caused by the target's radial velocity
    """
    if radial_velocity == 0:
        return(raw_wavelength)
    else:
        z = radial_velocity / (3 * 10**5)
        calibrated_wavelength = raw_wavelength / (1 + z)
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
    # Phase folds the times array
    phased_times = (times / period) % 1
    # sorts the magnitude and error data by phase-folded time instead of linear time
    phase_mag = magnitudes[np.argsort(phased_times)]
    phase_err = errors[np.argsort(phased_times)]
    top_sum = 0
    bottom_sum = 0
    w_i_sum = 0
    # This runs through the the sums in the Laffler-Kinnman equation
    for i in range(1, len(phase_mag)):
        top = phase_mag[i] - phase_mag[i - 1]
        bottom = phase_mag[i] - np.mean(phase_mag)
        w_i = (phase_err[i]**2 + phase_err[i - 1]**2)**(-1)
        top_sum += (top**2) * w_i
        bottom_sum += bottom**2
        w_i_sum += w_i
    # Computes final Theta
    Theta = top_sum / (bottom_sum * w_i_sum)
    return(Theta)


def hybrid_periodogram(times, magnitudes, errors, give_ls=False):
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
    # Lomb-Scargle
    LS = LombScargle(times, magnitudes)
    LS_frequency, LS_power = LS.autopower()
    # Converts frequency to period
    period = 1 / LS_frequency
    # Laffler-Kinman
    theta_list = []
    # for each period that was used in LombScargle, this computes Theta
    for p in period:
        theta_list.append(Theta(p, times, magnitudes, errors))
    theta_array = np.array(theta_list)
    # Finding Menorah (otherwise known as upper-case psi)
    menorah = (2 * LS_power) / theta_array
    if give_ls == True:
        return(period, menorah, LS)
    else:
        return period, menorah


def find_best_period(period, power):
    """
    Takes the results of a periodogram and searches it for the period with the
    highest power.

    Parameters
    ------------
    period : list or array
        The list of periods which have been input into the periodogram. These are
        potentially regular periods in the data which are evaluated.

    power : list or array
        The list of powers computed in the periodogram.

    Returns
    ------------
    best_period : float
        The period which had the highest power, the largest "peak" in the periodogram.
    """
    best_period_ind = np.argsort(power)
    best_period = (period[best_period_ind])[len(period) - 1]
    return(best_period)


def make_figure(x, y, err, title="", axis=["", ""]):
    """
    Creates a pretty errorbar figure

    Parameters
    ---------
    x : list
        data for the x axis of the plot

    y : list
        data for the y axis of the plot

    err : list
        list of errors for the data on the y axis

    title : string
        what you would like the title of the plot to be.

    """

    # Creates actual plot
    plt.errorbar(x, y, yerr=err)
    # Sets the title of the plot
    plt.title(title, size=23, fontname='Times New Roman')
    # Labels the x and y axis.
    plt.gca().set_xlabel(axis[0], size=17, fontname='Times New Roman')
    plt.ylabel(axis[1], size=17, fontname='Times New Roman')
    # Makes the axis look nice (Thank you Cayenne)
    plt.xticks(fontsize=15, fontname='Times New Roman', color='darkslategrey')
    plt.yticks(fontsize=15, fontname='Times New Roman', color='darkslategrey')
    # Creates a grid on the plot to make it more readable
    plt.grid()
    return True


def make_periodogram_figure(period, power):
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
    plt.figure(figsize=[10, 5]),
    # creates the actual figure
    plt.errorbar(period, power)
    # sets title
    plt.title('Periodogram', size=23, fontname='Times New Roman')
    # Labels the x-axis then the y-axis
    plt.gca().set_xlabel('Period [MJD]', size=17, fontname='Times New Roman')
    plt.ylabel('Peroidicity', size=17, fontname='Times New Roman')
    # Makes the axis look nice (Thank you Cayenne)
    plt.xticks(fontsize=15, fontname='Times New Roman', color='darkslategrey')
    plt.yticks(fontsize=15, fontname='Times New Roman', color='darkslategrey')
    # Sets up grid
    plt.grid()
    return True


def make_scatter(x, y, title=""):
    """
    Creates a pretty scatter plot

    Parameters
    ---------
    x : list
        data for the x axis of the plot

    y : list
        data for the y axis of the plot

    title : string
        what you would like the title of the plot to be.

    """
    # Crank up the size of the figure, for beauty purposes.
    plt.figure(figsize=[10, 5])
    # Creates actual plot
    plt.scatter(x, y)
    # Sets the title of the plot
    plt.title(title, size=23, fontname='Times New Roman')
    # Makes the axis look nice (Thank you Cayenne)
    plt.xticks(fontsize=15, fontname='Times New Roman', color='darkslategrey')
    plt.yticks(fontsize=15, fontname='Times New Roman', color='darkslategrey')
    # Creates a grid on the plot to make it more readable
    plt.grid()
    return True


def plot_phase_time(times, variable, period):
    """
    takes a set of time-series data, phase-folds it, then plots it as a scatter
    plot of the variable vs. phase.

    Parameters
    ---------
    times : list
        time data for the x axis of the plot

    variable : list
        data for the y axis of the plot

    period : float
        The period which you would like to fold the data by.
    """
    phase_times = (times / period) % 1
    make_scatter(phase_times, variable,
                 title='Time-Folded - Period=' + str(period))
    plt.show()
    return True


def polarization(Q, U):
    """
    computes % Polarization and Position angle based on lists of the stokes
    parameters Q and U.

    Parameters
    ---------
    Q : list
        A list of the stokes parameter Q.

    U : list
        A list of the stokes parameter U.


    Returns
    -------
    polarization : numpy array
        An array of % polarizations

    position_angle : numpy array
        An array of position angles
    """
    Q, U = np.array(Q), np.array(U)
    polarization = np.sqrt(Q**2 + U**2)
    position_angle = np.degrees(np.arctan2(U, Q))
    return(polarization, position_angle)


def easy_bin_mean(x, data, bins_num, type='mean'):
    """
    This takes data and bins it, this is literally just astropy's binned_statistic
    but in a format that is easier to understand and use.

    Parameters
    ---------
    x : list
        x axis to be binned.

    data : list
        A list of data to be grouped together and binned.

    bins_num : integer
        The number of bins that you wish to have.


    Returns
    -------
    bin_x : list
        the bin edges in the x data

    bin_data : list
        List of data which has been grouped into bins and averaged.
    """
    if type == 'mean':
        bin_data, bin_x, binnumber = stats.binned_statistic(
            x, data, statistic='mean', bins=bins_num)
    elif type == 'errors':
        bin_data, bin_x, binnumber = stats.binned_statistic(
            x, data, statistic=bin_error(data), bins=bins_num)
    else:
        print("ERROR: the function easy_bin_mean has not been passed a valid type.")
    return(bin_x, bin_data)


def match_rb_txt_files(rtxtfiles, btxtfiles):
    """
    Takes the HPOL data txt files and matches the red and blue observations that
    are from the same night.

    Parameters
    ------------
    rtxtfiles : list
        list of red observation txt files

    btxtfiles : list
        list of blue observation txt files

    Returns
    ------------
    matched_files : list
        A list where each entry is a list, [r, b], where is r is the red txt file
        name and b is the blue txt file name.
    """
    matched_files = []
    for r in rtxtfiles:
        for b in btxtfiles:
            if r[len(r) - 21:len(r) - 13] == b[len(b) - 21:len(b) - 13]:
                matched_files.append([r, b])
    return(matched_files)


def match_rb_fits_files(rfitsfiles, bfitsfiles):
    """
    Takes the HPOL data fits files and matches the red and blue observations that
    are from the same night.

    Parameters
    ------------
    rfitsfiles : list
        list of red observation fits files

    bfitsfiles : list
        list of blue observation fits files

    Returns
    ------------
    matched_files : list
        A list where each entry is a list, [r, b], where r is the red fits file
        name and b is the blue fits file name.
    """
    matched_files = []
    for r in rfitsfiles:
        for b in bfitsfiles:
            if r[len(r) - 17:len(r) - 9] == b[len(b) - 17:len(b) - 9]:
                matched_files.append([r, b])
    return(matched_files)


'''FLUX, POLARIZATION AND POSITION ANGLE FIGURES'''


def txt_pol_data(txt_file_name, bin_num, radial_velocity=0):
    """
    Extracts, calculates, and bins wavelength, flux, % polarization, position
    angle, and error from the HPOL txt files.

    Parameters
    ------------
    txt_file_name : string
        The string of the pathname to the txt file you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    Returns
    ------------
    wavelength_bin_edges[1:] : list
        The edges of wavelength bins starting at the second entry, this means that
        all wavelength dependent data will be asscosiated with the right edge of
        the bin they are in.

    bin_flux : list
        list of binned flux values

    pol : list
        list of % polarization at each wavelenght which have been calculated from
        binned values of Q and U stokes parameters

    pos : list
        list of postion angles at each wavelength which have been calculated form
        binned values of Q and U stokes parameters

    bin_err : list
        list of binned errors for the Q and U stokes parameters
    """
    if txt_file_name[len(txt_file_name) - 1] == 's':
        print("ERROR: fits and txt files have been switched! Check parameters.")
    # extract all info from txt file
    table = np.genfromtxt(txt_file_name, skip_header=1, names=[
                          'Wavelength', 'Flux', 'q', 'u', 'err'])
    wavelength = np.array(table["Wavelength"])
    Q = np.array(table["q"])
    U = np.array(table["u"])
    err = np.array(table["err"])
    flux = table["Flux"]
    # bin up all of the data
    wavelength_bin_edges, bin_Q = easy_bin_mean(wavelength, Q, bin_num)
    wavelength_bin_edges, bin_U = easy_bin_mean(wavelength, U, bin_num)
    wavelength_bin_edges, bin_flux = easy_bin_mean(wavelength, flux, bin_num)
    bin_err = bin_errors(err, wavelength, wavelength_bin_edges)
    # calculate polarization and position angle
    pol, pos = polarization(bin_Q, bin_U)
    # get rid of doppler shift
    wavelength = dedopler(wavelength, radial_velocity)
    return(wavelength_bin_edges[1:], bin_flux, pol, pos, bin_err)


def txt_QU_data(txt_file_name, bin_num, radial_velocity=0):
    """
    Extracts wavelength, flux, Q, U, and error from the given file

    Parameters
    ----------
    txt_file_name : string
        The string of the pathname to the txt file you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    Returns
    ----------
    wavelength_bin_edges[1:] : list
        The edges of wavelength bins starting at the second entry, this means that
        all wavelength dependent data will be asscosiated with the right edge of
        the bin they are in.

    bin_flux : list
        list of binned flux values

    bin_Q : list
        list of binned stokes Q parameter

    bin_U : list
        list of binned stokes U parameter

    bin_err : list
        list of binned errors for the Q and U stokes parameters
    """
    # extract all info from txt file
    table = np.genfromtxt(txt_file_name, skip_header=1, names=[
                          'Wavelength', 'Flux', 'q', 'u', 'err'])
    wavelength = np.array(table["Wavelength"])
    Q = np.array(table["q"])
    U = np.array(table["u"])
    err = np.array(table["err"])
    flux = table["Flux"]
    # bin up all of the data
    wavelength_bin_edges, bin_Q = easy_bin_mean(wavelength, Q, bin_num)
    wavelength_bin_edges, bin_U = easy_bin_mean(wavelength, U, bin_num)
    wavelength_bin_edges, bin_flux = easy_bin_mean(wavelength, flux, bin_num)
    bin_err = bin_errors(err, wavelength, wavelength_bin_edges)
    # get rid of doppler shift
    wavelength = dedopler(wavelength, radial_velocity)
    return(wavelength_bin_edges[1:], bin_flux, bin_Q, bin_U, bin_err)


def stack_txt_pol_data(txt_file_list, fits_file_list, bin_num, radial_velocity=0, window=[3000, 8000], vert_line=False):
    """
    Takes the files you give it and graphs flux, % pol, and position angle all
    vs. wavelength for each file. The files are stacked on eachother and color
    coded so that the lighter the line, the newer the observation

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    window : list
        A 2 piece list where the first is the lower x_lim and the second is the
        upper x_lim you want

    vert_line : Integer or Float
        If you would like a verticle line, then enter the x value you would like
        the line to appear at
    """
    all_wavelength = []
    all_flux = []
    all_pol = []
    all_pos = []
    all_MJD = []
    cmap = cm.get_cmap('inferno')
    plt.figure(figsize=[15, 20])
    for i in range(len(txt_file_list)):
        txt, fits = txt_file_list[i], fits_file_list[i]
        wavelength, flux, pol, pos, err = txt_pol_data(txt, bin_num)
        wavelength = dedopler(wavelength, radial_velocity)
        if fits[len(fits) - 26:len(fits) - 23] == 'ret':
            MJD, ave_pol, ave_err = get_ave_pol(fits)
        elif fits[len(fits) - 27:len(fits) - 24] == 'ccd':
            MJD, ave_pol, ave_err = get_ccd_ave_pol(fits)
        else:
            print("Sorry, I can't tell if this is a ccd image or a ret image.")
            break
        all_wavelength.append(wavelength)
        all_flux.append(flux)
        all_pol.append(pol)
        all_pos.append(pos)
        all_MJD.append(MJD)
    for wavelength, flux, pol, pos, MJD in zip(all_wavelength, all_flux, all_pol, all_pos, all_MJD):
        plt.subplot(3, 1, 1)
        plt.plot(wavelength, flux / np.median(flux), label=(txt[len(txt) - 20:len(
            txt) - 12]), c=cmap((MJD - min(all_MJD)) / (max(all_MJD) - min(all_MJD))))
        plt.plot([vert_line, vert_line], [0, 6])
        plt.title("Flux")
        plt.xlim(window[0], window[1])
        plt.subplot(3, 1, 2)
        plt.plot(wavelength, pol, label=(txt[len(txt) - 20:len(txt) - 12]), c=cmap(
            (MJD - min(all_MJD)) / (max(all_MJD) - min(all_MJD))))
        plt.plot([vert_line, vert_line], [min(pol), max(pol)])
        plt.title("% Polarization")
        plt.xlim(window[0], window[1])
        plt.subplot(3, 1, 3)
        plt.plot(wavelength, pos, label=(txt[len(txt) - 20:len(txt) - 12]), c=cmap(
            (MJD - min(all_MJD)) / (max(all_MJD) - min(all_MJD))))
        plt.plot([vert_line, vert_line], [min(pos), max(pos)])
        plt.title("Position Angle")
        plt.xlim(window[0], window[1])


def stack_ccd_txt_pol_data(txt_file_list, fits_file_list, bin_num, radial_velocity=0, window=[2800, 11000], vert_line=False):
    """
    Takes the files you give it and graphs flux, % pol, and position angle all
    vs. wavelength for each file. The files are stacked on eachother and color
    coded so that the lighter the line, the newer the observation. Same as above
    but for ccd files specifically

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    window : list
        A 2 piece list where the first is the lower x_lim and the second is the
        upper x_lim you want

    vert_line : Integer or Float
        If you would like a verticle line, then enter the x value you would like
        the line to appear at
    """
    all_wavelength = []
    all_flux = []
    all_pol = []
    all_pos = []
    all_MJD = []
    cmap = cm.get_cmap('viridis')
    plt.figure(figsize=[15, 20])
    for i in range(len(txt_file_list)):
        txt, fits = txt_file_list[i], fits_file_list[i]
        wavelength, flux, pol, pos, err = txt_pol_data(txt, bin_num)
        wavelength = dedopler(wavelength, radial_velocity)
        if fits[len(fits) - 26:len(fits) - 23] == 'ret':
            MJD, ave_pol, ave_err = get_ave_pol(fits)
        elif fits[len(fits) - 27:len(fits) - 24] == 'ccd':
            MJD, ave_pol, ave_err = get_ccd_ave_pol(fits)
        else:
            print("Sorry, I can't tell if this is a ccd image or a ret image.")
            break
        all_wavelength.append(wavelength)
        all_flux.append(flux)
        all_pol.append(pol)
        all_pos.append(pos)
        all_MJD.append(MJD)
    both_color_wave, r_color_wave, b_color_wave = [], [], []
    both_color_flux, r_color_flux, b_color_flux = [], [], []
    both_color_pol, r_color_pol, b_color_pol = [], [], []
    both_color_pos, r_color_pos, b_color_pos = [], [], []
    both_color_MJD, r_color_MJD, b_color_MJD = [], [], []

    for wavelength, flux, pol, pos, MJD in zip(all_wavelength, all_flux, all_pol, all_pos, all_MJD):
        plt.subplot(3, 1, 1)
        plt.plot(wavelength, flux, label=(txt[len(txt) - 20:len(txt) - 12]), c=cmap(
            (MJD - min(all_MJD)) / (max(all_MJD) - min(all_MJD))))
        plt.plot([vert_line, vert_line], [0, 6e-12])
        plt.title("Flux")
        plt.xlim(window[0], window[1])
        plt.subplot(3, 1, 2)
        plt.plot(wavelength, pol, label=(txt[len(txt) - 20:len(txt) - 12]), c=cmap(
            (MJD - min(all_MJD)) / (max(all_MJD) - min(all_MJD))))
        plt.plot([vert_line, vert_line], [min(pol), max(pol)])
        plt.title("% Polarization")
        plt.xlim(window[0], window[1])
        plt.subplot(3, 1, 3)
        plt.plot(wavelength, pos, label=(txt[len(txt) - 20:len(txt) - 12]), c=cmap(
            (MJD - min(all_MJD)) / (max(all_MJD) - min(all_MJD))))
        plt.plot([vert_line, vert_line], [min(pos), max(pos)])
        plt.title("Position Angle")
        plt.xlim(window[0], window[1])


def vert_stack_ccd_txt_pol_data(txt_file_list, fits_file_list, bin_num, radial_velocity=0, window=[2800, 11000], vert_line=False):
    """
    Takes the files you give it and graphs flux, % pol, and position angle all
    vs. wavelength for each file. The files are stacked one on top of eachother
    and color coded so that the lighter the line, the newer the observation

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    window : list
        A 2 piece list where the first is the lower x_lim and the second is the
        upper x_lim you want

    vert_line : Integer or Float
        If you would like a verticle line, then enter the x value you would like
        the line to appear at
    """
    all_wavelength = []
    all_flux = []
    all_pol = []
    all_pos = []
    all_MJD = []
    cmap = cm.get_cmap('viridis')

    for i in range(len(txt_file_list)):
        txt, fits = txt_file_list[i], fits_file_list[i]
        wavelength, flux, pol, pos, err = txt_pol_data(txt, bin_num)
        wavelength = dedopler(wavelength, radial_velocity)
        if fits[len(fits) - 26:len(fits) - 23] == 'ret':
            MJD, ave_pol, ave_err = get_ave_pol(fits)
        elif fits[len(fits) - 27:len(fits) - 24] == 'ccd':
            MJD, ave_pol, ave_err = get_ccd_ave_pol(fits)
        else:
            print("Sorry, I can't tell if this is a ccd image or a ret image.")
            break
        all_wavelength.append(wavelength)
        all_flux.append(flux)
        all_pol.append(pol)
        all_pos.append(pos)
        all_MJD.append(MJD)
    both_color_wave, r_color_wave, b_color_wave = [], [], []
    both_color_flux, r_color_flux, b_color_flux = [], [], []
    both_color_pol, r_color_pol, b_color_pol = [], [], []
    both_color_pos, r_color_pos, b_color_pos = [], [], []
    both_color_MJD, r_color_MJD, b_color_MJD = [], [], []
    i = 0
    plt.figure(figsize=[15, 20])
    for wavelength, flux, pol, pos, MJD in zip(all_wavelength, all_flux, all_pol, all_pos, all_MJD):
        plt.subplot(3, 1, 1)
        plt.plot(wavelength, flux + ((i * 10**-10) / 3), label=(txt[len(txt) - 20:len(txt) - 12]), c=cmap(
            (MJD - min(all_MJD)) / (max(all_MJD) - min(all_MJD))))
        plt.plot([vert_line, vert_line], [0, 6e-12])
        plt.title("Flux")
        plt.xlim(window[0], window[1])
        plt.subplot(3, 1, 2)
        plt.plot(wavelength, pol + (i / 10), label=(txt[len(txt) - 20:len(txt) - 12]), c=cmap(
            (MJD - min(all_MJD)) / (max(all_MJD) - min(all_MJD))))
        plt.plot([vert_line, vert_line], [min(pol), max(pol)])
        plt.title("% Polarization")
        plt.xlim(window[0], window[1])
        plt.subplot(3, 1, 3)
        plt.plot(wavelength, pos + (5 * i), label=(txt[len(txt) - 20:len(txt) - 12]), c=cmap(
            (MJD - min(all_MJD)) / (max(all_MJD) - min(all_MJD))))
        plt.plot([vert_line, vert_line], [min(pos), max(pos)])
        plt.title("Position Angle")
        plt.xlim(window[0], window[1])
        i += 1


'''AVERAGE POLARIZATION CURVES'''


def get_ave_pol(fits_file):
    """
    Takes the fits file and extracts time, average polarization and error

    Parameters
    ----------
    fits_file : string
        The pathname to the fits file you wish to extract from

    Returns
    ----------
    obs_time : Float
        MJD time of the observation in days

    ave_pol : Float
        Mean % polarization over all wavelengths. NOTE: I think this includes
        the messed up bits towards the edge, so don't trust this measurement too
        much

    ave_err : Float
        error on ave_pol
    """
    hdu = fits.open(fits_file)
    infohdu = hdu[0].header
    tablehdu = hdu[1].header
    ave_pol, ave_err = tablehdu["Ave-pol"], tablehdu["Ave-erro"]
    obs_time = infohdu["MJD_OBS"]
    return(obs_time, ave_pol, ave_err)


def get_fits_table(fits_file):
    """
    Takes the fits file and extracts wavelength, Q, U, and error

    Parameters
    ----------
    fits_file : string
        The pathname to the fits file you wish to extract from

    Returns
    ----------
    wavelength : list
        List of wavelength values for this observation

    Q : list
        The stokes Q parameter over the wavelength values above

    U : list
        The stokes U parameter over the wavelength values above

    err : list
        list of errors on the Q and U measurements
    """

    hdu = fits.open(fits_file)
    table = hdu[1].data
    wavelength = table['wavelength'][0]
    Q = table['q'][0]
    U = table['u'][0]
    err = table['error'][0]
    return(wavelength, Q, U, err)


def get_ccd_ave_pol(fits_file):
        """
        Takes the fits file and extracts time, average polarization and error
        now in ccd flavor!

        Parameters
        ----------
        fits_file : string
            The pathname to the ccd fits file you wish to extract from

        Returns
        ----------
        obs_time : Float
            MJD time of the observation in days

        ave_pol : Float
            Mean % polarization over all wavelengths. NOTE: I think this includes
            the messed up bits towards the edge, so don't trust this measurement too
            much

        ave_err : Float
            error on ave_pol
        """
    hdu = fits.open(fits_file)
    infohdu = hdu[0].header
    tablehdu = hdu[1].header
    ave_pol, ave_err = tablehdu["Ave-pol"], tablehdu["Ave-erro"]
    obs_time = infohdu["MJD-OBS"]
    return(obs_time, ave_pol, ave_err)


def ave_pol_curve(fits_file_list):
    """
    Takes the fits files and extracts time, average polarization and error

    Parameters
    ----------
    fits_file_list : list
        A list where each indice is the pathname to the ccd fits file you wish
        to extract from

    Returns
    ----------
    times_list : list
        list where each indice is MJD time of the observation in days

    average_pol_list : list
        list where each indice is Mean % polarization over all wavelengths.
        NOTE: I think this includes the messed up bits towards the edge, so don't
        trust this measurement too much

    average_err_list : list
        list where each indice is error on the average pol measurement
    """
    # Create empty lists for later
    times_list = []
    average_pol_list = []
    average_err_list = []
    # For every file given, find time,average polarization, and average error
    for files in fits_file_list:
        obs_time, ave_pol, ave_err = get_ave_pol(files)
        times_list.append(obs_time)
        average_pol_list.append(ave_pol)
        average_err_list.append(ave_err)
    # Sort all of these boys because SOME TIMES DON"T PLAY NICE
    times_list, average_pol_list, average_err_list = np.array(
        times_list), np.array(average_pol_list), np.array(average_err_list)
    times_sort = np.argsort(times_list)
    times_list = times_list[times_sort]
    average_pol_list = average_pol_list[times_sort]
    average_err_list = average_err_list[times_sort]
    return(times_list, average_pol_list, average_err_list)


def ccd_ave_pol_curve(fits_file_list):
    """
    Takes the fits files and extracts time, average polarization and error

    Parameters
    ----------
    fits_file_list : list
        A list where each indice is the pathname to the ccd fits file you wish
        to extract from

    Returns
    ----------
    times_list : list
        list where each indice is MJD time of the observation in days

    average_pol_list : list
        list where each indice is Mean % polarization over all wavelengths.
        NOTE: I think this includes the messed up bits towards the edge, so don't
        trust this measurement too much

    average_err_list : list
        list where each indice is error on the average pol measurement
    """
    # Create empty lists for later
    times_list = []
    average_pol_list = []
    average_err_list = []
    # For every file given, find time,average polarization, and average error
    for files in fits_file_list:
        obs_time, ave_pol, ave_err = get_ccd_ave_pol(files)
        times_list.append(obs_time)
        average_pol_list.append(ave_pol)
        average_err_list.append(ave_err)
    # Sort all of these boys because SOME TIMES DON"T PLAY NICE
    times_list, average_pol_list, average_err_list = np.array(
        times_list), np.array(average_pol_list), np.array(average_err_list)
    times_sort = np.argsort(times_list)
    times_list = times_list[times_sort]
    average_pol_list = average_pol_list[times_sort]
    average_err_list = average_err_list[times_sort]
    return(times_list, average_pol_list, average_err_list)


def calculate_average_polarization_ret(ret_fits_file_list):
    """
    A more robust and accurate way of measuring average pol and pos

    Parameters
    ----------
    ret_fits_file_list : list
        A list where each indice is the pathname to the ret fits file you wish
        to extract from

    Returns
    ----------
    times : list
        list where each indice is MJD time of the observation in days

    ave_pol_list : list
        list where each indice is Mean % polarization over all wavelengths.

    ave_pos_list : list
        list where each indice is Mean position angle over all wavelengths.

    ave_err_list : list
        list where each indice is error on the average pol measurement
    """
    times = []
    ave_pol_list = []
    ave_pos_list = []
    ave_err_list = []
    for ret in ret_fits_file_list:
        obs_time, ave_pol, ave_err = get_ave_pol(ret)
        wavelength, q, u, err = get_fits_table(ret)
        pol, pos = polarization(q, u)
        good_ind = []
        for i in wavelength:
            if i >= 4000 and i <= 7000:
                good_ind.append(True)
            else:
                good_ind.append(False)
        good_pol, good_pos = pol[good_ind], pos[good_ind]
        mean_pol = np.mean(good_pol)
        mean_pos = np.mean(good_pos)
        times.append(obs_time)
        ave_pol_list.append(mean_pol)
        ave_pos_list.append(mean_pos)
        ave_err_list.append(ave_err)
    times, ave_pol_list, ave_pos_list, ave_err_list = np.array(times), np.array(
        ave_pol_list), np.array(ave_pos_list), np.array(ave_err_list)
    times_sort = np.argsort(times)
    times = times[times_sort]
    ave_pos_list = ave_pos_list[times_sort]
    ave_pol_list = ave_pol_list[times_sort]
    ave_err_list = ave_err_list[times_sort]
    return(times, ave_pol_list, ave_pos_list, ave_err_list)


"""THIS NEEDS TO BE CHANGED, BOUNDS ARE WRONG"""


def calculate_average_polarization_ccd(ccd_fits_file_list):
    """
    A more robust and accurate way of measuring average pol and pos

    Parameters
    ----------
    ccd_fits_file_list : list
        A list where each indice is the pathname to the ccd fits file you wish
        to extract from

    Returns
    ----------
    times : list
        list where each indice is MJD time of the observation in days

    ave_pol_list : list
        list where each indice is Mean % polarization over all wavelengths.

    ave_pos_list : list
        list where each indice is Mean position angle over all wavelengths.

    ave_err_list : list
        list where each indice is error on the average pol measurement
    """
    times = []
    ave_pol_list = []
    ave_pos_list = []
    ave_err_list = []
    for ccd in ccd_fits_file_list:
        obs_time, ave_pol, ave_err = get_ccd_ave_pol(ccd)
        wavelength, q, u, err = get_fits_table(ccd)
        pol, pos = polarization(q, u)
        good_ind = []
        for i in wavelength:
            if i >= 4000 and i <= 7000:
                good_ind.append(True)
            else:
                good_ind.append(False)

        good_pol, good_pos = pol[good_ind], pos[good_ind]
        mean_pol = np.mean(good_pol)
        mean_pos = np.mean(good_pos)
        times.append(obs_time)
        ave_pol_list.append(mean_pol)
        ave_pos_list.append(mean_pos)
        ave_err_list.append(ave_err)
    times, ave_pol_list, ave_pos_list, ave_err_list = np.array(times), np.array(
        ave_pol_list), np.array(ave_pos_list), np.array(ave_err_list)
    times_sort = np.argsort(times)
    times = times[times_sort]
    ave_pos_list = ave_pos_list[times_sort]
    ave_pol_list = ave_pol_list[times_sort]
    ave_err_list = ave_err_list[times_sort]
    return(times, ave_pol_list, ave_pos_list, ave_err_list)


"""MEDIANING ALL OF THE THINGS"""


def get_all_fpp(txt_file_list, fits_file_list, bin_num, radial_velocity=0):
    """
    Takes the files you give it and extracts wavelength, flux, % Polarization,
    position angle, and error for each file. These values are stored in a large
    list where each indice is the wavelength or error or postion angle for one
    observation.

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    Returns
    ----------
    all_wavelength : List
        A list of wavelength arrays for each observation in the list of files
        given

    all_flux : List
        A list of flux arrays for each observation in the list of files
        given

    all_pol : List
        A list of % polarization arrays for each observation in the list of files
        given

    all_pos : List
        A list of postion angle arrays for each observation in the list of files
        given

    all_err : List
        A list of error arrays for each observation in the list of files
        given
    """
    all_wavelength = []
    all_flux = []
    all_pol = []
    all_pos = []
    all_err = []
    for txt, fits in zip(txt_file_list, fits_file_list):
        wavelength, flux, pol, pos, err = txt_pol_data(txt, bin_num)
        wavelength = dedopler(wavelength, radial_velocity)
        all_wavelength.append(wavelength)
        all_flux.append(flux)
        all_pol.append(pol)
        all_pos.append(pos)
        all_err.append(err)
    return(all_wavelength, all_flux, all_pol, all_pos, all_err)


def get_all_QU(txt_file_list, fits_file_list, bin_num, radial_velocity=0):
    """
    Takes the files you give it and extracts wavelength, flux, Stokes Q, Stokes
    U, and error for each file. These values are stored in a large list where
    each indice is the wavelength or error or postion angle for one observation.

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    Returns
    ----------
    all_wavelength : List
        A list of wavelength arrays for each observation in the list of files
        given

    all_flux : List
        A list of flux arrays for each observation in the list of files
        given

    all_Q : List
        A list of Stokes Q arrays for each observation in the list of files
        given

    all_U : List
        A list of Stokes U arrays for each observation in the list of files
        given

    all_err : List
        A list of error arrays for each observation in the list of files
        given
    """
    all_wavelength = []
    all_flux = []
    all_Q = []
    all_U = []
    all_err = []
    for txt, fits in zip(txt_file_list, fits_file_list):
        wavelength, flux, Q, U, err = txt_QU_data(txt, bin_num)
        wavelength = dedopler(wavelength, radial_velocity)
        all_wavelength.append(wavelength)
        all_flux.append(flux)
        all_Q.append(Q)
        all_U.append(U)
        all_err.append(err)
    return(all_wavelength, all_flux, all_Q, all_U, all_err)


def median_flux_pol_pos(txt_file_list, fits_file_list, bin_num, radial_velocity=0):
    """
    Takes the files you give it and extracts wavelength, flux, % Polarization,
    position angle, and error for each file. The median of each wavelength
    dependent value is then calculated for at each wavelength value across all
    the observations.

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    Returns
    ----------
    all_wavelength[0] : List
        An array of wavelengths across which the other values can be compared

    median_flux : List
        A list of median flux values across the values of wavelength

    median_pol : List
        A list of median % polarization values across the values of wavelength

    median_pos : List
        A list of median postion angle values across the values of wavelength

    median_err : List
        A list of median error values across the values of wavelength
    """

    all_wavelength, all_flux, all_pol, all_pos, all_err = get_all_fpp(
        txt_file_list, fits_file_list, bin_num, radial_velocity=radial_velocity)
    interp_flux = []
    interp_pol = []
    interp_pos = []
    interp_err = []
    for i in range(1, len(all_wavelength)):
        interpolated_flux = np.interp(
            all_wavelength[0], all_wavelength[i], all_flux[i])
        interp_flux.append(interpolated_flux)
        interpolated_pol = np.interp(
            all_wavelength[0], all_wavelength[i], all_pol[i])
        interp_pol.append(interpolated_pol)
        interpolated_pos = np.interp(
            all_wavelength[0], all_wavelength[i], all_pos[i])
        interp_pos.append(interpolated_pos)
        interpolated_err = np.interp(
            all_wavelength[0], all_wavelength[i], all_err[i])
        interp_err.append(interpolated_err)
    median_flux = np.median(np.array(interp_flux), axis=0)
    median_pol = np.median(np.array(interp_pol), axis=0)
    median_pos = np.median(np.array(interp_pos), axis=0)
    median_err = error_of_median(all_err)
    return(all_wavelength[0], median_flux, median_pol, median_pos, median_err)


def median_flux_Q_U(txt_file_list, fits_file_list, bin_num, radial_velocity=0):
    """
    Takes the files you give it and extracts wavelength, flux, Stokes Q, Stokes
    U, and error for each file. The median of each wavelength dependent value is
    then calculated for at each wavelength value across all the observations.

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    Returns
    ----------
    all_wavelength[0] : List
        An array of wavelengths across which the other values can be compared

    median_flux : List
        A list of median flux values across the values of wavelength

    median_Q : List
        A list of median Stokes Q values across the values of wavelength

    median_U : List
        A list of median Stokes U values across the values of wavelength

    median_err : List
        A list of median error values across the values of wavelength
    """
    all_wavelength, all_flux, all_Q, all_U, all_err = get_all_QU(
        txt_file_list, fits_file_list, bin_num, radial_velocity=radial_velocity)
    interp_flux = []
    interp_Q = []
    interp_U = []
    interp_err = []
    for i in range(1, len(all_wavelength)):
        interpolated_flux = np.interp(
            all_wavelength[0], all_wavelength[i], all_flux[i])
        interp_flux.append(interpolated_flux)
        interpolated_Q = np.interp(
            all_wavelength[0], all_wavelength[i], all_Q[i])
        interp_Q.append(interpolated_Q)
        interpolated_U = np.interp(
            all_wavelength[0], all_wavelength[i], all_U[i])
        interp_U.append(interpolated_U)
        interpolated_err = np.interp(
            all_wavelength[0], all_wavelength[i], all_err[i])
        interp_err.append(interpolated_err)
    median_flux = np.median(np.array(interp_flux), axis=0)
    median_Q = np.median(np.array(interp_Q), axis=0)
    median_U = np.median(np.array(interp_U), axis=0)
    median_err = error_of_median(all_err)
    return(all_wavelength[0], median_flux, median_Q, median_U, median_err)


def mean_flux_pol_pos(txt_file_list, fits_file_list, bin_num, radial_velocity=0):
    """
    Takes the files you give it and extracts wavelength, flux, % Polarization,
    position angle, and error for each file. The mean of each wavelength
    dependent value is then calculated for at each wavelength value across all
    the observations.

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    Returns
    ----------
    all_wavelength[0] : List
        An array of wavelengths across which the other values can be compared

    mean_flux : List
        A list of mean flux values across the values of wavelength

    mean_pol : List
        A list of mean % polarization values across the values of wavelength

    mean_pos : List
        A list of mean postion angle values across the values of wavelength

    mean_err : List
        A list of mean error values across the values of wavelength
    """
    all_wavelength, all_flux, all_pol, all_pos, all_err = get_all_fpp(
        txt_file_list, fits_file_list, bin_num, radial_velocity=radial_velocity)
    interp_flux = []
    interp_pol = []
    interp_pos = []
    interp_err = []
    for i in range(1, len(all_wavelength)):
        interpolated_flux = np.interp(
            all_wavelength[0], all_wavelength[i], all_flux[i])
        interp_flux.append(interpolated_flux)
        interpolated_pol = np.interp(
            all_wavelength[0], all_wavelength[i], all_pol[i])
        interp_pol.append(interpolated_pol)
        interpolated_pos = np.interp(
            all_wavelength[0], all_wavelength[i], all_pos[i])
        interp_pos.append(interpolated_pos)
        interpolated_err = np.interp(
            all_wavelength[0], all_wavelength[i], all_err[i])
        interp_err.append(interpolated_err)
    mean_flux = np.mean(np.array(interp_flux), axis=0)
    mean_pol = np.mean(np.array(interp_pol), axis=0)
    mean_pos = np.mean(np.array(interp_pos), axis=0)
    mean_err = np.sqrt(np.sum(np.array(all_err)**2, axis=0) / len(all_err))
    return(all_wavelength[0], mean_flux, mean_pol, mean_pos, mean_err)


def mean_flux_Q_U(txt_file_list, fits_file_list, bin_num, radial_velocity=0):
    """
    Takes the files you give it and extracts wavelength, flux, Stokes Q, Stokes
    U, and error for each file. The mean of each wavelength dependent value is
    then calculated for at each wavelength value across all the observations.

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    Returns
    ----------
    all_wavelength[0] : List
        An array of wavelengths across which the other values can be compared

    mean_flux : List
        A list of mean flux values across the values of wavelength

    mean_Q : List
        A list of mean Stokes Q values across the values of wavelength

    mean_U : List
        A list of mean Stokes U values across the values of wavelength

    mean_err : List
        A list of mean error values across the values of wavelength
    """
    all_wavelength, all_flux, all_Q, all_U, all_err = get_all_QU(
        txt_file_list, fits_file_list, bin_num, radial_velocity=radial_velocity)
    interp_flux = []
    interp_Q = []
    interp_U = []
    interp_err = []
    for i in range(1, len(all_wavelength)):
        interpolated_flux = np.interp(
            all_wavelength[0], all_wavelength[i], all_flux[i])
        interp_flux.append(interpolated_flux)
        interpolated_Q = np.interp(
            all_wavelength[0], all_wavelength[i], all_Q[i])
        interp_Q.append(interpolated_Q)
        interpolated_U = np.interp(
            all_wavelength[0], all_wavelength[i], all_U[i])
        interp_U.append(interpolated_U)
        interpolated_err = np.interp(
            all_wavelength[0], all_wavelength[i], all_err[i])
        interp_err.append(interpolated_err)
    mean_flux = np.mean(np.array(interp_flux), axis=0)
    mean_Q = np.mean(np.array(interp_Q), axis=0)
    mean_U = np.mean(np.array(interp_U), axis=0)
    mean_err = np.sqrt(np.sum(np.array(all_err)**2, axis=0) / len(all_err))
    return(all_wavelength[0], mean_flux, mean_Q, mean_U, mean_err)


def fpp_sub_mean(txt_file_list, fits_file_list, bin_num, radial_velocity=0):
    """
    Takes the files you give it and extracts wavelength, flux, % Polarization,
    position angle, and error for each file. The mean value of the Stokes Q and
    U parameters are subtracted from those values before calculating % Polarization
    and position angle to give a somewhat "normalized" version of these
    measurements. The median of each wavelength dependent value is then calculated
    for at each wavelength value across all the observations. 

    Parameters
    ----------
    txt_file_list : list or array
        A list where each indice is the string of the pathname to the txt file
        you wish to extract data from

    fits_file_list : list or array
        A list where each indice is the string of the pathname to the fits file
        you wish to extract data from

    bin_num : integer
        The number of bins you which to have. A higher number will mean smaller bins

    radial_velocity : float
        The radial velocity of the object you are observing. This is used to account
        for the doppler shift in wavelength dependent data.

    Returns
    ----------
    all_wavelength[0] : List
        An array of wavelengths across which the other values can be compared

    median_flux : List
        A list of median flux values across the values of wavelength

    median_pol : List
        A list of median % polarization values across the values of wavelength

    median_pos : List
        A list of median postion angle values across the values of wavelength

    median_err : List
        A list of median error values across the values of wavelength
    """
    all_wavelength, all_flux, all_Q, all_U, all_err = get_all_QU(
        txt_file_list, fits_file_list, bin_num, radial_velocity=radial_velocity)
    all_Q = sub_mean(all_Q)
    all_U = sub_mean(all_U)
    all_pol, all_pos = [], []
    for q, u in zip(all_Q, all_U):
        pol, pos = polarization(q, u)
        all_pol.append(pol)
        all_pos.append(pos)
    interp_flux = []
    interp_pol = []
    interp_pos = []
    interp_err = []
    for i in range(1, len(all_wavelength)):
        interpolated_flux = np.interp(
            all_wavelength[0], all_wavelength[i], all_flux[i])
        interp_flux.append(interpolated_flux)
        interpolated_pol = np.interp(
            all_wavelength[0], all_wavelength[i], all_pol[i])
        interp_pol.append(interpolated_pol)
        interpolated_pos = np.interp(
            all_wavelength[0], all_wavelength[i], all_pos[i])
        interp_pos.append(interpolated_pos)
        interpolated_err = np.interp(
            all_wavelength[0], all_wavelength[i], all_err[i])
        interp_err.append(interpolated_err)
    median_flux = np.median(np.array(interp_flux), axis=0)
    median_pol = np.median(np.array(interp_pol), axis=0)
    median_pos = np.median(np.array(interp_pos), axis=0)
    median_err = error_of_median(all_err)
    return(all_wavelength[0], median_flux, median_pol, median_pos, median_err)
