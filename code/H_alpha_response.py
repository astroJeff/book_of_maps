import sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


H_alpha_response = None



def get_power_law(x, y):
    """ Define power law function

    Inputs
    ------
    x : float
        Linear independent variable
    y : float
        Logarithmic dependent variable

    Returns
    -------
    m, b : floats
        slope and intercept of power law fit
     """

    log_x = np.log10(x)
    log_y = y

    N_points = len(x)

    m = (N_points * np.sum(log_x * log_y) - np.sum(log_x)*np.sum(log_y)) / (N_points * np.sum(log_x**2) - np.sum(log_x)**2)
    b = (np.sum(log_y) - m*np.sum(log_x))/N_points

    return m, b



def expand_arr(array):
    """
    Function to perform power law extrapolation

    Inputs
    ------
    array : np array
        Numpy array containing columns 'age' and 'log_N'

    Returns
    -------
    array : np array
        Numpy array of same dtype as input, but with values
        added from extrapolation to 1 Gyr

    """

    idx = np.where(array['age'] > 5.0)[0]
    m, b = get_power_law(array['age'][idx], array['log_N'][idx])

    # print(m, b)

    add_x = 1.0e3
    add_y = np.log10(add_x)*m + b

    # print(add_x, add_y)

    tmp_arr = np.ones(1, dtype=array.dtype)
    tmp_arr[0]['age'] = add_x
    tmp_arr[0]['log_N'] = add_y

    array = np.append(array, tmp_arr)

    return array



def load_response_function(Z=0.004):

    global H_alpha_response

    Z_allowed = [0.004, 0.008]

    if not Z in Z_allowed:
        print("You must provide a valid metallicity.")
        print("Allowed options:", Z_allowed)
        sys.exit()


    # Preparing to load up data
    dtype = [('age', 'f8'), ('log_N','f8')]

    if Z==0.004:
        filename = "../data/response_NLyc_z0.004.txt"
    elif Z==0.008:
        filename = "../data/response_NLyc_z0.008.txt"


    # Load up data from paper
    N_Lyc =  np.genfromtxt(filename, dtype=dtype, skip_header=9)

    # Expanding loaded up data
    N_Lyc = expand_arr(N_Lyc)


    # Transformation from Section 3.1 from Oti-Floranes & Mas Hesse (2010)
    # L(Hα) = 1.36x10^-12 N_Lyc erg s^-1
    H_alpha_response = interp1d(np.log10(N_Lyc['age']), np.log10(1.36e-12) + N_Lyc['log_N'],
                                bounds_error=False, fill_value=0.0)



def H_alpha_response_t(t, Z=0.004):
    """  Return the H_alpha response some time, t, after a star formation episode

    Inputs
    ------
    t : float
        time after star formation episode (Myr)
    Z : float
        Metallicity in solar units. Check load_response_function for options

    Returns
    -------
    L_Hα : float
        H_alpha luminosity at time t (erg/s)

    """

    global H_alpha_response

    if H_alpha_response is None: load_response_function(Z=Z)

    return 10**H_alpha_response(t)



def test_response(Z=0.004):
    """  Test H_alpha response function """


    global H_alpha_response


    load_response_function(Z=Z)



    tmp_x = np.linspace(1.0, 1.0e3-1, 10000)

    tmp_y = 10**H_alpha_response(np.log10(tmp_x))
    plt.plot(tmp_x, tmp_y)

    plt.xlim(0.1, 20)

    plt.yscale('log')
    plt.ylim(1.0e30, 1.0e36)
    plt.show()
