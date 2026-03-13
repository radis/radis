import numpy as np
from scipy.constants import c, h, k

from radis.lbl.envelope import compute_gaussian_envelope, compute_lorentzian_envelope

c_cm = 100 * c
c2 = h * c_cm / k


def init_L_params(na, gamma_arr, verbose=False):

    if verbose >= 2:
        print("Initializing Lorentzian parameters ")

    global _L_param_data
    _L_param_data = compute_lorentzian_envelope(na, gamma_arr)

    if verbose >= 2:
        print("done!")


def init_G_params(log_2vMm, verbose=False):

    if verbose >= 2:
        print("Initializing Gaussian parameters")

    global _G_param_data
    _G_param_data = compute_gaussian_envelope(log_2vMm)

    if verbose >= 2:
        print("done!")


def init_Q(Q_intp_list):
    global _Q_intp_list
    _Q_intp_list = Q_intp_list


def set_L_params(init_h, iter_h, epsilon=1e-4):

    global _L_param_data
    result = []
    for params in _L_param_data:
        A, B, X = params
        i = 0
        while X[i] < iter_h.log_rT:
            i += 1
        result.append(iter_h.log_rT * A[i] + B[i] + iter_h.log_2p)
    log_wL_min, log_wL_max = result

    log_wL_max += epsilon
    N = int(np.ceil((log_wL_max - log_wL_min) / init_h.dxL) + 1)

    iter_h.log_wL_min = log_wL_min
    iter_h.N_L = N


def set_G_params(init_h, iter_h, epsilon=1e-4):
    global _G_param_data

    log_2vMm_min, log_2vMm_max = _G_param_data
    log_wG_min = log_2vMm_min + iter_h.hlog_T
    log_wG_max = log_2vMm_max + iter_h.hlog_T
    log_wG_max += epsilon

    N = int(np.ceil((log_wG_max - log_wG_min) / init_h.dxG) + 1)

    iter_h.log_wG_min = log_wG_min
    iter_h.N_G = N


def set_pTQ(p, T, mole_fraction, iter_h, l=1.0, slit_FWHM=0.0):
    """


    Parameters
    ----------
    p : float
        pressure [bar].
    T : float
        temperature [K].
    mole_fraction : float
    iter_h : TYPE
        DESCRIPTION.
    l : TYPE, optional
        DESCRIPTION. The default is 1.0.
    slit_FWHM : TYPE, optional
        DESCRIPTION. The default is 0.0.

    Returns
    -------
    None.

    """
    global _Q_intp_list
    iter_h.p = p  # bar
    iter_h.log_2p = np.log(2 * p)
    iter_h.hlog_T = 0.5 * np.log(T)
    iter_h.log_rT = np.log(296.0 / T)
    iter_h.c2T = c2 / T
    iter_h.N = p * 1e5 / (1e6 * k * T)  # cm-3
    iter_h.x[0] = mole_fraction  # self
    iter_h.x[1] = 1 - mole_fraction  # air
    # iter_h.l = l //TODO: GPU calculation of absorbance not currently implemented
    # iter_h.slit_FWHM = slit_FWHM //TODO: GPU application of ILS not currently implemented

    for i in range(len(_Q_intp_list)):
        iter_h.Q[i] = _Q_intp_list[i](T)
