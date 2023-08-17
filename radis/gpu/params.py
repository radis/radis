import numpy as np
from scipy.constants import c, h, k

c_cm = 100 * c
c2 = h * c_cm / k

# TODO: the first two steps could be done on GPU if it makes it faster.
def init_L_params(na, gamma, verbose=False):

    if verbose >= 2:
        print("Initializing Lorentzian parameters ")

    try:
        from radis_cython_extensions import cy_init_L_params

        param_data = cy_init_L_params(na, gamma)

    except (ModuleNotFoundError):

        # Remove duplicates
        unique_lines = set([])
        for i in range(len(na)):
            unique_lines.add(str(na[i]) + " " + str(gamma[i]))

        # Only keep extremes
        max_dict = {}
        min_dict = {}
        for s in unique_lines:
            na_i, gamma_i = map(float, s.split())
            try:
                min_dict[na_i] = gamma_i if gamma_i < min_dict[na_i] else min_dict[na_i]
                max_dict[na_i] = gamma_i if gamma_i > max_dict[na_i] else max_dict[na_i]

            except (KeyError):
                min_dict[na_i] = gamma_i
                max_dict[na_i] = gamma_i

        # Check which ones are really at the top:
        result = []
        for test_dict in (min_dict, max_dict):

            keys = sorted(test_dict.keys(), reverse=(test_dict == min_dict))
            A = [keys[0]]
            B = [np.log(test_dict[keys[0]])]
            X = [-np.inf]

            for key in keys[1:]:
                for i in range(len(X)):
                    xi = (np.log(test_dict[key]) - B[i]) / (A[i] - key)
                    if xi >= X[i]:
                        if i < len(X) - 1:
                            if xi < X[i + 1]:
                                break
                        else:
                            break

                A = A[: i + 1] + [key]
                B = B[: i + 1] + [np.log(test_dict[key])]
                X = X[: i + 1] + [xi]

            X = X[1:] + [np.inf]
            result.append((A, B, X))

        param_data = tuple(result)

    if verbose >= 2:
        print("done!")

    global _L_param_data
    _L_param_data = param_data


def init_G_params(log_2vMm, verbose=False):

    if verbose >= 2:
        print("Initializing Gaussian parameters")

    try:
        from radis_cython_extensions import cy_init_G_params

        param_data = cy_init_G_params(log_2vMm)

    except (ModuleNotFoundError):

        param_data = (np.min(log_2vMm), np.max(log_2vMm))

    if verbose >= 2:
        print("done!")

    global _G_param_data
    _G_param_data = param_data


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
    iter_h.c2T = -c2 / T
    iter_h.N = mole_fraction * p * 1e5 / (1e6 * k * T)  # cm-3
    iter_h.l = l
    iter_h.slit_FWHM = slit_FWHM

    for i in range(len(_Q_intp_list)):
        iter_h.Q[i] = _Q_intp_list[i](T)
