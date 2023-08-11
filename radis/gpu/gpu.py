import ctypes

import numpy as np
from scipy.constants import N_A, c, h, k

c_cm = 100 * c
c2 = h * c_cm / k


class initData(ctypes.Structure):
    _fields_ = [
        ("v_min", ctypes.c_float),
        ("v_max", ctypes.c_float),
        ("dv", ctypes.c_float),
        ("N_v", ctypes.c_int),
        ("dxG", ctypes.c_float),
        ("dxL", ctypes.c_float),
        ("N_total", ctypes.c_int),
        ("Max_lines", ctypes.c_int),
        ("N_lines", ctypes.c_int),
        ("N_points_per_block", ctypes.c_int),
        ("N_threads_per_block", ctypes.c_int),
        ("N_blocks_per_grid", ctypes.c_int),
        ("N_points_per_thread", ctypes.c_int),
        ("N_iterations_per_thread", ctypes.c_int),
        ("shared_size_floats", ctypes.c_int),
        ("log_c2Mm", ctypes.c_float * 16),
    ]


class iterData(ctypes.Structure):
    _fields_ = [
        ("p", ctypes.c_float),
        ("log_2p", ctypes.c_float),
        ("hlog_T", ctypes.c_float),
        ("log_rT", ctypes.c_float),
        ("c2T", ctypes.c_float),
        ("N", ctypes.c_float),
        ("l", ctypes.c_float),
        ("slit_FWHM", ctypes.c_float),
        ("log_wG_min", ctypes.c_float),
        ("log_wL_min", ctypes.c_float),
        ("N_G", ctypes.c_int),
        ("N_L", ctypes.c_int),
        ("Q", ctypes.c_float * 16),
    ]


init_h = initData()
iter_h = iterData()


def py_calc_lorentzian_envelope_params(na, gamma, verbose=False):
    """


    Parameters
    ----------
    na : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION..

    """
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

    return tuple(result)


def py_calc_gaussian_envelope_params(log_2vMm, verbose=False):
    """


    Parameters
    ----------
    log_2vMm : TYPE
        DESCRIPTION.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return np.amin(log_2vMm), np.max(log_2vMm)


try:
    from radis_cython_extensions import (  # isort:skip
        calc_gaussian_envelope_params,
        calc_lorentzian_envelope_params,
    )
except (ModuleNotFoundError):
    calc_gaussian_envelope_params = py_calc_gaussian_envelope_params
    calc_lorentzian_envelope_params = py_calc_lorentzian_envelope_params


def init_gaussian_params(log_2vMm, verbose):
    """


    Parameters
    ----------
    log_2vMm : TYPE
        DESCRIPTION.
    verbose : TYPE
        DESCRIPTION.

    Returns
    -------
    param_data : TYPE
        DESCRIPTION.

    """

    if verbose >= 2:
        print("Initializing Gaussian parameters")

    ##fname = "Gaussian_minmax_" + str(len(log_2vMm)) + ".dat"
    ##    try:
    ##        param_data = pickle.load(open(fname, "rb"))
    ##        if verbose >= 2:
    ##            print(" (from cache)... ")
    ##
    ##    except (OSError, IOError):
    if True:
        if verbose >= 2:
            print("... ")

        param_data = calc_gaussian_envelope_params(log_2vMm, verbose)
    ##        pickle.dump(param_data, open(fname, "wb"))

    if verbose >= 2:
        print("done!")

    return param_data


def init_lorentzian_params(na, gamma, verbose):
    """


    Parameters
    ----------
    na : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    verbose : TYPE
        DESCRIPTION.

    Returns
    -------
    param_data : TYPE
        DESCRIPTION.

    """

    if verbose >= 2:
        print("Initializing Lorentzian parameters ")

    ##fname = "Lorenzian_minmax_" + str(len(gamma)) + ".dat"
    ##
    ##    try:
    ##        with open(fname, "rb") as f:
    ##            if verbose >= 2:
    ##                print(" (from cache)... ")
    ##            param_data = pickle.load(f)
    ##
    ##    except:
    if True:
        if verbose >= 2:
            print(" ... ")

        param_data = calc_lorentzian_envelope_params(na, gamma, verbose)
    ##        with open(fname, "wb") as f:
    ##            pickle.dump(param_data, f)

    if verbose >= 2:
        print("done!")

    return param_data


def calc_gaussian_params(gaussian_param_data, init_h, iter_h, epsilon=1e-4):
    """


    Parameters
    ----------
    gaussian_param_data : TYPE
        DESCRIPTION.
    init_h : TYPE
        DESCRIPTION.
    iter_h : TYPE
        DESCRIPTION.
    epsilon : TYPE, optional
        DESCRIPTION. The default is 1e-4.

    Returns
    -------
    None.

    """

    log_2vMm_min, log_2vMm_max = gaussian_param_data
    log_wG_min = log_2vMm_min + iter_h.hlog_T
    log_wG_max = log_2vMm_max + iter_h.hlog_T
    ##    print("wG:", log_wG_min, log_wG_max)
    log_wG_max += epsilon

    # dxG = (log_wG_max - log_wG_min) / (init_h.N_G - 1)
    N = int(np.ceil((log_wG_max - log_wG_min) / init_h.dxG) + 1)

    iter_h.log_wG_min = log_wG_min
    # iter_h.dxG = dxG
    iter_h.N_G = N


def calc_lorentzian_minmax(param_data, log_rT, log_2p):
    """


    Parameters
    ----------
    param_data : TYPE
        DESCRIPTION.
    log_rT : TYPE
        DESCRIPTION.
    log_2p : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    result = []
    for params in param_data:
        A, B, X = params
        i = 0
        while X[i] < log_rT:
            i += 1
        result.append(log_rT * A[i] + B[i] + log_2p)
    return tuple(result)


def calc_lorentzian_params(param_data, init_h, iter_h, epsilon=1e-4):
    """


    Parameters
    ----------
    param_data : TYPE
        DESCRIPTION.
    init_h : TYPE
        DESCRIPTION.
    iter_h : TYPE
        DESCRIPTION.
    epsilon : TYPE, optional
        DESCRIPTION. The default is 1e-4.

    Returns
    -------
    None.

    """

    log_wL_min, log_wL_max = calc_lorentzian_minmax(
        param_data, iter_h.log_rT, iter_h.log_2p
    )
    log_wL_max += epsilon
    # dxL = (log_wL_max - log_wL_min) / (init_h.N_L - 1)
    N = int(np.ceil((log_wL_max - log_wL_min) / init_h.dxL) + 1)

    iter_h.log_wL_min = log_wL_min
    iter_h.N_L = N


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

    iter_h.p = p  # bar
    iter_h.log_2p = np.log(2 * p)
    iter_h.hlog_T = 0.5 * np.log(T)
    iter_h.log_rT = np.log(296.0 / T)
    iter_h.c2T = -c2 / T
    iter_h.N = mole_fraction * p * 1e5 / (1e6 * k * T)  # cm-3
    iter_h.l = l
    iter_h.slit_FWHM = slit_FWHM

    for i in range(len(Q_interpolator_list)):
        iter_h.Q[i] = Q_interpolator_list[i](T)


def gpu_init(
    v_arr,
    dxG,
    dxL,
    iso,
    v0,
    da,
    gamma,
    na,
    S0,
    El,
    Mm_arr,
    Q_intp_list,
    verbose=True,
    gpu=False,
):
    """


    Parameters
    ----------
    v_arr : TYPE
        DESCRIPTION.
    dxG : TYPE
        DESCRIPTION.
    dxL : TYPE
        DESCRIPTION.
    iso : TYPE
        DESCRIPTION.
    v0 : TYPE
        DESCRIPTION.
    da : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    na : TYPE
        DESCRIPTION.
    S0 : TYPE
        DESCRIPTION.
    El : TYPE
        DESCRIPTION.
    Mm_arr  : TYPE
        DESCRIPTION.
    Q_intp_list: TYPE
        DESCRIPTION.
    verbose : TYPE, optional
        DESCRIPTION. The default is True.
    gpu : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """

    # ----------- setup global variables -----------------
    global init_h
    global block_preparation_step_size
    global iso_d, v0_d, da_d, S0_d, El_d, gamma_d, na_d
    global spectrum_in_d, transmittance_noslit_d, transmittance_FT_d
    global lorentzian_param_data, gaussian_param_data, Q_interpolator_list
    global ctx, cuda_module
    # -----------------------------------------------------

    from radis.gpu.cuda_driver import cuArray, cuContext, cuModule

    ctx = cuContext()
    cuda_module = cuModule(ctx, "kernels.ptx")

    init_h.v_min = np.min(v_arr)  # 2000.0
    init_h.v_max = np.max(v_arr)  # 2400.0
    init_h.dv = (v_arr[-1] - v_arr[0]) / (len(v_arr) - 1)  # 0.002
    init_h.N_v = len(v_arr)
    init_h.dxG = dxG
    init_h.dxL = dxL

    init_h.N_iterations_per_thread = 1024
    block_preparation_step_size = 128

    shared_size = 0x8000  # Bytes - Size of the shared memory
    init_h.shared_size_floats = shared_size // 4  # size of float

    # init_h.N_total = init_h.N_v * init_h.N_G * init_h.N_L
    # init_h.N_points_per_block = init_h.shared_size_floats // (init_h.N_v * init_h.N_G)

    init_h.N_threads_per_block = 1024
    init_h.N_blocks_per_grid = 4 * 256 * 256
    init_h.N_points_per_thread = init_h.N_points_per_block // init_h.N_threads_per_block

    if verbose >= 2:
        print()
        print("Spectral points per block  : {0}".format(init_h.N_points_per_block))
        print("Threads per block          : {0}".format(init_h.N_threads_per_block))
        print("Spectral points per thread : {0}".format(init_h.N_points_per_thread))
        print()

    init_h.N_lines = int(len(v0))

    if verbose == 1:
        print("Number of lines loaded: {0}".format(init_h.N_lines))
        print()

    if verbose >= 2:
        print("Allocating device memory and copying data...")

    # S_klm_d = zeros(
    #     (2 * init_h.N_v, init_h.N_G, init_h.N_L),
    #     order="C",
    #     dtype=float32,
    # )
    spectrum_in_d = cuArray(np.zeros(init_h.N_v + 1, dtype=np.complex64))
    transmittance_noslit_d = cuArray(np.zeros(init_h.N_v * 2, dtype=np.float32))
    transmittance_FT_d = cuArray(np.zeros(init_h.N_v + 1, dtype=np.complex64))

    # #Copy spectral data to device
    iso_d = cuArray(iso)
    v0_d = cuArray(v0)
    da_d = cuArray(da)
    S0_d = cuArray(S0)
    El_d = cuArray(El)
    gamma_d = cuArray(gamma)
    na_d = cuArray(na)

    if verbose >= 2:
        print("done!")

    log_c2Mm_arr = np.array(
        [0]
        + [
            0.5 * np.log(8 * k * np.log(2) / (c**2 * Mm * 1e-3 / N_A))
            for Mm in Mm_arr[1:]
        ]
    )
    for i in range(len(log_c2Mm_arr)):
        init_h.log_c2Mm[i] = log_c2Mm_arr[i]
    Q_interpolator_list = Q_intp_list

    log_2vMm = np.log(v0) + log_c2Mm_arr.take(iso)
    gaussian_param_data = init_gaussian_params(log_2vMm.astype(np.float32), verbose)
    lorentzian_param_data = init_lorentzian_params(na, gamma, verbose)

    if verbose >= 2:
        print("Copying initialization parameters to device memory...")
    cuda_module.setConstant("init_d", init_h)

    if verbose >= 2:
        print("done!")
        print("Copying spectral data to device memory...")


##    print("wG fast:",np.min(log_2vMm),np.max(log_2vMm))
##
##    log_wG_debug_h = np.log(v0) + 0.5*np.log(8*k*np.log(2)/(c**2*Mm_arr[iso]))
##    print("wG arr:",np.min(log_wG_debug_h),np.max(log_wG_debug_h))


def gpu_iterate(p, T, mole_fraction, l=1.0, slit_FWHM=0.0, verbose=0, gpu=False):
    """
    Parameters
    ----------
    p : float
        pressure [bar]
    T : float
        temperature [K]
    mole_fraction : float

    Other Parameters
    ----------------
    verbose : bool, optional
        The default is True.
    l : TYPE, optional
        DESCRIPTION. The default is 1.0.
    slit_FWHM : TYPE, optional
        DESCRIPTION. The default is 0.0.
    gpu : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    abscoeff_h : TYPE
        DESCRIPTION.
    transmittance_h : TYPE
        DESCRIPTION.

    """

    # ----------- setup global variables -----------------
    global init_h, iter_h
    global iso_d, v0_d, da_d, S0_d, El_d, gamma_d, na_d
    global spectrum_in_d, transmittance_noslit_d, transmittance_FT_d
    global ctx, cuda_module
    global start, stop, elapsedTime
    # ------------------------------------------------------

    from radis.gpu.cuda_driver import cuArray, cuFFT

    if verbose >= 2:
        print("Copying iteration parameters to device...")

    set_pTQ(p, T, mole_fraction, iter_h, l=l, slit_FWHM=slit_FWHM)

    calc_gaussian_params(
        gaussian_param_data,
        init_h,
        iter_h,
    )

    calc_lorentzian_params(
        lorentzian_param_data,
        init_h,
        iter_h,
    )

    cuda_module.setConstant("iter_d", iter_h)

    if verbose >= 2:
        print("done!")
        print("Filling LDM...")

    S_klm_d = cuArray(
        np.zeros((2 * init_h.N_v, iter_h.N_G, iter_h.N_L), order="C", dtype=np.float32)
    )

    spectrum_in_d[:] = 0.0
    transmittance_FT_d[:] = 0.0

    n_threads = init_h.N_threads_per_block
    n_blocks = init_h.N_lines // (n_threads * init_h.N_iterations_per_thread) + 1

    cuda_module.fillLDM(
        iso_d,
        v0_d,
        da_d,
        S0_d,
        El_d,
        gamma_d,
        na_d,
        S_klm_d,
        blocks=(n_blocks, 1, 1),
        threads=(n_threads, 1, 1),
    )

    ctx.synchronize()

    if verbose >= 2:
        print("Applying lineshapes...")

    S_klm_FT_d = cuArray(
        np.zeros((init_h.N_v, iter_h.N_G, iter_h.N_L), order="C", dtype=np.float32)
    )

    fft_fwd = cuFFT(S_klm_d, S_klm_FT_d, direction="fwd")
    fft_fwd.execute()

    ctx.synchronize()

    n_threads = init_h.N_threads_per_block
    n_blocks = (init_h.N_v + 1) // n_threads + 1
    cuda_module.applyLineshapes(
        S_klm_FT_d, spectrum_in_d, blocks=(n_blocks, 1, 1), threads=(n_threads, 1, 1)
    )

    ctx.synchronize()

    spectrum_out_d = cuArray(np.zeros(2 * init_h.N_v, dtype=np.float32))

    fft_bwd = cuFFT(spectrum_in_d, spectrum_out_d, direction="bwd")
    fft_bwd.execute()

    ctx.synchronize()

    if verbose >= 2:
        print("Done!")
        print("Calculating transmittance...")

    abscoeff_h = spectrum_out_d.d2h()[: init_h.N_v]

    ## Calc transmittance_noslit:
    n_threads = init_h.N_threads_per_block
    n_blocks = (2 * init_h.N_v) // n_threads + 1
    cuda_module.calcTransmittanceNoslit(
        spectrum_out_d,
        transmittance_noslit_d,
        blocks=(n_blocks, 1, 1),
        threads=(n_threads, 1, 1),
    )

    ctx.synchronize()

    if verbose >= 2:
        print("Done!")
        print("Applying slit function...")

    ## Apply slit function:
    transmittance_noslit_FT_d = cuArray(np.zeros(init_h.N_v, dtype=np.complex64))
    fft_fwd2 = cuFFT(transmittance_noslit_d, transmittance_noslit_FT_d, direction="fwd")
    fft_fwd2.execute()

    ctx.synchronize()

    n_threads = init_h.N_threads_per_block
    n_blocks = (init_h.N_v + 1) // n_threads + 1
    cuda_module.applyGaussianSlit(
        transmittance_noslit_FT_d,
        transmittance_FT_d,
        blocks=(n_blocks, 1, 1),
        threads=(n_threads, 1, 1),
    )

    ctx.synchronize()

    transmittance_d = cuArray(np.zeros(2 * init_h.N_v, dtype=np.float32))
    fft_bwd2 = cuFFT(transmittance_FT_d, transmittance_d, direction="bwd")
    fft_bwd2.execute()

    ctx.synchronize()

    transmittance_h = transmittance_d.d2h()[: init_h.N_v]

    if verbose >= 2:
        print("done!")

    if verbose == 1:
        print("Finished calculating spectrum!")

    ctx.destroy()
    return abscoeff_h, transmittance_h, iter_h
