import os.path

import numpy as np
from scipy.constants import N_A, c, h, k
from scipy.fft import next_fast_len

from radis.gpu.structs import initData, iterData
from radis.misc.utils import getProjectRoot

# import sys
# import matplotlib.pyplot as plt


c_cm = 100 * c
c2 = h * c_cm / k

_cuda_context_open = False


init_h = initData()
iter_h = iterData()

# TODO: the first two steps could be done on GPU if it makes it faster.
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
    global init_h, ctx, cu_mod, _cuda_context_open
    global lorentzian_param_data, gaussian_param_data, Q_interpolator_list
    # -----------------------------------------------------

    from radis.gpu.driver import CuArray, CuContext, CuFFT, CuModule

    if _cuda_context_open:
        # TODO: warn
        print("Only a single CUDA context allowed; please call gpu_exit() first.")
        return

    if verbose == 1:
        print("Number of lines loaded: {0}".format(len(v0)))
        print()

    ## First a CUDA context is created, then the .ptx file is read
    ## and made available as the CuModule object cu_mod

    ctx = CuContext()
    _cuda_context_open = True
    ptx_path = os.path.join(getProjectRoot(), "gpu", "kernels.ptx")
    if not os.path.exists(ptx_path):
        raise FileNotFoundError(ptx_path)
    cu_mod = CuModule(ctx, ptx_path)

    ## Next, the GPU is made aware of a number of parameters.
    ## Parameters that don't change during iteration are stored
    ## in iter_h. They are copied to the GPU through cu_mod.setConstant()

    if verbose >= 2:
        print("Copying initialization parameters to device memory...")

    init_h.v_min = np.min(v_arr)
    init_h.v_max = np.max(v_arr)
    init_h.dv = (v_arr[-1] - v_arr[0]) / (len(v_arr) - 1)  # TODO: get this from caller
    init_h.N_v = len(v_arr)
    init_h.N_v_FT = next_fast_len(2 * init_h.N_v)
    init_h.N_x_FT = init_h.N_v_FT // 2 + 1
    init_h.dxG = dxG
    init_h.dxL = dxL
    init_h.N_lines = int(len(v0))
    init_h.N_iterations_per_thread = 1024

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

    cu_mod.setConstant("init_d", init_h)

    gaussian_param_data = init_gaussian_params(log_2vMm.astype(np.float32), verbose)
    lorentzian_param_data = init_lorentzian_params(na, gamma, verbose)

    if verbose >= 2:
        print("done!")

    ## Next the block- and thread size of the CUDA kernels are set.
    ## This determines how the GPU internally divides up the work.

    if verbose >= 2:
        print("Allocating device memory and copying data...")

    NvFT = init_h.N_v_FT
    NxFT = NvFT // 2 + 1
    Ntpb = 1024  # threads per block #TODO: determine dynamically
    Nipt = init_h.N_iterations_per_thread
    Nli = init_h.N_lines
    threads = (Ntpb, 1, 1)

    cu_mod.fillLDM.setGrid((Nli // (Ntpb * Nipt) + 1, 1, 1), threads)
    cu_mod.applyLineshapes.setGrid((NxFT // Ntpb + 1, 1, 1), threads)
    cu_mod.calcTransmittanceNoslit.setGrid((NvFT // Ntpb + 1, 1, 1), threads)
    cu_mod.applyGaussianSlit.setGrid((NxFT // Ntpb + 1, 1, 1), threads)

    ## Next the variables are initialized on the GPU. Constant variables
    ## that don't change (i.e. pertaining to the database) are immediately
    ## copied to the GPU through CuArray.fromArray().
    ## Other variables are only allocated. S_klm_d and S_klm_FT_d are
    ## special cases because their shape changes during iteration.
    ## They are not allocated, only given a device pointer by which
    ## they can be referenced later.

    S_klm_d = CuArray(0, dtype=np.float32, init="defer")  # dev_ptr only
    S_klm_FT_d = CuArray(0, dtype=np.complex64, init="defer")  # dev_ptr only

    spectrum_in_d = CuArray(NxFT, dtype=np.complex64)  # malloc
    spectrum_out_d = CuArray(NvFT, dtype=np.float32)  # malloc

    transmittance_noslit_d = CuArray(NvFT, dtype=np.float32)  # malloc
    transmittance_noslit_FT_d = CuArray(NxFT, dtype=np.complex64)  # malloc

    transmittance_FT_d = CuArray(NxFT, dtype=np.complex64)  # malloc
    transmittance_d = CuArray(NvFT, dtype=np.float32)  # malloc

    cu_mod.fillLDM.setArgs(
        CuArray.fromArray(iso),
        CuArray.fromArray(v0),
        CuArray.fromArray(da),
        CuArray.fromArray(S0),
        CuArray.fromArray(El),
        CuArray.fromArray(gamma),
        CuArray.fromArray(na),
        S_klm_d,
    )
    cu_mod.applyLineshapes.setArgs(S_klm_FT_d, spectrum_in_d)
    cu_mod.calcTransmittanceNoslit.setArgs(spectrum_out_d, transmittance_noslit_d)
    cu_mod.applyGaussianSlit.setArgs(transmittance_noslit_FT_d, transmittance_FT_d)

    ## FFT's are performed through the CuFFT object. The required functions are internally
    ## loaded from the cufft library, not through the user kernels (.ptx files).

    cu_mod.fft_fwd = CuFFT(S_klm_d, S_klm_FT_d, direction="fwd", plan_fft=False)
    cu_mod.fft_rev = CuFFT(spectrum_in_d, spectrum_out_d, direction="rev")
    cu_mod.fft_fwd2 = CuFFT(
        transmittance_noslit_d, transmittance_noslit_FT_d, direction="fwd"
    )
    cu_mod.fft_rev2 = CuFFT(transmittance_FT_d, transmittance_d, direction="rev")

    if verbose >= 2:
        print("done!")


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
    global init_h, iter_h, cu_mod, _cuda_context_open
    # ------------------------------------------------------

    if not _cuda_context_open:
        # TODO: warn
        print("Must have an open CUDA context; please call gpu_init() first.")
        return

    if verbose >= 2:
        print("Copying iteration parameters to device...")

    ## First a number of parameters that change during iteration
    ## are computed and copied to the GPU.

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

    cu_mod.setConstant("iter_d", iter_h)

    ## Next the S_klm_d variable is reshaped to the correct shape,
    ## and filled with spectral data.

    if verbose >= 2:
        print("done!")
        print("Filling LDM...")

    S_klm_shape = (init_h.N_v_FT, iter_h.N_G, iter_h.N_L)
    cu_mod.fillLDM.args[-1].resize(S_klm_shape, init="zeros")
    cu_mod.fillLDM()

    ## Next the S_klm_FT_d is also reshaped, and the lineshapes are
    ## applied. This consists of an FT of the LDM, a product by the
    ## lineshape FTs & summing all G and L axes, and an inverse FT
    ## on the accumulated spectra.

    if verbose >= 2:
        print("done!")
        print("Applying lineshapes...")

    S_klm_FT_shape = (init_h.N_x_FT, iter_h.N_G, iter_h.N_L)
    cu_mod.fft_fwd.arr_out.resize(S_klm_FT_shape)
    cu_mod.fft_fwd.planMany()  # replan because size may have changed
    cu_mod.fft_fwd()
    cu_mod.applyLineshapes()
    cu_mod.fft_rev()

    if verbose >= 2:
        print("Done!")
        print("Calculating transmittance...")

    # ctx.synchronize()
    abscoeff_h = cu_mod.fft_rev.arr_out.getArray()[: init_h.N_v]

    ## To apply a slit function, first the transmittance is calculated.
    ## Then the convolution is applied by an FT, product with the
    ## instrument function's FT, followed by an inverse FT.

    if verbose >= 2:
        print("Done!")
        print("Applying slit function...")

    cu_mod.calcTransmittanceNoslit()
    cu_mod.fft_fwd2()
    cu_mod.applyGaussianSlit()
    cu_mod.fft_rev2()

    # ctx.synchronize()
    transmittance_h = cu_mod.fft_rev2.arr_out.getArray()[: init_h.N_v]

    if verbose >= 2:
        print("done!")

    if verbose == 1:
        print("Finished calculating spectrum!")

    return abscoeff_h, transmittance_h, iter_h


def gpu_exit():
    global ctx, _cuda_context_open
    ctx.destroy()
    _cuda_context_open = False
