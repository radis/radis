import os.path
from warnings import warn

import numpy as np
from scipy.constants import N_A, c, h, k
from scipy.fft import next_fast_len

from radis.gpu.params import (
    init_G_params,
    init_L_params,
    init_Q,
    set_G_params,
    set_L_params,
    set_pTQ,
)
from radis.gpu.structs import initData_t, iterData_t
from radis.misc.utils import getProjectRoot

c_cm = 100 * c
c2 = h * c_cm / k

cu_mod = None
init_h = initData_t()
iter_h = iterData_t()


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
    emulate=True,
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

    global cu_mod
    if cu_mod is not None:
        warn("Only a single CUDA context allowed; please call gpu_exit() first.")
        return

    if emulate:
        from radis.gpu.emulate import CuArray, CuContext, CuFFT, CuModule
    else:
        from radis.gpu.driver import CuArray, CuContext, CuFFT, CuModule

    ctx = CuContext()

    ##    if ctx is None:
    ##        warn(("Failed to load CUDA context, this happened either because"+
    ##              "CUDA is not installed properly, or you have no NVIDIA GPU."+
    ##              "Continuing with emulated GPU on CPU..."+
    ##              "This means *NO* GPU acceleration!"))
    ##
    ##        from radis.gpu.emulate import CuArray, CuContext, CuFFT, CuModule
    ##        ctx = CuContext()

    if verbose == 1:
        print("Number of lines loaded: {0}".format(len(v0)))
        print()

    ## First a CUDA context is created, then the .ptx file is read
    ## and made available as the CuModule object cu_mod

    ptx_path = os.path.join(getProjectRoot(), "gpu", "kernels.ptx")
    if not os.path.exists(ptx_path):
        raise FileNotFoundError(ptx_path)
    cu_mod = CuModule(ctx, ptx_path)

    ## Next, the GPU is made aware of a number of parameters.
    ## Parameters that don't change during iteration are stored
    ## in init_h. They are copied to the GPU through cu_mod.setConstant()

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

    init_Q(Q_intp_list)
    log_2vMm = np.log(v0) + log_c2Mm_arr.take(iso)

    cu_mod.setConstant("init_d", init_h)

    init_G_params(log_2vMm.astype(np.float32), verbose)
    init_L_params(na, gamma, verbose)

    if verbose >= 2:
        print("done!")

    ## Next the block- and thread size of the CUDA kernels are set.
    ## This determines how the GPU internally divides up the work.

    if verbose >= 2:
        print("Allocating device memory and copying data...")

    NvFT = init_h.N_v_FT
    NxFT = NvFT // 2 + 1
    Ntpb = ctx.getMaxThreadsPerBlock()
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

    if cu_mod is None:
        warn("Must have an open CUDA context; please call gpu_init() first.")
        return

    if verbose >= 2:
        print("Copying iteration parameters to device...")

    ## First a number of parameters that change during iteration
    ## are computed and copied to the GPU.
    set_pTQ(p, T, mole_fraction, iter_h, l=l, slit_FWHM=slit_FWHM)
    set_G_params(init_h, iter_h)
    set_L_params(init_h, iter_h)
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

    transmittance_h = cu_mod.fft_rev2.arr_out.getArray()[: init_h.N_v]

    if verbose >= 2:
        print("done!")

    if verbose == 1:
        print("Finished calculating spectrum!")

    return abscoeff_h, transmittance_h, iter_h


def gpu_exit(event=None):
    global cu_mod
    cu_mod.context.destroy()
    cu_mod = None
