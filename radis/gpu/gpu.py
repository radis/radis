import os.path
from warnings import warn

import numpy as np
from scipy.constants import N_A, c, k
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
from radis.misc.warning import NoGPUWarning

cu_mod = None
init_h = initData_t()
iter_h = iterData_t()


def gpu_init(
    vmin,
    Nv,
    dv,
    dxG,
    dxL,
    v0,
    da,
    na,
    S0,
    El,
    gamma,
    iso,
    Mm_arr,
    Q_intp_list,
    verbose=0,
    backend="gpu-cuda",
):
    """
    Initialize GPU-based calculation for emission and absorption spectra in spectroscopy.

    Parameters
    ----------
    vmin : float
        Minimum value frequency/wavenumber axis
    Nv : int
        Total number of frequency/wavenumber points.
    dv : float
        Stepsize of frequency/wavenumber axis (called wstep elsewhere in RADIS).
    dxG : float
        Relative grid spacing for Gaussian lineshapes.
    dxL : float
        Relative grid spacing for Lorentzian lineshapes.
    v0 : numpy.ndarray[np.float32]
        Array of line center frequencies (in cm-1).
    da : numpy.ndarray[np.float32]
        Pressure shift  (in cm-1.atm-1).
    na : numpy.ndarray[np.float32]
        Temperature dependency of Lorentzian widths
    S0 : numpy.ndarray[np.float32]
        Line intensity scaling factors.
    El : numpy.ndarray[np.float32]
        Lower level energy levels.
    gamma : numpy.ndarray[np.float32]
        Lorentzian width parameters.
    iso : numpy.ndarray[np.uint8]
        Index of isotopologue.
    Mm_arr : numpy.ndarray
        Molecular masses for all isotopologues in the database (Mm_arr[0] is always 0).
    Q_intp_list : list
        List of Q branch interpolators.
    verbose : bool, optional
        Print verbosity level. Default is 0.
    backend : str, optional
        Which backend to use currently only 'gpu-cuda' and 'cpu-cuda' available. Default is 'gpu-cuda'.

    Returns
    -------
    None.
    """

    global cu_mod
    if cu_mod is not None:
        warn("Only a single CUDA context allowed; please call gpu_exit() first.")
        return

    ## First a CUDA context is created, then the .ptx file is read
    ## and made available as the CuModule object cu_mod
    ## If this fails, None is returned and calculations are
    ## defaulted to CPU emulation

    if backend == "cpu-cuda":
        from radis.gpu.cuda.emulate import CuArray, CuContext, CuFFT, CuModule, CuTimer

        ctx = CuContext.Open(verbose=verbose)

    else:
        # Try to load GPU
        from radis.gpu.cuda.driver import CuContext

        ctx = CuContext.Open(
            verbose=verbose
        )  # Set verbosity to 2 or higher for printing driver output
        if ctx is None:
            warn(
                NoGPUWarning(
                    "Failed to load CUDA context, this happened either because"
                    + "CUDA is not installed properly, or you have no NVIDIA GPU. "
                    + "Continuing with emulated GPU on CPU..."
                    + "This means *NO* GPU acceleration!"
                )
            )

            # failed to init CUDA context, continue with CPU:
            from radis.gpu.cuda.emulate import (
                CuArray,
                CuContext,
                CuFFT,
                CuModule,
                CuTimer,
            )

            ctx = CuContext.Open(verbose=verbose)

        else:
            # successfully initialized CUDA context, continue with GPU:
            from radis.gpu.cuda.driver import CuArray, CuFFT, CuModule, CuTimer

    if verbose == 1:
        print("Number of lines loaded: {0}".format(len(v0)))
        print()

    ptx_path = os.path.join(getProjectRoot(), "gpu", "cuda", "build", "kernels.ptx")
    if not os.path.exists(ptx_path):
        raise FileNotFoundError(ptx_path)
    cu_mod = CuModule(ctx, ptx_path)  # gpu
    print("mode:", cu_mod.getMode())

    ## Next, the GPU is made aware of a number of parameters.
    ## Parameters that don't change during iteration are stored
    ## in init_h. They are copied to the GPU through cu_mod.setConstant()

    if verbose >= 2:
        print("Copying initialization parameters to device memory...")

    init_h.v_min = vmin
    init_h.dv = dv
    init_h.N_v = Nv
    init_h.N_v_FT = next_fast_len(2 * init_h.N_v)
    init_h.N_x_FT = init_h.N_v_FT // 2 + 1
    init_h.dxG = dxG
    init_h.dxL = dxL
    init_h.N_lines = int(len(v0))

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
    Nli = init_h.N_lines
    threads = (Ntpb, 1, 1)

    cu_mod.fillLDM.setGrid((Nli // Ntpb + 1, 1, 1), threads)
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

    S_klm_d = CuArray(0, dtype=np.float32, grow_only=True)
    S_klm_FT_d = CuArray(0, dtype=np.complex64, grow_only=True)

    spectrum_in_d = CuArray(NxFT, dtype=np.complex64)
    spectrum_out_d = CuArray(NvFT, dtype=np.float32)

    transmittance_noslit_d = CuArray(NvFT, dtype=np.float32)
    transmittance_noslit_FT_d = CuArray(NxFT, dtype=np.complex64)

    transmittance_FT_d = CuArray(NxFT, dtype=np.complex64)
    transmittance_d = CuArray(NvFT, dtype=np.float32)

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
    ## The FFT's need some memory as "work area". Because the different FFT's can
    ## reuse the work area, we make a CuArray at this scope that is passed to the
    ## CuFFT objects. The work area will be scaled according to needs by the CuFFT objects,
    ## so it can be initialized with a small value.

    workarea_d = CuArray(0, dtype=np.byte, grow_only=True)
    cu_mod.fft_fwd = CuFFT(S_klm_d, S_klm_FT_d, workarea=workarea_d, direction="fwd")
    cu_mod.fft_rev = CuFFT(
        spectrum_in_d, spectrum_out_d, workarea=workarea_d, direction="rev"
    )
    cu_mod.fft_fwd2 = CuFFT(
        transmittance_noslit_d,
        transmittance_noslit_FT_d,
        workarea=workarea_d,
        direction="fwd",
    )
    cu_mod.fft_rev2 = CuFFT(
        transmittance_FT_d, transmittance_d, workarea=workarea_d, direction="rev"
    )

    cu_mod.timer = CuTimer()

    if verbose >= 2:
        print("done!")

    return init_h


def gpu_iterate(p, T, mole_fraction, l=1.0, slit_FWHM=0.0, verbose=0):
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

    cu_mod.timer.reset()

    set_pTQ(p, T, mole_fraction, iter_h, l=l, slit_FWHM=slit_FWHM)
    set_G_params(init_h, iter_h)
    set_L_params(init_h, iter_h)
    cu_mod.setConstant("iter_d", iter_h)
    cu_mod.timer.lap("iter_params")

    ## Next the S_klm_d variable is reshaped to the correct shape,
    ## and filled with spectral data.

    if verbose >= 2:
        print("done!")
        print("Filling LDM...")

    S_klm_shape = (init_h.N_v_FT, iter_h.N_G, iter_h.N_L)

    cu_mod.fillLDM.args[-1].resize(S_klm_shape, init="zeros")
    cu_mod.fillLDM()
    cu_mod.timer.lap("fillLDM")

    ## Next the S_klm_FT_d is also reshaped, and the lineshapes are
    ## applied. This consists of an FT of the LDM, a product by the
    ## lineshape FTs & summing all G and L axes, and an inverse FT
    ## on the accumulated spectra.

    if verbose >= 2:
        print("done!")
        print("Applying lineshapes...")

    S_klm_FT_shape = (init_h.N_x_FT, iter_h.N_G, iter_h.N_L)
    cu_mod.fft_fwd.arr_out.resize(S_klm_FT_shape)
    cu_mod.fft_fwd()
    cu_mod.timer.lap("fft_fwd")

    cu_mod.applyLineshapes()
    cu_mod.timer.lap("applyLineshapes")

    cu_mod.fft_rev()
    cu_mod.timer.lap("fft_rev")

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
    cu_mod.timer.lap("calcTransmittanceNoslit")

    cu_mod.fft_fwd2()
    cu_mod.timer.lap("fft_fwd2")

    cu_mod.applyGaussianSlit()
    cu_mod.timer.lap("applyGaussianSlit")

    cu_mod.fft_rev2()
    cu_mod.timer.lap("fft_rev2")

    transmittance_h = cu_mod.fft_rev2.arr_out.getArray()[: init_h.N_v]

    if verbose >= 2:
        print("done!")

    if verbose == 1:
        print("Finished calculating spectrum!")

    cu_mod.timer.lap("total")
    times = cu_mod.timer.getTimes()

    return abscoeff_h, transmittance_h, iter_h, times


def gpu_exit(event=None):
    global cu_mod
    cu_mod.context.destroy()
    cu_mod = None