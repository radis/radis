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

gpu_mod = None
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
    gamma_arr,
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
    gamma_arr : numpy.ndarray[np.float32]
        (m,n) shaped array with Lorentzian width parameters, with n the number of lines in the database
        and m the number of collision partners included. This is usually at least two,
        with the first (m=0) always self broadening and the last (m=-1) always air broadening.
    iso : numpy.ndarray[np.uint8]
        Index of isotopologue.
    Mm_arr : numpy.ndarray
        Molecular masses for all isotopologues in the database (Mm_arr[0] is always 0).
    Q_intp_list : list
        List of Q branch interpolators.
    verbose : bool, optional
        Print verbosity level. Default is 0.
    backend :  ``'gpu-cuda'``, ``'cpu-cuda'``, optional
        Which backend to use; currently only CUDA backends (Nvidia) are supported. ``'cpu-cuda'`` runs the kernel on CPU. Default is ``'gpu-cuda'``.

    Returns
    -------
    init_h : radis.gpu.structs.initData_t
        structue with parameters used for GPU computation that are constant
        during iterations.
    """

    global gpu_mod
    if gpu_mod is not None:
        warn("Only a single GPU context allowed; please call gpu_exit() first.")
        return

    ## First a GPU context is created, then the .ptx file is read
    ## and made available as the GPUModule object gpu_mod
    ## If this fails, None is returned and calculations are
    ## defaulted to CPU emulation

    if backend == "cpu-cuda":
        from radis.gpu.cuda.emulate import CuContext as GPUContext

        ctx = GPUContext.Open(verbose=verbose)
        import radis.gpu.cuda.emulate as backend_module

    else:
        # Try to load GPU
        from radis.gpu.cuda.driver import CuContext as GPUContext

        ctx = GPUContext.Open(verbose=verbose)  # Set verbose to >=2 for comments
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
            from radis.gpu.cuda.emulate import CuContext as GPUContext

            ctx = GPUContext.Open(verbose=verbose)
            import radis.gpu.cuda.emulate as backend_module

        else:
            # successfully initialized CUDA context, continue with GPU:
            import radis.gpu.cuda.driver as backend_module

    GPUContext, GPUModule, GPUArray, GPUFFT, GPUTimer = backend_module.getClasses()

    if verbose:
        print("Number of lines loaded: {0}".format(len(v0)))
        print()

    ptx_path = os.path.join(getProjectRoot(), "gpu", "cuda", "build", "kernels.ptx")
    if not os.path.exists(ptx_path):
        raise FileNotFoundError(ptx_path)
    gpu_mod = GPUModule(ctx, ptx_path)  # gpu
    if verbose:
        print("mode:", gpu_mod.getMode())

    ## Next, the GPU is made aware of a number of parameters.
    ## Parameters that don't change during iteration are stored
    ## in init_h. They are copied to the GPU through gpu_mod.setConstant()

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
    init_h.N_collision_partners = gamma_arr.shape[0]

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

    gpu_mod.setConstant("init_d", init_h)

    init_G_params(log_2vMm.astype(np.float32), verbose)
    init_L_params(na, gamma_arr, verbose)

    if verbose >= 2:
        print("done!")

    ## Next the block- and thread size of the GPU kernels are set.
    ## This determines how the GPU internally divides up the work.

    if verbose >= 2:
        print("Allocating device memory and copying data...")

    NvFT = init_h.N_v_FT
    NxFT = NvFT // 2 + 1
    Ntpb = ctx.getMaxThreadsPerBlock()
    Nli = init_h.N_lines
    threads = (Ntpb, 1, 1)

    gpu_mod.fillLDM.setGrid((Nli // Ntpb + 1, 1, 1), threads)
    gpu_mod.applyLineshapes.setGrid((NxFT // Ntpb + 1, 1, 1), threads)
    gpu_mod.calcTransmittanceNoslit.setGrid((NvFT // Ntpb + 1, 1, 1), threads)
    gpu_mod.applyGaussianSlit.setGrid((NxFT // Ntpb + 1, 1, 1), threads)

    ## Next the variables are initialized on the GPU. Constant variables
    ## that don't change (i.e. pertaining to the database) are immediately
    ## copied to the GPU through GPUArray.fromArray().
    ## Other variables are only allocated. S_klm_d and S_klm_FT_d are
    ## special cases because their shape changes during iteration.
    ## They are not allocated, only given a device pointer by which
    ## they can be referenced later.

    S_klm_d = GPUArray(0, dtype=np.float32, grow_only=True)
    S_klm_FT_d = GPUArray(0, dtype=np.complex64, grow_only=True)

    spectrum_in_d = GPUArray(NxFT, dtype=np.complex64)
    spectrum_out_d = GPUArray(NvFT, dtype=np.float32)

    transmittance_noslit_d = GPUArray(NvFT, dtype=np.float32)
    transmittance_noslit_FT_d = GPUArray(NxFT, dtype=np.complex64)

    transmittance_FT_d = GPUArray(NxFT, dtype=np.complex64)
    transmittance_d = GPUArray(NvFT, dtype=np.float32)

    gpu_mod.fillLDM.setArgs(
        GPUArray.fromArray(iso),
        GPUArray.fromArray(v0),
        GPUArray.fromArray(da),
        GPUArray.fromArray(S0),
        GPUArray.fromArray(El),
        GPUArray.fromArray(gamma_arr),
        GPUArray.fromArray(na),
        S_klm_d,
    )
    gpu_mod.applyLineshapes.setArgs(S_klm_FT_d, spectrum_in_d)
    gpu_mod.calcTransmittanceNoslit.setArgs(spectrum_out_d, transmittance_noslit_d)
    gpu_mod.applyGaussianSlit.setArgs(transmittance_noslit_FT_d, transmittance_FT_d)

    ## FFT's are performed through the GPUFFT object. The required functions are internally
    ## loaded from the cufft library, not through the user kernels (.ptx files).
    ## The FFT's need some memory as "work area". Because the different FFT's can
    ## reuse the work area, we make a GPUArray at this scope that is passed to the
    ## GPUFFT objects. The work area will be scaled according to needs by the GPUFFT objects,
    ## so it can be initialized with a small value.

    workarea_d = GPUArray(0, dtype=np.byte, grow_only=True)
    gpu_mod.fft_fwd = GPUFFT(S_klm_d, S_klm_FT_d, workarea=workarea_d, direction="fwd")
    gpu_mod.fft_rev = GPUFFT(
        spectrum_in_d, spectrum_out_d, workarea=workarea_d, direction="rev"
    )
    gpu_mod.fft_fwd2 = GPUFFT(
        transmittance_noslit_d,
        transmittance_noslit_FT_d,
        workarea=workarea_d,
        direction="fwd",
    )
    gpu_mod.fft_rev2 = GPUFFT(
        transmittance_FT_d, transmittance_d, workarea=workarea_d, direction="rev"
    )

    gpu_mod.timer = GPUTimer()

    if verbose >= 2:
        print("done!")

    return init_h


def gpu_iterate(
    p,
    T,
    mole_fraction,
    verbose=0,
    # for GPU instrument functions (not currently supported):
    l=1.0,
    slit_FWHM=0.0,
):
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
    verbose : int, optional
        The default is 0.


    Returns
    -------
    abscoeff_h : numpy.ndarray[np.float32]
        array with absorbtion coefficients in (cm.-1)
    iter_h : radis.gpu.structs.iterData_t
        structue with parameters used for computation of abscoeff_h.
    times : dict
        dictionary with computation cumulative computation times for
        different stages of the GPU computation. The ``'total'`` key
        gives the total time.
    """

    if gpu_mod is None:
        warn("Must have an open GPU context; please call gpu_init() first.")
        return

    if verbose >= 2:
        print("Copying iteration parameters to device...")

    ## First a number of parameters that change during iteration
    ## are computed and copied to the GPU.

    gpu_mod.timer.reset()

    set_pTQ(p, T, mole_fraction, iter_h, l=l, slit_FWHM=slit_FWHM)
    set_G_params(init_h, iter_h)
    set_L_params(init_h, iter_h)
    gpu_mod.setConstant("iter_d", iter_h)
    gpu_mod.timer.lap("iter_params")

    ## Next the S_klm_d variable is reshaped to the correct shape,
    ## and filled with spectral data.

    if verbose >= 2:
        print("done!")
        print("Filling LDM...")

    S_klm_shape = (init_h.N_v_FT, iter_h.N_G, iter_h.N_L)

    gpu_mod.fillLDM.args[-1].resize(S_klm_shape, init="zeros")
    gpu_mod.fillLDM()
    gpu_mod.timer.lap("fillLDM")

    ## Next the S_klm_FT_d is also reshaped, and the lineshapes are
    ## applied. This consists of an FT of the LDM, a product by the
    ## lineshape FTs & summing all G and L axes, and an inverse FT
    ## on the accumulated spectra.

    if verbose >= 2:
        print("done!")
        print("Applying lineshapes...")

    S_klm_FT_shape = (init_h.N_x_FT, iter_h.N_G, iter_h.N_L)
    gpu_mod.fft_fwd.arr_out.resize(S_klm_FT_shape)
    gpu_mod.fft_fwd()
    gpu_mod.timer.lap("fft_fwd")

    gpu_mod.applyLineshapes()
    gpu_mod.timer.lap("applyLineshapes")

    gpu_mod.fft_rev()
    gpu_mod.timer.lap("fft_rev")

    if verbose >= 2:
        print("Done!")
        print("Calculating transmittance...")

    abscoeff_h = gpu_mod.fft_rev.arr_out.getArray()[: init_h.N_v]

    ## To apply a slit function, first the transmittance is calculated.
    ## Then the convolution is applied by an FT, product with the
    ## instrument function's FT, followed by an inverse FT.

    if verbose >= 2:
        print("Done!")

    ##    ##The code below is to process slits on the GPU, which is currently unsupported.
    ##
    ##        print("Applying slit function...")
    ##
    ##    gpu_mod.calcTransmittanceNoslit()
    ##    gpu_mod.timer.lap("calcTransmittanceNoslit")
    ##
    ##    gpu_mod.fft_fwd2()
    ##    gpu_mod.timer.lap("fft_fwd2")
    ##
    ##    gpu_mod.applyGaussianSlit()
    ##    gpu_mod.timer.lap("applyGaussianSlit")
    ##
    ##    gpu_mod.fft_rev2()
    ##    gpu_mod.timer.lap("fft_rev2")
    ##
    ##    transmittance_h = gpu_mod.fft_rev2.arr_out.getArray()[: init_h.N_v]
    ##
    ##    if verbose >= 2:
    ##        print("done!")

    if verbose == 1:
        print("Finished calculating spectrum!")

    gpu_mod.timer.lap("total")
    times = gpu_mod.timer.getTimes()

    ##    diffs = gpu_mod.timer.getDiffs()
    ##    print(diffs['fillLDM'], diffs['fft_fwd'])

    return abscoeff_h, iter_h, times


def gpu_exit(event=None):
    global gpu_mod
    gpu_mod.context.destroy()
    gpu_mod = None
