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

# from radis.misc.warning import NoGPUWarning


gpu_mod = None
init_h = initData_t()
iter_h = iterData_t()


def next_fast_len_even(n):
    n = next_fast_len(n)
    while n & 1:
        n = next_fast_len(n + 1)
    return n


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
    backend="gpu-vulkan",
    device_id=0,
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
    iso : numpy.ndarray[np.uint32]
        Index of isotopologue.
    Mm_arr : numpy.ndarray
        Molecular masses for all isotopologues in the database (Mm_arr[0] is always 0).
    Q_intp_list : list
        List of Q branch interpolators.
    verbose : bool, optional
        Print verbosity level. Default is 0.
    backend :  ``'gpu-cuda'``, ``'cpu-cuda'``, optional
        Which backend to use; currently only CUDA backends (Nvidia) are supported. ``'cpu-cuda'`` runs the kernel on CPU. Default is ``'gpu-cuda'``.
    device_id : int
        The id of the selected GPU. Check the console output for more details.
    Returns
    -------
    init_h : radis.gpu.structs.initData_t
        structue with parameters used for GPU computation that are constant
        during iterations.
    """

    global app

    from radis.gpu.vulkan.vulkan_compute_lib import GPUApplication, GPUArray, GPUStruct

    shader_path = os.path.join(getProjectRoot(), "gpu", "vulkan", "shaders")
    app = GPUApplication(deviceID=device_id, path=shader_path, verbose=verbose)

    ## Next, the GPU is made aware of a number of parameters.
    ## Parameters that don't change during iteration are stored
    ## in init_h. They are copied to the GPU through gpu_mod.setConstant()

    if verbose >= 2:
        print("Copying initialization parameters to device memory...")

    init_h.v_min = vmin
    init_h.dv = dv
    init_h.N_v = Nv
    init_h.N_v_FT = next_fast_len_even(2 * init_h.N_v)
    init_h.N_x_FT = init_h.N_v_FT // 2 + 1
    init_h.dxG = dxG
    init_h.dxL = dxL
    init_h.N_lines = int(len(v0))
    init_h.N_coll = gamma_arr.shape[0]

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

    # gpu_mod.setConstant("init_d", init_h)

    init_G_params(log_2vMm.astype(np.float32), verbose)
    init_L_params(na, gamma_arr, verbose)

    # Calulate params once to obtain N_G_max and N_L_max:
    p_max = 5.0  # bar #TODO: obtain this from defaults/keywords
    T_max = 3500.0  # K
    set_pTQ(p_max, T_max, 0.5, iter_h, l=1.0, slit_FWHM=0.0)
    set_G_params(init_h, iter_h)
    set_L_params(init_h, iter_h)
    N_G_max = iter_h.N_G
    N_L_max = iter_h.N_L

    if verbose >= 2:
        print("done!")

    ## Next the block- and thread size of the GPU kernels are set.
    ## This determines how the GPU internally divides up the work.

    if verbose >= 2:
        print("Allocating device memory and copying data...")

    ## Next the variables are initialized on the GPU. Constant variables
    ## that don't change (i.e. pertaining to the database) are immediately
    ## copied to the GPU through GPUArray.fromArray().
    ## Other variables are only allocated. S_klm_d and S_klm_FT_d are
    ## special cases because their shape changes during iteration.
    ## They are not allocated, only given a device pointer by which
    ## they can be referenced later.

    database_arrays = [iso, v0, da, S0, El, na, gamma_arr]
    N_db = np.sum(
        [np.sum(arr.shape[:-1]) if len(arr.shape) > 1 else 1 for arr in database_arrays]
    )

    app.init_d = GPUStruct.fromStruct(init_h, binding=0)
    app.iter_d = GPUStruct.fromStruct(iter_h, binding=1)

    app.database_d = GPUArray(
        (N_db, init_h.N_lines),
        np.float32,  # Not all arrays are np.float32, but it doesn't matter because all dtypes have size 4 (=sizeof(float32)), and the host never reads this buffer.
        binding=2,
    )

    byte_offset = 0
    for arr in database_arrays:
        byte_offset += app.database_d.setData(arr, byte_offset=byte_offset)

    app.S_klm_d = GPUArray((N_L_max, N_G_max, init_h.N_v_FT), np.float32, binding=3)
    app.S_klm_FT_d = GPUArray(
        (N_L_max, N_G_max, init_h.N_x_FT), np.complex64, binding=4
    )

    app.spectrum_FT_d = GPUArray((init_h.N_x_FT,), np.complex64, binding=5)
    app.spectrum_d = GPUArray((init_h.N_v_FT,), np.float32, binding=6)

    # Write command buffer:
    N_tpb = 128  # threads per block
    threads = (N_tpb, 1, 1)

    app.command_list = [
        app.addTimestamp("start"),
        app.clearBuffer(app.S_klm_d, timestamp=True),
        app.fillLDM((init_h.N_lines // N_tpb + 1, 1, 1), threads, timestamp=True),
        app.clearBuffer(app.S_klm_FT_d, timestamp=True),
        app.fft(app.S_klm_d, app.S_klm_FT_d, timestamp=True),
        # app.clearBuffer(app.spectrum_FT_d, timestamp=True),
        app.applyLineshapes(
            (init_h.N_x_FT // N_tpb + 1, 1, 1), threads, timestamp=True
        ),
        app.clearBuffer(app.spectrum_d, timestamp=True),
        app.ifft(app.spectrum_FT_d, app.spectrum_d, timestamp=True),
    ]

    app.writeCommandBuffer()

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
    if app is None:
        warn("No GPUApplication initialized; please call gpu_init() first.")
        return

    if verbose >= 2:
        print("Copying iteration parameters to device...")

    set_pTQ(p, T, mole_fraction, iter_h, l=l, slit_FWHM=slit_FWHM)
    set_G_params(init_h, iter_h)
    set_L_params(init_h, iter_h)
    app.iter_d.setData(iter_h)

    if verbose >= 2:
        print("Running compute pipelines...")

    app.run()
    gpu_times = app.get_timestamps()

    abscoeff_h = np.copy(app.spectrum_d.getData()[: init_h.N_v])

    if verbose == 1:
        print("Finished calculating spectrum!")

    return abscoeff_h, iter_h, gpu_times


def gpu_exit(event=None):
    global app
    # TODO: free(), del, *and* =None might be redundant..
    if app:
        app.free()
        del app
    app = None
