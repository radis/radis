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
from radis.gpu.structs import initData_t, iterData_t, workGroupSizeArray_t
from radis.misc.utils import getProjectRoot
from ctypes import sizeof

# from radis.misc.warning import NoGPUWarning




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
    T_max_parsum=None,
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
        The id of the selected GPU. Check the console output for more details. ``'cpu-cuda'`` runs the kernel on CPU. Default is ``'gpu-cuda'``.
    T_max_parsum : int
        The maximum temperature at which the partition function of all the isotopologues is calculated (at higher T, the partition function should raise an error).
    Returns
    -------
    init_h : radis.gpu.structs.initData_t
        structue with parameters used for GPU computation that are constant
        during iterations.
    """

    global app, init_h, iter_h

    from radis.gpu.vulkan.vulkan_compute_lib import GPUApplication, GPUBuffer

    shader_path = os.path.join(getProjectRoot(), "gpu", "vulkan", "shaders")
    app = GPUApplication(deviceID=device_id, path=shader_path, verbose=verbose)

    ## Next, the GPU is made aware of a number of parameters.
    ## Parameters that don't change during iteration are stored
    ## in init_h. They are copied to the GPU through gpu_mod.setConstant()

    if verbose >= 2:
        print("Copying initialization parameters to device memory...")
    
    gpu_mod = None

    app.init_d = GPUBuffer(sizeof(initData_t), usage='uniform', binding=0)  #host_visible, device_local, (host_cached)
    app.iter_d = GPUBuffer(sizeof(iterData_t), usage='uniform', binding=1)  #host_visible, device_local, (host_cached) 

    # init_h = initData_t()
    # iter_h = iterData_t()
    
    # app.init_d = GPUStruct.fromStruct(init_h, binding=0)
    # app.iter_d = GPUStruct.fromStruct(iter_h, binding=1)



    init_h = app.init_d.getHostStructPtr(initData_t)
    iter_h = app.iter_d.getHostStructPtr(iterData_t)

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
    if T_max_parsum is not None:
        T_max = T_max_parsum  # K
    else:
        T_max = 3500.0  # default value for now

        
    set_pTQ(p_max, T_max, 0.5, iter_h, l=1.0, slit_FWHM=0.0)
    set_G_params(init_h, iter_h)
    set_L_params(init_h, iter_h)

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



    
    # app.database_d = GPUArray(
    #     (N_db, init_h.N_lines),
    #     np.float32,  # Not all arrays are np.float32, but it doesn't matter because all dtypes have size 4 (=sizeof(float32)), and the host never reads this buffer.
    #     binding=2,
    # )
    app.database_d = GPUBuffer(N_db * v0.nbytes, binding=2)

    byte_offset = 0
    for arr in database_arrays:
        #byte_offset += app.database_d.setData(arr, byte_offset=byte_offset)
        byte_offset += app.database_d.copyToBuffer(arr, device_offset=byte_offset, chunksize=32*1024*1024)
        
    app.S_kl_d = GPUBuffer(fftSize=init_h.N_v_FT, binding=3) #req: large, device_local
    app.spectrum_d = GPUBuffer(fftSize=init_h.N_v_FT, binding=4) #req: device_local, host_visible, host_cached, (large)

    app.indirect_d = GPUBuffer(sizeof(workGroupSizeArray_t), usage='indirect', binding=10) #host_visible, device_local, (host_cached)
    indirect_h = app.setIndirectBuffer(app.indirect_d, workGroupSizeArray_t) #TODO: set fwd/inv


    # app.S_klm_d = GPUArray((N_L_max, N_G_max, init_h.N_v_FT), np.float32, binding=3)
    # app.S_klm_FT_d = GPUArray(
    #     (N_L_max, N_G_max, init_h.N_x_FT), np.complex64, binding=4
    # )

    # app.spectrum_FT_d = GPUArray((init_h.N_x_FT,), np.complex64, binding=5)
    # app.spectrum_d = GPUArray((init_h.N_v_FT,), np.float32, binding=6)

    # Write command buffer:
    N_tpb = 128  # threads per block
    threads = (N_tpb, 1, 1)

    # app.command_list = [
    #     app.cmdAddTimestamp("start"),
    #     app.cmdClearBuffer(app.S_klm_d, timestamp=True),
    #     app.cmdFillLDM((init_h.N_lines // N_tpb + 1, 1, 1), threads, timestamp=True),
    #     app.cmdClearBuffer(app.S_klm_FT_d, timestamp=True),
    #     app.cmdFFT(app.S_klm_d, app.S_klm_FT_d, timestamp=True),
    #     # app.cmdClearBuffer(app.spectrum_FT_d, timestamp=True),
    #     app.cmdApplyLineshapes(
    #         (init_h.N_x_FT // N_tpb + 1, 1, 1), threads, timestamp=True
    #     ),
    #     app.cmdClearBuffer(app.spectrum_d, timestamp=True),
    #     app.cmdIFFT(app.spectrum_FT_d, app.spectrum_d, timestamp=True),
    # ]

    app.appendCommands([
        app.cmdAddTimestamp('Initialization'),
        app.indirect_d.cmdTransferStagingBuffer('H2D'),
        app.iter_d.cmdTransferStagingBuffer('H2D'),
        app.cmdAddTimestamp('Zeroing buffers'),
        app.S_kl_d.cmdClearBuffer(),
        app.spectrum_d.cmdClearBuffer(),
        app.cmdAddTimestamp('Line addition'),
        app.cmdScheduleShader('cmdFillLDM.spv', (init_h.N_lines // N_tpb + 1, 1, 1), threads),
        app.cmdAddTimestamp('FFT fwd'),
        app.cmdFFT(app.S_kl_d, name='FFTa'),
        app.cmdAddTimestamp('Apply conv.'),
        app.cmdScheduleShader('cmdApplyLineshapes.spv', (init_h.N_x_FT // N_tpb + 1, 1, 1), threads),
        #app.cmdScheduleShader('cmdTestApplyLineshapesP.spv', (init_h.N_x_FT // N_tpb + 1, N_G*N_L, 1), threads),
        app.cmdAddTimestamp('FFT inv'),
        app.cmdIFFT(app.spectrum_d, name='FFTb'), 
        app.spectrum_d.cmdTransferStagingBuffer('D2H'),
        app.cmdAddTimestamp('End'),
    ])

    app.writeCommandBuffer()

    app.registerBatchSizeUpdateFunction(app.S_kl_d.setBatchSize)
    app.registerBatchSizeUpdateFunction(app.setFwdFFTWorkGroupSize)


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
    #app.iter_d.setData(iter_h)

    if verbose >= 2:
        print("Running compute pipelines...")

    app.setBatchSize(iter_h.N_G * iter_h.N_L)
    app.run()
    gpu_times = app.get_timestamps()

    #abscoeff_h = np.copy(app.spectrum_d.getData()[: init_h.N_v])
    abscoeff_h = np.zeros(init_h.N_v, dtype=np.float32)
    app.spectrum_d.toArray(abscoeff_h)

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
