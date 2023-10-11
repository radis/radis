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

    Returns
    -------
    init_h : radis.gpu.structs.initData_t
        structue with parameters used for GPU computation that are constant
        during iterations.
    """

    global app
    
    from radis.gpu.vulkan.vulkan_compute_lib import ComputeApplication, StructBuffer,  ArrayBuffer    
    app = ComputeApplication(deviceID=0)

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

    # gpu_mod.setConstant("init_d", init_h)
    app.init_d = StructBuffer.fromStruct(init_h, binding=0, app=app)
    app.iter_d = StructBuffer.fromStruct(iter_h, binding=1, app=app)

    init_G_params(log_2vMm.astype(np.float32), verbose)
    init_L_params(na, gamma_arr, verbose)
    
    #Calulate params once to obtain N_G_max and N_L_max:
    p_max = 5.0 #bar #TODO: obtain this from defaults/keywords
    T_max = 3500.0 #K
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

    NvFT = init_h.N_v_FT
    NxFT = NvFT // 2 + 1
    Ntpb = 1024 #TODO: Get these through vulkan
    Nli = init_h.N_lines
    threads = (Ntpb, 1, 1, Nli)

    # gpu_mod.fillLDM.setGrid((Nli // Ntpb + 1, 1, 1), threads)
    # gpu_mod.applyLineshapes.setGrid((NxFT // Ntpb + 1, 1, 1), threads)
    #gpu_mod.calcTransmittanceNoslit.setGrid((NvFT // Ntpb + 1, 1, 1), threads)
    #gpu_mod.applyGaussianSlit.setGrid((NxFT // Ntpb + 1, 1, 1), threads)

    ## Next the variables are initialized on the GPU. Constant variables
    ## that don't change (i.e. pertaining to the database) are immediately
    ## copied to the GPU through GPUArray.fromArray().
    ## Other variables are only allocated. S_klm_d and S_klm_FT_d are
    ## special cases because their shape changes during iteration.
    ## They are not allocated, only given a device pointer by which
    ## they can be referenced later.


    database_arrays = [iso, v0, da, S0, El, na, gamma_arr]
    database_length = np.sum([np.sum(arr.shape[:-1]) if len(arr.shape)>1 else 1 for arr in database_arrays])

    database_SSBO = ArrayBuffer((database_length, Nli),
                                    np.float32, #TODO: Not all arrays are np.float32, but it doesn't matter because all dtypes have size 4 (=sizeof(float32)), and we never read this buffer.
                                    binding=2, app=app)
                                    
    byte_offset = 0
    for arr in database_arrays:
        byte_offset += database_SSBO.setData(arr, byte_offset=byte_offset)
        
    S_klm_d = ArrayBuffer((N_L_max, N_G_max, 2*NxFT), np.float32, binding=3, app=app)
    S_klm_FT_d = ArrayBuffer((N_L_max, N_G_max, NxFT), np.complex64, binding=4, app=app)

    I_k_FT_d = ArrayBuffer((NxFT,), np.complex64, binding=5, app=app)
    app.I_k_d = ArrayBuffer((2*NxFT,), np.float32, binding=6, app=app)
    
    from radis.gpu.vulkan.pyvkfft_vulkan import prepare_fft
    app.fft_LDM = prepare_fft(S_klm_d, S_klm_FT_d, compute_app=app)
    app.fft_spec = prepare_fft(app.I_k_d, I_k_FT_d, compute_app=app)
    
    # Write command buffer:
    shader_path = os.path.join(getProjectRoot(), "gpu", "vulkan")

    app.schedule_shader(os.path.join(shader_path,'fillLDM.spv'), 
                        (Nli // Ntpb + 1, 1, 1), (Ntpb, 1, 1, Nli))
    app.fft_LDM.fft(app._commandBuffer, S_klm_d._buffer, S_klm_FT_d._buffer)
    app.schedule_shader(os.path.join(shader_path, 'applyLineshapes.spv'), 
                    (NxFT // Ntpb + 1, 1, 1), (Ntpb, 1, 1))
    app.fft_spec.ifft(app._commandBuffer, I_k_FT_d._buffer, app.I_k_d._buffer)

    app.endCommandBuffer()

    
    
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

    # if gpu_mod is None:
        # warn("Must have an open GPU context; please call gpu_init() first.")
        # return

    if verbose >= 2:
        print("Copying iteration parameters to device...")

    set_pTQ(p, T, mole_fraction, iter_h, l=l, slit_FWHM=slit_FWHM)
    set_G_params(init_h, iter_h)
    set_L_params(init_h, iter_h)    
    app.iter_d.setData(iter_h)
    
    if verbose >= 2:
        print("Running compute pipelines...")
    
    app.run()
    abscoeff_h = app.I_k_d.getData()[: init_h.N_v]
    
    if verbose == 1:
        print("Finished calculating spectrum!")

    times = [0, 0, 0, 0]

    return abscoeff_h, iter_h, times


def gpu_exit(event=None):
    global app
    del app
    