import ctypes
import pickle
from os.path import join

import cupy as cp
import numpy as np

from radis.misc.utils import getProjectRoot


class initData(ctypes.Structure):
    _fields_ = [
        ("v_min", ctypes.c_float),
        ("v_max", ctypes.c_float),
        ("dv", ctypes.c_float),
        ("N_v", ctypes.c_int),
        ("N_wG", ctypes.c_int),
        ("N_wL", ctypes.c_int),
        ("N_wG_x_N_wL", ctypes.c_int),
        ("N_total", ctypes.c_int),
        ("Max_lines", ctypes.c_int),
        ("N_lines", ctypes.c_int),
        ("N_points_per_block", ctypes.c_int),
        ("N_threads_per_block", ctypes.c_int),
        ("N_blocks_per_grid", ctypes.c_int),
        ("N_points_per_thread", ctypes.c_int),
        ("Max_iterations_per_thread", ctypes.c_int),
        ("shared_size_floats", ctypes.c_int),
    ]


try:
    from radis_cython_extensions import blockData

except (ModuleNotFoundError):

    class blockData(ctypes.Structure):
        _fields_ = [("line_offset", ctypes.c_int), ("iv_offset", ctypes.c_int)]


class iterData(ctypes.Structure):
    _fields_ = [
        ("p", ctypes.c_float),
        ("log_p", ctypes.c_float),
        ("hlog_T", ctypes.c_float),
        ("log_rT", ctypes.c_float),
        ("c2T", ctypes.c_float),
        ("N", ctypes.c_float),
        ("log_wG_min", ctypes.c_float),
        ("log_wL_min", ctypes.c_float),
        ("log_dwG", ctypes.c_float),
        ("log_dwL", ctypes.c_float),
        ("blocks", blockData * 4096),
    ]


init_params_h = initData()
iter_params_h = iterData()

host_params_h_start = cp.cuda.Event()
host_params_h_stop = cp.cuda.Event()
host_params_h_start_DLM = cp.cuda.Event()
host_params_h_stop_DLM = cp.cuda.Event()
host_params_h_data_start = cp.cuda.Event()
host_params_h_data_stop = cp.cuda.Event()

# TO-DO: read and compile CUDA code at install time, then pickle the cuda object
cuda_fname = join(getProjectRoot(), "lbl", "gpu.cu")
with open(cuda_fname, "rb") as f:
    cuda_code = f.read().decode()

cuda_module = cp.RawModule(code=cuda_code)
fillDLM = cuda_module.get_function("fillDLM")
applyLineshapes = cuda_module.get_function("applyLineshapes")


def py_calc_lorentzian_envelope_params(na, log_2gs, verbose=False):
    # Remove duplicates
    unique_lines = set([])
    for i in range(len(na)):
        unique_lines.add(str(na[i]) + " " + str(log_2gs[i]))

    # Only keep extremes
    max_dict = {}
    min_dict = {}
    for s in unique_lines:
        na_i, log_2gs_i = map(float, s.split())
        try:
            max_dict[na_i] = log_2gs_i if log_2gs_i > max_dict[na_i] else max_dict[na_i]
        except (KeyError):
            max_dict[na_i] = log_2gs_i

        try:
            min_dict[na_i] = log_2gs_i if log_2gs_i < min_dict[na_i] else min_dict[na_i]
        except (KeyError):
            min_dict[na_i] = log_2gs_i

    # Check which ones are really at the top:
    keys = sorted(max_dict.keys())

    topA = [keys[0]]
    topB = [max_dict[keys[0]]]
    topX = [-np.inf]

    for key in keys[1:]:
        for i in range(len(topX)):
            xi = (max_dict[key] - topB[i]) / (topA[i] - key)
            if xi >= topX[i]:
                if i < len(topX) - 1:
                    if xi < topX[i + 1]:
                        break
                else:
                    break

        topA = topA[: i + 1] + [key]
        topB = topB[: i + 1] + [max_dict[key]]
        topX = topX[: i + 1] + [xi]

    topX = topX[1:] + [np.inf]

    # Check which ones are really at the bottom:
    keys = sorted(min_dict.keys())[::-1]

    bottomA = [keys[0]]
    bottomB = [min_dict[keys[0]]]
    bottomX = [-np.inf]

    for key in keys[1:]:
        for i in range(len(bottomX)):
            xi = (min_dict[key] - bottomB[i]) / (bottomA[i] - key)
            if xi >= bottomX[i]:
                if i < len(bottomX) - 1:
                    if xi < bottomX[i + 1]:
                        break
                else:
                    break

        bottomA = bottomA[: i + 1] + [key]
        bottomB = bottomB[: i + 1] + [min_dict[key]]
        bottomX = bottomX[: i + 1] + [xi]

    bottomX = bottomX[1:] + [np.inf]

    return (topA, topB, topX, bottomA, bottomB, bottomX)


def py_calc_gaussian_envelope_params(log_2vMm, verbose=False):
    return np.amin(log_2vMm), np.max(log_2vMm)


try:
    from radis_cython_extensions import (
        calc_gaussian_envelope_params,
        calc_lorentzian_envelope_params,
        prepare_blocks,
    )
except (ModuleNotFoundError):
    calc_gaussian_envelope_params = py_calc_gaussian_envelope_params
    calc_lorentzian_envelope_params = py_calc_lorentzian_envelope_params
    # TO-DO: currently we don't have a purely Python implementation for prepare_blocks
    # prepare_blocks = py_prepare_blocks


def init_gaussian_params(log_2vMm, verbose_gpu):

    if verbose_gpu >= 2:
        print("Initializing Gaussian parameters")

    fname = "Gaussian_minmax_" + str(len(log_2vMm)) + ".dat"
    try:
        param_data = pickle.load(open(fname, "rb"))
        if verbose_gpu >= 2:
            print(" (from cache)... ")

    except (OSError, IOError):
        if verbose_gpu >= 2:
            print("... ")

        param_data = calc_gaussian_envelope_params(log_2vMm, verbose_gpu)
        pickle.dump(param_data, open(fname, "wb"))

    if verbose_gpu >= 2:
        print("done!")

    return param_data


def init_lorentzian_params(log_2gs, na, verbose_gpu):

    if verbose_gpu >= 2:
        print("Initializing Lorentzian parameters ")

    fname = "Lorenzian_minmax_" + str(len(log_2gs)) + ".dat"

    try:
        with open(fname, "rb") as f:
            if verbose_gpu >= 2:
                print(" (from cache)... ")
            param_data = pickle.load(f)

    except:
        if verbose_gpu >= 2:
            print(" ... ")

        param_data = calc_lorentzian_envelope_params(log_2gs, na, verbose_gpu)
        with open(fname, "wb") as f:
            pickle.dump(param_data, f)

    if verbose_gpu >= 2:
        print("done!")

    return param_data


def calc_gaussian_params(
    gaussian_param_data, init_params_h, iter_params_h, epsilon=1e-4
):

    host_params_h_log_2vMm_min, host_params_h_log_2vMm_max = gaussian_param_data
    log_wG_min = host_params_h_log_2vMm_min + iter_params_h.hlog_T
    log_wG_max = host_params_h_log_2vMm_max + iter_params_h.hlog_T + epsilon
    log_dwG = (log_wG_max - log_wG_min) / (init_params_h.N_wG - 1)

    iter_params_h.log_wG_min = log_wG_min
    iter_params_h.log_dwG = log_dwG


def calc_lorentzian_params(param_data, init_params_h, iter_params_h, epsilon=1e-4):

    (
        top_x,
        top_a,
        top_b,
        bottom_x,
        bottom_a,
        bottom_b,
    ) = param_data

    for i in range(bottom_x.size):
        if iter_params_h.log_rT < bottom_x[i]:
            log_wL_min = (
                iter_params_h.log_rT * bottom_a[i] + bottom_b[i] + iter_params_h.log_p
            )
            break

    for i in range(top_x.size):
        if iter_params_h.log_rT < top_x[i]:
            log_wL_max = (
                iter_params_h.log_rT * top_a[i]
                + top_b[i]
                + iter_params_h.log_p
                + epsilon
            )
            break

    log_dwL = (log_wL_max - log_wL_min) / (init_params_h.N_wL - 1)
    log_wL_min, log_dwL = -3.8011555671691895, 0.07515250146389008
    iter_params_h.log_wL_min = log_wL_min
    iter_params_h.log_dwL = log_dwL


def set_pT(p, T, mole_fraction, iter_params_h):

    c2 = 1.4387773538277204  # K.cm
    k = 1.38064852e-23  # J.K-1
    iter_params_h.p = p  # bar
    iter_params_h.log_p = np.log(p)
    iter_params_h.hlog_T = 0.5 * np.log(T)
    iter_params_h.log_rT = np.log(296.0 / T)
    iter_params_h.c2T = -c2 / T
    iter_params_h.N = mole_fraction * p * 1e5 / (1e6 * k * T)  # cm-3

    ## TO-DO: These are molecule/isotopologue specific params and should not be compiled
    # cdef float B  = <float>     0.3902 #cm-1
    # cdef float w1 = <float>  1354.31 #cm-1
    # cdef float w2 = <float>   672.85 #cm-1
    # cdef float w3 = <float>  2396.32 #cm-1

    # cdef int d1 = 1
    # cdef int d2 = 2
    # cdef int d3 = 1
    # cdef float gr = 0.5

    # cdef float Trot = T
    # cdef float Tv12 = T
    # cdef float Tv3  = T

    # cdef float Qr = gr * Trot/(c2 * B)*np.exp(c2*B/(<float>3.0*Trot)) #McDowell 1978
    # cdef float Qv1 = 1 / np.power(1 - np.exp(-c2 * w1 / Tv12), d1)
    # cdef float Qv2 = 1 / np.power(1 - np.exp(-c2 * w2 / Tv12), d2)
    # cdef float Qv3 = 1 / np.power(1 - np.exp(-c2 * w3 / Tv3 ), d3)
    # cdef float Q = Qr * Qv1 * Qv2 * Qv3

    # iter_params_h.Q = Q


def set_constant_memory(var_str):
    var_h = globals()[var_str + "_h"]
    memptr_d = cuda_module.get_global(var_str + "_d")
    ptr = ctypes.cast(ctypes.pointer(var_h), ctypes.c_void_p)
    struct_size = ctypes.sizeof(var_h)
    memptr_d.copy_from_host(ptr, struct_size)


def gpu_init(
    v_arr, N_wG, N_wL, iso, v0, da, log_2gs, na, log_2vMm, S0, El, Q, verbose_gpu
):

    # ----------- setup global variables -----------------
    global init_params_h
    global host_params_h_dec_size
    global host_params_h_block_preparation_step_size
    global host_params_h_iso_d
    global host_params_h_v0_d
    global host_params_h_v0_dec
    global host_params_h_da_d
    global host_params_h_da_dec
    global host_params_h_S0_d
    global host_params_h_El_d
    global host_params_h_Q_d
    global host_params_h_log_2gs_d
    global host_params_h_na_d
    global host_params_h_log_2vMm_d
    global host_params_h_DLM_d_in
    global host_params_h_spectrum_d_in
    global host_params_h_I_add
    global host_params_h_data_start
    global host_params_h_data_stop
    global host_params_h_elapsedTimeData

    global lorentzian_param_data
    global gaussian_param_data

    global cuda_module
    global database_path
    global N_lines_to_load
    # -----------------------------------------------------

    init_params_h.v_min = np.min(v_arr)  # 2000.0
    init_params_h.v_max = np.max(v_arr)  # 2400.0
    init_params_h.dv = (v_arr[-1] - v_arr[0]) / (len(v_arr) - 1)  # 0.002
    init_params_h.N_v = len(
        v_arr
    )  # int((init_params_h.v_max - init_params_h.v_min)/init_params_h.dv)

    init_params_h.N_wG = N_wG
    init_params_h.N_wL = N_wL
    ##    spectrum_h = np.zeros(init_params_h.N_v, dtype=np.float32)

    init_params_h.Max_iterations_per_thread = 1024
    host_params_h_block_preparation_step_size = 128

    host_params_h_shared_size = 0x8000  # Bytes - Size of the shared memory
    ##    host_params_h_Min_threads_per_block = 128   # Ensures a full warp from each of the 4 processors
    ##    host_params_h_Max_threads_per_block = 1024  # Maximum determined by device parameters
    init_params_h.shared_size_floats = host_params_h_shared_size // 4  # size of float

    init_params_h.N_wG_x_N_wL = init_params_h.N_wG * init_params_h.N_wL
    init_params_h.N_total = init_params_h.N_wG_x_N_wL * init_params_h.N_v
    init_params_h.N_points_per_block = (
        init_params_h.shared_size_floats // init_params_h.N_wG_x_N_wL
    )

    init_params_h.N_threads_per_block = 1024
    init_params_h.N_blocks_per_grid = 4 * 256 * 256
    init_params_h.N_points_per_thread = (
        init_params_h.N_points_per_block // init_params_h.N_threads_per_block
    )

    if verbose_gpu >= 2:
        print()
        print(
            "Spectral points per block  : {0}".format(init_params_h.N_points_per_block)
        )
        print(
            "Threads per block          : {0}".format(init_params_h.N_threads_per_block)
        )
        print(
            "Spectral points per thread : {0}".format(init_params_h.N_points_per_thread)
        )
        print()

    host_params_h_v0_dec = np.zeros(
        len(v0) // init_params_h.N_threads_per_block, dtype=np.float32
    )
    for i in range(0, len(v0) // init_params_h.N_threads_per_block):
        host_params_h_v0_dec[i] = v0[i * init_params_h.N_threads_per_block]
    host_params_h_dec_size = host_params_h_v0_dec.size
    host_params_h_da_dec = np.zeros(
        len(v0) // init_params_h.N_threads_per_block, dtype=np.float32
    )
    for i in range(0, len(v0) // init_params_h.N_threads_per_block):
        host_params_h_da_dec[i] = da[i * init_params_h.N_threads_per_block]

    lorentzian_param_data = init_lorentzian_params(log_2gs, na, verbose_gpu)
    gaussian_param_data = init_gaussian_params(log_2vMm, verbose_gpu)
    init_params_h.N_lines = int(len(v0))

    if verbose_gpu == 1:
        print("Number of lines loaded: {0}".format(init_params_h.N_lines))
        print()

    if verbose_gpu >= 2:
        print("Allocating device memory and copying data...")

    host_params_h_data_start.record()

    host_params_h_DLM_d_in = cp.zeros(
        (2 * init_params_h.N_v, init_params_h.N_wG, init_params_h.N_wL),
        order="C",
        dtype=cp.float32,
    )
    host_params_h_spectrum_d_in = cp.zeros(init_params_h.N_v + 1, dtype=cp.complex64)

    host_params_h_I_add = cp.zeros(init_params_h.N_lines, dtype=cp.float32)

    if verbose_gpu >= 2:
        print("Copying initialization parameters to device memory...")
    set_constant_memory("init_params")

    if verbose_gpu >= 2:
        print("done!")
        print("Copying spectral data to device memory...")

    # #Copy spectral data to device
    host_params_h_iso_d = cp.array(iso)
    host_params_h_v0_d = cp.array(v0)
    host_params_h_da_d = cp.array(da)
    host_params_h_S0_d = cp.array(S0)
    host_params_h_El_d = cp.array(El)
    host_params_h_log_2gs_d = cp.array(log_2gs)
    host_params_h_na_d = cp.array(na)
    host_params_h_log_2vMm_d = cp.array(log_2vMm)
    host_params_h_Q_d = cp.array(Q)

    host_params_h_data_stop.record()
    host_params_h_data_stop.synchronize()
    host_params_h_elapsedTimeData = cp.cuda.get_elapsed_time(
        host_params_h_data_start, host_params_h_data_stop
    )

    if verbose_gpu >= 2:
        print("done!")

    if verbose_gpu >= 2:
        print(
            "Time to copy data from host to device = {0} ms".format(
                host_params_h_elapsedTimeData
            )
        )


def gpu_iterate(p, T, mole_fraction, Ia_arr, molarmass_arr, verbose_gpu):

    # ----------- setup global variables -----------------

    global host_params_h_start

    global init_params_h, iter_params_h
    global host_params_h_iso
    global host_params_h_v0_d
    global host_params_h_da_d
    global host_params_h_S0_d
    global host_params_h_El_d
    global host_params_h_Q_d
    global host_params_h_log_2gs_d
    global host_params_h_na_d
    global host_params_h_log_2vMm_d
    global host_params_h_stop
    global host_params_h_elapsedTime

    global host_params_h_DLM_d_in
    global host_params_h_spectrum_d_in

    global host_params_h_I_add

    global cuda_module
    global host_params_h_v0_dec
    global host_params_h_da_dec
    global DLM
    global I_ADD
    # ------------------------------------------------------

    host_params_h_start.record()
    # test comment
    ##    cdef int n_blocks
    set_pT(p, T, mole_fraction, iter_params_h)

    calc_gaussian_params(
        gaussian_param_data,
        init_params_h,
        iter_params_h,
    )

    calc_lorentzian_params(
        lorentzian_param_data,
        init_params_h,
        iter_params_h,
    )

    n_blocks = prepare_blocks(
        host_params_h_v0_dec,
        host_params_h_da_dec,
        host_params_h_dec_size,
        host_params_h_block_preparation_step_size,
        iter_params_h,
        init_params_h,
    )

    if verbose_gpu >= 2:
        print("Copying iteration parameters to device...")
    set_constant_memory("iter_params")

    if verbose_gpu >= 2:
        print("done!")

    # Zero DLM:
    host_params_h_DLM_d_in.fill(0)
    host_params_h_spectrum_d_in.fill(0)

    host_params_h_I_add.fill(0)

    if verbose_gpu >= 2:
        print("Getting ready...")

    host_params_h_start_DLM.record()

    fillDLM(
        (n_blocks,),
        (init_params_h.N_threads_per_block,),
        (
            host_params_h_iso_d,
            host_params_h_v0_d,
            host_params_h_da_d,
            host_params_h_S0_d,
            host_params_h_El_d,
            host_params_h_log_2gs_d,
            host_params_h_na_d,
            host_params_h_log_2vMm_d,
            host_params_h_DLM_d_in,
            host_params_h_Q_d,
            host_params_h_I_add,
        ),
    )

    cp.cuda.runtime.deviceSynchronize()

    # This makes the DLM array available in the calling module
    DLM = cp.asnumpy(host_params_h_DLM_d_in)
    I_ADD = cp.asnumpy(host_params_h_I_add)
    # DLM /= 0.8816117 #difference between current code and benchmark.
    host_params_h_stop_DLM.record()
    host_params_h_stop_DLM.synchronize()
    host_params_h_elapsedTimeDLM = cp.cuda.get_elapsed_time(
        host_params_h_start_DLM, host_params_h_stop_DLM
    )

    if verbose_gpu >= 2:
        print("<<<LAUNCHED>>> ")

    cp.cuda.runtime.deviceSynchronize()

    if verbose_gpu >= 2:
        # FFT
        print("Performing Fourier transform...", end=" ")
    host_params_h_DLM_d_out = cp.fft.rfft(host_params_h_DLM_d_in, axis=0)

    if verbose_gpu >= 2:
        print("done!")

    cp.cuda.runtime.deviceSynchronize()
    n_threads = 1024
    n_blocks = (init_params_h.N_v + 1) // n_threads + 1

    if verbose_gpu >= 2:
        print("Applying lineshapes...")

    applyLineshapes(
        (n_blocks,),
        (n_threads,),
        (
            host_params_h_DLM_d_out,
            host_params_h_spectrum_d_in,
        ),
    )

    cp.cuda.runtime.deviceSynchronize()

    if verbose_gpu >= 2:
        print("done!")

        # inverse FFT
        print("Performing inverse Fourier transform...")

    host_params_h_spectrum_d_out = cp.fft.irfft(host_params_h_spectrum_d_in)
    cp.cuda.runtime.deviceSynchronize()

    if verbose_gpu >= 2:
        print("done!")
    spectrum_h = (
        host_params_h_spectrum_d_out.get()[: init_params_h.N_v] / init_params_h.dv
    )

    host_params_h_stop.record()
    host_params_h_stop.synchronize()
    host_params_h_elapsedTime = cp.cuda.get_elapsed_time(
        host_params_h_start, host_params_h_stop
    )

    if verbose_gpu == 1:
        print("[rG = {0}%".format((np.exp(iter_params_h.log_dwG) - 1) * 100))
        print("rL = {0}%]".format((np.exp(iter_params_h.log_dwL) - 1) * 100))
        print("Runtime: {0}".format(host_params_h_elapsedTimeDLM))
        print(" + {0}".format(host_params_h_elapsedTime - host_params_h_elapsedTimeDLM))
        print(" = {0} ms".format(host_params_h_elapsedTime))
        print("Finished calculating spectrum!")

    return spectrum_h
