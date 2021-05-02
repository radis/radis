import ctypes

import cupy as cp
import numpy as np

from radis_cython_extensions import (
    calc_gaussian_params,
    calc_lorentzian_params,
    init_gaussian_params,
    init_lorentzian_params,
    prepare_blocks,
)

host_params_h_start = cp.cuda.Event()
host_params_h_stop = cp.cuda.Event()
host_params_h_start_DLM = cp.cuda.Event()
host_params_h_stop_DLM = cp.cuda.Event()
host_params_h_data_start = cp.cuda.Event()
host_params_h_data_stop = cp.cuda.Event()

cuda_code = r"""
#include<cupy/complex.cuh>
extern "C"{

struct initData {
	float v_min;
	float v_max;
	float dv;
	int N_v;
	int N_wG;
	int N_wL;
	int N_wG_x_N_wL;
	int N_total;
	int Max_lines;
	int N_lines;
	int N_points_per_block;
	int N_threads_per_block;
	int N_blocks_per_grid;
	int N_points_per_thread;
	int	Max_iterations_per_thread;
	int shared_size_floats;
};

struct blockData {
	int line_offset;
	int iv_offset;
};

struct iterData {
	float p;
	float log_p;
	float hlog_T;
	float log_rT;
	float c2T;
    float N;
	float log_wG_min;
	float log_wL_min;
	float log_dwG;
	float log_dwL;
	blockData blocks[4096];
};

__device__ __constant__ initData init_params_d;
__device__ __constant__ iterData iter_params_d;

__global__ void fillDLM(
    int* iso,
	float* v0,
	float* da,
	float* S0,
	float* El,
	float* log_2gs,
	float* na,
	float* log_2vMm,
	float* global_DLM,
    float* Q,               // Q is an array of size max(isotopes_id) + 1
    float* I_add_arr) {

	// Some overhead for "efficient" block allocation:
	blockData block = iter_params_d.blocks[blockIdx.x + gridDim.x * blockIdx.y];
	int block_id = blockIdx.x + gridDim.x * blockIdx.y;
	int N_iterations = (iter_params_d.blocks[block_id + 1].line_offset - iter_params_d.blocks[block_id].line_offset) / init_params_d.N_threads_per_block;
	int DLM_offset = iter_params_d.blocks[block_id].iv_offset * init_params_d.N_wG_x_N_wL;
	int iv_offset = iter_params_d.blocks[block_id].iv_offset;

	int NwL = init_params_d.N_wL;
	int NwGxNwL = init_params_d.N_wG_x_N_wL;

	////Allocate and zero the Shared memory
	//extern __shared__ float shared_DLM[];

	float* DLM = global_DLM;

	for (int n = 0; n < N_iterations; n++) { // eliminate for-loop

		// >>: Process from left to right edge:
		int i = iter_params_d.blocks[block_id].line_offset + threadIdx.x + n * blockDim.x;

		if (i < init_params_d.N_lines) {
			//Calc v
			float v_dat = v0[i] + iter_params_d.p * da[i];
			float iv = (v_dat - init_params_d.v_min) / init_params_d.dv; //- iv_offset;
			int iv0 = (int)iv;
			int iv1 = iv0 + 1  ;

            //arr_idx[i] = iv0;

			//^4

			if ((iv0 >= 0) && (iv1 < init_params_d.N_v)) {

				//Calc wG
				float log_wG_dat = log_2vMm[i] + iter_params_d.hlog_T;
				float iwG = (log_wG_dat - iter_params_d.log_wG_min) / iter_params_d.log_dwG;
				int iwG0 = (int)iwG;
				int iwG1 = iwG0 + 1;
				//^8

                //arr_idx[i] = iwG0;

				//Calc wL
				float log_wL_dat = log_2gs[i] + iter_params_d.log_p + na[i] * iter_params_d.log_rT;
				float iwL = (log_wL_dat - iter_params_d.log_wL_min) / iter_params_d.log_dwL;
				int iwL0 = (int)iwL;
				int iwL1 = iwL0 + 1;
				//^12

                //arr_idx[i] = iwL0;

				//Calc I
				float I_add = iter_params_d.N * S0[i] * (expf(iter_params_d.c2T * El[i]) - expf(iter_params_d.c2T * (El[i] + v0[i]))) / Q[iso[i]];

                I_add_arr[i] = I_add;

				float av = iv - iv0;
				float awG = (iwG - iwG0) * expf((iwG1 - iwG) * iter_params_d.log_dwG);
				float awL = (iwL - iwL0) * expf((iwL1 - iwL) * iter_params_d.log_dwL);

				float aV00 = (1 - awG) * (1 - awL);
				float aV01 = (1 - awG) * awL;
				float aV10 = awG * (1 - awL);
				float aV11 = awG * awL;

				float Iv0 = I_add * (1 - av);
				float Iv1 = I_add * av;

				atomicAdd(&DLM[iwL0 + iwG0 * NwL + iv0 * NwGxNwL], aV00 * Iv0);
				atomicAdd(&DLM[iwL0 + iwG0 * NwL + iv1 * NwGxNwL], aV00 * Iv1);
				atomicAdd(&DLM[iwL0 + iwG1 * NwL + iv0 * NwGxNwL], aV01 * Iv0);
				atomicAdd(&DLM[iwL0 + iwG1 * NwL + iv1 * NwGxNwL], aV01 * Iv1);
				atomicAdd(&DLM[iwL1 + iwG0 * NwL + iv0 * NwGxNwL], aV10 * Iv0);
				atomicAdd(&DLM[iwL1 + iwG0 * NwL + iv1 * NwGxNwL], aV10 * Iv1);
				atomicAdd(&DLM[iwL1 + iwG1 * NwL + iv0 * NwGxNwL], aV11 * Iv0);
				atomicAdd(&DLM[iwL1 + iwG1 * NwL + iv1 * NwGxNwL], aV11 * Iv1);
			}
		}
	}
}

__global__ void applyLineshapes(complex<float>* DLM, complex<float>* spectrum) {

	const float pi = 3.141592653589793f;
	const float r4log2 = 0.36067376022224085f; // = 1 / (4 * ln(2))
	int iv = threadIdx.x + blockDim.x * blockIdx.x;

	if (iv < init_params_d.N_v + 1) {

		float x = iv / (2 * init_params_d.N_v * init_params_d.dv);
		float mul = 0.0;
        complex<float> out_complex = 0;
        // float out_re = 0.0;
		// float out_im = 0.0;
		float wG, wL;
		int index;

		for (int iwG = 0; iwG < init_params_d.N_wG; iwG++) {
			wG = expf(iter_params_d.log_wG_min + iwG * iter_params_d.log_dwG);
			for (int iwL = 0; iwL < init_params_d.N_wL; iwL++) {
				//index = iwG + iwL * init_params_d.N_wG + iv * init_params_d.N_wG_x_N_wL;
                index = iv + iwG * (init_params_d.N_v+1) + iwL * init_params_d.N_wG * (init_params_d.N_v+1);
				wL = expf(iter_params_d.log_wL_min + iwL * iter_params_d.log_dwL);
				mul = expf(-r4log2 * powf(pi * x * wG, 2) - pi * x * wL);
                out_complex += mul* DLM[index];
                //out_complex += DLM[index];
			}
		}
        complex<float> temp(out_complex.real(), out_complex.imag());
		spectrum[iv] = temp;
	}
}
}"""


def set_pT(p, T, mole_fraction):

    # ----------- setup global variables -----------------
    global iter_params_h
    # ------------------------------------------------------

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


cuda_module = cp.RawModule(code=cuda_code)
fillDLM = cuda_module.get_function("fillDLM")
applyLineshapes = cuda_module.get_function("applyLineshapes")


def init(v_arr, N_wG, N_wL, iso, v0, da, log_2gs, na, log_2vMm, S0, El, Q, verbose_gpu):

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

    ##    cdef np.ndarray[dtype=np.int32_t, ndim=1] spec_h_iso = iso
    ##    cdef np.ndarray[dtype=np.float32_t, ndim=1] spec_h_v0 = v0
    ##    cdef np.ndarray[dtype=np.float32_t, ndim=1] spec_h_da = da
    ##    cdef np.ndarray[dtype=np.float32_t, ndim=1] spec_h_log_2gs = log_2gs
    ##    cdef np.ndarray[dtype=np.float32_t, ndim=1] spec_h_na = na
    ##    cdef np.ndarray[dtype=np.float32_t, ndim=1] spec_h_log_2vMm = log_2vMm
    ##    cdef np.ndarray[dtype=np.float32_t, ndim=1] spec_h_S0 = S0
    ##    cdef np.ndarray[dtype=np.float32_t, ndim=1] spec_h_El = El
    ##    cdef np.ndarray[dtype=np.float32_t, ndim=1] spec_h_Q = Q

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

    init_lorentzian_params(log_2gs, na, verbose_gpu)
    init_gaussian_params(log_2vMm, verbose_gpu)
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

    # CTYPES #2
    if verbose_gpu >= 2:
        print("Copying initialization parameters to device memory...")  # , end = " ")
    memptr_init_params_d = cuda_module.get_global("init_params_d")
    init_params_ptr = ctypes.cast(ctypes.pointer(init_params_h), ctypes.c_void_p)
    init_params_size = ctypes.sizeof(init_params_h)
    memptr_init_params_d.copy_from_host(init_params_ptr, init_params_size)

    if verbose_gpu >= 2:
        print("done!")
        print("Copying spectral data to device memory...")  # , end = " ")

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


def iterate(p, T, mole_fraction, Ia_arr, molarmass_arr, verbose_gpu):

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
    set_pT(p, T, mole_fraction)
    calc_gaussian_params()
    calc_lorentzian_params()
    n_blocks = prepare_blocks()

    if verbose_gpu >= 2:
        print("Copying iteration parameters to device...")  # , end = " ")

    ## CTYPES #1
    memptr_iter_params_d = cuda_module.get_global("iter_params_d")
    iter_params_ptr = ctypes.cast(ctypes.pointer(iter_params_h), ctypes.c_void_p)
    struct_size = ctypes.sizeof(iter_params_h)
    memptr_iter_params_d.copy_from_host(iter_params_ptr, struct_size)

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
        print("<<<LAUNCHED>>> ")  # , end = " ")

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
        print("Applying lineshapes...")  # , end = " ")

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
        print("Performing inverse Fourier transform...")  # , end = " ")

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
        print(
            "[rG = {0}%".format((np.exp(iter_params_h.log_dwG) - 1) * 100)
        )  # , end = " ")
        print(
            "rL = {0}%]".format((np.exp(iter_params_h.log_dwL) - 1) * 100)
        )  # , end = " ")
        print("Runtime: {0}".format(host_params_h_elapsedTimeDLM))  # , end = "")
        print(
            " + {0}".format(host_params_h_elapsedTime - host_params_h_elapsedTimeDLM)
        )  # , end = "")
        print(" = {0} ms".format(host_params_h_elapsedTime))
        print("Finished calculating spectrum!")

    return spectrum_h
