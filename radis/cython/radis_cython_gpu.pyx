# distutils: language=c++

import pickle
import ctypes

from libcpp.map cimport map as mapcpp
from libcpp.set cimport set
from libcpp.utility cimport pair
from libcpp.vector cimport vector

import cython
from cpython cimport array
from cython.operator import dereference

import cupy as cp
import cupyx
import cupyx.scipy.fftpack

import numpy as np
cimport numpy as np


cdef float epsilon   = <float> 0.0001
cdef float FLOAT_MAX = <float> 1e30
cdef float FLOAT_MIN = <float>-1e30


cdef  int host_params_h_block_preparation_step_size
cdef  int host_params_h_Min_threads_per_block
cdef  int host_params_h_Max_threads_per_block

cdef  float host_params_h_log_2vMm_min
cdef  float host_params_h_log_2vMm_max

cdef  vector[float] host_params_h_top_x
cdef  vector[float] host_params_h_top_a
cdef  vector[float] host_params_h_top_b
cdef  vector[float] host_params_h_bottom_x
cdef  vector[float] host_params_h_bottom_a
cdef  vector[float] host_params_h_bottom_b

cdef  size_t host_params_h_dec_size
cdef  int host_params_h_shared_size

host_params_h_start = cp.cuda.Event()
host_params_h_stop = cp.cuda.Event()
host_params_h_start_DLM = cp.cuda.Event()
host_params_h_stop_DLM = cp.cuda.Event()
host_params_h_data_start = cp.cuda.Event()
host_params_h_data_stop = cp.cuda.Event()

cdef float host_params_h_elapsedTime
cdef float host_params_h_elapsedTimeDLM
cdef float host_params_h_elapsedTimeData

host_params_h_iso_d = None
host_params_h_v0_d = None
host_params_h_da_d = None
host_params_h_S0_d = None
host_params_h_El_d = None
host_params_h_log_2gs_d = None
host_params_h_na_d = None
host_params_h_log_2vMm_d = None
host_params_h_DLM_d = None
host_params_h_spectrum_d = None
host_params_h_Q_d = None
host_params_h_I_add = None

# defined in 'iterate'
host_params_h_DLM_d_in = None
host_params_h_DLM_d_out = None

host_params_h_spectrum_d_in = None
host_params_h_spectrum_d_out = None

#set default path:
database_path = '/home/pankaj/radis-lab/data-2000-2400/'

#Default number of lines to load
#This is outside of the struct to allow it to be both None and int type;
#Setting N_lines_to_load to None loads all lines in the file.
N_lines_to_load = None


class initData(ctypes.Structure):
    _fields_= [
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

        ("shared_size_floats", ctypes.c_int)
    ]

class blockData(ctypes.Structure):
    _fields_=[
        ("line_offset", ctypes.c_int),
        ("iv_offset", ctypes.c_int)
    ]

class iterData(ctypes.Structure):
    _fields_=[
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

        ("blocks", blockData * 4096)
    ]


init_params_h = initData()
iter_params_h = iterData()

cuda_code = r'''
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
}'''

cuda_module = cp.RawModule(code=cuda_code)
fillDLM = cuda_module.get_function('fillDLM')
applyLineshapes = cuda_module.get_function('applyLineshapes')


####################################


cdef void set_pT(float p, float T, float mole_fraction):

    # ----------- setup global variables -----------------
    global iter_params_h
    #------------------------------------------------------

    cdef float c2 = <float>1.4387773538277204 #K.cm
    cdef float k = <float>1.38064852e-23 #J.K-1
    iter_params_h.p = p #bar
    iter_params_h.log_p = np.log(p)
    iter_params_h.hlog_T = 0.5 * np.log(T)
    iter_params_h.log_rT = np.log(296.0/T)
    iter_params_h.c2T = -c2/T
    iter_params_h.N = mole_fraction * p*1e5 / (1e6 * k * T) #cm-3

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

    #iter_params_h.Q = Q

def read_npy(fname, arr):
    print("Loading {0}...".format(fname))
    arr = np.load(fname)
    print("Done!")


def set_path(path):
    global database_path
    database_path = path

def set_N_lines(int N):
    global N_lines_to_load
    N_lines_to_load = N

# CUSTOM COMPARATOR to sort map keys in non increasing order
cdef extern from *:
    """
    struct greater {
        bool operator () (const float x, const float y) const {return x > y;}
    };
    """
    ctypedef struct greater:
        float a
        float b


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

cdef void init_lorentzian_params(np.ndarray[dtype=np.float32_t, ndim=1] log_2gs, np.ndarray[dtype=np.float32_t, ndim=1] na, verbose_gpu):

    # ----------- setup global variables -----------------
    global host_params_h_top_a
    global host_params_h_top_b
    global host_params_h_top_x
    global host_params_h_bottom_a
    global host_params_h_bottom_b
    global host_params_h_bottom_x
    #------------------------------------------------------

    cdef set[pair[float,float]] unique_set
    cdef float float_pair[2]


    cdef vector[pair[float,float]] duplicates_removed
    cdef vector[float] na_short
    cdef vector[float] log_2gs_short
    cdef mapcpp[float, float, greater] bottom_envelope_map
    cdef mapcpp[float, float] top_envelope_map

    cdef vector[float] top_a
    cdef vector[float] top_b
    cdef vector[float] top_x

    cdef vector[float] bottom_a
    cdef vector[float] bottom_b
    cdef vector[float] bottom_x

    if verbose_gpu >= 2:
        print("Initializing Lorentzian parameters ", end = "")

    cdef size_t top_size = 0
    cdef size_t bottom_size = 0

    fname = "Lorenzian_minmax_" + str(len(log_2gs)) + ".dat"

    cdef unsigned int i,na_len
    cdef float na_i, log_2gs_i

    try:
        with open(fname, 'rb') as f:

            if verbose_gpu >= 2:
                print(" (from cache)... ", end=" ")

            lt = pickle.load(f)

            top_size = lt[0]
            host_params_h_top_a.resize(top_size)
            host_params_h_top_b.resize(top_size)
            host_params_h_top_x.resize(top_size)

            # now read top_size bits 3 times to fill the above 3 vectors
            host_params_h_top_a = lt[1]
            host_params_h_top_b = lt[2]
            host_params_h_top_x = lt[3]

            bottom_size = lt[4]
            host_params_h_bottom_a.resize(bottom_size)
            host_params_h_bottom_b.resize(bottom_size)
            host_params_h_bottom_x.resize(bottom_size)

            host_params_h_bottom_a = lt[5]
            host_params_h_bottom_b = lt[6]
            host_params_h_bottom_x = lt[7]

    except:
        if verbose_gpu >= 2:
            print(" ... ", end = " ")

        na_len = na.size
        for i in range(na_len):
            # Somewhat of a hack; all of the structs I tried show Python interaction
            # when storing values. float[2] didn't have this problem, so I fill the
            # float[2] array, and then convert the float* pointer to a const
            # pair[float,float]* pointer, which is then dereferenced by the [0].
            # well it gets the job done I suppose, there are no yellow lines inside
            # of this loop anymore.

            float_pair[0] = na[i]
            float_pair[1] = log_2gs[i]
            unique_set.insert((<const pair[float,float]*>float_pair)[0])

        duplicates_removed.assign(unique_set.begin(), unique_set.end())


        for na_i, log_2gs_i in duplicates_removed:
            na_short.push_back(na_i)
            log_2gs_short.push_back(log_2gs_i)

        for i in range(len(na_short)):
            na_i = na_short[i]
            log_2gs_i = log_2gs_short[i]

            if bottom_envelope_map.count(na_i):
                if log_2gs_i < bottom_envelope_map.at(na_i):
                    bottom_envelope_map[na_i] = log_2gs_i
            else:
                bottom_envelope_map.insert({na_i, log_2gs_i})

            if top_envelope_map.count(na_i):
                if log_2gs_i > top_envelope_map.at(na_i):
                    top_envelope_map[na_i] = log_2gs_i
            else:
                top_envelope_map.insert({na_i, log_2gs_i})

        top_a = { dereference(top_envelope_map.begin()).first }
        top_b = { dereference(top_envelope_map.begin()).second }
        top_x = { FLOAT_MIN }

        idx = 0
        for first_el, second_el in top_envelope_map:
            if idx != 0:
                for i in range(len(top_x)):
                    x_ij = (second_el - top_b[i]) / (top_a[i] - first_el)
                    if x_ij >= top_x[i]:
                        if i < top_x.size() - 1:
                            if x_ij < top_x[i+1]:
                                break;
                        else:
                            break

                top_a.resize(i+1)
                top_b.resize(i+1)
                top_x.resize(i+1)

                top_a.push_back(first_el)
                top_b.push_back(second_el)
                top_x.push_back(x_ij)

            idx+=1

        top_x.erase(top_x.begin())
        top_x.push_back(FLOAT_MAX)

        host_params_h_top_a = top_a
        host_params_h_top_b = top_b
        host_params_h_top_x = top_x
        top_size = top_x.size()

        bottom_a = { dereference(bottom_envelope_map.begin()).first }
        bottom_b = { dereference(bottom_envelope_map.begin()).second }
        bottom_x = { FLOAT_MIN }

        idx = 0

        for first_el, second_el in bottom_envelope_map:
            if idx != 0:
                for i in range(len(bottom_x)):
                    x_ij = (second_el - bottom_b[i]) / (bottom_a[i] - first_el)
                    if x_ij >= bottom_x[i]:
                        if i < bottom_x.size() - 1:
                            if x_ij < bottom_x[i+1]:
                                break
                        else:
                            break

                bottom_a.resize(i + 1)
                bottom_b.resize(i + 1)
                bottom_x.resize(i + 1)

                bottom_a.push_back(first_el)
                bottom_b.push_back(second_el)
                bottom_x.push_back(x_ij)

            idx+=1

        bottom_x.erase(bottom_x.begin())
        bottom_x.push_back(FLOAT_MAX)

        host_params_h_bottom_a = bottom_a
        host_params_h_bottom_b = bottom_b
        host_params_h_bottom_x = bottom_x
        bottom_size = bottom_x.size()

        lt = [top_size,
            host_params_h_top_a,
            host_params_h_top_b,
            host_params_h_top_x,
            bottom_size,
            host_params_h_bottom_a,
            host_params_h_bottom_b,
            host_params_h_bottom_x]

        with open(fname, 'wb') as f:
            pickle.dump(lt, f)

    if verbose_gpu >= 2:
        print("done!")
    return

cdef void calc_lorentzian_params():

    # ----------- setup global variables -----------------
    global host_params_h_top_x
    global host_params_h_bottom_x
    global host_params_h_top_a
    global host_params_h_bottom_a
    global host_params_h_top_b
    global host_params_h_bottom_b
    global iter_params_h
    global epsilon
    #------------------------------------------------------

    cdef float log_wL_min
    cdef float log_wL_max

    for i in range(host_params_h_bottom_x.size()):
        if iter_params_h.log_rT < host_params_h_bottom_x[i]:
            log_wL_min = iter_params_h.log_rT * host_params_h_bottom_a[i] + host_params_h_bottom_b[i]  + iter_params_h.log_p
            break

    for i in range(host_params_h_top_x.size()):
        if iter_params_h.log_rT < host_params_h_top_x[i]:
            log_wL_max = iter_params_h.log_rT * host_params_h_top_a[i] + host_params_h_top_b[i]  + iter_params_h.log_p + epsilon
            break

    cdef float log_dwL = (log_wL_max - log_wL_min) / (init_params_h.N_wL - 1)

    iter_params_h.log_wL_min = log_wL_min
    iter_params_h.log_dwL = log_dwL
    return


cdef void init_gaussian_params(np.ndarray[dtype=np.float32_t, ndim=1] log_2vMm, verbose_gpu):

    # ----------- setup global variables -----------------
    global host_params_h_log_2vMm_min
    global host_params_h_log_2vMm_max
    #------------------------------------------------------

    cdef float log_2vMm_min
    cdef float log_2vMm_max
    if verbose_gpu >= 2:
        print("Initializing Gaussian parameters", end="")

    fname = "Gaussian_minmax_" + str(len(log_2vMm)) + ".dat"
    try:
        lt = pickle.load(open(fname, "rb"))
        if verbose_gpu >= 2:
            print(" (from cache)... ", end=" ")
        log_2vMm_min = lt[0]
        log_2vMm_max = lt[1]
    except (OSError, IOError) as e:
        if verbose_gpu >= 2:
            print("... ", end=" ")
        log_2vMm_min = np.amin(log_2vMm)
        log_2vMm_max = np.amax(log_2vMm)
        lt = [log_2vMm_min, log_2vMm_max]
        pickle.dump(lt, open(fname, "wb"))

    host_params_h_log_2vMm_min = log_2vMm_min
    host_params_h_log_2vMm_max = log_2vMm_max

    if verbose_gpu >= 2:
        print("done!")

    return


cdef void calc_gaussian_params():

    # ----------- setup global variables -----------------
    global host_params_h_log_2vMm_min
    global host_params_h_log_2vMm_max
    global init_params_h, iter_params_h
    global epsilon
    #------------------------------------------------------

    cdef float log_wG_min = host_params_h_log_2vMm_min + iter_params_h.hlog_T
    cdef float log_wG_max = host_params_h_log_2vMm_max + iter_params_h.hlog_T + epsilon
    cdef float log_dwG = (log_wG_max - log_wG_min) / (init_params_h.N_wG - 1)

    iter_params_h.log_wG_min = log_wG_min
    iter_params_h.log_dwG = log_dwG

    return

cdef int prepare_blocks():

    # ----------- setup global variables -----------------
    global host_params_h_v0_dec
    global host_params_h_da_dec
    global host_params_h_dec_size
    global host_params_h_block_preparation_step_size

    global iter_params_h, init_params_h
    #------------------------------------------------------

    cdef np.ndarray[dtype=np.float32_t, ndim=1] v0 = host_params_h_v0_dec
    cdef np.ndarray[dtype=np.float32_t, ndim=1] da = host_params_h_da_dec


    cdef float v_prev
    cdef float dvdi
    cdef int i = 0
    cdef int n = 0
    cdef int step = host_params_h_block_preparation_step_size

    new_block = blockData()

    cdef float v_cur = v0[0] + iter_params_h.p * da[0]
    cdef float v_max = v_cur + init_params_h.N_points_per_block * init_params_h.dv
    cdef int i_max = init_params_h.Max_iterations_per_thread

    new_block.line_offset = 0
    new_block.iv_offset = int(((v_cur - init_params_h.v_min) / init_params_h.dv))
    while True:
        i += step
        if i > host_params_h_dec_size:
            iter_params_h.blocks[n] = new_block

            n+=1
            new_block.line_offset = i * init_params_h.N_threads_per_block

            iter_params_h.blocks[n] = new_block
            break

        v_prev = v_cur
        v_cur = v0[i] + iter_params_h.p * da[i]
        if ((v_cur > v_max) or (i >= i_max)) :
            # if (v_cur > v_max) :
            #     dvdi = (v_cur - v_prev) / <float>step
            #     i -= int(((v_cur - v_max) / dvdi)) + 1
            #     v_cur = v0[i] + iter_params_h.p * da[i]
            iter_params_h.blocks[n] = new_block
            n+=1
            new_block.iv_offset = int(((v_cur - init_params_h.v_min) / init_params_h.dv))
            new_block.line_offset = i * init_params_h.N_threads_per_block
            v_max = v_cur + (init_params_h.N_points_per_block) * init_params_h.dv
            i_max = i + init_params_h.Max_iterations_per_thread

    return n


