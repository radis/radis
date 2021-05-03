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

#TO-DO: should this be in gpu.py or in radis_cython_gpu.pyx?
class blockData(ctypes.Structure):
    _fields_ = [("line_offset", ctypes.c_int), ("iv_offset", ctypes.c_int)]


cdef float epsilon   = <float> 0.0001
cdef float FLOAT_MAX = <float> 1e30
cdef float FLOAT_MIN = <float>-1e30


cdef  int host_params_h_block_preparation_step_size
cdef  int host_params_h_Min_threads_per_block
cdef  int host_params_h_Max_threads_per_block

cdef  float host_params_h_log_2vMm_min
cdef  float host_params_h_log_2vMm_max

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





cdef  vector[float] host_params_h_top_x
cdef  vector[float] host_params_h_top_a
cdef  vector[float] host_params_h_top_b
cdef  vector[float] host_params_h_bottom_x
cdef  vector[float] host_params_h_bottom_a
cdef  vector[float] host_params_h_bottom_b

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)




def calc_lorentzian_envelope_params(
    np.ndarray[dtype=np.float32_t, ndim=1] log_2gs,
    np.ndarray[dtype=np.float32_t, ndim=1] na,
    verbose_gpu):

    cdef set[pair[float,float]] unique_set
    cdef float float_pair[2]

    cdef vector[pair[float,float]] duplicates_removed
    cdef vector[float] na_short
    cdef vector[float] log_2gs_short
    cdef mapcpp[float, float, greater] bottom_envelope_map
    cdef mapcpp[float, float] top_envelope_map

    if verbose_gpu >= 2:
        print("Initializing Lorentzian parameters ")

    cdef size_t top_size = 0
    cdef size_t bottom_size = 0

    fname = "Lorenzian_minmax_" + str(len(log_2gs)) + ".dat"

    cdef unsigned int i,na_len
    cdef float na_i, log_2gs_i

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

    # make two new vectors where all duplicates are removed:
    for na_i, log_2gs_i in duplicates_removed:
        na_short.push_back(na_i)
        log_2gs_short.push_back(log_2gs_i)

    # identify candidates that might be part of the envelope:
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

    # For all candidates check which ones are actually part of the envelope:
    # First for the top:
    top_a = [dereference(top_envelope_map.begin()).first]
    top_b = [dereference(top_envelope_map.begin()).second]
    top_x = [FLOAT_MIN]

    idx = 0
    for first_el, second_el in top_envelope_map:
        if idx != 0:
            for i in range(len(top_x)):
                x_ij = (second_el - top_b[i]) / (top_a[i] - first_el)
                if x_ij >= top_x[i]:
                    if i < len(top_x) - 1:
                        if x_ij < top_x[i+1]:
                            break;
                    else:
                        break

            top_a.append(first_el)
            top_b.append(second_el)
            top_x.append(x_ij)

        idx+=1

    top_x = top_x[1:] + [FLOAT_MAX]

    #Then for the bottom:
    bottom_a = [dereference(bottom_envelope_map.begin()).first]
    bottom_b = [dereference(bottom_envelope_map.begin()).second]
    bottom_x = [FLOAT_MIN]

    idx = 0
    for first_el, second_el in bottom_envelope_map:
        if idx != 0:
            for i in range(len(bottom_x)):
                x_ij = (second_el - bottom_b[i]) / (bottom_a[i] - first_el)
                if x_ij >= bottom_x[i]:
                    if i < len(bottom_x) - 1:
                        if x_ij < bottom_x[i+1]:
                            break
                    else:
                        break

            bottom_a.append(first_el)
            bottom_b.append(second_el)
            bottom_x.append(x_ij)

        idx+=1

    bottom_x = bottom_x[1:] + [FLOAT_MAX]

    return (np.array(top_a),
            np.array(top_b),
            np.array(top_x),
            np.array(bottom_a),
            np.array(bottom_b),
            np.array(bottom_x))


def calc_gaussian_envelope_params(
    np.ndarray[dtype=np.float32_t, ndim=1] log_2vMm,
    verbose_gpu):

    log_2vMm_min = np.amin(log_2vMm)
    log_2vMm_max = np.amax(log_2vMm)

    return log_2vMm_min, log_2vMm_max


def prepare_blocks(
    host_params_h_v0_dec,
    host_params_h_da_dec,
    host_params_h_dec_size,
    host_params_h_block_preparation_step_size,
    iter_params_h,
    init_params_h):

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


