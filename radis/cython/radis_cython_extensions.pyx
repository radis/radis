#cython: language_level=3
import numpy as np
cimport numpy as np

import cython
from cpython cimport array
from cython.operator import dereference

import ctypes
from libcpp cimport bool
from libcpp.map cimport map as mapcpp
from libcpp.set cimport set
from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libc.string cimport memcpy
from libc.math cimport log
#from libcpp.functional cimport greater

# std::greater is currently not implemented in cython (but will be in the future)
# for now we have to import it in this way:
#cdef extern from "<functional>" namespace "std" nogil:
#    # Comparisons
#    cdef cppclass greater[T=*]:
#        greater() except +
#        bool operator()(const T& lhs, const T& rhs) except +

cdef extern from *:
    """
    struct greatercpp {
        bool operator () (const float x, const float y) const {return x > y;}
    };
    """
    ctypedef struct greatercpp:
        float a
        float b


#This can be modified to see if the compiled/built/installed version is current:
__version__ = 0.021

# Fast compiled version of add_at():
def add_at(arr, k,l,m, I):
    if arr.dtype != I.dtype:
        I = I.astype(arr.dtype)
    if arr.dtype == np.float32:
        return add_at_32(arr, k, l, m, I)
    else:
        return add_at_64(arr, k, l, m, I)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def add_at_64(np.ndarray[np.float64_t, ndim=3] arr,
              np.ndarray[np.int32_t, ndim=1] k_arr,
              np.ndarray[np.int32_t, ndim=1] l_arr,
              np.ndarray[np.int32_t, ndim=1] m_arr,
              np.ndarray[np.float64_t, ndim=1] values):

    cdef unsigned int N = values.size
    cdef unsigned int i

    for i in range(N):
        arr[k_arr[i],l_arr[i],m_arr[i]] += values[i]

    return arr


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def add_at_32(np.ndarray[np.float32_t, ndim=3] arr,
              np.ndarray[np.int32_t, ndim=1] k_arr,
              np.ndarray[np.int32_t, ndim=1] l_arr,
              np.ndarray[np.int32_t, ndim=1] m_arr,
              np.ndarray[np.float32_t, ndim=1] values):

    cdef unsigned int N = values.size
    cdef unsigned int i

    for i in range(N):
        arr[k_arr[i],l_arr[i],m_arr[i]] += values[i]

    return arr



# Fast compiled functions for determining the min/max broadening envelope.


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def calc_lorentzian_envelope_params(
    np.ndarray[dtype=np.float32_t, ndim=1] na,
    np.ndarray[dtype=np.float32_t, ndim=1] gamma,
    verbose_gpu = False):

    cdef set[pair[float,float]] unique_set
    cdef float float_pair[2]

    cdef vector[pair[float,float]] duplicates_removed
    cdef vector[float] na_short
    cdef vector[float] gamma_short
    cdef mapcpp[float, float, greatercpp] bottom_envelope_map
    cdef mapcpp[float, float] top_envelope_map


    cdef float FLOAT_MAX = <float> 1e30
    cdef float FLOAT_MIN = <float>-1e30

    if verbose_gpu >= 2:
        print("Initializing Lorentzian parameters ")

    cdef size_t top_size = 0
    cdef size_t bottom_size = 0

    fname = "Lorenzian_minmax_" + str(len(gamma)) + ".dat"

    cdef unsigned int i,na_len
    cdef float na_i, gamma_i

    na_len = na.size
    for i in range(na_len):
        # Somewhat of a hack; all of the structs I tried show Python interaction
        # when storing values. float[2] didn't have this problem, so I fill the
        # float[2] array, and then convert the float* pointer to a const
        # pair[float,float]* pointer, which is then dereferenced by the [0].
        # well it gets the job done I suppose, there are no yellow lines inside
        # of this loop anymore.

        float_pair[0] = na[i]
        float_pair[1] = gamma[i]
        unique_set.insert((<const pair[float,float]*>float_pair)[0])

    duplicates_removed.assign(unique_set.begin(), unique_set.end())

    # make two new vectors where all duplicates are removed:
    for na_i, gamma_i in duplicates_removed:
        na_short.push_back(na_i)
        gamma_short.push_back(gamma_i)

    # identify candidates that might be part of the envelope:
    for i in range(len(na_short)):
        na_i = na_short[i]
        gamma_i = gamma_short[i]

        if bottom_envelope_map.count(na_i):
            if gamma_i < bottom_envelope_map.at(na_i):
                bottom_envelope_map[na_i] = gamma_i
        else:
            bottom_envelope_map.insert({na_i, gamma_i})

        if top_envelope_map.count(na_i):
            if gamma_i > top_envelope_map.at(na_i):
                top_envelope_map[na_i] = gamma_i
        else:
            top_envelope_map.insert({na_i, gamma_i})

    # For all candidates check which ones are actually part of the envelope:
    # First for the top:
    top_a = [dereference(top_envelope_map.begin()).first]
    top_b = [log(dereference(top_envelope_map.begin()).second)]
    top_x = [FLOAT_MIN]

    idx = 0
    for first_el, second_el in top_envelope_map:
        if idx != 0:
            for i in range(len(top_x)):
                x_ij = (log(second_el) - top_b[i]) / (top_a[i] - first_el)
                if x_ij >= top_x[i]:
                    if i < len(top_x) - 1:
                        if x_ij < top_x[i+1]:
                            break;
                    else:
                        break

            top_a = top_a[:i+1] + [first_el]
            top_b = top_b[:i+1] + [log(second_el)]
            top_x = top_x[:i+1] + [x_ij]

        idx+=1

    top_x = top_x[1:] + [FLOAT_MAX]

    #Then for the bottom:
    bottom_a = [dereference(bottom_envelope_map.begin()).first]
    bottom_b = [log(dereference(bottom_envelope_map.begin()).second)]
    bottom_x = [FLOAT_MIN]

    idx = 0
    for first_el, second_el in bottom_envelope_map:
        if idx != 0:
            for i in range(len(bottom_x)):
                x_ij = (log(second_el) - bottom_b[i]) / (bottom_a[i] - first_el)
                if x_ij >= bottom_x[i]:
                    if i < len(bottom_x) - 1:
                        if x_ij < bottom_x[i+1]:
                            break
                    else:
                        break

            bottom_a = bottom_a[:i+1] + [first_el]
            bottom_b = bottom_b[:i+1] + [log(second_el)]
            bottom_x = bottom_x[:i+1] + [x_ij]

        idx+=1

    bottom_x = bottom_x[1:] + [FLOAT_MAX]

    return ((bottom_a, bottom_b, bottom_x),
            (top_a, top_b, top_x))


def calc_gaussian_envelope_params(
    np.ndarray[dtype=np.float32_t, ndim=1] log_2vMm,
    verbose_gpu=False):

    log_2vMm_min = np.amin(log_2vMm)
    log_2vMm_max = np.amax(log_2vMm)

    return log_2vMm_min, log_2vMm_max

# Below are CPU functions extracted from the GPU module

cimport cpu_gpu_agnostic as cga
#from radis.lbl.gpu import initData, iterData

cpdef void fillLDM(gridDim, blockDim, args):

   cdef np.ndarray[np.uint8_t, ndim=1] iso
   cdef np.ndarray[np.float32_t, ndim=1] v0
   cdef np.ndarray[np.float32_t, ndim=1] da
   cdef np.ndarray[np.float32_t, ndim=1] S0
   cdef np.ndarray[np.float32_t, ndim=1] El
   cdef np.ndarray[np.float32_t, ndim=1] gamma
   cdef np.ndarray[np.float32_t, ndim=1] na
   cdef np.ndarray[np.float32_t, ndim=3] S_klm

   iso, v0, da, S0, El, gamma, na, S_klm = args

   cga.set_dims(blockDim[0],gridDim[0])
   cga.fillLDM(&iso[0], &v0[0], &da[0], &S0[0], &El[0], &gamma[0], &na[0], &S_klm[0,0,0])


cpdef void applyLineshapes(gridDim, blockDim, args):

    cdef np.ndarray[np.complex64_t, ndim=3] S_klm_FT
    cdef np.ndarray[np.complex64_t, ndim=1] abscoeff
    S_klm_FT, abscoeff = args

    cga.set_dims(blockDim[0], gridDim[0])
    cga.applyLineshapes(&S_klm_FT[0,0,0], &abscoeff[0])


cpdef void calcTransmittanceNoslit(gridDim, blockDim, args):

    cdef np.ndarray[np.float32_t, ndim=1] abscoeff,
    cdef np.ndarray[np.float32_t, ndim=1] transmittance_noslit
    abscoeff, transmittance_noslit = args

    cga.set_dims(blockDim[0], gridDim[0])
    cga.calcTransmittanceNoslit(&abscoeff[0], &transmittance_noslit[0])


cpdef void applyGaussianSlit(gridDim, blockDim, args):

    cdef np.ndarray[np.complex64_t, ndim=1] transmittance_noslit_FT,
    cdef np.ndarray[np.complex64_t, ndim=1] transmittance_FT
    transmittance_noslit_FT, transmittance_FT = args

    cga.set_dims(blockDim[0],gridDim[0])
    cga.applyGaussianSlit(&transmittance_noslit_FT[0], &transmittance_FT[0])


cpdef set_init_params(params):
    cdef void* cpp_ptr = cga.get_init_ptr()
    cdef size_t size = ctypes.sizeof(params)
    cdef char[:] struct_memview = bytearray(params)
    memcpy(cpp_ptr, &struct_memview[0], size)


cpdef set_iter_params(params):
    cdef void* cpp_ptr = cga.get_iter_ptr()
    cdef size_t size = ctypes.sizeof(params)
    cdef char[:] struct_memview = bytearray(params)
    memcpy(cpp_ptr, &struct_memview[0], size)

cpdef get_T():
    return cga.get_T()
