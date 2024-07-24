#cython: language_level=3
import cython
import numpy as np

cimport numpy as np

#This can be modified to see if the compiled/built/installed version is current:
__version__ = 0.021

# Fast compiled version of add_at():
ctypedef fused float_t:
    np.float32_t
    np.float64_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def add_at(float_t[:,:,::1] arr,
           np.int32_t[::1] k_arr,
           np.int32_t[::1] l_arr,
           np.int32_t[::1] m_arr,
           float_t[::1] values):

    cdef unsigned int N = values.size
    cdef unsigned int i

    for i in range(N):
        arr[k_arr[i],l_arr[i],m_arr[i]] += values[i]

    return arr
