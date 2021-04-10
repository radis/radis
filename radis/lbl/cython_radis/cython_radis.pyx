import numpy as np
cimport numpy as np
import cython

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