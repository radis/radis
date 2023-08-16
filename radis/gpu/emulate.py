import os.path
from ctypes import (
    POINTER,
    byref,
    c_char,
    c_int,
    c_longlong,
    c_short,
    c_void_p,
    cast,
    cdll,
    memmove,
    sizeof,
    windll,
)
from os import name as os_name

import numpy as np
from scipy.fft import irfft, rfft

from radis.gpu.structs import blockDim_t, gridDim_t
from radis.misc.utils import getProjectRoot

# from os.path import dirname


CUFFT_R2C = 0x2A
CUFFT_C2R = 0x2C

# This number does not have a meaningful interpretation in a CPU,
# so we just default to a typical number for GPU.
_max_threads_per_block = 1024


class CuContext:
    def __init__(self, device_id=0, flags=0):

        # private:
        self._context = c_void_p(123)
        self._device = c_void_p(456)

    @staticmethod
    def getDeviceList():
        import platform

        return [platform.processor()]

    def printDeviceCapabilities(self):
        print("> GPU emulated by CPU")
        print("  Total amount of global memory:   0 bytes")

    def getMaxThreadsPerBlock(self):
        return _max_threads_per_block

    def setCurrent(self):
        pass

    def synchronize(self):
        pass

    def destroy(self):
        pass


class CuModule:
    def __init__(self, context, module_name):
        # public:
        lib_ext = ".dll" if os_name == "nt" else ".so"
        self.module_name = os.path.splitext(module_name)[0] + lib_ext
        self.context = context

        # private:
        # radis_path = dirname(dirname(__file__))
        radis_path = getProjectRoot()

        lib_obj = windll if os_name == "nt" else cdll
        self._module = lib_obj.LoadLibrary(
            os.path.join(radis_path, "gpu", self.module_name)
        )
        self._func_dict = {}
        self._global_dict = {}

    def __getattr__(self, attr):
        try:
            self._func_dict[attr]
            return self._func_dict[attr]

        except (KeyError):
            _function = getattr(self._module, attr)
            self._func_dict[attr] = CuFunction(_function)
            self._func_dict[attr].module = self
            return self._func_dict[attr]

    def setConstant(self, name, c_val):
        try:
            _var, _size = self._global_dict[name]

        except (KeyError):
            _var = type(c_val).in_dll(self._module, name)
            _size = sizeof(c_val)
            self._global_dict[name] = (_var, _size)

        memmove(byref(_var), byref(c_val), _size)

    def getConstant(self, name):
        try:
            _var, _size = self._global_dict[name]

        except (KeyError):
            pass  # TODO: catch exception

        return _var


class CuFunction:
    def __init__(self, _function):

        # public:
        self.module = None
        self.args = []
        self.blocks = (1, 1, 1)
        self.threads = (1, 1, 1)
        self.sync = False

        # private:
        self._function = _function

    def setGrid(self, blocks=(1, 1, 1), threads=(1, 1, 1)):
        self.blocks = blocks
        self.threads = threads

    def setArgs(self, *vargs):
        self.args = vargs

    def __call__(self, *vargs, blocks=None, threads=None, sync=None):

        self.args = self.args if not len(vargs) else vargs
        self.blocks = self.blocks if blocks is None else blocks
        self.threads = self.threads if threads is None else threads
        self.sync = self.sync if sync is None else sync

        self.module.setConstant("blockDim", blockDim_t(*self.threads))
        self.module.setConstant("gridDim", gridDim_t(*self.blocks))

        c_args = [arr._ptr for arr in self.args]
        self._function(*c_args)


class CuArray:
    def __init__(self, shape, dtype=np.float32, init="empty"):
        self._ptr = c_void_p()
        self.resize(shape, dtype, init)

    def resize(self, shape=None, dtype=None, init="empty"):
        shape_tuple = shape if isinstance(shape, tuple) else (shape,)
        self.shape = self.shape if shape is None else shape_tuple
        self.dtype = self.dtype if dtype is None else np.dtype(dtype)
        self.size = int(np.prod(shape))
        self.itemsize = self.dtype.itemsize
        self.nbytes = self.itemsize * self.size

        if init not in ("zeros", "empty"):
            return

        elif init == "empty":
            self._arr = np.empty(self.shape, dtype=self.dtype)

        elif init == "zeros":
            self._arr = np.zeros(self.shape, dtype=self.dtype)

        self._ptr = c_void_p(self._arr.ctypes.data)

    @staticmethod
    def fromArray(arr):

        obj = CuArray(arr.shape, arr.dtype, init=None)
        obj.setArray(arr)

        return obj

    def zeroFill(self):
        self.getArray()[:] = 0

    def setArray(self, arr):

        params_changed = arr.shape != self.shape or arr.dtype != self.dtype
        uninitialized_memory = self._ptr.value is None
        if params_changed or uninitialized_memory:
            self.resize(arr.shape, arr.dtype, None)
        self._ptr = c_void_p(arr.ctypes.data)
        self._arr = arr

    def getArray(self):
        my_ctype = {1: c_char, 2: c_short, 4: c_int, 8: c_longlong}[self.itemsize]
        my_cptr = cast(self._ptr.value, POINTER(my_ctype))
        arr = np.ctypeslib.as_array(my_cptr, self.shape).view(self.dtype)
        return arr


class CuFFT:
    def __init__(self, arr_in, arr_out, direction="fwd", plan_fft=True):

        # public:
        self.arr_in = arr_in
        self.arr_out = arr_out

        # private:
        self._direction = direction
        self._arr = arr_in if direction == "fwd" else arr_out
        self._fft_type = CUFFT_R2C if direction == "fwd" else CUFFT_C2R
        self._plan = c_void_p(789)

    def planMany(self):
        pass

    def __call__(self):

        np_arr_in = self.arr_in.getArray()
        np_arr_out = self.arr_out.getArray()

        if self._direction == "fwd":
            np_arr_out[:] = rfft(np_arr_in, axis=0, norm="backward")
        else:
            n = np_arr_out.shape[0]
            np_arr_out[:] = irfft(np_arr_in, n=n, axis=0, norm="forward")

    def destroy(self):
        pass
