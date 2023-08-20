# import sys
from ctypes import (
    byref,
    c_char_p,
    c_int,
    c_long,
    c_longlong,
    c_float,
    c_size_t,
    c_void_p,
    cast,
    cdll,
    create_string_buffer,
    windll,
)
from os import name as os_name

import numpy as np
from nvidia.cufft import __path__ as cufft_path

verbose = False


def cu_print(*vargs):
    global verbose
    if verbose:
        print(*vargs)


# TODO: establish cuda version requirement


CUDA_SUCCESS = 0
CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 1
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76

CUFFT_R2C = 0x2A
CUFFT_C2R = 0x2C

lib = None
lib_cufft = None
LoadLibrary = windll.LoadLibrary if os_name == 'nt' else cdll.LoadLibrary


class CuContext:
    def __init__(self, device, context):

        # private:
        self._device = device
        self._context = context
        

    @staticmethod
    def Open(device_id=0, flags=0):
        global lib
        
        # private:
        _context = c_void_p(0)
        _device = c_void_p(0)

        cuda_name = 'nvcuda.dll' if os_name == 'nt' else 'libcuda.so'
        try:
            lib = LoadLibrary(cuda_name)
        except(FileNotFoundError):
            print("Can't find {:s}...".format(cuda_name))
            return None
            
        err = lib.cuInit(0)

        deviceCount = c_int(0)
        lib.cuDeviceGetCount(byref(deviceCount))

        if deviceCount == 0:
            print("Error: no devices supporting CUDA\n")
            return None
        
        lib.cuDeviceGet(byref(_device), device_id)

        err = lib.cuCtxCreate_v2(byref(_context), flags, _device)
        if err != CUDA_SUCCESS:
            print("Error initializing the CUDA context.")
            lib.cuCtxDestroy_v2(self._context)
            return None

        return CuContext(_device, _context)


    @staticmethod
    def getDeviceList():
        dev_list = []
        _deviceCount = c_int(0)

        lib.cuInit(0)
        lib.cuDeviceGetCount(byref(_deviceCount))

        for i in range(_deviceCount.value):

            _device = c_long(0)
            lib.cuDeviceGet(byref(_device), i)

            _name = create_string_buffer(100)
            lib.cuDeviceGetName(_name, 100, _device)
            dev_list.append(_name.value.decode())

        return dev_list

    def printDeviceCapabilities(self):

        _major = c_long(0)
        _minor = c_long(0)
        lib.cuDeviceGetAttribute(
            byref(_major), CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, self._device
        )
        lib.cuDeviceGetAttribute(
            byref(_minor), CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, self._device
        )
        print(
            "> GPU Device has SM {:d}.{:d} compute capability".format(
                _major.value, _minor.value
            )
        )

        _totalGlobalMem = c_size_t(0)
        lib.cuDeviceTotalMem_v2(byref(_totalGlobalMem), self._device)
        print(
            "  Total amount of global memory:   {:d} bytes".format(
                _totalGlobalMem.value
            )
        )
        print(
            "  64-bit Memory Address:           {:s}".format(
                "YES" if _totalGlobalMem.value > (2 << 31) else "NO"
            )
        )

    def getMaxThreadsPerBlock(self):
        _Ntpb = c_long(0)
        lib.cuDeviceGetAttribute(
            byref(_Ntpb), CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, self._device
        )
        return _Ntpb.value

    def setCurrent(self):
        cu_print(lib.cuCtxSetCurrent(self._context), "ctx.setcurrent")

    def synchronize(self):
        cu_print(lib.cuCtxSynchronize(), "ctx.sync")

    def destroy(self):
        try:
            cu_print(lib.cuCtxDestroy_v2(self._context), "ctx.destroy")
        except (AttributeError):
            pass

    def __del__(self):
        self.destroy()


class CuModule:
    def __init__(self, context, module_name):

        # public:
        self.module_name = module_name
        self.context = context
        self.mode = 'GPU'

        # private:
        self._module = c_void_p(0)
        self._func_dict = {}
        self._global_dict = {}

        _module_name = c_char_p(self.module_name.encode())
        err = lib.cuModuleLoad(byref(self._module), _module_name)
        if err != CUDA_SUCCESS:
            if err == 222:
                print("* CUDA Driver too old, please update driver\n")
            print(
                "* Error loading the module {:s}\n".format(_module_name.value.decode())
            )
            self.context.destroy()
            # sys.exit()

    def __getattr__(self, attr):
        try:
            self._func_dict[attr]
            return self._func_dict[attr]

        except (KeyError):
            _function = c_void_p(0)
            _kernel_name = c_char_p(attr.encode())
            err = lib.cuModuleGetFunction(byref(_function), self._module, _kernel_name)
            if err != CUDA_SUCCESS:
                print(
                    "* Error getting kernel function {:s}".format(
                        _kernel_name.value.decode()
                    )
                )
                self.context.destroy()
                # sys.exit()

            self._func_dict[attr] = CuFunction(_function)
            self._func_dict[attr].module = self
            return self._func_dict[attr]

    def _getGlobal(self, name, ctype=None):
        try:
            return self._global_dict[name]

        except (KeyError):
            if ctype is not None:
                _var = c_void_p()
                _size = c_size_t()
                _type = ctype
                cu_print(
                    lib.cuModuleGetGlobal_v2(
                        byref(_var), byref(_size), self._module, c_char_p(name.encode())
                    ),
                    "mod.global",
                )
                self._global_dict[name] = (_var, _size, _type)
                return self._global_dict[name]

            else:
                print('Constant type must be passed first time it is called!')
                return
    def getMode(self):
        return self.mode
            
    def setConstant(self, name, cval):
        _var, _size, _type = self._getGlobal(name, type(cval))
        cu_print(lib.cuMemcpyHtoD_v2(_var, byref(cval), _size), "mod.HtoD")
        

    def getConstant(self, name, ctype=None):
        _var, _size, _type = self._getGlobal(name, ctype)
        cvar = _type()    
        cu_print(lib.cuMemcpyDtoH_v2(byref(cvar), _var,  _size), "mod.DtoH")
        return cvar

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

        voidPtrArr = len(self.args) * c_void_p
        c_args = voidPtrArr(*[cast(byref(arr._ptr), c_void_p) for arr in self.args])

        cu_print(
            lib.cuLaunchKernel(
                self._function, *self.blocks, *self.threads, 0, 0, c_args, 0
            ),
            "func.kernel",
        )

        if self.sync:
            self.module.context.synchronize()


class CuArray:
    def __init__(self, shape, dtype=np.float32, init="empty"):
        self._ptr = c_void_p()
        self.resize(shape, dtype, init)

    def resize(self, shape=None, dtype=None, init="empty"):
        shape_tuple = shape if isinstance(shape, tuple) else (shape,)
        self.shape = self.shape if shape is None else shape_tuple
        self.dtype = self.dtype if dtype is None else np.dtype(dtype)
        self.size = int(np.prod(self.shape))
        self.itemsize = self.dtype.itemsize
        self.nbytes = self.itemsize * self.size

        if init not in ("zeros", "empty"):
            return

        if self._ptr.value is not None:
            cu_print(lib.cuMemFree_v2(self._ptr), "arr.free")

        cu_print(lib.cuMemAlloc_v2(byref(self._ptr), self.nbytes), "arr.alloc")

        if init == "zeros":
            self.zeroFill()

    @staticmethod
    def fromArray(arr):

        obj = CuArray(arr.shape, arr.dtype, init=None)
        obj.setArray(arr)

        return obj

    def zeroFill(self):
        cu_print(lib.cuMemsetD8_v2(self._ptr, 0, self.nbytes), "arr.zeros")

    def setArray(self, arr):

        params_changed = arr.shape != self.shape or arr.dtype != self.dtype
        uninitialized_memory = self._ptr.value is None
        if params_changed or uninitialized_memory:
            self.resize(arr.shape, arr.dtype, "empty")

        cu_print(
            lib.cuMemcpyHtoD_v2(self._ptr, c_void_p(arr.ctypes.data), self.nbytes),
            "arr.HtoD",
        )

    def getArray(self):

        arr = np.empty(self.shape, dtype=self.dtype)
        cu_print(
            lib.cuMemcpyDtoH_v2(c_void_p(arr.ctypes.data), self._ptr, self.nbytes),
            "arr.DtoH",
        )

        return arr

    def free(self):
        cu_print(lib.cuMemFree_v2(self._ptr), "arr.free")


    def __del__(self):

        try:
            self.free()
        except (AttributeError):
            pass


class CuFFT:
    def __init__(self, arr_in, arr_out, direction="fwd", plan_fft=True):
        global lib_cufft
        # public:
        self.arr_in = arr_in
        self.arr_out = arr_out

        # private:
        self._direction = direction
        self._arr = arr_in if direction == "fwd" else arr_out
        self._fft_type = CUFFT_R2C if direction == "fwd" else CUFFT_C2R
        self._plan = c_void_p(0)


        cufft_name = '\\bin\\cufft64_10.dll' if os_name == 'nt' else '\\lib\\libcufft.so.11'
        try:
            if lib_cufft is None:
                lib_cufft = LoadLibrary(cufft_path[0] + cufft_name)
        except(FileNotFoundError):
            print("Can't find {:s}...".format(cufft_name))
            return None

        cu_print(lib_cufft.cufftCreate(byref(self._plan)), "fft.plan create")

        if plan_fft:
            self.planMany()

    def planMany(self):

        batch = int(np.prod(self._arr.shape[1:]))
        oneInt = 1 * c_int
        _n = oneInt(self._arr.shape[0])
        stride = batch

        cu_print(
            lib_cufft.cufftPlanMany(
                byref(self._plan),
                1,  # rank,
                _n,  # n
                oneInt(0),  # inembed,
                stride,  # istride
                1,  # idist,
                oneInt(0),  # onembed,
                stride,  # ostride
                1,  # odist,
                self._fft_type,
                batch,
            ),
            "fft.plan many",
        )

    def __call__(self):
        if self._direction == "fwd":
            cu_print(
                lib_cufft.cufftExecR2C(self._plan, self.arr_in._ptr, self.arr_out._ptr),
                "fft.fwd",
            )
        else:
            cu_print(
                lib_cufft.cufftExecC2R(self._plan, self.arr_in._ptr, self.arr_out._ptr),
                "fft.rev",
            )

    def destroy(self):
        try:
            cu_print(lib_cufft.cufftDestroy(self._plan), "fft.destroy")
        except (AttributeError):
            pass

    def __del__(self):
        try:
            if self._plan is not c_void_p(0):
                self.destroy()
        except(OSError):
            pass


class CuTimer:
    def __init__(self, stream=0, flags=0):
        #private:
        self._stream = c_int(stream)
        self._start = c_void_p(0)
        self._stop = c_void_p(0)

        #public:
        self.times = {}

        cu_print(lib.cuEventCreate(byref(self._start), flags|1), 'time.create')
        cu_print(lib.cuEventCreate(byref(self._stop), flags|1), 'time.create')
        self.reset()
                     
    def reset(self):
        cu_print(lib.cuEventRecord(self._start, self._stream), 'time.record')

    def lap(self, name=None):
        if name is None:
            name = 'event{:d}'.format(len(self.times.keys()))
        self.times[name] = self()    

    def __call__(self):
        _time = c_float(0) #elapsed time in ms
        cu_print(lib.cuEventRecord(self._stop, self._stream), 'time.record')
        cu_print(lib.cuEventSynchronize(self._stop), 'time.sync')
        cu_print(lib.cuEventElapsedTime(byref(_time), self._start, self._stop), 'time.time')
        return _time.value

    def getTimes(self):
        return self.times

    def getDiffs(self):
        vals = [*self.times.values()]
        diffs = [vals[0]] + [vals[i] - vals[i-1] for i in range(1, len(vals))]
        return dict(zip(self.times.keys(), diffs))

    
