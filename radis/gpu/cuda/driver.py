# import sys
from ctypes import (
    byref,
    c_char_p,
    c_float,
    c_int,
    c_long,
    c_size_t,
    c_void_p,
    cast,
    create_string_buffer,
)
from os import name as os_name
from warnings import warn

if os_name == "nt":
    from ctypes import windll as dllobj
else:
    from ctypes import cdll as dllobj


import numpy as np
from nvidia.cufft import __path__ as cufft_path

_verbose = False


def getClasses():
    return CuContext, CuModule, CuArray, CuFFT, CuTimer


def getCUDAVersion(ptx_file):
    # Reads the version of the CUDA Toolkit that was used to compile PTX file
    with open(ptx_file, "rb") as f:
        for i in range(5):
            line = f.readline()
        version_string = line.decode().split(",")[2]
        version_string = version_string[version_string.find("V") + 1 :]
        major, minor, patch = map(int, version_string.split("."))
        return major, minor, patch


def getRequiredDriverVersion(ptx_file):
    # Returns the minimal required CUDA driver version
    major, minor, patch = getCUDAVersion(ptx_file)
    if major == 11:
        if minor == 0:
            drv_ver = "v451.22" if os_name == "nt" else "v450.36.06"
        else:
            drv_ver = "v452.39" if os_name == "nt" else "v450.80.02"

    elif major == 12:
        drv_ver = "v527.41" if os_name == "nt" else "v525.60.13"

    else:
        drv_ver = "[latest version]"

    return drv_ver


def cu_print(*vargs):
    global _verbose
    if _verbose:
        print(*vargs)


CUDA_SUCCESS = 0
CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 1
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76

CUFFT_R2C = 0x2A
CUFFT_C2R = 0x2C

lib = None
lib_cufft = None

_arrays = []
_plans = []
_modules = []


class CuContext:
    def __init__(self, device, context, stream):

        # private:
        self._device = device
        self._context = context
        self._stream = stream

    @staticmethod
    def setVerbosity(level):
        global _verbose
        if level >= 2:
            _verbose = level

    @staticmethod
    def Open(device_id=0, flags=0, verbose=None):
        global lib

        if verbose is not None:
            CuContext.setVerbosity(verbose)

        # private:
        _context = c_void_p(0)
        _device = c_void_p(0)
        _stream = c_void_p(0)

        cuda_name = "nvcuda.dll" if os_name == "nt" else "libcuda.so"
        try:
            lib = dllobj.LoadLibrary(cuda_name)

        except Exception:
            if os_name == "nt":
                print("Can't find {:s}...".format(cuda_name))
                return None

            else:
                cuda_name = "libcuda.so.1"
                try:
                    lib = dllobj.LoadLibrary(cuda_name)
                except (OSError):
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
            lib.cuCtxDestroy_v2(_context)
            return None

        # cu_print(lib.cuStreamCreate(byref(_stream), 0), "ctx.create stream")

        return CuContext(_device, _context, _stream)

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

    def syncStream(self):
        cu_print(lib.cuStreamSynchronize(0), "ctx.sync stream")  # _ptsz

    def destroy(self):
        global _arrays, _plans, _modules

        while len(_arrays):
            arr = _arrays.pop(0)
            arr.free()

        while len(_plans):
            plan = _plans.pop(0)
            plan.destroy()

        while len(_modules):
            mod = _modules.pop(0)
            mod.unload()

        if self._stream.value is not None:
            cu_print(lib.cuStreamDestroy_v2(self._stream), "ctx.destroy stream")

        try:
            cu_print(lib.cuCtxDestroy_v2(self._context), "ctx.destroy")
        except (AttributeError):
            pass


class CuModule:
    def __init__(self, context, module_name):

        # public:
        self.module_name = module_name
        self.context = context
        self.mode = "GPU"
        global _modules
        _modules.append(self)

        # private:
        self._module = c_void_p(0)
        self._func_dict = {}
        self._global_dict = {}

        _module_name = c_char_p(self.module_name.encode())
        err = lib.cuModuleLoad(byref(self._module), _module_name)
        if err != CUDA_SUCCESS:
            if err == 222:
                drv_ver = getRequiredDriverVersion(module_name)
                print(
                    "\n\n*** CUDA Driver too old!***\nMinimally required driver version is {:s}\n".format(
                        drv_ver
                    )
                    + "Please update driver by downloading it at:\n-> www.nvidia.com/download/index.aspx \n"
                )

            warn(
                "* Error loading the module {:s}\n".format(_module_name.value.decode())
            )
            self.context.destroy()
            return

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
                cu_print("size: ", _size.value)
                return self._global_dict[name]

            else:
                print("Constant type must be passed first time it is called!")
                return

    def getMode(self):
        return self.mode

    def setConstant(self, name, cval):
        _var, _size, _type = self._getGlobal(name, type(cval))
        cu_print(lib.cuMemcpyHtoD_v2(_var, byref(cval), _size), "mod.HtoD")

    def getConstant(self, name, ctype=None):
        _var, _size, _type = self._getGlobal(name, ctype)
        cvar = _type()
        cu_print(lib.cuMemcpyDtoH_v2(byref(cvar), _var, _size), "mod.DtoH")
        return cvar

    def unload(self):
        cu_print(lib.cuModuleUnload(self._module), "mod.unload")


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

        if len(vargs):
            self.args = vargs
        if blocks is not None:
            self.blocks = blocks
        if threads is not None:
            self.threads = threads
        if sync is not None:
            self.sync = sync

        _stream = c_void_p(0)  # default stream; else: self.module.context._stream
        voidPtrArr = len(self.args) * c_void_p
        _c_args = voidPtrArr(*[cast(byref(arr._ptr), c_void_p) for arr in self.args])
        _extra_args = (1 * c_void_p)(0)

        cu_print(
            lib.cuLaunchKernel(  # _ptsz
                self._function,
                *self.blocks,
                *self.threads,
                0,  # shared memory size #TODO: get capability; set max size
                _stream,
                _c_args,
                _extra_args,
            ),
            "func.kernel",
        )

        if self.sync:
            self.module.context.synchronize()


class CuArray:
    def __init__(self, shape, dtype=np.float32, init="empty", grow_only=False):
        self._ptr = c_void_p()
        self._nbytes_alloc = 0
        self.grow_only = grow_only
        self.resize(shape, dtype, init)
        global _arrays
        _arrays.append(self)

    def resize(self, shape=None, dtype=None, init="empty"):

        shape_tuple = shape if isinstance(shape, tuple) else (shape,)
        self.shape = self.shape if shape is None else shape_tuple
        self.dtype = self.dtype if dtype is None else np.dtype(dtype)
        self.size = int(np.prod(self.shape))
        self.itemsize = self.dtype.itemsize
        self.nbytes = self.itemsize * self.size

        if init not in ("zeros", "empty"):
            return

        if self.nbytes > self._nbytes_alloc or (
            (self.nbytes < self._nbytes_alloc) and (not self.grow_only)
        ):
            if self._ptr.value is not None:
                self.free()
                self._ptr = c_void_p(0)

            cu_print(lib.cuMemAlloc_v2(byref(self._ptr), self.nbytes), "arr.alloc")
            self._nbytes_alloc = self.nbytes

            if self._ptr.value is not None:
                cu_print(hex(self._ptr.value), self.nbytes)

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


class CuFFT:
    def __init__(self, arr_in, arr_out, workarea=None, direction="fwd"):
        global lib_cufft, _plans
        _plans.append(self)

        # public:
        self.arr_in = arr_in
        self.arr_out = arr_out
        self.workarea = (
            CuArray(0, dtype=np.byte, grow_only=True) if workarea is None else workarea
        )

        # private:
        self._direction = direction
        self._fft_type = CUFFT_R2C if direction == "fwd" else CUFFT_C2R
        self._arr = arr_in if direction == "fwd" else arr_out
        self._plans = {}

        cufft_name = (
            "\\bin\\cufft64_10.dll" if os_name == "nt" else "/lib/libcufft.so.10"
        )
        try:
            if lib_cufft is None:
                lib_cufft = dllobj.LoadLibrary(cufft_path[0] + cufft_name)
        except (FileNotFoundError):
            print("Can't find {:s}...".format(cufft_name))
            return None

    def _getPlan(self):
        try:
            return self._plans[self._arr.shape]

        except (KeyError):

            _plan = c_void_p(0)
            cu_print(lib_cufft.cufftCreate(byref(_plan)), "fft.plan create")
            cu_print(lib_cufft.cufftSetAutoAllocation(_plan, c_int(0)), "fft.setalloc")

            batch = int(np.prod(self._arr.shape[1:]))
            oneInt = 1 * c_int
            _n = oneInt(self._arr.shape[0])
            stride = batch
            oneSizeT = 1 * c_size_t
            _worksize = oneSizeT(0)

            cu_print(
                lib_cufft.cufftMakePlanMany(
                    _plan,
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
                    _worksize,
                ),
                "fft.plan many",
            )

            self.workarea.resize((_worksize[0],), dtype=np.byte)
            cu_print(
                lib_cufft.cufftSetWorkArea(_plan, self.workarea._ptr),
                "fft.set workarea",
            )

            self._plans[self._arr.shape] = [_plan, int(self.workarea._ptr.value)]
            return self._plans[self._arr.shape]

    def __call__(self):

        _plan, ptrval = self._getPlan()

        if ptrval != self.workarea._ptr.value:
            cu_print(
                lib_cufft.cufftSetWorkArea(_plan, self.workarea._ptr),
                "fft.set workarea",
            )
            self._plans[self._arr.shape][1] = int(self.workarea._ptr.value)

        if self._direction == "fwd":
            cu_print(
                lib_cufft.cufftExecR2C(_plan, self.arr_in._ptr, self.arr_out._ptr),
                "fft.fwd",
            )
        else:
            cu_print(
                lib_cufft.cufftExecC2R(_plan, self.arr_in._ptr, self.arr_out._ptr),
                "fft.rev",
            )

    def destroy(self):

        for key in self._plans:
            _plan, workarea = self._plans[key]
            cu_print(lib_cufft.cufftDestroy(_plan), "fft.destroy")


class CuTimer:
    def __init__(self, stream=0, flags=0):
        # private:
        self._stream = c_int(stream)
        self._start = c_void_p(0)
        self._stop = c_void_p(0)

        # public:
        self.times = {}

        cu_print(lib.cuEventCreate(byref(self._start), flags), "time.create")
        cu_print(lib.cuEventCreate(byref(self._stop), flags), "time.create")
        self.reset()

    def reset(self):
        cu_print(lib.cuEventRecord(self._start, self._stream), "time.record")  # _ptsz

    def lap(self, name=None):
        if name is None:
            name = "event{:d}".format(len(self.times.keys()))
        self.times[name] = self()

    def __call__(self):
        _time = c_float(0)  # elapsed time in ms
        cu_print(lib.cuEventRecord(self._stop, self._stream), "time.record")  # _ptsz
        cu_print(lib.cuEventSynchronize(self._stop), "time.sync")
        cu_print(
            lib.cuEventElapsedTime(byref(_time), self._start, self._stop), "time.time"
        )
        return _time.value

    def getTimes(self):
        return self.times

    def getDiffs(self):
        vals = [*self.times.values()]
        diffs = [vals[0]] + [vals[i] - vals[i - 1] for i in range(1, len(vals))]
        return dict(zip(self.times.keys(), diffs))
