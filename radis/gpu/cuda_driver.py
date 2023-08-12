# import sys
from ctypes import (
    byref,
    c_char_p,
    c_int,
    c_long,
    c_longlong,
    c_size_t,
    c_void_p,
    cast,
    create_string_buffer,
    windll,
)

import numpy as np
from nvidia.cufft import __path__ as cufft_path

lib = windll.LoadLibrary("nvcuda.dll")
lib_cufft = windll.LoadLibrary(cufft_path[0] + "\\bin\\cufft64_11.dll")

verbose = True


def cu_print(*vargs):
    global verbose
    if verbose:
        print(*vargs)


# TODO: make compatible with linux/mac
# TODO: implement cpu compatibility mode
# TODO: establish cuda version requirement


CUDA_SUCCESS = 0
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76

CUFFT_R2C = 0x2A
CUFFT_C2R = 0x2C


class CuContext:
    def __init__(self, device_id=0, flags=0):
        err = lib.cuInit(0)

        deviceCount = c_int(0)
        lib.cuDeviceGetCount(byref(deviceCount))

        if deviceCount == 0:
            print("Error: no devices supporting CUDA\n")
            # sys.exit()

        self.device = c_long(0)
        lib.cuDeviceGet(byref(self.device), device_id)

        self.context = c_longlong(0)
        err = lib.cuCtxCreate_v2(byref(self.context), flags, self.device)
        if err != CUDA_SUCCESS:
            print("* Error initializing the CUDA context.")
            lib.cuCtxDestroy_v2(self.context)
            # sys.exit()

    @staticmethod
    def getDeviceList():
        dev_list = []
        deviceCount = c_int(0)

        lib.cuInit(0)
        lib.cuDeviceGetCount(byref(deviceCount))

        for i in range(deviceCount.value):

            device = c_long(0)
            lib.cuDeviceGet(byref(device), i)

            name = create_string_buffer(100)
            lib.cuDeviceGetName(name, 100, device)
            dev_list.append(name.value.decode())

        return dev_list

    def printDeviceCapabilities(self):

        major = c_long(0)
        minor = c_long(0)
        lib.cuDeviceGetAttribute(
            byref(major), CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, self.device
        )
        lib.cuDeviceGetAttribute(
            byref(minor), CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, self.device
        )
        print(
            "> GPU Device has SM {:d}.{:d} compute capability".format(
                major.value, minor.value
            )
        )

        totalGlobalMem = c_size_t(0)
        lib.cuDeviceTotalMem_v2(byref(totalGlobalMem), self.device)
        print(
            "  Total amount of global memory:   {:d} bytes".format(totalGlobalMem.value)
        )
        print(
            "  64-bit Memory Address:           {:s}".format(
                "YES" if totalGlobalMem.value > (2 << 31) else "NO"
            )
        )

    def setCurrent(self):
        cu_print(lib.cuCtxSetCurrent(self.context), "ctx.setcurrent")

    def synchronize(self):
        cu_print(lib.cuCtxSynchronize(), "ctx.sync")

    def destroy(self):
        try:
            cu_print(lib.cuCtxDestroy_v2(self.context), "ctx.destroy")
        except (AttributeError):
            pass

    def __del__(self):
        self.destroy()


class CuArray:
    def __init__(self, shape, dtype=np.float32, init="empty"):
        self.dev_ptr = c_void_p()
        self.resize(shape, dtype, init)

    def resize(self, shape=None, dtype=None, init="empty"):
        shape_tuple = shape if isinstance(shape, tuple) else (shape,)
        self.shape = self.shape if shape is None else shape_tuple
        self.dtype = self.dtype if dtype is None else np.dtype(dtype)
        self.size = int(np.prod(shape))
        self.itemsize = self.dtype.itemsize
        self.bytesize = self.itemsize * self.size

        if init not in ("zeros", "empty"):
            return

        if self.dev_ptr.value is not None:
            cu_print(lib.cuMemFree_v2(self.dev_ptr), "arr.free")

        cu_print(lib.cuMemAlloc_v2(byref(self.dev_ptr), self.bytesize), "arr.alloc")

        if init == "zeros":
            self.zeroFill()

    @staticmethod
    def fromArray(arr):

        obj = CuArray(arr.shape, arr.dtype, init=None)
        obj.setArray(arr)

        return obj

    def zeroFill(self):
        cu_print(lib.cuMemsetD8_v2(self.dev_ptr, 0, self.bytesize), "arr.zeros")

    def setArray(self, arr):

        params_changed = arr.shape != self.shape or arr.dtype != self.dtype
        uninitialized_memory = self.dev_ptr.value is None
        if params_changed or uninitialized_memory:
            self.resize(arr.shape, arr.dtype, "empty")

        cu_print(
            lib.cuMemcpyHtoD_v2(self.dev_ptr, c_void_p(arr.ctypes.data), self.bytesize),
            "arr.HtoD",
        )

    def getArray(self):

        arr = np.empty(self.shape, dtype=self.dtype)
        cu_print(
            lib.cuMemcpyDtoH_v2(c_void_p(arr.ctypes.data), self.dev_ptr, self.bytesize),
            "arr.DtoH",
        )

        return arr

    def __del__(self):

        try:
            cu_print(lib.cuMemFree_v2(self.dev_ptr), "arr.free")
        except (AttributeError):
            pass


class CuFunction:
    def __init__(self, fptr):
        self.fptr = fptr
        self.retvars = None
        self.blocks = (1, 1, 1)
        self.threads = (1, 1, 1)
        self.sync = False
        self.context_obj = None

    def setGrid(self, blocks=(1, 1, 1), threads=(1, 1, 1)):
        self.blocks = blocks
        self.threads = threads

    def __call__(self, *vargs, blocks=None, threads=None, sync=None):

        self.blocks = self.blocks if blocks is None else blocks
        self.threads = self.threads if threads is None else threads
        self.sync = self.sync if sync is None else sync

        voidPtrArr = len(vargs) * c_void_p
        cargs = voidPtrArr(*[cast(byref(arr.dev_ptr), c_void_p) for arr in vargs])

        cu_print(
            lib.cuLaunchKernel(self.fptr, *self.blocks, *self.threads, 0, 0, cargs, 0),
            "func.kernel",
        )

        if self.sync:
            self.context_obj.synchronize()


class CuModule:
    def __init__(self, context_obj, module_name):
        self.module_name = module_name
        module_file = c_char_p(self.module_name.encode())
        self.context_obj = context_obj
        self.module = c_longlong(0)
        err = lib.cuModuleLoad(byref(self.module), module_file)
        if err != CUDA_SUCCESS:
            print(err)
            print(
                "* Error loading the module {:s}\n".format(module_file.value.decode())
            )
            self.context_obj.destroy()
            # sys.exit()

        self.func_dict = {}
        self.global_dict = {}

    def __getattr__(self, attr):
        try:
            self.func_dict[attr]
            return self.func_dict[attr]

        except (KeyError):
            function = c_longlong(0)
            kernel_name = c_char_p(attr.encode())
            err = lib.cuModuleGetFunction(byref(function), self.module, kernel_name)
            if err != CUDA_SUCCESS:
                print(
                    "* Error getting kernel function {:s}".format(
                        kernel_name.value.decode()
                    )
                )
                self.context_obj.destroy()
                # sys.exit()

            self.func_dict[attr] = CuFunction(function)
            self.func_dict[attr].context_obj = self.context_obj
            return self.func_dict[attr]

    def setConstant(self, name, c_val):
        try:
            var, size = self.global_dict[name]

        except (KeyError):
            var = c_void_p()
            size = c_long()
            cu_print(
                lib.cuModuleGetGlobal_v2(
                    byref(var), byref(size), self.module, c_char_p(name.encode())
                ),
                "mod.global",
            )
            self.global_dict[name] = (var, size)

        cu_print(lib.cuMemcpyHtoD_v2(var, byref(c_val), size), "mod.HtoD")


class CuFFT:
    def __init__(self, arr_in, arr_out, direction="fwd", plan_fft=True):
        self.arr_in = arr_in
        self.arr_out = arr_out
        self._arr = arr_in if direction == "fwd" else arr_out
        self.direction = direction
        self._cufftType = CUFFT_R2C if direction == "fwd" else CUFFT_C2R

        self.plan_ptr = c_longlong(0)
        cu_print(lib_cufft.cufftCreate(byref(self.plan_ptr)), "fft.plan create")

        if plan_fft:
            self.planMany()

    def planMany(self):

        batch = int(np.prod(self._arr.shape[1:]))
        oneInt = 1 * c_int
        n = oneInt(self._arr.shape[0])
        stride = batch

        cu_print(
            lib_cufft.cufftPlanMany(
                byref(self.plan_ptr),
                1,  # rank,
                n,  # n
                oneInt(0),  # inembed,
                stride,  # istride
                1,  # idist,
                oneInt(0),  # onembed,
                stride,  # ostride
                1,  # odist,
                self._cufftType,
                batch,
            ),
            "fft.plan many",
        )

    def execute(self):
        if self.direction == "fwd":
            cu_print(
                lib_cufft.cufftExecR2C(
                    self.plan_ptr, self.arr_in.dev_ptr, self.arr_out.dev_ptr
                ),
                "fft.fwd",
            )
        else:
            cu_print(
                lib_cufft.cufftExecC2R(
                    self.plan_ptr, self.arr_in.dev_ptr, self.arr_out.dev_ptr
                ),
                "fft.rev",
            )

    def destroy(self):
        try:
            cu_print(lib_cufft.cufftDestroy(self.plan_ptr), "fft.destroy")
        except (AttributeError):
            pass

    def __del__(self):
        self.destroy()
