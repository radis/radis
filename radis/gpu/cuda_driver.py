import sys
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

from nvidia.cufft import __path__ as cufft_path

from radis.gpu import __path__ as gpu_path

lib = windll.LoadLibrary("nvcuda.dll")
lib_cufft = windll.LoadLibrary(cufft_path[0] + "\\bin\\cufft64_11.dll")


# TODO: make compatible with linux/mac
# TODO: implement cpu compatibility mode
# TODO: establish cuda version requirement


CUDA_SUCCESS = 0
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR = 75
CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR = 76

CUFFT_R2C = 0x2A
CUFFT_C2R = 0x2C


class cuContext:
    def __init__(self, device_id=0, flags=0):
        err = lib.cuInit(0)

        deviceCount = c_int(0)
        lib.cuDeviceGetCount(byref(deviceCount))

        if deviceCount == 0:
            print("Error: no devices supporting CUDA\n")
            sys.exit()

        self.device = c_long(0)
        lib.cuDeviceGet(byref(self.device), device_id)

        self.context = c_longlong(0)
        err = lib.cuCtxCreate_v2(byref(self.context), flags, self.device)
        if err != CUDA_SUCCESS:
            print("* Error initializing the CUDA context.")
            lib.cuCtxDestroy_v2(self.context)
            sys.exit()

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

    def synchronize(self):
        lib.cuCtxSynchronize()

    def destroy(self):
        lib.cuCtxDestroy_v2(self.context)

    def __del__(self):
        self.destroy()


class cuArray:
    def __init__(self, arr):
        self.arr_h = arr
        self.dtype = arr.dtype
        self.size = arr.size
        self.itemsize = arr.dtype.itemsize
        self.bytesize = arr.itemsize * arr.size
        self.shape = arr.shape

        self.arr_d = c_void_p()
        lib.cuMemAlloc_v2(byref(self.arr_d), self.bytesize)
        self.h2d()
        self.is_uptodate = True

    def __getitem__(self, i):
        if not self.is_uptodate:
            self.d2h()
            self.is_uptodate = True
        return self.arr_h[i]

    def __setitem__(self, i, a):

        self.arr_h[i] = a
        self.h2d()

    def h2d(self):

        lib.cuMemcpyHtoD_v2(self.arr_d, c_void_p(self.arr_h.ctypes.data), self.bytesize)

    def d2h(self):

        lib.cuMemcpyDtoH_v2(c_void_p(self.arr_h.ctypes.data), self.arr_d, self.bytesize)
        return self.arr_h

    # def __del__(self):

    # lib.cuMemFree_v2(self.arr_d)


class cuFunction:
    def __init__(self, fptr):
        self.fptr = fptr
        self.retvars = None

    def set_grid(self, blocks=(1, 1, 1), threads=(1, 1, 1)):
        self.blocks = blocks
        self.threads = threads

    def set_retvars(self, retvars):
        self.retvars = retvars

    def __call__(self, *vargs, **kwargs):

        try:
            self.blocks = kwargs["blocks"]
        except (KeyError):
            pass

        try:
            self.threads = kwargs["threads"]
        except (KeyError):
            pass

        try:
            self.retvars = kwargs["retvars"]
        except (KeyError):
            pass

        if self.retvars is None:
            self.retvars = len(vargs) * [False]
            self.retvars[-1] = True

        voidPtrArr = len(vargs) * c_void_p
        cargs = voidPtrArr(*[cast(byref(arr.arr_d), c_void_p) for arr in vargs])

        lib.cuLaunchKernel(self.fptr, *self.blocks, *self.threads, 0, 0, cargs, 0)

        for arr, is_retvar in zip(vargs, self.retvars):
            if is_retvar:
                arr.is_uptodate = False


class cuModule:
    def __init__(self, context, module_name):
        self.module_name = gpu_path[0] + "\\" + module_name
        module_file = c_char_p(self.module_name.encode())
        self.context_obj = context
        self.module = c_longlong(0)
        err = lib.cuModuleLoad(byref(self.module), module_file)
        if err != CUDA_SUCCESS:
            print(err)
            print(
                "* Error loading the module {:s}\n".format(module_file.value.decode())
            )
            self.context_obj.destroy()
            sys.exit()

        self.func_dict = {}

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
                sys.exit()

            self.func_dict[attr] = cuFunction(function)
            return self.func_dict[attr]

    def setConstant(self, name, c_val):
        var = c_void_p()
        size = c_long()

        lib.cuModuleGetGlobal_v2(
            byref(var), byref(size), self.module, c_char_p(name.encode())
        )
        lib.cuMemcpyHtoD_v2(var, byref(c_val), size)


class cuFFT:
    def __init__(self, arr_in, arr_out, direction="fwd"):
        self.arr_in = arr_in
        self.arr_out = arr_out
        self.direction = direction
        self.plan = c_longlong(0)

        lib_cufft.cufftCreate(byref(self.plan))
        lib_cufft.cufftPlan1d(byref(self.plan), arr_in.size, CUFFT_R2C, 1)

        arr = arr_in if direction == "fwd" else arr_out

        n = (1 * c_int)(arr.shape[0])
        istride = arr.shape[1] if len(arr.shape) > 1 else 0
        ostride = istride
        idist = 1
        odist = 1
        inembed = (1 * c_int)(0)
        onembed = (1 * c_int)(0)

        cufftType = CUFFT_R2C if direction == "fwd" else CUFFT_C2R
        batch = arr.shape[1] if len(arr.shape) > 1 else 1

        lib_cufft.cufftPlanMany(
            byref(self.plan),
            1,
            n,
            inembed,
            istride,
            idist,
            onembed,
            ostride,
            odist,
            cufftType,
            batch,
        )

    def execute(self):
        if self.direction == "fwd":
            lib_cufft.cufftExecR2C(self.plan, self.arr_in.arr_d, self.arr_out.arr_d)
        else:
            lib_cufft.cufftExecC2R(self.plan, self.arr_in.arr_d, self.arr_out.arr_d)

    def destroy(self):
        lib_cufft.cufftDestroy(self.plan)

    def __del__(self):
        self.destroy()
