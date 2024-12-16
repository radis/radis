# -*- coding: utf-8 -*-

# PyVkFFT
#   (c) 2021- : ESRF-European Synchrotron Radiation Facility
#       authors:
#         Vincent Favre-Nicolin, favre@esrf.fr

import ctypes
import os
from os import name as os_name

import numpy as np
from vulkan import ffi

# from pyvkfft.tune import tune_vkfft
##import platform
##import sysconfig


def getVulkanPtr(obj):
    return int(ffi.cast("unsigned long long", obj))


# Type to pass array size, omit and batch as int arrays
if os_name == "nt":  # windows
    libname = "vkfft_vulkan.dll"
    vkfft_long_type = np.int32
else:  # unix
    libname = "vkfft_vulkan.so"
    vkfft_long_type = np.int64

# Regrettably, VkFFT was compiled using "long" types, which have different sizes on Windows/Linux; we patch this here.
# TODO: recompile VkFFT replacing "long" with "int" (always 32 bit) or "long long' (always 64 bit).
ctype_int_size_p = np.ctypeslib.ndpointer(
    dtype=vkfft_long_type, ndim=1, flags="C_CONTIGUOUS"
)


from radis.gpu.vulkan.pyvkfft_base import (
    VkFFTApp as VkFFTAppBase,  # , check_vkfft_result
)

# from pyvkfft_base import VkFFTApp as VkFFTAppBase  # , check_vkfft_result


##def load_library(basename):
##    if platform.system() == "Windows":
##        # We patched build_ext so the module is a .so and not a dll
##        ext = ".so"
##    else:
##        ext = sysconfig.get_config_var("EXT_SUFFIX")
##    return ctypes.cdll.LoadLibrary(
##        os.path.join(os.path.dirname(__file__) or os.path.curdir, basename + ext)
##    )


vkfft_path = os.path.join(os.path.dirname(__file__), "bin", libname)
_vkfft_vulkan = ctypes.cdll.LoadLibrary(vkfft_path)


##stdout = os.dup(1)
##silent = os.open(os.devnull, os.O_WRONLY)


def prepare_fft(arr_in, arr_out=None, ndim=1, norm=1, compute_app=None, tune=False):

    tune_config = {"backend": "pycuda"} if tune else None

    if arr_out is None:
        arr_out = arr_in
        inplace = True
    else:
        inplace = False

    return VkFFTApp(
        arr_in.shape,
        arr_in.dtype,
        buffer_src=arr_in._buffer,
        buffer_dst=arr_out._buffer,
        physical_device=compute_app._physicalDevice,
        device=compute_app._device,
        queue=compute_app._queue,
        command_pool=compute_app._commandPool,
        fence=compute_app._fence,
        ndim=ndim,
        inplace=inplace,
        norm=norm,
        r2c=True,
        strides=arr_in.strides,
        tune_config=tune_config,
    )


class _types:
    """Aliases"""

    VkFFTConfiguration = ctypes.c_void_p
    VkBuffer = ctypes.c_void_p
    VkPhysicalDevice = ctypes.c_void_p
    VkDevice = ctypes.c_void_p
    VkQueue = ctypes.c_void_p
    VkCommandPool = ctypes.c_void_p
    VkFence = ctypes.c_void_p
    VkCommandBuffer = ctypes.c_void_p
    vkfft_app = ctypes.c_void_p


_vkfft_vulkan.get_dev_props.restype = ctypes.c_int
_vkfft_vulkan.get_dev_props.argtypes = [ctypes.c_void_p, ctypes.c_char_p]

_vkfft_vulkan.get_dev_props2.restype = ctypes.c_int
_vkfft_vulkan.get_dev_props2.argtypes = [_types.VkFFTConfiguration, ctypes.c_char_p]

_vkfft_vulkan.sync_app.restype = ctypes.c_int
_vkfft_vulkan.sync_app.argtypes = [_types.vkfft_app]


_vkfft_vulkan.make_config.restype = _types.VkFFTConfiguration
_vkfft_vulkan.make_config.argtypes = [
    ctype_int_size_p,
    ctypes.c_size_t,
    _types.VkBuffer,
    _types.VkBuffer,
    ctypes.POINTER(_types.VkPhysicalDevice),
    ctypes.POINTER(_types.VkDevice),
    ctypes.POINTER(_types.VkQueue),
    ctypes.POINTER(_types.VkCommandPool),
    ctypes.POINTER(_types.VkFence),
    ctypes.c_ulonglong,
    ctypes.c_int,
    ctypes.c_size_t,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_size_t,
    ctype_int_size_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctype_int_size_p,
]

_vkfft_vulkan.init_app.restype = _types.vkfft_app
_vkfft_vulkan.init_app.argtypes = [
    _types.VkFFTConfiguration,
    ctypes.POINTER(ctypes.c_int),
]

_vkfft_vulkan.fft.restype = ctypes.c_int
_vkfft_vulkan.fft.argtypes = [_types.vkfft_app, ctypes.c_void_p, ctypes.c_void_p]

_vkfft_vulkan.ifft.restype = ctypes.c_int
_vkfft_vulkan.ifft.argtypes = [_types.vkfft_app, ctypes.c_void_p, ctypes.c_void_p]

_vkfft_vulkan.free_app.restype = None
_vkfft_vulkan.free_app.argtypes = [_types.vkfft_app]

_vkfft_vulkan.free_config.restype = None
_vkfft_vulkan.free_config.argtypes = [_types.VkFFTConfiguration]

_vkfft_vulkan.vkfft_max_fft_dimensions.restype = ctypes.c_uint32
_vkfft_vulkan.vkfft_max_fft_dimensions.argtypes = None


class VkFFTApp(VkFFTAppBase):
    """
    VkFFT application interface, similar to a cuFFT plan.
    """

    def __init__(
        self,
        shape_in,
        dtype: type,
        buffer_src,
        buffer_dst,
        physical_device,
        device,
        queue,
        command_pool,
        fence,
        ndim=None,
        inplace=True,
        norm=1,
        r2c=False,
        dct=False,
        axes=None,
        strides=None,
        tune_config=None,
        **kwargs,
    ):

        """

        :param shape: the shape of the array to be transformed. The number
            of dimensions of the array can be larger than the FFT dimensions,
            but only for 1D and 2D transforms. 3D FFT transforms can only
            be done on 3D arrays.
        :param dtype: the numpy dtype of the source array (can be complex64 or complex128)
        :param ndim: the number of dimensions to use for the FFT. By default,
            uses the array dimensions. Can be smaller, e.g. ndim=2 for a 3D
            array to perform a batched 3D FFT on all the layers. The FFT
            is always performed along the last axes if the array's number
            of dimension is larger than ndim, i.e. on the x-axis for ndim=1,
            on the x and y axes for ndim=2.
        :param inplace: if True (the default), performs an inplace transform and
            the destination array should not be given in fft() and ifft().

        ~~~~~
        :param stream: the pycuda.driver.Stream or cupy.cuda.Stream to use
            for the transform. This can also be the pointer/handle (int) to the
            cuda stream object. If None, the default stream will be used.
        ~~~~~

        :param norm: if 0 (unnormalised), every transform multiplies the L2
            norm of the array by its size (or the size of the transformed
            array if ndim<d.ndim).
            if 1 (the default) or "backward", the inverse transform divides
            the L2 norm by the array size, so FFT+iFFT will keep the array norm.
            if "ortho", each transform will keep the L2 norm, but that will
            involve an extra read & write operation.
        :param r2c: if True, will perform a real->complex transform, where the
            complex destination is a half-hermitian array.
            For an inplace transform, if the input data shape is (...,nx), the input
            float array should have a shape of (..., nx+2), the last two columns
            being ignored in the input data, and the resulting
            complex array (using pycuda's GPUArray.view(dtype=np.complex64) to
            reinterpret the type) will have a shape (..., nx//2 + 1).
            For an out-of-place transform, if the input (real) shape is (..., nx),
            the output (complex) shape should be (..., nx//2+1).
            Note that for C2R transforms with ndim>=2, the source (complex) array
            is modified.
        :param dct: used to perform a Direct Cosine Transform (DCT) aka a R2R transform.
            An integer can be given to specify the type of DCT (1, 2, 3 or 4).
            if dct=True, the DCT type 2 will be performed, following scipy's convention.
        :param axes: a list or tuple of axes along which the transform should be made.
            if None, the transform is done along the ndim fastest axes, or all
            axes if ndim is None. Not allowed for R2C transforms
        :param strides: the array strides - needed if not C-ordered.
        :param tune_config: this can be used to automatically generate an
            optimised set of VkFFT parameters by testing various configurations
            and measuring the FFT speed, in a manner similar to fftw's FFTW_MEASURE.
            This should be a dictionary including the backend used and the parameter
            values which will be tested.
            This is EXPERIMENTAL, as wrong parameters may lead to crashes.
            Note that this will allocate temporary GPU arrays, unless the arrays
            to used have been passed as parameters ('dest' and 'src').
            Examples:
            tune={'backend':'cupy} - minimal example, will automatically test a small
            set of parameters (4 to 10 tests). Recommended !
            tune={'backend':'cupy, 'warpSize':[8,16,32,64,128]}: this will test
            5 possible values for the warpSize.
            tune={'backend':'cupy, 'groupedBatch':[[-1,-1,-1],[8,8,8], [4,16,16}:
            this will test 3 possible values for groupedBatch. This one is more
            tricky to use.
            tune={'backend':'cupy, 'warpSize':[8,16,32,64,128], 'src':a}: this
            will test 5 possible values for the warpSize, with a given source GPU
            array. This would only be valid for an inplace transform as no
            destination array is given.

        :raises RuntimeError: if the initialisation fails, e.g. if the CUDA
            driver has not been properly initialised, or if the transform dimensions
            are not allowed by VkFFT.
        """
        # if tune_config is not None:
        #     kwargs = tune_vkfft(tune_config, shape=shape, dtype=dtype, ndim=ndim, inplace=inplace, stream=stream,
        #                         norm=norm, r2c=r2c, dct=dct, axes=axes, strides=strides, verbose=False,
        #                         **kwargs)[0]
        super().__init__(
            shape_in,
            dtype,
            ndim=ndim,
            inplace=inplace,
            norm=norm,
            r2c=r2c,
            dct=dct,
            axes=axes,
            strides=strides,
            **kwargs,
        )

        self.bufferSrc = _types.VkBuffer(getVulkanPtr(buffer_src))
        self.bufferDest = _types.VkBuffer(getVulkanPtr(buffer_dst))

        self.physicalDevice = _types.VkPhysicalDevice(getVulkanPtr(physical_device))
        self.device = _types.VkDevice(getVulkanPtr(device))
        self.queue = _types.VkQueue(getVulkanPtr(queue))
        self.commandPool = _types.VkCommandPool(getVulkanPtr(command_pool))
        self.fence = _types.VkFence(getVulkanPtr(fence))

        # buf = ctypes.create_string_buffer(256)

        self.config = self._make_config()

        if self.config is None:
            raise RuntimeError(
                "Error creating VkFFTConfiguration. Was the CUDA context properly initialised ?"
            )
        res = ctypes.c_int(0)

        self.app = _vkfft_vulkan.init_app(self.config, ctypes.byref(res))

        #!!!
        # check_vkfft_result(res, shape, dtype, ndim, inplace, norm, r2c, dct, axes, "cuda")
        if self.app is None:
            raise RuntimeError(
                "Error {:d}  creating VkFFTApplication. Was the Vulkan driver initialised ?".format(
                    res.value
                )
            )

    def __del__(self):
        """Takes care of deleting allocated memory in the underlying
        VkFFTApplication and VkFFTConfiguration.
        """

        # print('VkFFT.app.__del__... ')
        # print('free_app @',hex(self.app), end='... ')
        if self.app is not None:
            _vkfft_vulkan.free_app(self.app)
            self.app = None
        # print('Done!')

        # print('free_config @',hex(self.config), end='... ')
        if self.config is not None:
            _vkfft_vulkan.free_config(self.config)
            self.config = None
        # print('Done!')

        # print('VkFFT.app.__del__ done!')

    def _make_config(self):
        """Create a vkfft configuration for a FFT transform"""
        if len(self.shape) > vkfft_max_fft_dimensions():
            raise RuntimeError(
                f"Too many FFT dimensions after collapsing non-transform axes: "
                f"{len(self.shape)}>{vkfft_max_fft_dimensions()}"
            )

        shape = np.ones(vkfft_max_fft_dimensions(), dtype=vkfft_long_type)
        # shape[:len(self.shape)] = self.shape

        skip = np.zeros(vkfft_max_fft_dimensions(), dtype=vkfft_long_type)
        # skip[:len(self.skip_axis)] = self.skip_axis

        shape[: len(self.shape)] = self.shape[::-1]
        skip[1 : len(self.shape)] = 1
        FFTdim = len(self.shape)
        n_batch = 1

        grouped_batch = np.empty(vkfft_max_fft_dimensions(), dtype=vkfft_long_type)
        grouped_batch.fill(-1)
        grouped_batch[: len(self.groupedBatch)] = self.groupedBatch

        if self.r2c and self.inplace:
            # the last two columns are ignored in the R array, and will be used
            # in the C array with a size nx//2+1
            shape[0] -= 2

        # s = 0
        # if self.stream is not None:
        # # if isinstance(self.stream, cp.cuda.Stream):
        # # s = self.stream.ptr
        # if s == 0 and isinstance(self.stream, int):
        # # Assume the ptr or handle was passed
        # s = self.stream

        if self.norm == "ortho":
            norm = 0
        else:
            norm = self.norm

        # We pass fake buffer pointer addresses to VkFFT. The real ones will be
        # given when performing the actual FFT.
        # dest_gpudata = 2
        # if self.inplace:
        # dest_gpudata = 0

        # print('physicalDevice:', '0x'+hex(self.compute_app.getVulkanPtr('_physicalDevice'))[2:].upper())

        # ptr = ctypes.c_void_p(self.compute_app.getVulkanPtr('_physicalDevice'))
        # _vkfft_vulkan.get_dev_props(ctypes.byref(ptr), buf)
        return _vkfft_vulkan.make_config(
            shape,
            FFTdim,
            self.bufferSrc,
            self.bufferDest,
            # ctypes.c_void_p(0), ctypes.c_void_p(0),
            ctypes.byref(self.physicalDevice),
            ctypes.byref(self.device),
            ctypes.byref(self.queue),
            ctypes.byref(self.commandPool),
            ctypes.byref(self.fence),
            0,
            norm,
            self.precision,
            int(self.r2c),
            int(self.dct),
            int(self.disableReorderFourStep),
            int(self.registerBoost),
            int(self.use_lut),
            int(self.keepShaderCode),
            n_batch,
            skip,
            int(self.coalescedMemory),
            int(self.numSharedBanks),
            int(self.aimThreads),
            int(self.performBandwidthBoost),
            int(self.registerBoostNonPow2),
            int(self.registerBoost4Step),
            int(self.warpSize),
            grouped_batch,
        )

    def getBufSize(self, src):
        self.bufferSrc = _types.VkBuffer(getVulkanPtr(src))
        m = _vkfft_vulkan.get_buf_size(self.bufferSrc, ctypes.byref(self.device))
        print("Returned buffer size: ", m)

    # def sync(self):
    # res = _vkfft_vulkan.sync_app(self.app)

    def fft(self, cmd_buf, src, dest=None):
        """
        Compute the forward FFT

        :param src: the source pycuda.gpuarray.GPUArray or cupy.ndarray
        :param dest: the destination GPU array. Should be None for an inplace transform
        :raises RuntimeError: in case of a GPU kernel launch error
        :return: the transformed array. For a R2C inplace transform, the complex view of the
            array is returned.
        """

        if dest is None:
            assert self.inplace
            dest = src

        self.bufferSrc = _types.VkBuffer(getVulkanPtr(src))
        self.bufferDest = _types.VkBuffer(getVulkanPtr(dest))
        self.commandBufferFwd = _types.VkCommandBuffer(getVulkanPtr(cmd_buf))

        _vkfft_vulkan.fft(
            self.app,
            ctypes.byref(self.commandBufferFwd),
            ctypes.byref(self.bufferSrc),
            ctypes.byref(self.bufferDest),
        )

    def ifft(self, cmd_buf, src, dest=None):
        """
        Compute the backward FFT

        :param src: the source pycuda.gpuarray.GPUArray or cupy.ndarray
        :param dest: the destination GPU array. Should be None for an inplace transform
        :raises RuntimeError: in case of a GPU kernel launch error
        :return: the transformed array. For a C2R inplace transform, the float view of the
            array is returned.
        """
        if dest is None:
            assert self.inplace
            dest = src

        self.bufferSrc = _types.VkBuffer(getVulkanPtr(src))
        self.bufferDest = _types.VkBuffer(getVulkanPtr(dest))
        self.commandBufferRev = _types.VkCommandBuffer(getVulkanPtr(cmd_buf))

        _vkfft_vulkan.ifft(
            self.app,
            ctypes.byref(self.commandBufferRev),
            ctypes.byref(self.bufferDest),
            ctypes.byref(self.bufferSrc),
        )


def vkfft_version():
    """
    Get VkFFT version
    :return: version as X.Y.Z
    """
    int_ver = _vkfft_vulkan.vkfft_version()
    return "%d.%d.%d" % (int_ver // 10000, (int_ver % 10000) // 100, int_ver % 100)


def vkfft_max_fft_dimensions():
    """
    Get the maximum number of dimensions VkFFT can handle. This is
    set at compile time. VkFFT default is 4, pyvkfft sets this to 8.
    Note that consecutive non-transformed are collapsed into a single
    axis, reducing the effective number of dimensions.

    :return: VKFFT_MAX_FFT_DIMENSIONS
    """
    return _vkfft_vulkan.vkfft_max_fft_dimensions()


##def cuda_runtime_version(raw=False):
##    """
##    Get CUDA runtime version
##
##    :param raw: if True, return the version as X*1000+Y*10+Z
##    :return: version as X.Y.Z
##    """
##    int_ver = _vkfft_cuda.cuda_runtime_version()
##    if raw:
##        return raw
##    return "%d.%d.%d" % (int_ver // 1000, (int_ver % 1000) // 10, int_ver % 10)
##
##
##def cuda_driver_version(raw=False):
##    """
##    Get CUDA driver version
##
##    :param raw: if True, return the version as X*1000+Y*10+Z
##    :return: version as X.Y.Z
##    """
##    int_ver = _vkfft_cuda.cuda_driver_version()
##    if raw:
##        return raw
##    return "%d.%d.%d" % (int_ver // 1000, (int_ver % 1000) // 10, int_ver % 10)
##
##
##def cuda_compile_version(raw=False):
##    """
##    Get CUDA version against which pyvkfft was compiled
##
##    :param raw: if True, return the version as X*1000+Y*10+Z
##    :return: version as X.Y.Z
##    """
##    if raw:
##        return raw
##    int_ver = _vkfft_cuda.cuda_compile_version()
##    return "%d.%d.%d" % (int_ver // 1000, (int_ver % 1000) // 10, int_ver % 10)
