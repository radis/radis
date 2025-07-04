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
# vkfft_path = os.path.join(os.path.dirname(__file__), '..', 'src', 'vulkan')
vkfft_path = os.path.join(os.path.dirname(__file__), "bin")
if os_name == "nt":  # windows
    vkfft_path = os.path.join(vkfft_path, "vkfft_vulkan.dll")
    vkfft_long_type = np.int32
else:  # unix
    vkfft_path = os.path.join(vkfft_path, "vkfft_vulkan.so")
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

_vkfft_vulkan = ctypes.cdll.LoadLibrary(vkfft_path)


##stdout = os.dup(1)
##silent = os.open(os.devnull, os.O_WRONLY)


def prepare_fft(
    buffer, name="", norm=1, compute_app=None, indirectOffset=0, exclusivePlan=0
):

    return VkFFTApp(
        buffer=buffer,
        indirectOffset=indirectOffset,
        indirectBuffer=compute_app.indirect_d,
        indirectHost=compute_app._indirect_h,
        physical_device=compute_app._physicalDevice,
        device=compute_app._device,
        queue=compute_app._queue,
        command_pool=compute_app._commandPool,
        fence=compute_app._fence,
        ndim=1,
        inplace=True,
        norm=norm,
        r2c=True,
        tune_config=None,
        name=name,
        exclusivePlan=exclusivePlan,
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
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_size_t,
    _types.VkBuffer,
    _types.VkBuffer,
    # ctypes.c_int,
    # _types.VkBuffer,
    # ctypes.c_int,
    ctypes.c_int,
    _types.VkBuffer,
    ctypes.c_int,
    ctypes.c_void_p,
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
    ctypes.c_int,
    ctype_int_size_p,
    ctypes.c_char_p,
    ctypes.c_int,
]

_vkfft_vulkan.init_app.restype = _types.vkfft_app
_vkfft_vulkan.init_app.argtypes = [
    _types.VkFFTConfiguration,
    ctypes.POINTER(ctypes.c_int),
]

_vkfft_vulkan.fft.restype = ctypes.c_int
_vkfft_vulkan.fft.argtypes = [
    _types.vkfft_app,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
]

_vkfft_vulkan.ffto.restype = ctypes.c_int
_vkfft_vulkan.ffto.argtypes = [
    _types.vkfft_app,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_int,
    ctypes.c_int,
]

_vkfft_vulkan.ifft.restype = ctypes.c_int
_vkfft_vulkan.ifft.argtypes = [
    _types.vkfft_app,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
]

_vkfft_vulkan.iffto.restype = ctypes.c_int
_vkfft_vulkan.iffto.argtypes = [
    _types.vkfft_app,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_int,
    ctypes.c_int,
]

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
        buffer,
        # shape_in,
        # dtype: type,
        # buffer_size,
        # buffer_src,
        # buffer_dst,
        # currentBatchUBO,
        # currentBatchUBOOffset,
        indirectOffset,
        indirectBuffer,
        indirectHost,
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
        # strides=None,
        tune_config=None,
        name="",
        exclusivePlan=0,
        **kwargs,
    ):

        shape_in = (
            (buffer._fftSize,)
            if buffer._batchSize == 1
            else (buffer._batchSize, buffer._fftSize)
        )
        dtype = buffer._dtype
        super().__init__(
            shape_in,
            dtype,
            ndim=ndim,
            inplace=inplace,
            norm=norm,
            r2c=r2c,
            dct=dct,
            axes=axes,
            **kwargs,
        )

        self.bufferSize = buffer._bufferSize
        self.bufferSrc = _types.VkBuffer(getVulkanPtr(buffer._buffer))
        self.bufferDest = _types.VkBuffer(getVulkanPtr(buffer._buffer))
        self.indirectOffset = indirectOffset
        self.indirectBuffer = _types.VkBuffer(getVulkanPtr(indirectBuffer._buffer))
        self.indirectHost = indirectBuffer._structPtr.contents
        self.physicalDevice = _types.VkPhysicalDevice(getVulkanPtr(physical_device))
        self.device = _types.VkDevice(getVulkanPtr(device))
        self.queue = _types.VkQueue(getVulkanPtr(queue))
        self.commandPool = _types.VkCommandPool(getVulkanPtr(command_pool))
        self.fence = _types.VkFence(getVulkanPtr(fence))
        self.name = name.encode()
        self.exclusivePlan = exclusivePlan
        # buf = ctypes.create_string_buffer(256)

        self.config = self._make_config()

        if self.config is None:
            raise RuntimeError(
                "Error creating VkFFTConfiguration. Was the Vulkan context properly initialised ?"
            )
        res = ctypes.c_int(0)

        self.app = _vkfft_vulkan.init_app(self.config, ctypes.byref(res))

        if self.app is None:
            raise RuntimeError(
                f"Error {res.value:d}  creating VkFFTApplication. Was the Vulkan driver initialised ?"
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
        skip = np.zeros(vkfft_max_fft_dimensions(), dtype=vkfft_long_type)
        grouped_batch = -np.ones(vkfft_max_fft_dimensions(), dtype=vkfft_long_type)

        shape[0] = self.shape[-1]
        # skip[1 : len(self.shape)] = 1
        # FFTdim = 1
        n_batch = 1 if len(self.shape) == 1 else self.shape[-2]

        indirectDispatch = 1 if self.exclusivePlan == 1 else 0

        return _vkfft_vulkan.make_config(
            shape,
            self.bufferSize,
            0,  # self.bufferSize,
            1,  # FFTdim,
            self.bufferSrc,
            None,  # self.bufferDest,
            indirectDispatch,
            self.indirectBuffer,
            self.indirectOffset,
            ctypes.byref(self.indirectHost),
            ctypes.byref(self.physicalDevice),
            ctypes.byref(self.device),
            ctypes.byref(self.queue),
            ctypes.byref(self.commandPool),
            ctypes.byref(self.fence),
            0,  # isCompilerInitialized
            1,  # norm,
            self.precision,
            1,  # int(self.r2c),
            0,  # int(self.dct),
            0,  # int(self.disableReorderFourStep),
            int(self.registerBoost),
            int(self.use_lut),
            int(self.keepShaderCode),
            # min(6,n_batch),
            n_batch,
            skip,
            int(self.coalescedMemory),
            int(self.numSharedBanks),
            int(self.aimThreads),
            int(self.performBandwidthBoost),
            int(self.registerBoostNonPow2),
            int(self.registerBoost4Step),
            int(self.warpSize),
            0,  # 1,#int(1), #specify offset at launch
            grouped_batch,
            self.name,
            self.exclusivePlan,
            0,  # debugEnable
        )

    def getBufSize(self, src):
        self.bufferSrc = _types.VkBuffer(getVulkanPtr(src))
        m = _vkfft_vulkan.get_buf_size(self.bufferSrc, ctypes.byref(self.device))
        print("Returned buffer size: ", m)

    # def sync(self):
    # res = _vkfft_vulkan.sync_app(self.app)

    def fft(self, cmd_buf, src):
        """
        Compute the forward FFT

        :return: the transformed array. For a R2C inplace transform, the complex view of the
            array is returned.
        """

        self.bufferSrc = _types.VkBuffer(getVulkanPtr(src))
        self.commandBufferFwd = _types.VkCommandBuffer(getVulkanPtr(cmd_buf))

        _vkfft_vulkan.fft(
            self.app,
            ctypes.byref(self.commandBufferFwd),
            ctypes.byref(self.bufferSrc),
            None,
        )

    def ifft(self, cmd_buf, src):
        """
        Compute the backward FFT

        :return: the transformed array. For a C2R inplace transform, the float view of the
            array is returned.
        """

        self.bufferSrc = _types.VkBuffer(getVulkanPtr(src))
        self.commandBufferRev = _types.VkCommandBuffer(getVulkanPtr(cmd_buf))

        _vkfft_vulkan.ifft(
            self.app,
            ctypes.byref(self.commandBufferRev),
            None,
            ctypes.byref(self.bufferSrc),
        )

    def setFFTWorkGroupSize(self, N):
        for wg in self.indirectHost:
            # print(wg.id, ':', wg.x, wg.y, wg.z, '->', end=' ')
            if wg.id > 0:
                try:
                    ax = "xyz"[wg.id]
                    setattr(wg, ax, N)
                    # print(wg.id)
                except IndexError:
                    # the wg.id is set by VkFFT to either 0, 1, or 2,
                    # to specify which axis contains the batch number.
                    # this exception is to catch issues with the wg.id being any other number.
                    # However, this issue was likely solved by zeroing the values in
                    # vulkan_compute_lib.setIndirectBuffer().
                    print("Vulkan Workgroup ID error:", wg.id)
            # print(wg.x, wg.y, wg.z)


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
