import os.path
from ctypes import sizeof

import numpy as np
from scipy.constants import N_A, c, k
from scipy.fft import next_fast_len

from radis.gpu.params import (
    init_G_params,
    init_L_params,
    init_Q,
    set_G_params,
    set_L_params,
    set_pTQ,
)
from radis.gpu.structs import initData_t, iterData_t, workGroupSizeArray_t
from radis.misc.utils import getProjectRoot

shader_path = os.path.join(getProjectRoot(), "gpu", "vulkan", "shaders")


def next_fast_len_even(n):
    # Find an FFT length that can be processed quickly.
    # The length needs to be divisible by 2 because of the
    # length doubling for prevention of circular convolution.
    n = next_fast_len(n)
    while n & 1:
        n = next_fast_len(n + 1)
    return n


from radis.gpu.vulkan.vulkan_compute_lib import GPUApplication, GPUBuffer


class gpuApp(GPUApplication):
    def __init__(
        self,
        vmin,
        Nv,
        dv,
        dxG,
        dxL,
        v0,
        da,
        na,
        S0,
        El,
        gamma_arr,
        iso,
        Mm_arr,
        Q_intp_list,
        verbose=0,
        device_id=0,
    ):
        """
        Initialize gpuApp object required for GPU calculations

        Parameters
        ----------
        vmin : float
            Minimum value frequency/wavenumber axis
        Nv : int
            Total number of frequency/wavenumber points.
        dv : float
            Stepsize of frequency/wavenumber axis (called wstep elsewhere in RADIS).
        dxG : float
            Relative grid spacing for Gaussian lineshapes.
        dxL : float
            Relative grid spacing for Lorentzian lineshapes.
        v0 : numpy.ndarray[np.float32]
            Array of line center frequencies (in cm-1).
        da : numpy.ndarray[np.float32]
            Pressure shift  (in cm-1.atm-1).
        na : numpy.ndarray[np.float32]
            Temperature dependency of Lorentzian widths
        S0 : numpy.ndarray[np.float32]
            Line intensity scaling factors.
        El : numpy.ndarray[np.float32]
            Lower level energy levels.
        gamma_arr : numpy.ndarray[np.float32]
            (m,n) shaped array with Lorentzian width parameters, with n the number of lines in the database
            and m the number of collision partners included. This is usually at least two,
            with the first (m=0) always self broadening and the last (m=-1) always air broadening.
        iso : numpy.ndarray[np.uint32]
            Index of isotopologue.
        Mm_arr : numpy.ndarray
            Molecular masses for all isotopologues in the database (Mm_arr[0] is always 0).
        Q_intp_list : list
            List of Q branch interpolators.
        verbose : bool, optional
            Print verbosity level. Default is 0.
        device_id: int, str
            Select the GPU device. If ``int``, specifies the device index, which is printed for convenience during GPU initialization with backend='vulkan' (default).
            If ``str``, return the first device that includes the specified string (case in-sesitive). If not found, return the device at index 0.
            default = 0

        """

        self.backend = "vulkan"

        # Initialize the vulkan GPUApplication:
        super().__init__(deviceID=device_id, path=shader_path, verbose=verbose)

        ## Next, the GPU is made aware of a number of parameters.
        ## Parameters that don't change during iteration are stored
        ## in init_h.

        ## When using the GPU, we distinguish objects for host (=CPU side) use (_h)
        ## and device (=GPU side) use (_d).

        ## The init_d/init_h pair is created by first making a GPUBuffer object
        ## init_d. Then that object makes the init_h object, which guarantees
        ## that data written to init_h can be transferred to init_d.

        ## init_d is a Uniform buffer, which means that the contents don't change
        ## during execution of the shader (=GPU code). This means that iter_d
        ## is also a Uniform buffer, because it only changes before/after execution
        ## of the shader.

        ## Radis configures Uniform buffers as using memory that is visible to both
        ## host (CPU) and device (GPU). For this reason, no explicit transfers between
        ## host and device memory are required.

        if verbose >= 2:
            print("Copying initialization parameters to device memory...")

        self.init_d = GPUBuffer(
            sizeof(initData_t), usage="uniform", binding=0
        )  # host_visible, device_local, (host_cached)

        self.init_h = self.init_d.getHostStructPtr(initData_t)

        self.init_h.v_min = vmin
        self.init_h.dv = dv
        self.init_h.N_v = Nv
        self.init_h.N_v_FT = next_fast_len_even(2 * self.init_h.N_v)
        self.init_h.N_x_FT = self.init_h.N_v_FT // 2 + 1
        self.init_h.dxG = dxG
        self.init_h.dxL = dxL
        self.init_h.N_lines = int(len(v0))
        self.init_h.N_coll = gamma_arr.shape[0]

        log_c2Mm_arr = np.array(
            [0]
            + [
                0.5 * np.log(8 * k * np.log(2) / (c**2 * Mm * 1e-3 / N_A))
                for Mm in Mm_arr[1:]
            ]
        )
        for i in range(len(log_c2Mm_arr)):
            self.init_h.log_c2Mm[i] = log_c2Mm_arr[i]

        init_Q(Q_intp_list)
        log_2vMm = np.log(v0) + log_c2Mm_arr.take(iso)

        init_G_params(log_2vMm.astype(np.float32), verbose)
        init_L_params(na, gamma_arr, verbose)  # TODO: do this on GPU?

        # Parameters that *do* change are stored in iter_h.
        # We only assign values to them in the iterate() method,
        # but already allocate the buffers here.

        self.iter_d = GPUBuffer(
            sizeof(iterData_t), usage="uniform", binding=1
        )  # host_visible, device_local, (host_cached)
        self.iter_h = self.iter_d.getHostStructPtr(iterData_t)

        if verbose >= 2:
            print("done!")

        if verbose >= 2:
            print("Allocating device memory and copying data...")

        ## Next the database memory is initialized. The database is potentially the largest buffer,
        ## Which could take up gigabytes of device memory. Radis tells the Vulkan API to allocate
        ## this buffer as 'storage' buffer. In dedicated GPU's, this is the largest memory heap,
        ## but usually not directly visible to the host. In order to copy memory from device to the host,
        ## internally a staging buffer is created. the staging buffer is small host memory where data is temporarily copied to,
        ## which can be used to transfer data from host to device. This transfer thus consists of two separate transfers:
        ## 1. transfer data from RAM to staging buffer, 2. transfer data from staging buffer to device memory.
        ## The staging buffer is typically smaller than the data to-be-transferred, so this process repeats until all data is copied.
        ## The staging buffer size is given by the chunksize keyword, which is only used the first time when the buffer is initialized, and ignored
        ## during consecutive calls.
        ## In the code below we copy multiple arrays to the device, so for each array the staging buffer transfers are repeated internally until each array is fully copied.

        ## By using database_d.copyToBuffer(), the transfer from staging buffer to device buffer happens immediately.
        ## Alternatively, one could split up the transfer by first using .fromArray() to copy data to the staging buffer, and then
        ## call .transferStagingBuffer() separately. We will use this later for buffers that are copied every iteration.

        database_arrays = [iso, v0, da, S0, El, na, gamma_arr]
        N_db = np.sum(
            [
                np.sum(arr.shape[:-1]) if len(arr.shape) > 1 else 1
                for arr in database_arrays
            ]
        )
        self.database_d = GPUBuffer(N_db * v0.nbytes, usage="storage", binding=2)

        byte_offset = 0
        for arr in database_arrays:
            byte_offset += self.database_d.copyToBuffer(
                arr, device_offset=byte_offset, chunksize=32 * 1024 * 1024
            )

        ## The S_klm_d buffer is purely a device buffer, which stores the result of the lineshape-distribution algorithm
        ## before the FT is applied. Due to changing N_G and N_L per iteration, this buffer may also need to change size.
        ## Because allocating new buffers is costly and may result in fragmentation of device memory, we only change the buffer size
        ## if it needs to be increased. When we do, we increase it by some margin (default is 50% extra), so as to not have to re-allocate everytime
        ## the N_G*N_L size increases. Increasing the buffer size will happen during the iteration step,
        ## but to allow for dynamic resizing we specify it as a FFT buffer by using fftSize. fftSize specifies the size of the FT, which may be
        ## smaller than the actual buffer size. In fact, in order to use in-place transforms with VkFFT, the length of the buffer is that of the reverse FT,
        ## i.e. N_v_FT//2+1 of type complex64 (8 byte), instead of N_v_FT of type float32 (4 byte). This means that the buffer size is 2 elements larger per
        ## nu-axis.

        ## The S_klm_d buffer is (implicitly) created with batchSize = 1. During iteration the batchSize will be set to N_G * N_L.

        self.S_klm_d = GPUBuffer(
            fftSize=self.init_h.N_v_FT, usage="storage", binding=3
        )  # req: large, device_local

        ## The spectrum_d buffer receives the resulting spectrum. This is a storage buffer like S_klm_d, but does not have to change size dynamically.
        ## It *does* have to transfer data from device buffers to host buffers, and since it's contents change during shader execution, we cannot use Uniform for this.
        ## To make this transfer as efficient as possible, the device-to-host (staging buffe) transfer is recorded in the command buffer (see ahead),
        ## while the staging buffer to RAM is done in python code after the comman buffer is done executing.

        self.spectrum_d = GPUBuffer(
            fftSize=self.init_h.N_v_FT, usage="storage", binding=4
        )  # req: device_local, host_visible, host_cached, (large)

        ## The indirect_d buffer is a special buffer that is used to specify the number of workgroups (i.e number of instances of a shader thread)
        ## dynamically. This is needed because the forward FFT will change size as N_G or N_L change. To address this, the FFT plan is created with
        ## the batchSize number of the highest N_G*N_L*margin(=1.5), but less threads are launched depending on the current value of N_G*N_L.
        ## We tell the gpuApp specifically that this is our indirect buffer, so that VkFFT can fill it with the workgroup sizes it would need for the batchSize
        ## it was created for. We then every iteration set the workgroup sizes to the right number, so that we don't perform superfluous FFT's.

        self.indirect_d = GPUBuffer(
            sizeof(workGroupSizeArray_t), usage="indirect", binding=10
        )  # host_visible, device_local, (host_cached)
        self.setIndirectBuffer(self.indirect_d, workGroupSizeArray_t)

        # Write command buffer:
        N_tpb = 128  # threads per block
        threads = (N_tpb, 1, 1)

        # Finally it is time to write the command buffer. The command buffer is a list of commands that the GPU needs to execute.
        # Timestamp labels can be added as desired, to help with benchmarking GPU performance. The command buffer is written at this step,
        # but will not be ran until the iterate() stage.

        # The command buffer starts with copying the new iter_d and indirect_d data from the staging buffer to the device.
        # This is however not needed explicitly, because we took care to use memory that is already device_local.
        # Then the S_klm_d and spectrum_d buffers are zeroed. The "cmdFillLDM.spv" shader is executed, which is a separate file in the shader folder.
        # This shader populates the S_klm_d array. Next the forward FFT is performed, followed by cmdApplyLineshapes.spv (another user defined shader),
        # and finally the contents of spectrum_d are transferred to the staging buffer on the host side again.

        self.appendCommands(
            [
                # self.cmdAddTimestamp("Transfer staging buffers"),
                # self.indirect_d.cmdTransferStagingBuffer("H2D"),
                # self.iter_d.cmdTransferStagingBuffer("H2D"),
                self.cmdAddTimestamp("Zeroing buffers"),
                self.S_klm_d.cmdClearBuffer(),
                self.spectrum_d.cmdClearBuffer(),
                self.cmdAddTimestamp("Line addition"),
                self.cmdScheduleShader(
                    "cmdFillLDM.spv", (self.init_h.N_lines // N_tpb + 1, 1, 1), threads
                ),
                self.cmdAddTimestamp("FFT fwd"),
                self.cmdFFT(self.S_klm_d, name="FFTa"),
                self.cmdAddTimestamp("selfly conv."),
                self.cmdScheduleShader(
                    "cmdApplyLineshapes.spv",
                    (self.init_h.N_x_FT // N_tpb + 1, 1, 1),
                    threads,
                ),
                # self.cmdScheduleShader('cmdTestApplyLineshapesP.spv', (self.init_h.N_x_FT // N_tpb + 1, N_G*N_L, 1), threads),
                self.cmdAddTimestamp("FFT inv"),
                self.cmdIFFT(self.spectrum_d, name="FFTb"),
                self.cmdAddTimestamp("Tranfer staging buffer (result)"),
                self.spectrum_d.cmdTransferStagingBuffer("D2H"),
                self.cmdAddTimestamp("End"),
            ]
        )

        self.writeCommandBuffer()

        ## During iteration, we only run .setBatchSize() to update the Vulkan buffers.
        ## Here we specify what needs to happen when .setBatchSize() is called.
        ## that is two things: the size of S_klm_d needs to be updated, and the workgroup size (indirect_d buffer contents)
        ## must be updated.

        ## If S_klm_d.setBatchSize() needs to increase the buffer size, the internal forward FFT plan (stored in ._fftAppFwd)
        ## will be removed such that it will be rewritten with the new maximum sizes. Note that this happens much less often
        ## than a change in N_G*N_L, because of the 50% extra margin. When it does happen, the performance of that particular iteration
        ## will drop, but this is only once every so often.

        self.registerBatchSizeUpdateFunction(self.S_klm_d.setBatchSize)
        self.registerBatchSizeUpdateFunction(self.setFwdFFTWorkGroupSize)

        if verbose >= 2:
            print("done!")

    def iterate(
        self,
        p,
        T,
        mole_fraction,
        verbose=0,
        # TODO: for GPU instrument functions (not currently supported):
        l=1.0,
        slit_FWHM=0.0,
    ):
        """
        Parameters
        ----------
        p : float
            pressure [bar]
        T : float
            temperature [K]
        mole_fraction : float

        Other Parameters
        ----------------
        verbose : int, optional
            The default is 0.


        Returns
        -------
        abscoeff_h : numpy.ndarray[np.float32]
            array with absorbtion coefficients in (cm.-1)
        iter_h : radis.gpu.structs.iterData_t
            structue with parameters used for computation of abscoeff_h.
        times : dict
            dictionary with computation cumulative computation times for
            different stages of the GPU computation. The ``'total'`` key
            gives the total time.
        """

        if verbose >= 2:
            print("Copying iteration parameters to device...")

        ## Update the contents of iter_h, which as a device_local buffer also
        ## immediately updates the GPU side buffer

        set_pTQ(p, T, mole_fraction, self.iter_h, l=l, slit_FWHM=slit_FWHM)
        set_G_params(self.init_h, self.iter_h)
        set_L_params(self.init_h, self.iter_h)

        if verbose >= 2:
            print("Running compute pipelines...")

        ## Update the batch size, which internally re-allocates buffers and re-plan the FFT if required,
        ## as well as update the workgroup size in the indirect_d buffer
        self.setBatchSize(self.iter_h.N_G * self.iter_h.N_L)

        ## RUN THE COMMAND BUFFER
        self.run()

        ## get the timestamp values
        gpu_times = self.get_timestamps()

        ## Copy the results into a numpy array. Remember that the device-to-host transfer is recorded in the command buffer,
        ## here we only have to copy the contents of the staging buffer to a numpy array using spectrum_d.toArray().
        abscoeff_h = np.zeros(self.init_h.N_v, dtype=np.float32)
        self.spectrum_d.toArray(abscoeff_h)

        if verbose == 1:
            print("Finished calculating spectrum!")

        return abscoeff_h, gpu_times

    def get_griddims(self):
        # use separate getter so as not to increase the app.iter_h reference count
        return self.iter_h.N_L, self.iter_h.N_G

    def __del__(self):
        # print('>>> Deleting gpuApp...')
        self.init_h = None
        self.iter_h = None
        self.free()
