# modified from https://github.com/Erkaman/vulkan_minimal_compute

# Patching vulkan lib:
import collections

collections.Iterable = collections.abc.Iterable
import ctypes
import os
from functools import partial

import numpy as np
import vulkan as vk

from radis.gpu.vulkan.pyvkfft_vulkan import prepare_fft

QUERY_POOL_SIZE = 32  # Max number of queries (=timestamps)


class GPUApplication(object):
    def __init__(self, deviceID=0, path="./", verbose=True):
        # In order to use Vulkan, you must create an instance
        self._instance = None
        self._debugReportCallback = None
        self._physicalDevice = None
        self._device = None
        self._nextDstBinding = 0
        self._timestampLabels = []
        self._bufferObjects = []
        self.verbose = verbose
        # self._deviceID = deviceID
        self._shaderPath = path
        self._fftApps = {}

        self._computeShaderModules = []
        self._descriptorPools = []
        self._descriptorSetLayouts = []
        self._pipelineLayouts = []
        self._pipelines = []

        self._commandPool = None
        self._commandBuffer = None
        self.command_list = []
        self._queryPool = None

        self._enabledLayers = []
        self._enabledExtensions = []

        self._queue = None
        self._queueFamilyIndex = -1

        # Initialize vulkan
        self.createInstance()
        self.selectPhysicalDevice(deviceID, verbose=verbose)
        self.createDevice()
        self.init_shaders()

        self.createCommandBuffer()

    def free(self):

        self.command_list = []

        fft_keys = [*self._fftApps.keys()]
        for key in fft_keys:
            del self._fftApps[key]

        for bufferObject in self._bufferObjects:
            bufferObject.free()

        for computeShaderModule in self._computeShaderModules:
            vk.vkDestroyShaderModule(self._device, computeShaderModule, None)
        for descriptorPool in self._descriptorPools:
            vk.vkDestroyDescriptorPool(self._device, descriptorPool, None)
        for descriptorSetLayout in self._descriptorSetLayouts:
            vk.vkDestroyDescriptorSetLayout(self._device, descriptorSetLayout, None)
        for pipelineLayout in self._pipelineLayouts:
            vk.vkDestroyPipelineLayout(self._device, pipelineLayout, None)
        for pipeline in self._pipelines:
            vk.vkDestroyPipeline(self._device, pipeline, None)

        # TODO: [developers-note]: VkFFT used to free some resources, which has now all moved to here. The two freeings below have been removed from VkFFT and not been replaced here:
        # if (config->physicalDevice) free(config->physicalDevice);
        # if (config->queue) free(config->queue);

        if self._queryPool:
            vk.vkDestroyQueryPool(self._device, self._queryPool, None)
        if self._commandPool:
            vk.vkDestroyCommandPool(self._device, self._commandPool, None)
        if self._fence:
            vk.vkDestroyFence(self._device, self._fence, None)
        if self._device:
            vk.vkDestroyDevice(self._device, None)
        if self._instance:
            vk.vkDestroyInstance(self._instance, None)

    def __setattr__(self, name, val):
        if isinstance(val, GPUObject):
            val.app = self
            val._init_buffer()
            val._delayedSetData()
        self.__dict__[name] = val

    def init_shaders(self):
        shader_fnames = [f for f in os.listdir(self._shaderPath) if f[-3:] == "spv"]
        for shader_fname in shader_fnames:
            fun_name = shader_fname.split(".")[0]  # TODO: do this with os.path.basename
            named_shader = partial(self.schedule_shader, shader_fname)
            setattr(self.__class__, fun_name, named_shader)

    def schedule_shader(
        self,
        shader_fname=None,
        global_workgroup=(1, 1, 1),
        local_workgroup=(1, 1, 1),
        sync=True,
        timestamp=False,
    ):
        def func(
            self,
            shader_fname=None,
            global_workgroup=(1, 1, 1),
            local_workgroup=(1, 1, 1),
            sync=True,
            timestamp=False,
        ):

            shader_fpath = os.path.join(self._shaderPath, shader_fname)

            # Descriptor set:
            descriptorSetLayout = self.createDescriptorSetLayout()
            descriptorPool, descriptorSet = self.createDescriptorSet(
                descriptorSetLayout
            )

            # Compute Pipeline:
            pipeline, pipelineLayout, computeShaderModule = self.createComputePipeline(
                shader_fpath, local_workgroup, descriptorSetLayout
            )
            self.bindAndDispatch(
                global_workgroup, descriptorSet, pipeline, pipelineLayout
            )
            if sync:
                self.sync()

            self._computeShaderModules.append(computeShaderModule)
            self._descriptorPools.append(descriptorPool)
            self._descriptorSetLayouts.append(descriptorSetLayout)
            self._pipelineLayouts.append(pipelineLayout)
            self._pipelines.append(pipeline)

            if timestamp == True:
                self.cmdAddTimestamp(shader_fname.split(".")[0]).writeCommand()
            elif isinstance(timestamp, str):
                self.cmdAddTimestamp(timestamp).writeCommand()
            else:
                pass

        return GPUCommand(
            func,
            [
                self,
            ],
            {
                "shader_fname": shader_fname,
                "global_workgroup": global_workgroup,
                "local_workgroup": local_workgroup,
                "sync": sync,
                "timestamp": timestamp,
            },
        )

    def cmdFFT(self, buf, buf_FT, timestamp=False):
        def func(self, buf, buf_FT, timestamp=False):
            key = (id(buf), id(buf_FT))

            try:
                fft_app = self._fftApps[key]
            except (KeyError):
                fft_app = prepare_fft(buf, buf_FT, compute_app=self)
                self._fftApps[key] = fft_app

            fft_app.fft(self._commandBuffer, buf._buffer, buf_FT._buffer)

            if timestamp == True:
                self.cmdAddTimestamp("fft").writeCommand()
            elif isinstance(timestamp, str):
                self.cmdAddTimestamp(timestamp).writeCommand()
            else:
                pass

        return GPUCommand(func, [self, buf, buf_FT], {"timestamp": timestamp})

    def cmdIFFT(self, buf_FT, buf, timestamp=False):
        def func(self, buf_FT, buf, timestamp=False):
            key = (id(buf), id(buf_FT))

            try:
                fft_app = self._fftApps[key]
            except (KeyError):
                fft_app = prepare_fft(buf, buf_FT, compute_app=self)
                self._fftApps[key] = fft_app

            fft_app.ifft(self._commandBuffer, buf_FT._buffer, buf._buffer)

            if timestamp == True:
                self.cmdAddTimestamp("fft").writeCommand()
            elif isinstance(timestamp, str):
                self.cmdAddTimestamp(timestamp).writeCommand()
            else:
                pass

        return GPUCommand(func, [self, buf_FT, buf], {"timestamp": timestamp})

    def run(self):
        self.runCommandBuffer()

    @staticmethod
    def debugReportCallbackFn(*args):
        print("Debug Report: {} {}".format(args[5], args[6]))
        return 0

    def createInstance(self):
        # Next, we actually create the instance.

        # Contains application info. This is actually not that important.
        # The only real important field is apiVersion.
        applicationInfo = vk.VkApplicationInfo(
            sType=vk.VK_STRUCTURE_TYPE_APPLICATION_INFO,
            pApplicationName="RADIS GPU app",
            applicationVersion=0,
            pEngineName="RADIS Vulkan engine",
            engineVersion=0,
            apiVersion=vk.VK_API_VERSION_1_0,
        )

        createInfo = vk.VkInstanceCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
            flags=0,
            pApplicationInfo=applicationInfo,
            # Give our desired layers and extensions to vulkan.
            enabledLayerCount=len(self._enabledLayers),
            ppEnabledLayerNames=self._enabledLayers,
            enabledExtensionCount=len(self._enabledExtensions),
            ppEnabledExtensionNames=self._enabledExtensions,
        )

        # Actually create the instance.
        # Having created the instance, we can actually start using vulkan.
        self._instance = vk.vkCreateInstance(createInfo, None)

    def selectPhysicalDevice(self, deviceID, verbose=True):
        devices = vk.vkEnumeratePhysicalDevices(self._instance)

        if isinstance(deviceID, str):
            test_str = deviceID.upper()
            for deviceID, device in enumerate(devices):
                props = vk.vkGetPhysicalDeviceProperties(device)
                devname = vk.ffi.string(props.obj.deviceName).decode()
                if test_str in devname.upper():
                    break
            else:
                deviceID = 0
        else:
            if deviceID < 0:
                deviceID += len(devices)
            if not 0 <= deviceID < len(devices):
                deviceID = 0

        self._deviceID = deviceID
        self._physicalDevice = devices[self._deviceID]

        if verbose:
            print("Vulkan version: ", vk.__version__)
            print("Selected card (deviceID={:d}):".format(self._deviceID))
            for i, device in enumerate(devices):
                props = vk.vkGetPhysicalDeviceProperties(device)
                devname = vk.ffi.string(props.obj.deviceName).decode()
                print(
                    "[{:s}] {:d}: {:s}".format(
                        "X" if i == self._deviceID else " ", i, devname
                    )
                )
            print("")

    # Returns the index of a queue family that supports compute operations.
    def getComputeQueueFamilyIndex(self):
        # Retrieve all queue families.
        queueFamilies = vk.vkGetPhysicalDeviceQueueFamilyProperties(
            self._physicalDevice
        )

        # Now find a family that supports compute.
        for i, props in enumerate(queueFamilies):
            if props.queueCount > 0 and props.queueFlags & vk.VK_QUEUE_COMPUTE_BIT:
                # found a queue with compute. We're done!
                return i

        return -1

    def createDevice(self):
        # We create the logical device in this function.

        self._queueFamilyIndex = self.getComputeQueueFamilyIndex()
        # When creating the device, we also specify what queues it has.
        queueCreateInfo = vk.VkDeviceQueueCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO,
            queueFamilyIndex=self._queueFamilyIndex,  # find queue family with compute capability.
            queueCount=1,  # create one queue in this family. We don't need more.
            pQueuePriorities=[
                1.0
            ],  # we only have one queue, so this is not that imporant.
        )

        # Now we create the logical device. The logical device allows us to interact with the physical device.
        # Specify any desired device features here. We do not need any for this application, though.
        deviceFeatures = vk.VkPhysicalDeviceFeatures()
        deviceCreateInfo = vk.VkDeviceCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO,
            enabledLayerCount=len(self._enabledLayers),
            ppEnabledLayerNames=self._enabledLayers,
            pQueueCreateInfos=queueCreateInfo,
            queueCreateInfoCount=1,
            pEnabledFeatures=deviceFeatures,
        )

        self._device = vk.vkCreateDevice(self._physicalDevice, deviceCreateInfo, None)
        self._queue = vk.vkGetDeviceQueue(self._device, self._queueFamilyIndex, 0)

        # We create a fence.
        fenceCreateInfo = vk.VkFenceCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_FENCE_CREATE_INFO, flags=0
        )
        self._fence = vk.vkCreateFence(self._device, fenceCreateInfo, None)
        self._memoryBarrier = vk.VkMemoryBarrier(
            sType=vk.VK_STRUCTURE_TYPE_MEMORY_BARRIER,
            pNext=0,
            srcAccessMask=vk.VK_ACCESS_SHADER_WRITE_BIT,
            dstAccessMask=vk.VK_ACCESS_SHADER_READ_BIT,
        )

        queryPoolCreateInfo = vk.VkQueryPoolCreateInfo(
            vk.VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO,
            None,
            0,
            vk.VK_QUERY_TYPE_TIMESTAMP,
            QUERY_POOL_SIZE,
            0,
        )

        self._queryPool = vk.vkCreateQueryPool(self._device, queryPoolCreateInfo, None)

    def createCommandBuffer(self):
        # We are getting closer to the end. In order to send commands to the device(GPU),
        # we must first record commands into a command buffer.
        # To allocate a command buffer, we must first create a command pool. So let us do that.
        commandPoolCreateInfo = vk.VkCommandPoolCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,
            flags=0,
            # the queue family of this command pool. All command buffers allocated from this command pool,
            # must be submitted to queues of this family ONLY.
            queueFamilyIndex=self._queueFamilyIndex,
        )

        self._commandPool = vk.vkCreateCommandPool(
            self._device, commandPoolCreateInfo, None
        )

        # Now allocate a command buffer from the command pool.
        commandBufferAllocateInfo = vk.VkCommandBufferAllocateInfo(
            sType=vk.VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
            commandPool=self._commandPool,
            # if the command buffer is primary, it can be directly submitted to queues.
            # A secondary buffer has to be called from some primary command buffer, and cannot be directly
            # submitted to a queue. To keep things simple, we use a primary command buffer.
            level=vk.VK_COMMAND_BUFFER_LEVEL_PRIMARY,
            commandBufferCount=1,
        )

        self._commandBuffer = vk.vkAllocateCommandBuffers(
            self._device, commandBufferAllocateInfo
        )[0]

    def writeCommandBuffer(self):

        # Now we shall start recording commands into the newly allocated command buffer.
        beginInfo = vk.VkCommandBufferBeginInfo(
            sType=vk.VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
            flags=0,  # VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT #VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT
        )
        vk.vkBeginCommandBuffer(self._commandBuffer, beginInfo)
        for obj in self.command_list:
            obj.writeCommand()

        vk.vkEndCommandBuffer(self._commandBuffer)

    # find memory type with desired properties.
    def findMemoryType(self, memoryTypeBits, properties):
        memoryProperties = vk.vkGetPhysicalDeviceMemoryProperties(self._physicalDevice)

        # How does this search work?
        # See the documentation of VkPhysicalDeviceMemoryProperties for a detailed description.
        for i, mt in enumerate(memoryProperties.memoryTypes):
            if (
                memoryTypeBits & (1 << i)
                and (mt.propertyFlags & properties) == properties
            ):
                return i

        return -1

    def createDescriptorSetLayout(self):
        # Here we specify a descriptor set layout. This allows us to bind our descriptors to
        # resources in the shader.

        # Here we specify a binding of type VK_DESCRIPTOR_TYPE_STORAGE_BUFFER to the binding point
        # 0. This binds to
        #   layout(std140, binding = 0) buffer buf
        # in the compute shader.
        descriptorSetLayoutBindings = [
            vk.VkDescriptorSetLayoutBinding(
                binding=bufferObject._dstBinding,
                descriptorType=bufferObject._descriptorType,
                descriptorCount=1,
                stageFlags=vk.VK_SHADER_STAGE_COMPUTE_BIT,
            )
            for bufferObject in self._bufferObjects
        ]

        descriptorSetLayoutCreateInfo = vk.VkDescriptorSetLayoutCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
            bindingCount=len(descriptorSetLayoutBindings),
            pBindings=descriptorSetLayoutBindings,
        )

        # Create the descriptor set layout.
        return vk.vkCreateDescriptorSetLayout(
            self._device, descriptorSetLayoutCreateInfo, None
        )

    def createDescriptorSet(self, descriptorSetLayout):
        # So we will allocate a descriptor set here.
        # But we need to first create a descriptor pool to do that.

        # Our descriptor pool can only allocate a single storage buffer.
        descriptorPoolSize = vk.VkDescriptorPoolSize(
            type=vk.VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            descriptorCount=len(self._bufferObjects),
        )

        descriptorPoolCreateInfo = vk.VkDescriptorPoolCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
            maxSets=1,  # we only need to allocate one descriptor set from the pool.
            poolSizeCount=1,
            pPoolSizes=descriptorPoolSize,
        )

        # create descriptor pool.
        descriptorPool = vk.vkCreateDescriptorPool(
            self._device, descriptorPoolCreateInfo, None
        )

        # With the pool allocated, we can now allocate the descriptor set.
        descriptorSetAllocateInfo = vk.VkDescriptorSetAllocateInfo(
            sType=vk.VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
            descriptorPool=descriptorPool,
            descriptorSetCount=1,
            pSetLayouts=[descriptorSetLayout],
        )

        # allocate descriptor set.
        descriptorSet = vk.vkAllocateDescriptorSets(
            self._device, descriptorSetAllocateInfo
        )[0]

        # Next, we need to connect our actual storage buffer with the descrptor.
        # We use vkUpdateDescriptorSets() to update the descriptor set.

        writeDescriptorSets = [
            vk.VkWriteDescriptorSet(
                sType=vk.VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
                dstSet=descriptorSet,
                dstBinding=bufferObject._dstBinding,
                descriptorCount=1,
                descriptorType=bufferObject._descriptorType,
                pBufferInfo=vk.VkDescriptorBufferInfo(
                    buffer=bufferObject._buffer,
                    offset=0,
                    range=bufferObject._bufferSize,
                ),
            )
            for bufferObject in self._bufferObjects
        ]

        # perform the update of the descriptor set.
        vk.vkUpdateDescriptorSets(
            self._device, len(writeDescriptorSets), writeDescriptorSets, 0, None
        )
        return descriptorPool, descriptorSet

    def createComputePipeline(
        self, shaderFileName, localWorkGroup, descriptorSetLayout
    ):
        # We create a compute pipeline here.

        # Create a shader module. A shader module basically just encapsulates some shader code.
        with open(shaderFileName, "rb") as comp:
            code = comp.read()

            createInfo = vk.VkShaderModuleCreateInfo(
                sType=vk.VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
                codeSize=len(code),
                pCode=code,
            )

            computeShaderModule = vk.vkCreateShaderModule(
                self._device, createInfo, None
            )

        # specialize local workgroup size:
        dtype = "int"
        buffer = vk.ffi.new(dtype + "[]", localWorkGroup)
        pData = vk.ffi.cast("void*", buffer)

        dsize = vk.ffi.sizeof(dtype)
        entries = [
            vk.VkSpecializationMapEntry(constantID=i, offset=i * dsize, size=dsize)
            for i in range(len(buffer))
        ]

        specializationInfo = vk.VkSpecializationInfo(
            mapEntryCount=len(entries),
            pMapEntries=entries,
            dataSize=vk.ffi.sizeof(buffer),
            pData=pData,
        )

        # Now let us actually create the compute pipeline.
        # A compute pipeline is very simple compared to a graphics pipeline.
        # It only consists of a single stage with a compute shader.
        # So first we specify the compute shader stage, and it's entry point(main).
        shaderStageCreateInfo = vk.VkPipelineShaderStageCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
            stage=vk.VK_SHADER_STAGE_COMPUTE_BIT,
            module=computeShaderModule,
            pSpecializationInfo=specializationInfo,
            pName="main",
        )

        # The pipeline layout allows the pipeline to access descriptor sets.
        # So we just specify the descriptor set layout we created earlier.
        pipelineLayoutCreateInfo = vk.VkPipelineLayoutCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
            setLayoutCount=1,
            pSetLayouts=[descriptorSetLayout],
        )
        pipelineLayout = vk.vkCreatePipelineLayout(
            self._device, pipelineLayoutCreateInfo, None
        )

        pipelineCreateInfo = vk.VkComputePipelineCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO,
            stage=shaderStageCreateInfo,
            layout=pipelineLayout,
        )

        # Now, we finally create the compute pipeline.
        pipelines = vk.vkCreateComputePipelines(
            self._device, vk.VK_NULL_HANDLE, 1, pipelineCreateInfo, None
        )
        # print(self._device, pipelines)
        if len(pipelines) == 1:
            pipeline = pipelines[0]
        else:
            raise Exception("Could not create compute pipeline")

        return pipeline, pipelineLayout, computeShaderModule

    def bindAndDispatch(self, globalWorkGroup, descriptorSet, pipeline, pipelineLayout):
        # We need to bind a pipeline, AND a descriptor set before we dispatch.
        # The validation layer will NOT give warnings if you forget these, so be very careful not to forget them.
        vk.vkCmdBindPipeline(
            self._commandBuffer, vk.VK_PIPELINE_BIND_POINT_COMPUTE, pipeline
        )
        vk.vkCmdBindDescriptorSets(
            self._commandBuffer,
            vk.VK_PIPELINE_BIND_POINT_COMPUTE,
            pipelineLayout,
            0,
            1,
            [descriptorSet],
            0,
            None,
        )

        # Calling vkCmdDispatch basically starts the compute pipeline, and executes the compute shader.
        # The number of workgroups is specified in the arguments.
        # If you are already familiar with compute shaders from OpenGL, this should be nothing new to you.
        vk.vkCmdDispatch(self._commandBuffer, *globalWorkGroup)

    # def endCommandBuffer(self):
    # vk.vkEndCommandBuffer(self._commandBuffer)

    def runCommandBuffer(self):
        # Now we shall finally submit the recorded command buffer to a queue.
        submitInfo = vk.VkSubmitInfo(
            sType=vk.VK_STRUCTURE_TYPE_SUBMIT_INFO,
            commandBufferCount=1,  # submit a single command buffer
            pCommandBuffers=[self._commandBuffer],  # the command buffer to submit.
        )

        # We submit the command buffer on the queue, at the same time giving a fence.
        vk.vkQueueSubmit(self._queue, 1, submitInfo, self._fence)

        # The command will not have finished executing until the fence is signalled.
        # So we wait here.
        # We will directly after this read our buffer from the GPU,
        # and we will not be sure that the command has finished executing unless we wait for the fence.
        # Hence, we use a fence here.
        vk.vkWaitForFences(self._device, 1, [self._fence], vk.VK_TRUE, 100000000000)
        vk.vkResetFences(self._device, 1, [self._fence])

    def cmdClearBuffer(self, buffer_obj, timestamp=False):
        def func(self, buffer_obj, timestamp=False):
            vk.vkCmdFillBuffer(
                self._commandBuffer, buffer_obj._buffer, 0, buffer_obj.nbytes, 0
            )
            if timestamp == True:
                self.cmdAddTimestamp("clearBuffer").writeCommand()
            elif isinstance(timestamp, str):
                self.cmdAddTimestamp(timestamp).writeCommand()
            else:
                pass

        return GPUCommand(func, [self, buffer_obj], {"timestamp": timestamp})

    def sync(self):
        vk.vkCmdPipelineBarrier(
            self._commandBuffer,
            vk.VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            vk.VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            0,
            1,
            self._memoryBarrier,
            0,
            0,
            0,
            0,
        )

    def cmdAddTimestamp(self, label=None):
        def func(self, label=None):
            query = len(self._timestampLabels)
            if label is None:
                label = "timestamp_{:d}".format(query)
            # TODO: deal with duplicates
            self._timestampLabels.append(label)

            if query == 0:
                vk.vkCmdResetQueryPool(
                    self._commandBuffer, self._queryPool, 0, QUERY_POOL_SIZE
                )

            vk.vkCmdWriteTimestamp(
                self._commandBuffer,
                vk.VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                self._queryPool,
                query,
            )

            if query == QUERY_POOL_SIZE + 1:
                print("Max pool size reached!!!")
                # TODO: handle exception

        return GPUCommand(
            func,
            [
                self,
            ],
            {"label": label},
        )

    def get_timestamps(self):

        dtype = "unsigned int"
        dsize = vk.ffi.sizeof(dtype)
        length = len(self._timestampLabels)
        queryResult = vk.ffi.new(dtype + "[{:d}]".format(length))
        pData = vk.ffi.cast("void*", queryResult)

        vk.vkGetQueryPoolResults(
            self._device, self._queryPool, 0, length, dsize * length, pData, dsize, 0
        )

        result_arr = np.array([*queryResult]) * 1e-6  # in ms
        result_arr -= result_arr[0]
        result_dict = dict(zip(self._timestampLabels, result_arr))
        result_dict["total"] = result_arr[-1]

        return result_dict


class GPUCommand:
    def __init__(self, func, vargs, kwargs):
        self.func = func
        self.vargs = vargs
        self.kwargs = kwargs

    def writeCommand(self):
        self.func(*self.vargs, **self.kwargs)


class GPUObject:
    def __init__(self, bufferSize=0, uniform=False, binding=None, app=None):
        # customization:

        if uniform:
            self._usage = vk.VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT
            self._descriptorType = vk.VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER
        else:
            self._usage = vk.VK_BUFFER_USAGE_STORAGE_BUFFER_BIT
            self._descriptorType = vk.VK_DESCRIPTOR_TYPE_STORAGE_BUFFER

        self._bufferSize = bufferSize
        self._dstBinding = binding

        self._isInitialized = False
        self._delayedSetDataList = []

        self.app = app
        if self.app is not None:
            self._init_buffer()

    def _init_buffer(self):
        self._device = self.app._device

        if self._dstBinding is None:
            self._dstBinding = self.app._nextDstBinding
        self.app._nextDstBinding = self._dstBinding + 1

        self._buffer = None
        self._bufferMemory = None
        self._pmappedMemory = None

        # create buffer:
        bufferCreateInfo = vk.VkBufferCreateInfo(
            sType=vk.VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO,
            size=self._bufferSize,  # buffer size in bytes.
            usage=self._usage,
            sharingMode=vk.VK_SHARING_MODE_EXCLUSIVE,  # buffer is exclusive to a single queue family at a time.
        )

        self._buffer = vk.vkCreateBuffer(self._device, bufferCreateInfo, None)
        memoryRequirements = vk.vkGetBufferMemoryRequirements(
            self._device, self._buffer
        )

        index = self.app.findMemoryType(
            memoryRequirements.memoryTypeBits,
            vk.VK_MEMORY_PROPERTY_HOST_COHERENT_BIT
            | vk.VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT
            | vk.VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT,
        )

        allocateInfo = vk.VkMemoryAllocateInfo(
            sType=vk.VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,
            allocationSize=memoryRequirements.size,
            memoryTypeIndex=index,
        )

        self._bufferMemory = vk.vkAllocateMemory(self._device, allocateInfo, None)
        vk.vkBindBufferMemory(self._device, self._buffer, self._bufferMemory, 0)

        self._pmappedMemory = vk.vkMapMemory(
            self._device, self._bufferMemory, 0, self._bufferSize, 0
        )
        self.app._bufferObjects.append(self)

        self._isInitialized = True

    def _delayedSetData(self):
        while len(self._delayedSetDataList):
            vargs = self._delayedSetDataList.pop(0)
            self.setData(*vargs)

    def free(self):
        if self._pmappedMemory:
            vk.vkUnmapMemory(self._device, self._bufferMemory)
            self._pmappedMemory = None
        if self._bufferMemory:
            vk.vkFreeMemory(self._device, self._bufferMemory, None)
        if self._buffer:
            vk.vkDestroyBuffer(self._device, self._buffer, None)


class GPUArray(GPUObject):
    def __init__(
        self, shape=(1,), dtype=np.int32, strides=None, binding=None, app=None
    ):
        self.shape = shape
        self.dtype = np.dtype(dtype)
        self.itemsize = self.dtype.itemsize
        self.size = int(np.prod(self.shape))
        self.nbytes = self.size * self.itemsize

        if strides is None:
            self._calcStrides()
        else:
            self.strides = strides

        super().__init__(
            bufferSize=self.size * self.itemsize,
            uniform=False,
            binding=binding,
            app=app,
        )

        self._arr = None
        if self._isInitialized:
            self._arr = np.frombuffer(self._pmappedMemory, self.dtype).reshape(
                self.shape
            )

    def _init_buffer(self):
        super()._init_buffer()
        if self._arr is None:
            self._arr = np.frombuffer(self._pmappedMemory, self.dtype).reshape(
                self.shape
            )

    @staticmethod
    def fromArr(arr, binding=None, app=None):
        bufferObject = GPUArray(
            arr.shape, arr.dtype, arr.strides, binding=binding, app=app
        )
        bufferObject.setData(arr)
        return bufferObject

    def _calcStrides(self, order="c"):

        sshape = np.zeros_like(self.shape)
        sshape[
            0
        ] = 1  # TODO: Check this, should maybe be sshape[-1]=1; shape[:-1] = self.shape[:-1]??
        sshape[1:] = self.shape[:-1]
        sshape *= self.itemsize

        if order == "c":
            self.strides = np.multiply.accumulate(sshape[::-1])[::-1]
        else:
            self.strides = np.multiply.accumulate(sshape)

    def setData(self, arr, byte_offset=0):
        if self._isInitialized:
            ctypes.memmove(
                self._arr.ctypes.data + byte_offset, arr.ctypes.data, arr.nbytes
            )
        else:
            self._delayedSetDataList.append((arr, byte_offset))

        return arr.nbytes

    def getData(self):
        return self._arr


class GPUStruct(GPUObject):
    def __init__(self, bufferSize=0, binding=None, app=None):
        super().__init__(bufferSize=bufferSize, uniform=True, binding=binding, app=app)

    def fromStruct(struct, binding=None, app=None):
        bufferSize = ctypes.sizeof(struct)
        bufferObject = GPUStruct(bufferSize=bufferSize, binding=binding, app=app)
        bufferObject.setData(struct)
        return bufferObject

    # def fromStruct(structPtr, binding=None, app=None):
    # bufferSize = ffi.sizeof(structPtr[0])
    # bufferObject = GPUStruct(bufferSize=bufferSize, binding=binding, app=app)
    # bufferObject.setData(structPtr)
    # return bufferObject

    # def setData(self, structPtr):
    # ffi.memmove(self._pmappedMemory, structPtr, self._bufferSize)

    def setData(self, struct):
        if self._isInitialized:
            ptr = ctypes.c_void_p(
                int(
                    vk.ffi.cast(
                        "unsigned long long", vk.ffi.from_buffer(self._pmappedMemory)
                    )
                )
            )
            ctypes.memmove(ptr, ctypes.byref(struct), self._bufferSize)
        else:
            self._delayedSetDataList.append((struct,))

    def getData(self):
        # not implemented
        pass
