
import numpy as np
import os.path

from radis.misc.utils import getProjectRoot
from radis.gpu.vulkan.vulkan_compute_lib import GPUApplication, GPUBuffer

data_str = """ 55    2.2487642.700E-164 8.704E-07.07970.08664763.01990.76-.000265             41             41                    R  0      457665 5 8 2 2 1 7     6.0    2.0
 55    2.2797856.503E-162 8.685E-07.07970.08663631.47680.76-.000265             40             40                    R  0      457665 5 8 2 2 1 7     6.0    2.0
 56    2.2880133.667E-166 9.362E-07.07970.08665317.82150.76-.000265             41             41                    R  0      457665 5 8 2 2 1 7    36.0   12.0
 55    2.3106831.737E-159 8.638E-07.07970.08662478.03720.76-.000265             39             39                    R  0      457665 5 8 2 2 1 7     6.0    2.0
 56    2.3202598.961E-164 9.357E-07.07970.08664183.55380.76-.000265             40             40                    R  0      457665 5 8 2 2 1 7    36.0   12.0
 53    2.3253851.473E-164 1.002E-06.07970.08665843.01970.76-.000265             41             41                    R  0      457665 5 8 2 2 1 7     3.0    1.0
 52    2.3316165.471E-164 1.014E-06.07970.08665930.07650.76-.000265             41             41                    R  0      457665 5 8 2 2 1 7     6.0    2.0
 55    2.3414745.145E-157 8.563E-07.07970.08661302.71350.76-.000265             38             38                    R  0      457665 5 8 2 2 1 7     6.0    2.0
 56    2.3523682.435E-161 9.321E-07.07970.08663026.83920.76-.000265             39             39                    R  0      457665 5 8 2 2 1 7    36.0   12.0
 53    2.3588383.645E-162 1.003E-06.07970.08664706.40710.76-.000265             40             40                    R  0      457665 5 8 2 2 1 7     3.0    1.0
"""

fname = '05_HITEMP2019.par'

with open(fname, 'rb') as f:
    buf = f.read()

arr = np.frombuffer(buf, dtype=np.ubyte)


shader_path = os.path.join(getProjectRoot(), "gpu", "vulkan", "shaders")
app = GPUApplication(deviceID=0, path=shader_path, verbose=True)

app.buf_d = GPUBuffer(bufferSize=arr.nbytes, binding=0)
app.buf_d.fromArray(arr)

app.output_d = GPUBuffer(128*4*4, binding=1)


N_tpb = 128  # threads per block
threads = (N_tpb, 1, 1)

app.appendCommands(
    [
        app.cmdAddTimestamp("Parsing"),
        app.cmdScheduleShader(
            "cmdFillLDM.spv", (4, 1, 1), threads
        ),
        app.cmdAddTimestamp("End"),
    ]
)

app.writeCommandBuffer()
app.run()



