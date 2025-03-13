
import numpy as np
import os.path
import sys

from radis.misc.utils import getProjectRoot
from radis.gpu.vulkan.vulkan_compute_lib import GPUApplication, GPUBuffer

data_str = """
 62    0.031765 5.400E-34 5.865E-17.05130.072  690.05010.750.000000    0 0 0 0 1A1    0 0 0 0 1A1   11F1  3        11F2  3     5323331111 4 1 1 1   138.0  138.0
 62    0.034993 4.310E-34 1.558E-16.04850.071  950.53240.750.000000    0 0 0 0 1A1    0 0 0 0 1A1   13F1  4        13F2  3     5323331111 4 1 1 1   162.0  162.0
 62    0.041964 1.790E-33 6.428E-17.05510.078  470.85300.750.000000    0 0 0 0 1A1    0 0 0 0 1A1    9A2  1         9A1  1     5323331111 4 1 1 1   190.0  190.0
 62    0.054618 2.520E-33 4.772E-16.04940.077  815.12760.750.000000    0 0 0 0 1A1    0 0 0 0 1A1   12A1  2        12A2  1     5323331111 4 1 1 1   250.0  250.0
 62    0.059566 1.150E-33 1.602E-15.04180.075 1252.06160.750.000000    0 0 0 0 1A1    0 0 0 0 1A1   15A2  2        15A1  1     5323331111 4 1 1 1   310.0  310.0
 62    0.085666 6.600E-34 8.557E-17.04940.070  376.75320.750.000000    0 0 0 0 1A1    0 0 0 0 1A1    8E  2          8E   1     5323331111 4 1 1 1    68.0   68.0
 62    0.092576 5.290E-34 4.942E-17.05510.074  376.75130.750.000000    0 0 0 0 1A1    0 0 0 0 1A1    8F2  2         8F1  1     5323331111 4 1 1 1   102.0  102.0
 62    0.101145 1.460E-33 3.165E-16.05130.071  575.21130.750.000000    0 0 0 0 1A1    0 0 0 0 1A1   10F2  3        10F1  1     5323331111 4 1 1 1   126.0  126.0
 62    0.112260 5.160E-34 1.130E-15.04460.066 1095.92110.750.000000    0 0 0 0 1A1    0 0 0 0 1A1   14F2  3        14F1  2     5323331111 4 1 1 1   174.0  174.0
 61    0.117500 1.555E-40 0.000E+00.05190.068 6553.60760.66-.000000                                                  el        003333715947141313     0.0    0.0
"""



# str_list = [
#     '0.031765',
#     '0.034993',
#     '0.041964',
#     '0.054618',
#     '0.059566',
#     '0.085666',
#     ]


# for line in str_list:
#     res = np.zeros(6, dtype=np.float32)
#     for i, s in enumerate(line[2:]):
#         res[i] = np.float32(ord(s) - ord('0'))*np.float32(10**-i)
#     print(np.sum(res))

# sys.exit()

line_length = 161

# fname = '05_HITEMP2019.par'
fname = '06_HITEMP2020.par'
# fname = '06_HITEMP2020_small.par'

# sname = '06_HITEMP2020_small.par'
# with open(fname, 'rb') as f, open(sname, 'wb') as fw:
#     fw.write(f.read(10_000*line_length))


print('Initializing GPU... ',end='')
shader_path = os.path.join(getProjectRoot(), "gpu", "vulkan", "shaders")
app = GPUApplication(deviceID=0, path=shader_path, verbose=True)
print('Done!')


N_lines_read = 10_000_000
nbytes = N_lines_read * line_length
print('Loading limited to {:d} MB'.format(nbytes//1024//1024))

print('Loading data... ', end='')
arr = np.fromfile(fname, dtype=np.ubyte, count=N_lines_read*line_length)
N_lines = arr.nbytes // line_length

app.buf_d = GPUBuffer(bufferSize=arr.nbytes, binding=0)
app.buf_d.fromArray(arr)


app.output_d = GPUBuffer(bufferSize=N_lines*4, binding=1)

print('Done!')



N_tpb = 128  # threads per block
threads = (N_tpb, 1, 1)


print('Writing command buffer... ',end='')
app.appendCommands(
    [
        app.cmdAddTimestamp("Parsing"),
        app.cmdScheduleShader(
            "cmdParseHitemp.spv", (N_lines // N_tpb + 1, 1, 1), threads
        ),
        app.cmdAddTimestamp("End"),
    ]
)
app.writeCommandBuffer()
print('Done!')

print('Running app... ', end='')
app.run()
times = app.get_timestamps()
print('Done! (in {:.3f} ms)'.format(times['Total']))

arr_out = np.zeros(N_lines, dtype=np.float32)

app.output_d.toArray(arr_out)

app.free()
app = None

for a in arr_out[:5]: 
    print(a)
print('...')
for a in arr_out[-5:]: 
    print(a)







