from ctypes import c_float, c_int

import matplotlib.pyplot as plt
import numpy as np

# from driver import CuArray, CuContext, CuFFT, CuModule
from emulate import CuArray, CuContext, CuFFT, CuModule
from numpy.fft import rfftfreq
from numpy.random import rand, randint, seed
from scipy.fft import next_fast_len

L = lambda t, w: 2 / (w * np.pi) * 1 / (1 + 4 * (t / w) ** 2)
L_FT = lambda f, w: np.exp(-np.pi * np.abs(f) * w)

t_max = 100.0
Nt = 100005
Nt = next_fast_len(Nt)
print("Nt = {:d}".format(Nt))
t_arr = np.linspace(0, t_max, Nt)
dt = t_arr[1]
f_arr = rfftfreq(Nt, dt)
Nf = len(f_arr)
Ntpb = 1020  # threads per block
wL = 0.2
seed(0)
Nl = 20


def mock_spectrum(Nt, Nl):
    line_index = randint(0, Nt, Nl)
    line_strength = rand(Nl)
    I_arr = np.zeros(Nt, dtype=np.float32)
    np.add.at(I_arr, line_index, line_strength)
    return I_arr


I_arr = mock_spectrum(Nt, Nl)


## CUDA application:


print("Init....... ", end="")
ctx = CuContext()
mod = CuModule(ctx, "test_kernel.ptx")

mod.setConstant("Nt", c_int(Nt))
mod.setConstant("Nf", c_int(Nf))
mod.setConstant("dt", c_float(dt))

data_in_d = CuArray.fromArray(I_arr)
data_FT_d = CuArray(Nf, np.complex64)
data_out_d = CuArray(Nt, np.float32)

fft_fwd = CuFFT(data_in_d, data_FT_d, direction="fwd")
mod.applyLineshapes.setArgs(data_FT_d)
mod.applyLineshapes.setGrid((Nf // Ntpb + 1, 1, 1), (Ntpb, 1, 1))
fft_rev = CuFFT(data_FT_d, data_out_d, direction="rev")
print("Done!")


print("w = {:4.2f}... ".format(wL), end="\n")
mod.setConstant("wL", c_float(wL))
fft_fwd()
mod.applyLineshapes()
fft_rev()
I_out = fft_rev.arr_out.getArray() * wL / Nt
print("Done!")

plt.plot(I_out)
plt.show()
