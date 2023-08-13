from ctypes import c_float, c_longlong

import matplotlib.pyplot as plt
import numpy as np
from driver import CuArray, CuContext, CuFFT, CuModule
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
Ntpb = 1024  # threads per block
w0 = 0.2
seed(0)
Nl = 20


def mock_spectrum(Nt, Nl):
    line_index = randint(0, Nt, Nl)
    line_strength = rand(Nl)
    I_arr = np.zeros(Nt, dtype=np.float32)
    np.add.at(I_arr, line_index, line_strength)
    return I_arr


I_arr = mock_spectrum(Nt, Nl)
I_arr2 = mock_spectrum(Nt, Nl)


## CUDA application:


def init():
    print("Init....... ", end="")
    ctx = CuContext()
    mod = CuModule(ctx, "test_kernel.ptx")

    data_in1_d = CuArray.fromArray(I_arr)
    data_in2_d = CuArray.fromArray(I_arr2)

    data_FT_d = CuArray(Nf, np.complex64)
    data_out_d = CuArray(Nt, np.float32)

    mod.applyLineshapes.setArgs(data_FT_d)
    mod.applyLineshapes.setGrid((Nf // Ntpb + 1, 1, 1), (Ntpb, 1, 1))

    mod.setConstant("Nt", c_longlong(Nt))
    mod.setConstant("Nf", c_longlong(Nf))
    mod.setConstant("dt", c_float(dt))

    fft_fwd = CuFFT(data_in1_d, data_FT_d, direction="fwd")
    fft_rev = CuFFT(data_FT_d, data_out_d, direction="rev")
    print("Done!")
    return mod, (data_in1_d, data_in2_d), (fft_fwd, fft_rev)


def iterate(mod, I_arr_d, ffts, wL):

    print("w = {:4.2f}... ".format(wL), end="\n")
    fft_fwd, fft_rev = ffts
    fft_fwd.arr_in = I_arr_d
    mod.setConstant("wL", c_float(wL))
    fft_fwd()
    mod.applyLineshapes()
    fft_rev()
    I_out = fft_rev.arr_out.getArray() * wL / Nt

    print("Done!")

    return I_out


mod, data_d, ffts = init()
plt.plot(t_arr, iterate(mod, data_d[0], ffts, w0))
plt.plot(t_arr, iterate(mod, data_d[1], ffts, w0))

##for w in [0.02, 0.05, 0.1, 0.2, 0.5]:
##    I_out = iterate(mod, data_d, w)
##    plt.plot(t_arr,I_out)

plt.show()
