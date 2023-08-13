from ctypes import c_float, c_longlong

import matplotlib.pyplot as plt
import numpy as np
from driver import CuArray, CuContext, CuFFT, CuModule
from matplotlib.widgets import Slider
from numpy.fft import rfftfreq
from numpy.random import rand, randint, seed
from scipy.fft import next_fast_len

L = lambda t, w: 2 / (w * np.pi) * 1 / (1 + 4 * (t / w) ** 2)
L_FT = lambda f, w: np.exp(-np.pi * np.abs(f) * w)

t_max = 100.0
Nt = 300005
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


## CUDA application:


def init():
    print("Init... ", end="")
    ctx = CuContext()
    cm = CuModule(ctx, "test_kernel.ptx")

    data_in_d = CuArray.fromArray(I_arr)
    data_FT_d = CuArray(Nf, np.complex64)
    data_out_d = CuArray(Nt, np.float32)

    cm.applyLineshapes.setArgs(data_FT_d)
    cm.applyLineshapes.setGrid((Nf // Ntpb + 1, 1, 1), (Ntpb, 1, 1))

    cm.setConstant("Nt", c_longlong(Nt))
    cm.setConstant("Nf", c_longlong(Nf))
    cm.setConstant("dt", c_float(dt))

    cm.fft_fwd = CuFFT(data_in_d, data_FT_d, direction="fwd")
    cm.fft_rev = CuFFT(data_FT_d, data_out_d, direction="rev")
    print("Done!")
    return cm


def iterate(cm, wL):

    cm.setConstant("wL", c_float(wL))
    cm.fft_fwd()
    cm.applyLineshapes()
    cm.fft_rev()

    return cm.fft_rev.arr_out.getArray() / Nt


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)

cm = init()
plt.axhline(0, c="k", lw=1, alpha=0.5)
(l,) = plt.plot(t_arr, iterate(cm, w0), lw=0.5)

axw = plt.axes([0.25, 0.1, 0.65, 0.03])
sw = Slider(axw, "Width", 0.0, 2.0, valinit=w0)


def update(val):
    w = sw.val
    I = iterate(cm, w)
    l.set_ydata(I)
    fig.canvas.draw_idle()


sw.on_changed(update)
plt.show()
