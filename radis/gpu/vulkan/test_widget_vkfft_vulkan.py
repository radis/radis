import time
from ctypes import Structure, c_float, c_int

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
from numpy.fft import rfftfreq
from numpy.random import rand, randint, seed
from pyvkfft_vulkan import prepare_fft
from scipy.fft import next_fast_len
from vulkan_compute_lib import ArrayBuffer, ComputeApplication, StructBuffer

# from driver import CuArray, CuContext, CuModule


L = lambda t, w: 2 / (w * np.pi) * 1 / (1 + 4 * (t / w) ** 2)
L_FT = lambda f, w: np.exp(-np.pi * np.abs(f) * w)


def next_fast_len_even(n):
    n = next_fast_len(n)
    while n & 1:
        n = next_fast_len(n + 1)
    return n


t_max = 100.0
Nt = 300005
Nt = next_fast_len_even(Nt)
print("Nt = {:d}".format(Nt))
t_arr = np.linspace(0, t_max, Nt)
dt = t_arr[1]
f_arr = rfftfreq(Nt, dt)
Nf = len(f_arr)
Ntpb = 1024  # threads per block
w0 = 0.2
seed(0)
Nl = 200

Nb = 12  # batch


def mock_spectrum(Nt, Nl, m=1):

    I_arr = np.zeros((m, Nt + 2), dtype=np.float32)

    for k in range(m):
        line_index = randint(0, Nt, Nl)
        line_strength = rand(Nl)
        np.add.at(I_arr, (k, line_index), line_strength)
    return I_arr


I_arr = mock_spectrum(Nt, Nl, m=Nb)

## Vulkan application:


class initData(Structure):
    _fields_ = [
        ("Nt", c_int),
        ("Nf", c_int),
        ("Nb", c_int),
        ("dt", c_float),
        ("wL", c_float),
    ]


def initialize():
    global app, data_in_d, data_FT_d, data_out_d, init_h
    print("Init... ", end="")

    app = ComputeApplication(deviceID=0)

    # Buffers:
    app.init_h = initData()
    app.init_d = StructBuffer.fromStruct(app.init_h, app=app)

    app.data_in_d = ArrayBuffer.fromArr(I_arr, binding=2, app=app)
    app.data_FT_d = ArrayBuffer.fromArr(np.zeros((Nb, Nf), np.complex64), app=app)
    app.data_out_d = ArrayBuffer.fromArr(np.zeros((Nb, Nt + 2), np.float32), app=app)

    # Init FFT's:
    # TODO: init FFT without buffers
    app.fft_fwd = prepare_fft(app.data_in_d, app.data_FT_d, compute_app=app)
    app.fft_inv = prepare_fft(app.data_out_d, app.data_FT_d, compute_app=app)

    # Shaders:
    app.fft_fwd.fft(app._commandBuffer, app.data_in_d._buffer, app.data_FT_d._buffer)
    app.init_shader("test_shader.spv", (Nf // Ntpb + 1, Nb, 1), (Ntpb, 1, 1))
    app.fft_inv.ifft(app._commandBuffer, app.data_FT_d._buffer, app.data_out_d._buffer)
    app.init_shader("test_shader2.spv", (Nt // Ntpb + 1, 1, 1), (Ntpb, 1, 1))
    app.endCommandBuffer()

    # del app.fft_fwd
    # del app.fft_inv

    print("Done!")
    return app


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)

app = initialize()

app.init_h.Nt = Nt
app.init_h.Nf = Nf
app.init_h.Nb = Nb
app.init_h.dt = dt
app.init_h.wL = w0
app.init_d.setData(app.init_h)

app.run()

plt.axhline(1, c="k", lw=1, alpha=0.5)

res = app.data_out_d.getData()
lines = plt.plot(t_arr, res.T[:Nt, :1], lw=0.5)

axw = plt.axes([0.25, 0.1, 0.65, 0.03])
sw = Slider(axw, "Width", 0.0, 2.0, valinit=w0)


def update(val):
    t0 = time.perf_counter()
    app.init_h.wL = sw.val
    app.init_d.setData(app.init_h)
    app.run()
    res = app.data_out_d.getData()
    for i in range(len(lines)):
        lines[i].set_ydata(res[i, :Nt])
    fig.canvas.draw_idle()
    t1 = time.perf_counter()
    ax.set_title("Time: {:6.1f} ms".format((t1 - t0) * 1e3))


sw.on_changed(update)
plt.show()
