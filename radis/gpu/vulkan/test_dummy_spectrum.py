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


t_min = 0.0
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
Nl = 30


def mock_spectrum(Nt, Nl, m=1):

    I_arr = np.zeros((m, Nt + 2), dtype=np.float32)

    for k in range(m):
        line_index = randint(0, Nt, Nl)
        line_strength = rand(Nl)
        np.add.at(I_arr, (k, line_index), line_strength)
    return I_arr


def init_w_axis(dx, log_wi):
    log_w_min = np.min(log_wi)
    log_w_max = np.max(log_wi) + 1e-4
    N = int(np.ceil((log_w_max - log_w_min) / dx)) + 1
    return log_w_min, log_w_max, N


dxL = 0.2

I_data = rand(Nl).astype(np.float32)
log_w_data = np.log(0.1 + 1.9 * rand(Nl).astype(np.float32))
t0_data = ((t_max - t_min) * rand(Nl) + t_min).astype(np.float32)


log_w_min, log_w_max, Nw = init_w_axis(dxL, log_w_data)
print(Nw)


def calc_LDM(t0_data, log_w_data, I_data):
    LDM = np.zeros((Nw, 2 * Nf), dtype=np.float32)

    ki = (t0_data - t_min) / dt
    k0i = ki.astype(np.int32)
    k1i = k0i + 1
    avi = ki - k0i

    li = (log_w_data - log_w_min) / dxL
    l0i = li.astype(np.int32)
    l1i = l0i + 1
    aLi = li - l0i

    np.add.at(LDM, (l0i, k0i), (1 - aLi) * (1 - avi) * I_data)
    np.add.at(LDM, (l0i, k1i), (1 - aLi) * avi * I_data)
    np.add.at(LDM, (l1i, k0i), aLi * (1 - avi) * I_data)
    np.add.at(LDM, (l1i, k1i), aLi * avi * I_data)

    return LDM


# I_arr = mock_spectrum(Nt, Nl, m=Nb)


LDM = calc_LDM(t0_data, log_w_data, I_data)


##def calc_spectrum(t_arr, t0_data, log_w_data, I_data):
##    I_arr = np.zeros(len(t_arr), dtype=np.float32)


## Vulkan application:

c_float_arr_N = Nw * c_float


class initData(Structure):
    _fields_ = [
        ("Nt", c_int),
        ("Nf", c_int),
        ("Nb", c_int),
        ("Nw", c_int),
        ("dt", c_float),
        ("wL", c_float_arr_N),
    ]


def initialize():
    global app, data_in_d, data_FT_d, data_out_d, init_h
    print("Init... ", end="")

    app = ComputeApplication(deviceID=0)

    # Buffers:
    app.init_h = initData()
    app.init_d = StructBuffer.fromStruct(app.init_h, app=app)

    ##    app.data_in_d = ArrayBuffer((Nb, Nt+2), np.float32, binding=2, app=app)
    ##    app.data_in_d.setData(I_arr[0,:], byte_offset = 0)
    ##    app.data_in_d.setData(I_arr[1:,:], byte_offset = I_arr[0].nbytes)

    app.data_in_d = ArrayBuffer.fromArr(LDM, binding=2, app=app)
    app.data_FT_d = ArrayBuffer((Nw, Nf), np.complex64, binding=3, app=app)
    app.data_FT2_d = ArrayBuffer((Nf,), np.complex64, binding=4, app=app)
    app.data_out_d = ArrayBuffer((Nt + 2,), np.float32, binding=5, app=app)

    # Init FFT's:
    # TODO: init FFT without buffers
    app.fft_fwd = prepare_fft(app.data_in_d, app.data_FT_d, compute_app=app)
    app.fft_inv = prepare_fft(app.data_out_d, app.data_FT2_d, compute_app=app)

    # Shaders:
    app.fft_fwd.fft(app._commandBuffer, app.data_in_d._buffer, app.data_FT_d._buffer)
    app.schedule_shader("test_shader2.spv", (Nf // Ntpb + 1, 1, 1), (Ntpb, 1, 1))
    app.fft_inv.ifft(app._commandBuffer, app.data_FT2_d._buffer, app.data_out_d._buffer)
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
##app.init_h.Nb = Nb
app.init_h.Nw = Nw
app.init_h.dt = dt
app.init_h.wL = c_float_arr_N(*[np.exp(log_w_min + i * dxL) for i in range(Nw)])
app.init_d.setData(app.init_h)

app.run()

plt.axhline(0, c="k", lw=1, alpha=0.5)

res = app.data_out_d.getData()
# res = LDM
lines = plt.plot(t_arr, res.T[:Nt], lw=0.5)

axw = plt.axes([0.25, 0.1, 0.65, 0.03])
sw = Slider(axw, "Width", 0.0, 2.0, valinit=w0)


def update(val):
    t0 = time.perf_counter()
    app.init_h.wL[:] = [np.exp(log_w_min + i * dxL) for i in range(Nw)]
    app.init_d.setData(app.init_h)
    app.run()
    res = app.data_out_d.getData()
    for i in range(len(lines)):
        lines[i].set_ydata(res[:Nt])
    fig.canvas.draw_idle()
    t1 = time.perf_counter()
    ax.set_title("Time: {:6.1f} ms".format((t1 - t0) * 1e3))


sw.on_changed(update)
plt.show()
