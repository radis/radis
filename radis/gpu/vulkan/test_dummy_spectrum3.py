import time
from ctypes import Structure, c_float, c_int

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
from numpy.fft import rfftfreq
from numpy.random import rand, seed
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


v_min = 0.0
v_max = 100.0
N_v = 300005
N_v_FT = next_fast_len_even(N_v)
print("N_v = {:d}".format(N_v))
v_arr = np.linspace(0, v_max, N_v_FT)
dv = v_arr[1]
f_arr = rfftfreq(N_v_FT, dv)
N_x_FT = len(f_arr)
Ntpb = 1024  # threads per block
w0 = 0.2
seed(0)
N_lines = 30
T0 = 1000.0


def init_w_axis(dx, log_wi):
    log_w_min = np.min(log_wi)
    log_w_max = np.max(log_wi) + 1e-4
    N = int(np.ceil((log_w_max - log_w_min) / dx)) + 1
    return log_w_min, log_w_max, N


dxL = 0.2

E_data = 1000.0 * rand(N_lines).astype(np.float32)
log_w_data = np.log(0.1 + 0.5 * rand(N_lines).astype(np.float32))
v0_data = ((v_max - v_min) * rand(N_lines) + v_min).astype(np.float32)


log_w_min, log_w_max, N_L = init_w_axis(dxL, log_w_data)
print(N_L)


def gL(t, t0, wL):
    gamma = wL / 2
    return gamma / (np.pi * ((t - t0) ** 2 + gamma**2))


def calc_spectrum(t_arr, v0_data, w_data, E_data, T):
    I_data = np.exp(-E_data / T)
    I_arr = np.zeros(len(t_arr), dtype=np.float32)
    for i in range(len(v0_data)):
        I_arr += gL(t_arr, v0_data[i], w_data[i]) * I_data[i]

    return I_arr


## Vulkan application:

c_float_arr_N = N_L * c_float


class initData(Structure):
    _fields_ = [
        ("N_v", c_int),
        ("N_v_FT", c_int),
        ("N_x_FT", c_int),
        # ("N_L", c_int),
        ("N_lines", c_int),
        ("dv", c_float),
        ("v_min", c_float),
        ("dxL", c_float),
        # ("log_w_min", c_float),
        # ("T", c_float),
    ]


class iterData(Structure):
    _fields_ = [
        # ("N_v", c_int),
        # ("N_v_FT", c_int),
        # ("N_x_FT", c_int),
        ("N_L", c_int),
        # ("N_lines", c_int),
        # ("dv", c_float),
        # ("v_min", c_float),
        # ("dxL", c_float),
        ("log_w_min", c_float),
        ("T", c_float),
    ]


def initialize():
    global app, data_in_d, data_FT_d, data_out_d, init_h
    print("Init... ", end="")

    app = ComputeApplication(deviceID=0)

    # Buffers:
    app.init_h = initData()
    app.init_d = StructBuffer.fromStruct(app.init_h, app=app)

    app.iter_h = iterData()
    app.iter_d = StructBuffer.fromStruct(app.iter_h, app=app)

    database_arrays = [v0_data, log_w_data, E_data]

    app.database_SSBO_d = ArrayBuffer(
        (len(database_arrays), N_lines), np.float32, binding=2, app=app
    )
    byte_offset = 0
    for arr in database_arrays:
        byte_offset += app.database_SSBO_d.setData(arr, byte_offset=byte_offset)

    app.data_LDM_d = ArrayBuffer((N_L, N_v_FT), np.float32, binding=3, app=app)
    # app.data_LDM_d.setData(LDM)

    app.data_LDM_FT_d = ArrayBuffer((N_L, N_x_FT), np.complex64, binding=4, app=app)
    app.data_spectrum_FT_d = ArrayBuffer((N_x_FT,), np.complex64, binding=5, app=app)
    app.data_spectrum_d = ArrayBuffer((N_v_FT,), np.float32, binding=6, app=app)

    # Init FFT's:
    # TODO: init FFT without buffers
    app.fft_fwd = prepare_fft(app.data_LDM_d, app.data_LDM_FT_d, compute_app=app)
    app.fft_inv = prepare_fft(
        app.data_spectrum_d, app.data_spectrum_FT_d, compute_app=app
    )

    # Shaders:
    app.schedule_shader("fillLDM.spv", (N_lines // Ntpb + 1, 1, 1), (Ntpb, 1, 1))
    app.fft_fwd.fft(
        app._commandBuffer, app.data_LDM_d._buffer, app.data_LDM_FT_d._buffer
    )
    app.schedule_shader("applyLineshapes.spv", (N_x_FT // Ntpb + 1, 1, 1), (Ntpb, 1, 1))
    app.fft_inv.ifft(
        app._commandBuffer, app.data_spectrum_FT_d._buffer, app.data_spectrum_d._buffer
    )
    app.endCommandBuffer()

    # del app.fft_fwd
    # del app.fft_inv

    print("Done!")
    return app


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)

app = initialize()

app.init_h.N_v = N_v
app.init_h.N_v_FT = N_v_FT
app.init_h.N_x_FT = N_x_FT

app.init_h.N_lines = N_lines

app.init_h.dv = dv
app.init_h.v_min = v_min
app.init_h.dxL = dxL
app.init_d.setData(app.init_h)

app.iter_h.N_L = N_L
app.iter_h.log_w_min = log_w_min
app.iter_h.T = T0  # K
app.iter_d.setData(app.iter_h)


app.run()

plt.axhline(0, c="k", lw=1, alpha=0.5)

res = app.data_spectrum_d.getData()
lines = plt.plot(v_arr[:N_v], res.T[:N_v], lw=1)

I_ref = calc_spectrum(v_arr, v0_data, np.exp(log_w_data), E_data, T0)
plt.plot(v_arr[:N_v], I_ref[:N_v], "k--", lw=1)


axw = plt.axes([0.25, 0.1, 0.65, 0.03])
sw = Slider(axw, "T(K)", 300.0, 3000.0, valinit=T0)


def update(val):

    app.iter_h.T = sw.val
    app.iter_d.setData(app.iter_h)
    app.data_LDM_d.setData(np.zeros((N_L, N_v_FT), dtype=np.float32))
    t0 = time.perf_counter()
    app.run()
    t1 = time.perf_counter()
    res = app.data_spectrum_d.getData()
    lines[0].set_ydata(res[:N_v])
    fig.canvas.draw_idle()
    ax.set_title("Time: {:6.1f} ms".format((t1 - t0) * 1e3))


sw.on_changed(update)
plt.show()
