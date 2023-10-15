import time
from ctypes import Structure, c_float, c_int

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
from numpy.fft import rfftfreq
from pyvkfft_vulkan import prepare_fft
from scipy.constants import N_A, c, h, k
from scipy.fft import next_fast_len
from vulkan_compute_lib import ArrayBuffer, ComputeApplication, StructBuffer

c_cm = 100 * c
c2 = h * c_cm / k
c_float_arr_16 = c_float * 16


L = lambda t, w: 2 / (w * np.pi) * 1 / (1 + 4 * (t / w) ** 2)
L_FT = lambda f, w: np.exp(-np.pi * np.abs(f) * w)


def next_fast_len_even(n):
    n = next_fast_len(n)
    while n & 1:
        n = next_fast_len(n + 1)
    return n


v_min = 2150.0
v_max = 2450.0
dv = 0.002
N_v = 150001
N_v_FT = next_fast_len_even(2 * N_v)
print("N_v    = {:d}".format(N_v))
print("N_v_FT = {:d}".format(N_v_FT))
v_arr = np.arange(N_v) * dv + v_min
x_arr = rfftfreq(N_v_FT, dv)
N_x_FT = len(x_arr)
Ntpb = 1024  # threads per block

T0 = 3000.0
p0 = 1.0  # bar
x = 0.8

log_2p = np.log(2 * p0)
hlog_T = 0.5 * np.log(T0)
log_rT = np.log(296.0 / T0)


def init_w_axis(dx, log_wi):
    log_w_min = np.min(log_wi)
    log_w_max = np.max(log_wi) + 1e-4
    N = int(np.ceil((log_w_max - log_w_min) / dx)) + 1
    return log_w_min, log_w_max, N


# {0:iso, 1:v0, 2:da, 3:S0, 4:El, 5:na, 6:gamma_arr}
database_arrays = [np.load("data_arr_{:d}.npy".format(i)) for i in range(7)]
N_lines = database_arrays[0].size

dxG = 0.13753508031368256
dxL = 0.20180289447307587

# Calc G params:
Mm_list = [4.3989830e01, 4.4993183e01, 4.5994076e01]
iso = database_arrays[0]
v0 = database_arrays[1]
log_c2Mm_arr = c_float_arr_16(
    0, *[0.5 * np.log(8 * k * np.log(2) / (c**2 * Mm * 1e-3 / N_A)) for Mm in Mm_list]
)
log_wG_data0 = np.log(v0) + np.array(log_c2Mm_arr)[iso]
log_wG_data = log_wG_data0 + hlog_T
log_wG_min, log_wG_max, N_G = init_w_axis(dxG, log_wG_data)


# Calc L params:
na = database_arrays[-2]
gamma = database_arrays[-1][0]
log_wL_data0 = np.log(gamma)
log_wL_data = log_wL_data0 + log_2p + na * log_rT
log_wL_min, log_wL_max, N_L = init_w_axis(dxL, log_wL_data)
N_coll = 1


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
        ("v_min", c_float),
        ("dv", c_float),
        ("N_v", c_int),
        ("N_v_FT", c_int),
        ("N_x_FT", c_int),
        ("dxG", c_float),
        ("dxL", c_float),
        ("N_lines", c_int),
        ("N_coll", c_int),
        ("log_c2Mm", c_float_arr_16),
    ]


class iterData(Structure):
    _fields_ = [
        ("p", c_float),
        ("log_2p", c_float),
        ("hlog_T", c_float),
        ("log_rT", c_float),
        ("c2T", c_float),
        ("N", c_float),
        ("x", c_float_arr_16),
        # ("l", c_float),
        # ("slit_FWHM", c_float),
        ("log_wG_min", c_float),
        ("log_wL_min", c_float),
        ("N_G", c_int),
        ("N_L", c_int),
        ("Q", c_float_arr_16),
    ]


def initialize():
    global app, data_in_d, data_FT_d, data_out_d, init_h
    print("Init... ", end="")

    app = ComputeApplication(deviceID=1)

    # Buffers:
    app.init_h = initData()
    app.init_d = StructBuffer.fromStruct(app.init_h, app=app)

    app.iter_h = iterData()
    app.iter_d = StructBuffer.fromStruct(app.iter_h, app=app)

    # database_arrays = [v0_data, np.exp(log_w_data), E_data]
    database_length = np.sum(
        [np.sum(arr.shape[:-1]) if len(arr.shape) > 1 else 1 for arr in database_arrays]
    )
    app.database_SSBO_d = ArrayBuffer(
        (database_length, N_lines), np.float32, binding=2, app=app
    )

    byte_offset = 0
    for arr in database_arrays:
        byte_offset += app.database_SSBO_d.setData(arr, byte_offset=byte_offset)

    app.data_LDM_d = ArrayBuffer((N_L, N_G, N_v_FT), np.float32, binding=3, app=app)
    # app.data_LDM_d.setData(LDM)

    app.data_LDM_FT_d = ArrayBuffer(
        (N_L, N_G, N_x_FT), np.complex64, binding=4, app=app
    )
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
app.init_h.dxG = dxG
app.init_h.dxL = dxL
app.init_h.N_coll = N_coll
app.init_h.log_c2Mm = log_c2Mm_arr

app.init_d.setData(app.init_h)


app.iter_h.p = p0
app.iter_h.hlog_T = hlog_T
app.iter_h.log_2p = np.log(2 * p0)
app.iter_h.log_rT = np.log(296.0 / T0)
app.iter_h.c2T = c2 / T0  # K
app.iter_h.N = p0 * 1e5 / (1e6 * k * T0)  # cm-3
app.iter_h.x = c_float_arr_16(x, 1 - x)
app.iter_h.log_wG_min = log_wG_min
app.iter_h.log_wL_min = log_wL_min
app.iter_h.N_G = N_G
app.iter_h.N_L = N_L
app.iter_h.Q = c_float_arr_16(*(16 * [1]))
app.iter_d.setData(app.iter_h)


app.run()

plt.axhline(0, c="k", lw=1, alpha=0.5)
res = app.data_spectrum_d.getData()
lines = plt.plot(v_arr[:N_v], res.T[:N_v], lw=0.5)
plt.xlim(v_max, v_min)

axw = plt.axes([0.25, 0.1, 0.65, 0.03])
sw = Slider(axw, "T(K)", 300.0, 3000.0, valinit=T0)


def update(val):
    T = sw.val

    app.iter_h.log_2p = np.log(2 * p0)
    app.iter_h.hlog_T = 0.5 * np.log(T)
    app.iter_h.log_rT = np.log(296.0 / T)
    app.iter_h.c2T = c2 / T
    app.iter_h.N = p0 * 1e5 / (1e6 * k * T)  # cm-3

    log_wG_data = log_wG_data0 + app.iter_h.hlog_T
    log_wL_data = log_wL_data0 + app.iter_h.log_2p + na * app.iter_h.log_rT

    log_wG_min, log_wG_max, N_G = init_w_axis(dxG, log_wG_data)
    log_wL_min, log_wL_max, N_L = init_w_axis(dxL, log_wL_data)

    app.iter_h.log_wG_min = log_wG_min
    app.iter_h.log_wL_min = log_wL_min
    app.iter_h.N_G = N_G
    app.iter_h.N_L = N_L

    app.iter_d.setData(app.iter_h)
    app.data_LDM_d.setData(np.zeros((N_L, N_G, N_v_FT), dtype=np.float32))
    t0 = time.perf_counter()
    app.run()
    t1 = time.perf_counter()
    res = app.data_spectrum_d.getData()
    lines[0].set_ydata(res[:N_v])
    fig.canvas.draw_idle()
    ax.set_title("Time: {:6.1f} ms".format((t1 - t0) * 1e3))


sw.on_changed(update)
plt.show()
