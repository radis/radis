import matplotlib.pyplot as plt
import numpy as np
from gpu import gpu_exit, gpu_init, gpu_iterate
from matplotlib.widgets import Slider

cdsd_path = "C:/CDSD4000/npy/"

print("Loading database... ", end="")
N = int(2e7)
iso = np.array(np.load(cdsd_path + "iso.npy", mmap_mode="r")[-N:])  # iso 1, 2, 3
v0 = np.array(np.load(cdsd_path + "v0.npy", mmap_mode="r")[-N:])
da = np.array(np.load(cdsd_path + "da.npy", mmap_mode="r")[-N:])
gs = np.array(np.load(cdsd_path + "gs.npy", mmap_mode="r")[-N:])
na = np.array(np.load(cdsd_path + "na.npy", mmap_mode="r")[-N:])
S0 = np.array(np.load(cdsd_path + "S0.npy", mmap_mode="r")[-N:])
El = np.array(np.load(cdsd_path + "El.npy", mmap_mode="r")[-N:])
print("Done!")

Q_intp_list = 4 * [lambda T: 1.0]
Mm_arr = np.array([0.0, 43.98983, 44.993183, 45.994076], dtype=np.float32)

v_min = 2150.0
v_max = 2450.0
dv = 0.002
v_arr = np.arange(v_min, v_max, dv, dtype=np.float32)

dxG = 0.1
dxL = 0.2

gpu_init(
    v_arr,
    dxG,
    dxL,
    iso,
    v0,
    da,
    gs,
    na,
    S0,
    El,
    Mm_arr,
    Q_intp_list,
    verbose=True,
    emulate=False,
)

p = 1.0  # bar
T = 1500.0  # K
x = 0.1
abscoeff, transmittance, iter_h, times = gpu_iterate(p, T, x)
abscoeff /= np.max(abscoeff)
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.20)
ax.set_xlim(v_max, v_min)
(line,) = plt.plot(v_arr, abscoeff, lw=0.5)

axs = fig.add_axes([0.2, 0.05, 0.6, 0.03])
slider = Slider(
    ax=axs,
    label="Temperature",
    valmin=300.0,
    valmax=4000.0,
    valinit=T,
)


def update(val):
    T = slider.val
    abscoeff, transmittance, iter_h, times = gpu_iterate(p, T, x)
    abscoeff /= np.max(abscoeff)
    line.set_ydata(abscoeff)
    ax.set_title("{:8.2f} ms".format(times["fft_rev2"]))
    fig.canvas.draw_idle()


slider.on_changed(update)


plt.show()
gpu_exit()
