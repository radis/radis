# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 08:12:00 2025

@author: Nicolas Minesi
"""

from radis import SpectrumFactory  # , plot_diff

conditions = {
    # 'wavenum_min': 500,
    # 'wavenum_max': 10000,
    # "wavenum_min": 2007,
    # "wavenum_max": 2009.5,
    # "molecule": "CO",
    # "wavenum_min": 2000,
    # "wavenum_max": 2400,
    # "molecule": "CO2",
    "wavenum_min": 4028,
    "wavenum_max": 4031,
    "molecule": "H2O",
    "path_length": 0.1,
    "mole_fraction": 0.2,
    "isotope": 1,  # or "1,2" or "all"
    "pressure": 1e-5,
    "wstep": 0.002,
    # also measure interpolation time
    "save_memory": False,
    "truncation": 5,  # cm-1; Default value is 50 cm-1
}
Tgas = 700
databank = "hitemp"
database = "2010"  # "default"
#%% ############## No optimization used ##############
# conditions["optimization"] = None

# ## Using a convolution of a Gaussian and Lorentzian, no optimization
# conditions["broadening_method"] = "convolve"  # Voigt = numpy.convolve("Gaussian", "Lorentzian")
# sf = SpectrumFactory(**conditions)
# sf.fetch_databank(databank)
# s_conv = sf.eq_spectrum(Tgas=Tgas)
# s_conv.name = f"Convolution : {s_conv.c['calculation_time']:.1f}s"

# # Using a polynomial approximation, no optimization
# conditions["broadening_method"] = "voigt"  # Voigt = polynomial approximation derived from Whithing
# sf = SpectrumFactory(**conditions)
# sf.fetch_databank(databank, database=database)
# s_poly = sf.eq_spectrum(Tgas=Tgas)
# s_poly.name = f"Poly. : {s_poly.c['calculation_time']:.1f}s"

# # # # Compare the spectra
# plot_diff(s_conv, s_poly, "absorbance", yscale="log")

#%% ############## Using LDM optimization ##############
# LDM_dic = {}
# for LDM in ["simple", "min-RMS"]:
#     conditions["optimization"] = LDM
#     for method in ["convolve", "voigt", "fft"]:
#         conditions["broadening_method"] = method
#         if method == "fft":
#             conditions["truncation"] = None
#         else:
#             conditions["truncation"] = 1 #cm-1 (Default value = 50)
#         sf = SpectrumFactory(**conditions)
#         sf.fetch_databank(databank)
#         s = sf.eq_spectrum(Tgas=Tgas)
#         s.name = f"{method} - {LDM}.: {1000*s.c['calculation_time']:.1f} ms"
#         LDM_dic[f"{LDM}_{method}"] = s

#         # Compare the spectrum to the reference
#         plot_diff(s_conv, s, "absorbance", yscale="log")

# plot_diff(s_conv, LDM_dic["linear_convolve"], "absorbance", yscale="log")
# plot_diff(s_conv, LDM_dic["linear_voigt"], "absorbance", yscale="log")
# plot_diff(s_conv, LDM_dic["linear_fft"], "absorbance", yscale="log")

import matplotlib.pyplot as plt
from tqdm import tqdm

#%% Compute reference spectrum
from radis import get_residual

# Prepare the parameter grid for progress bar
# cutoff_values = [1e-27, 1e-28, 1e-29, 1e-35]
trunc_values0 = [0.05, 0.1, 0.5, 1, 5]
# LDM_values = [None, "simple"]
LDM_values = [None, "simple", "min-RMS"]
# method_values = ["convolve", "voigt", "fft"]
method_values = ["convolve", "voigt"]

# Calculate reference
conditions["optimization"] = None
conditions["broadening_method"] = "convolve"
# conditions["cutoff"] = min(cutoff_values)


from radis import load_spec

s_load = load_spec("ref_spectrum.spec")
if s_load.conditions["truncation"] >= max(trunc_values0):
    print("*** loading reference ***")
    s_ref = s_load
else:
    print(f"*** calculating reference with truncation = {max(trunc_values0)} ***")
    assert conditions["broadening_method"] == "convolve"
    assert conditions["optimization"] == None

    conditions["verbose"] = True
    conditions["truncation"] = max(trunc_values0)

    sf = SpectrumFactory(**conditions)
    sf.fetch_databank(databank)
    s_ref = sf.eq_spectrum(Tgas=Tgas)
    s_ref.store("ref_spectrum.spec", if_exists_then="replace")

#%%
conditions["verbose"] = False
total_iterations = sum(
    1 for trunc in trunc_values0 for LDM in LDM_values for method in method_values
)
progress = tqdm(total=total_iterations, desc="Processing")

results = {}
for LDM in LDM_values:
    if LDM is None:
        trunc_values = [0.1, 0.5, 1]
    else:
        trunc_values = trunc_values0
    for method in method_values:
        res_list, time_list = [], []
        for trunc in trunc_values:
            conditions["truncation"] = trunc
            conditions["optimization"] = LDM
            conditions["broadening_method"] = method

            sf = SpectrumFactory(**conditions)
            sf.fetch_databank(databank)
            s = sf.eq_spectrum(Tgas=Tgas)

            res_list.append(get_residual(s_ref, s, var="abscoeff"))
            time_list.append(s.c["calculation_time"])
            progress.update(1)
        results[f"{LDM} - {method}"] = {
            "trunc": trunc_values.copy(),
            "residual": res_list,
            "time": time_list,
        }
progress.close()

# for LDM in LDM_values:
#     res_list, time_list = [], []
#     for trunc in trunc_values:
#         conditions["truncation"] = trunc
#         sf = SpectrumFactory(**conditions)
#         sf.fetch_databank(databank)
#         s = sf.eq_spectrum(Tgas=Tgas)
#         res_list.append(get_residual(s_ref, s, var='abscoeff'))
#         time_list.append(s.c['calculation_time'])
#         progress.update(1)

#     label = f"{LDM} - {conditions['broadening_method']}"
#     plt.figure("Residual")
#     plt.plot(trunc_values, res_list, label=label)
#     plt.figure("Time")
#     plt.plot(trunc_values, time_list, label=label)

# progress.close()

#%%
linestyles = {None: "-", "simple": "--", "min-RMS": ":"}

plt.figure("Time")
for LDM in LDM_values:
    for method in method_values:
        tag = f"{LDM} - {method}"
        plt.semilogy(
            results[tag]["trunc"],
            results[tag]["time"],
            label=tag,
            linestyle=linestyles[LDM],
        )
plt.xlabel("Truncation")
plt.ylabel("Computation time (s)")
plt.legend()

plt.figure("Residual")
for LDM in LDM_values:
    for method in method_values:
        tag = f"{LDM} - {method}"
        plt.semilogy(
            results[tag]["trunc"],
            results[tag]["residual"],
            label=tag,
            linestyle=linestyles[LDM],
        )
plt.xlabel("Truncation")
plt.ylabel("Residual (abscoeff)")
plt.legend()
