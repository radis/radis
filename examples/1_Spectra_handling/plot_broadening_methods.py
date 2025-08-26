# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 08:12:00 2025

@author: Nicolas Minesi
"""
import os

from radis import load_spec, plot_diff


def save_if_better(s, path):
    if not os.path.exists(path):
        s_ref.store("ref_spectrum.spec", if_exists_then="replace")
    else:

        s_load = load_spec(path)
        loaded_cdt = s_load.conditions
        new_cdt = s.conditions

        save_new = True
        for c in new_cdt.keys():
            if c == "truncation":
                if loaded_cdt[c] >= loaded_cdt[c]:
                    save_new = False
                    print(
                        f"Loaded truncation is better '{c}': {loaded_cdt[c]} >= {new_cdt[c]}"
                    )
            if new_cdt[c] != new_cdt[c]:
                save_new = False
                print("Not the same conditions '{c}': {loaded_cdt[c]} != {new_cdt[c]}")

        if save_new:
            print("Saving new s at {path}")
            return s_ref.store(path, if_exists_then="replace", compress=True)


from radis import SpectrumFactory  # , plot_diff

trunc_ref = 5
conditions = {
    # 'wavenum_min': 500,
    # 'wavenum_max': 10000,
    # "wavenum_min": 2007,
    # "wavenum_max": 2009.5,
    # "molecule": "CO",
    # "wavenum_min": 2000,
    # "wavenum_max": 2400,
    # "molecule": "CO2",
    "wavelength_min": 2400,
    "wavelength_max": 3000,
    "molecule": "H2O",
    "path_length": 0.1,
    "mole_fraction": 0.2,
    "isotope": 1,  # or "1,2" or "all"
    "pressure": 1,
    "wstep": 0.002,
    # also measure interpolation time
    "truncation": trunc_ref,  # cm-1; Default value is 50 cm-1
    "cutoff": 1e-27,  # cm-1/(#.cm-2); Default value ie 1e-27
}
Tgas = 400
databank = "hitemp"
database = "2010"  # "default"
#%% ############## No optimization used ##############
# conditions["optimization"] = None

# ## Using a convolution of a Gaussian and Lorentzian, no optimization
# conditions["broadening_method"] = "convolve"  # Voigt = numpy.convolve("Gaussian", "Lorentzian")
# sf = SpectrumFactory(**conditions)
# sf.fetch_databank(databank)
# s_conv0 = sf.eq_spectrum(Tgas=Tgas)

# s_conv0.name = f"Convolution : {s_conv0.c['calculation_time']:.1f}s"
# save_if_better(s_conv0, "H20_fundamental\\ref_spectrum.spec")

# # Using a polynomial approximation, no optimization
# conditions["broadening_method"] = "voigt"  # Voigt = polynomial approximation derived from Whithing
# sf = SpectrumFactory(**conditions)
# sf.fetch_databank(databank, database=database)
# s_poly = sf.eq_spectrum(Tgas=Tgas)
# s_poly.name = f"Poly. : {s_poly.c['calculation_time']:.1f}s"

# # # # Compare the spectra
# plot_diff(s_conv0, s_poly, "absorbance",
#           yscale="log",
#           # yscale="linear",
#           )

#%% ############## Using LDM optimization ##############
"""
The objective of this section is to show how fast are the LDM optimizations.
The calculated spectra are compared to a reference (no optimization) where the
Voigts are calculated via numpy.convolve of a Gaussian and a Lorentzian.
The reference is loaded to save time.
**Conclusions:**
It takes significantly less time to compute the spectrum, whatever the broadening
method. The relative difference is alsways smaller than diff/ref < 1e-6. This
difference can be improved by changing from 'LDM = linear' to 'LDM = min-RMS'
for a small cost in time.
"""
s_conv = load_spec("H20_fundamental\\ref_spectrum.spec")
LDM_dic = {}
for method in ["convolve", "voigt", "fft"]:
    for LDM in ["simple", "min-RMS"]:
        conditions["optimization"] = LDM
        conditions["broadening_method"] = method
        if method == "fft":
            conditions["truncation"] = None
        else:
            conditions["truncation"] = s_conv.conditions[
                "truncation"
            ]  # cm-1 (Default value = 50)
        sf = SpectrumFactory(**conditions)
        sf.fetch_databank(databank)
        s = sf.eq_spectrum(Tgas=Tgas)
        s.name = f"{method} - {LDM}.: {s.c['calculation_time']:.1f} s"
        LDM_dic[f"{LDM}_{method}"] = s

        # Compare the spectrum to the reference
        plot_diff(s_conv, s, "absorbance", yscale="log")

#%% Compute reference spectrum
s_ref = load_spec("H20_fundamental\\ref_spectrum.spec")

from radis import get_residual

# Prepare the parameter grid for progress bar
# cutoff_values = [1e-27, 1e-28, 1e-29, 1e-35]
trunc_values0 = [0.05, 0.1, 0.5, 1, 5]
LDM_values = ["simple", "min-RMS"]
# LDM_values = [None, "simple", "min-RMS"]
# method_values = ["convolve", "voigt", "fft"]
method_values = ["convolve", "voigt"]

# Calculate reference
# conditions["optimization"] = None
# conditions["broadening_method"] = "convolve"
# conditions["cutoff"] = min(cutoff_values)


# s_load = load_spec("ref_spectrum.spec")
# if s_load.conditions["truncation"] >= max(trunc_values0) and False:
#     print("*** loading reference ***")
#     s_ref = s_load
# else:
#     print(f"*** calculating reference with truncation = {max(trunc_values0)} ***")
#     assert conditions["broadening_method"] == "convolve"
#     assert conditions["optimization"] == None

#     conditions["verbose"] = True
#     conditions["truncation"] = max(trunc_values0)

#     sf = SpectrumFactory(**conditions)
#     sf.fetch_databank(databank)
#     s_ref = sf.eq_spectrum(Tgas=Tgas)
#     s_ref.store("ref_spectrum.spec", if_exists_then="replace")


import matplotlib.pyplot as plt
from tqdm import tqdm

conditions["verbose"] = False
total_iterations = sum(
    1 for trunc in trunc_values0 for LDM in LDM_values for method in method_values
)
progress = tqdm(total=total_iterations, desc="Processing")

results = {}
for LDM in LDM_values:
    if LDM is None:
        trunc_values = [0.05, 0.1, 0.5]
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

## Add fft method
for LDM in LDM_values:
    conditions["truncation"] = None
    conditions["optimization"] = LDM
    conditions["broadening_method"] = "fft"

    sf = SpectrumFactory(**conditions)
    sf.fetch_databank(databank)
    s_fft = sf.eq_spectrum(Tgas=Tgas)

    results[f"{LDM} - fft"] = {
        "residual": get_residual(s_ref, s_fft, var="abscoeff"),
        "time": s_fft.c["calculation_time"],
    }

linestyles = {None: ":", "simple": "-", "min-RMS": "--"}

plt.figure("Time")
for LDM in LDM_values:
    for method in method_values:
        tag = f"{LDM} - {method}"
        plt.plot(
            results[tag]["trunc"],
            results[tag]["time"],
            label=tag,
            linestyle=linestyles[LDM],
            # fmt='o',
        )
    tag = f"{LDM} - fft"
    plt.hlines(
        results[tag]["time"],
        label="fft",
        linestyle=linestyles[LDM],
        xmin=min(trunc_values0),
        xmax=max(trunc_values0),
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
    plt.hlines(
        get_residual(s_fft, s, var="abscoeff"),
        label="fft",
        xmin=min(trunc_values0),
        xmax=max(trunc_values0),
    )
plt.xlabel("Truncation")
plt.ylabel("Residual (abscoeff)")
plt.legend()
