# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 08:12:00 2025

@author: Nicolas Minesi
"""
import os

from radis import load_spec, plot_diff


def save_if_better(s, path):
    if not os.path.exists(path):
        s.store("ref_spectrum.spec", if_exists_then="replace")
    else:

        s_load = load_spec(path)
        loaded_cdt = s_load.conditions
        new_cdt = s.conditions

        save_new = True
        for c in new_cdt.keys():
            if c == "truncation":
                if loaded_cdt[c] >= new_cdt[c]:
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

trunc_ref = 8
conditions = {
    # 'wavenum_min': 500,
    # 'wavenum_max': 10000,
    # "wavenum_min": 2007,
    # "wavenum_max": 2009.5,
    "molecule": "CO",
    "wavenum_min": 1850,
    "wavenum_max": 2325,
    # "molecule": "CO2",
    # "wavelength_min": 2400,
    # "wavelength_max": 3000,
    # "molecule": "H2O",
    "path_length": 0.1,
    "mole_fraction": 0.2,
    "isotope": 1,  # or "1,2" or "all"
    "pressure": 1,
    "wstep": 0.002,
    # also measure interpolation time
    "truncation": trunc_ref,  # cm-1; Default value is 50 cm-1
    "cutoff": 1e-27,  # cm-1/(#.cm-2); Default value ie 1e-27 but the value is increased in this example to accelerate the computation
}
Tgas = 400
databank = "hitran"
database = "default"  # "2019", "2010"
#%% ############## No optimization used ##############
conditions["optimization"] = None

## Using a convolution of a Gaussian and Lorentzian, no optimization
conditions[
    "broadening_method"
] = "convolve"  # Voigt = numpy.convolve("Gaussian", "Lorentzian")
sf = SpectrumFactory(**conditions)
sf.fetch_databank(databank)
s_conv0 = sf.eq_spectrum(Tgas=Tgas)

s_conv0.name = f"Convolution : {s_conv0.c['calculation_time']:.1f}s"
# save_if_better(s_conv0, "H20_fundamental\\ref_spectrum.spec") # TODO decide if we keep or not

# Using a polynomial approximation, no optimization
conditions[
    "broadening_method"
] = "voigt"  # Voigt = polynomial approximation derived from Whithing
sf = SpectrumFactory(**conditions)
sf.fetch_databank(databank, database=database)
s_poly = sf.eq_spectrum(Tgas=Tgas)
s_poly.name = f"Poly. : {s_poly.c['calculation_time']:.1f}s"

# # # Compare the spectra
plot_diff(
    s_conv0,
    s_poly,
    "absorbance",
    yscale="log",
    # yscale="linear",
)
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
# s_conv = load_spec("H20_fundamental\\ref_spectrum.spec")
s_conv = s_conv0
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
"""
The spectrum computation can be accelerated by truncation the lineshape. The
impact of the truncation on the performance is demonstrated here for several
LDM optizations and broadening methods. For all methods and optimizations,
the computation time increases linearly with truncation width. The residual
between the truncation spectrum and the reference one decreases exponentially.
A decent residual of 1e-6 can be achieved with both LDM methods in X seconds.

Observation 1: Residual does not improve for truncation higher than 5 cm-1 (for
a reference truncation at 10 cm-1).
Observation 2: The fft method is the longuest. However, any broadening method
with LDM is several orders of magnitude faster than without LDM optimization.

"""
# s_ref = load_spec("H20_fundamental\\ref_spectrum.spec")
s_ref = s_conv0

from radis import get_residual

# Prepare the parameter grid for progress bar
# cutoff_values = [1e-27, 1e-28, 1e-29, 1e-35]
trunc_values0 = [0.05, 0.1, 0.5, 1, 2, 3, 4, 5, 7]
# LDM_values = ["simple", "min-RMS"]
LDM_values = [None, "simple", "min-RMS"]
# LDM_values = [None, "simple"]
# method_values = ["convolve", "voigt", "fft"]
method_values = ["convolve", "voigt"]

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
        trunc_values = trunc_values0[:-1]  # to avoid the longest
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

## Add fft method (no truncation)
for LDM in ["simple", "min-RMS"]:
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
    if LDM is not None:
        tag = f"{LDM} - fft"
        plt.plot(
            max(trunc_values0) * 1.01,
            results[tag]["time"],
            "*",
            label=tag,
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
    if LDM is not None:
        tag = f"{LDM} - fft"
        plt.plot(
            max(trunc_values0) * 1.01,
            results[tag]["residual"],
            "*",
            label=tag,
        )
plt.xlabel("Truncation")
plt.ylabel("Residual (abscoeff)")
plt.legend()
