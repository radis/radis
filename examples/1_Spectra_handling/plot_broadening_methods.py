# -*- coding: utf-8 -*-
"""
=========================================================
Performance increase of broadening_methods and LDM optimizations.
=========================================================

Auto-download, calculate, and compare CO spectra calculated with several
numerical approximation of the Voigt functions. We also show that the LDM
optimization ("simple" or "min-RMS") leads to a significantly faster
computation for almost no impact on the Voigt computation.

"""
from radis import SpectrumFactory, plot_diff

trunc_ref = 8
conditions = {
    "molecule": "CO",
    "wavenum_min": 1850,
    "wavenum_max": 2325,
    # "molecule": "H2O",        #un-comment these lines to try with H2O
    # "wavelength_min": 2400,
    # "wavelength_max": 3000,
    "path_length": 0.1,
    "mole_fraction": 0.2,
    "isotope": 1,  # or "1,2" or "all"
    "pressure": 1,
    "wstep": 0.002,
    # also measure interpolation time
    "truncation": trunc_ref,  # cm-1; Default value is 50 cm-1 but the value is decreased in this example to accelerate the computation
    "cutoff": 1e-27,  # cm-1/(#.cm-2); Default value ie 1e-27 but the value is increased in this example to accelerate the computation
}
Tgas = 400
databank, database = "hitran", "default"
# databank, database =  = "hitemp", "2010" #for H2O, on your machine (will take longer to compute the reference spectrum)

#%% ############## Reference spectrum - No optimization used ##############
conditions["optimization"] = None

## Using a convolution of a Gaussian and Lorentzian, no optimization
conditions[
    "broadening_method"
] = "convolve"  # Voigt = numpy.convolve("Gaussian", "Lorentzian")
sf = SpectrumFactory(**conditions)
sf.fetch_databank(databank)
s_conv = sf.eq_spectrum(Tgas=Tgas)

s_conv.name = f"Convolution : {s_conv.c['calculation_time']:.1f}s"

# Using a polynomial approximation, no optimization
conditions[
    "broadening_method"
] = "voigt"  # Voigt = polynomial approximation derived from Whithing
sf = SpectrumFactory(**conditions)
sf.fetch_databank(databank, database=database)
s_poly = sf.eq_spectrum(Tgas=Tgas)
s_poly.name = f"Poly. : {s_poly.c['calculation_time']:.1f}s"

# Compare the spectra
plot_diff(
    s_conv,
    s_poly,
    var="absorbance",
    yscale="log",
    # yscale="linear",
)
#%% ############## Using LDM optimization ##############
"""
The objective of this section is to show how fast are the LDM optimizations.
The calculated spectra are compared to a reference (no optimization) where the
Voigts are calculated via numpy.convolve of a Gaussian and a Lorentzian.

**Conclusions:**
Using the 'linear' LDM optimization, it takes significantly less time to compute
the spectrum, whatever the broadening method employed. The relative difference is
always smaller than diff/ref < 1e-6. This difference can be improved by
changing from 'LDM = linear' to 'LDM = min-RMS' for a small cost in time.
"""
s_ref = s_conv
LDM_dic = {}
for method in ["convolve", "voigt", "fft"]:
    for LDM in ["simple", "min-RMS"]:
        conditions["optimization"] = LDM
        conditions["broadening_method"] = method
        if method == "fft":
            conditions["truncation"] = None
        else:
            conditions["truncation"] = s_ref.conditions[
                "truncation"
            ]  # cm-1 (Default value = 50)
        sf = SpectrumFactory(**conditions)
        sf.fetch_databank(databank)
        s = sf.eq_spectrum(Tgas=Tgas)
        s.name = f"{method} - {LDM}.: {s.c['calculation_time']:.1f} s"
        LDM_dic[f"{LDM}_{method}"] = s

        # Compare the spectrum to the reference
        plot_diff(s_ref, s, "absorbance", yscale="log")

#%% Compute reference spectrum
"""
The spectrum computation can be accelerated by truncation the lineshape. The
impact of the truncation on the performance is demonstrated here for several
LDM optizations and broadening methods. For all methods and optimizations,
the computation time increases linearly with truncation width. The residual
between the truncation spectrum and the reference one decreases exponentially.

For any broadening method, a decent residual of 1e-6 can be achieved with both
LDM methods and the computation takes from 0.1 to 2 seconds. This has to be
compared to the combination 'LDM = None' and 'broadening_method = convolve'
which takes approximatelly 20 seconds.
"""
s_ref = s_conv

import matplotlib.pyplot as plt
from tqdm import tqdm

from radis import get_residual

# Prepare the parameter grid for progress bar
trunc_values0 = [0.05, 0.1, 0.5, 1, 2, 3, 4, 5, 7]
LDM_values = [None, "simple", "min-RMS"]
# LDM_values = [None, "simple"]
method_values = ["convolve", "voigt"]

conditions["verbose"] = False
total_iterations = (
    sum(1 for trunc in trunc_values0 for LDM in LDM_values for method in method_values)
    + 2
)
progress = tqdm(total=total_iterations, desc="Processing")

results = {}
for LDM in LDM_values:
    for method in method_values:
        res_list, time_list = [], []

        trunc_values = trunc_values0
        if LDM is None and method == "convolve":
            trunc_values = trunc_values0[:-1]  # to avoid the longest computation

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
    progress.update(1)

linestyles = {None: ":", "simple": "-", "min-RMS": "--"}
colors = {"convolve": "k", "voigt": "r", "fft": "b"}

plt.figure("Time")
for LDM in LDM_values:
    for method in method_values:
        tag = f"{LDM} - {method}"
        plt.plot(
            results[tag]["trunc"],
            results[tag]["time"],
            f"{colors[method]}o",
            linestyle=linestyles[LDM],
            label=tag,
        )
    if LDM is not None:
        tag = f"{LDM} - fft"
        plt.plot(
            max(trunc_values0) * 1.01,
            results[tag]["time"],
            f"{colors['fft']}*",
            label=tag,
        )
plt.xlabel("Truncation")
plt.ylabel("Computation time (s)")
plt.legend()

plt.figure("Residual")
for LDM in LDM_values:
    for method, c in zip(method_values, colors):
        tag = f"{LDM} - {method}"
        plt.semilogy(
            results[tag]["trunc"],
            results[tag]["residual"],
            f"{colors[method]}o",
            linestyle=linestyles[LDM],
            label=tag,
        )
    if LDM is not None:
        tag = f"{LDM} - fft"
        plt.plot(
            max(trunc_values0) * 1.01,
            results[tag]["residual"],
            f"{colors['fft']}*",
            label=tag,
        )
plt.xlabel("Truncation")
plt.ylabel("Residual (abscoeff)")
plt.legend()
