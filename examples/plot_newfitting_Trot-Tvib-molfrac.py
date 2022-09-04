# -*- coding: utf-8 -*-
"""

================================================================================
Fit a non-LTE spectrum with multiple fit parameters
================================================================================

With the new fitting module introduced in :py:func:`~radis.tools.new_fitting.fit_spectrum` function,
and in the example of Tgas fitting using new fitting module, we can see its 1-temperature fitting
performance for equilibrium condition.

This example features how new fitting module can fit non-equilibrium spectrum, with multiple fit
parameters, such as vibrational/rotational temperatures, or mole fraction, etc.

This is a real fitting case introduced by Mr. Corentin Grimaldi, featuring CO molecule emitting on
a wide range of spectrum. This case also features a user-defined trapezoid slit function, which is
accounted by the distorion (dispersion) of the slit, as a result of the influences from experimental
parameters and spectrometer dispersion parameters during the experiment.

"""

import numpy as np

from radis import load_spec
from radis.test.utils import getTestFile
from radis.tools.new_fitting import fit_spectrum

# The following script is for Mr. Correnti Grimaldi's cases, used for his own slit settings.
# This is not native to current fitting modules.

# ----------------------- Begin of Mr. Grimaldi's script ----------------------- #


def slit_dispersion(w):
    phi = -6.33
    f = 750
    gr = 300
    m = 1
    phi *= -2 * np.pi / 360
    d = 1e-3 / gr
    disp = (
        w
        / (2 * f)
        * (-np.tan(phi) + np.sqrt((2 * d / m / (w * 1e-9) * np.cos(phi)) ** 2 - 1))
    )
    return disp  # nm/mm


slit = 1500  # µm
pitch = 20  # µm
top_slit_um = slit - pitch  # µm
base_slit_um = slit + pitch  # µm
center_slit = 5090
dispersion = slit_dispersion(center_slit)
top_slit_nm = top_slit_um * 1e-3 * dispersion
base_slit_nm = base_slit_um * 1e-3 * dispersion * 1.33

# ------------------------ End of Mr. Grimaldi's script ------------------------ #


# Load experimental spectrum.
specName = (
    "Corentin_0_100cm_DownSampled_20cm_10pctCO2_1-wc-gw450-gr300-sl1500-acc5000-.spec"
)
s_experimental = load_spec(getTestFile(specName)).offset(-0.6, "nm")

# Experimental conditions which will be used for spectrum modeling. Basically, these are known ground-truths.
experimental_conditions = {
    "molecule": "CO",  # Molecule ID
    "isotope": "1,2,3",  # Isotope ID, can have multiple at once
    "wmin": 2270,  # Starting wavelength/wavenumber to be cropped out from the original experimental spectrum.
    "wmax": 2400,  # Ending wavelength/wavenumber for the cropping range.
    "wunit": "nm",  # Accompanying unit of those 2 wavelengths/wavenumbers above.
    "pressure": 1.01325,  # Partial pressure of gas, in "bar" unit.
    "path_length": 1 / 195,  # Experimental path length, in "cm" unit.
    "slit": {  # Experimental slit. In simple form: "[value] [unit]", i.e. "-0.2 nm". In complex form: a dict with parameters of apply_slit()
        "slit_function": (top_slit_nm, base_slit_nm),
        "unit": "nm",
        "shape": "trapezoidal",
        "center_wavespace": center_slit,
        "slit_dispersion": slit_dispersion,
        "inplace": False,
    },
    "wstep": 0.001,  # Resolution of wavenumber grid, in cm-1.
    "cutoff": 0,  # (RADIS native) Discard linestrengths that are lower that this to reduce calculation time, in cm-1.
    "databank": "hitemp",  # Databank used for calculation. Must be stated.
}

# List of parameters to be fitted, accompanied by its initial value
fit_parameters = {
    "Tvib": 6000,  # Vibrational temperature, in K.
    "Trot": 4000,  # Rotational temperature, in K.
    "mole_fraction": 0.05,  # Species mole fraction, from 0 to 1.
}

# List of bounding ranges applied for those fit parameters above.
# You can skip this step and let it use default bounding ranges, but this is not recommended.
# Bounding range must be at format [<lower bound>, <upper bound>].
bounding_ranges = {
    "Tvib": [2000, 7000],
    "Trot": [2000, 7000],
    "mole_fraction": [0, 0.1],
}

# Fitting pipeline setups.
fit_properties = {
    "method": "lbfgsb",  # Preferred fitting method from the 17 confirmed methods of LMFIT stated in week 4 blog. By default, "leastsq".
    "fit_var": "radiance",  # Spectral quantity to be extracted for fitting process, such as "radiance", "absorbance", etc.
    "normalize": False,  # Either applying normalization on both spectra or not.
    "max_loop": 300,  # Max number of loops allowed. By default, 200.
    "tol": 1e-20,  # Fitting tolerance, only applicable for "lbfgsb" method.
}


# Conduct the fitting process!
s_best, result, log = fit_spectrum(
    s_exp=s_experimental,  # Experimental spectrum.
    fit_params=fit_parameters,  # Fit parameters.
    bounds=bounding_ranges,  # Bounding ranges for those fit parameters.
    model=experimental_conditions,  # Experimental ground-truths conditions.
    pipeline=fit_properties,  # # Fitting pipeline references.
)


# Now investigate the result logs for additional information about what's going on during the fitting process

print("\nResidual history: \n")
print(log["residual"])

print("\nFitted values history: \n")
for fit_val in log["fit_vals"]:
    print(fit_val)

print("\nTotal fitting time: ")
print(log["time_fitting"], end=" s\n")
