# -*- coding: utf-8 -*-
"""

================================================================================
Example #2: multiparameters (T, x_CO, offset)
================================================================================

With the new fitting module introduced in :py:func:`~radis.tools.new_fitting.fit_spectrum` function,
and in the example of Tgas fitting using new fitting module, we can see its 1-temperature fitting
performance for equilibrium condition.

This example features how new fitting module can fit an equilibrium spectrum, with multiple fit
parameters, including gas temperature, mole fraction, and wavelength offset.

This example presents the fit of a real experiment (CO in argon) from Minesi et al. (2022) - doi:10.1007/s00340-022-07931-7
The data are stored in a MatLab file. It is worth noticing that, the result differ
slightly from ground-truth (7515 K for index = 9), due to the fact that RADIS uses HITRAN air broadening parameters for
calculation, while this experiment was originally conducted in argon. Future updates on other molecules'
broadening coefficients will increase the accuracy of these cases with non-air diluents.

"""

import astropy.units as u
import scipy.io

from radis import Spectrum
from radis.test.utils import getTestFile
from radis.tools.new_fitting import fit_spectrum

# -------------------- Step 1. Load experimental spectrum -------------------- #

data_file = "trimmed_1857_VoigtCO_Minesi.mat"
data = scipy.io.loadmat(getTestFile(data_file), simplify_cells=True)["CO_resu_Voigt"]
index = (
    9  # Ground truth is 7515 K - Minesi et al. (2022) - doi:10.1007/s00340-022-07931-7
)
s_experimental = Spectrum.from_array(
    data["nu"], data["A_exp"][:, index], "absorbance", wunit="cm-1", unit=""
)  # adimensioned


# -------------------- Step 2. Fill ground-truths and data -------------------- #


# Experimental conditions which will be used for spectrum modeling. Basically, these are known ground-truths.
experimental_conditions = {
    "molecule": "CO",  # Molecule ID
    "isotope": "1",  # Isotope ID, can have multiple at once
    "wmin": 2010.6,  # Starting wavelength/wavenumber to be cropped out from the original experimental spectrum.
    "wmax": 2011.6,  # Ending wavelength/wavenumber for the cropping range.
    "wunit": "cm-1",  # Accompanying unit of those 2 wavelengths/wavenumbers above.
    "pressure": 1
    * u.bar,  # Total pressure of gas, in "bar" unit by default, but you can use Astropy units too.
    "path_length": 10
    * u.cm,  # Experimental path length, in "cm" unit by default, but you can use Astropy units too.
    "databank": "exomol",  # hitemp is recommended
    "wstep": "auto",
}

# List of parameters to be fitted.
fit_parameters = {
    "Tgas": 5000,  # Fit parameter, accompanied by its initial value.
    "mole_fraction": 0.05,  # Species mole fraction, from 0 to 1.
    "offset": "0 cm-1",  # Experimental offset, must be a blank space separating offset amount and unit.
    "pressure": 0.8,
}

# List of bounding ranges applied for those fit parameters above.
bounding_ranges = {
    "Tgas": [2000, 9000],
    "mole_fraction": [0, 1],
    "offset": [-0.03, 0.02],
    "pressure": [0.1, 2],
}

# Fitting pipeline setups.
fit_properties = {
    "method": "lbfgsb",  # Preferred fitting method. By default, "leastsq".
    "fit_var": "absorbance",  # Spectral quantity to be extracted for fitting process, such as "radiance", "absorbance", etc.
    "normalize": False,  # Either applying normalization on both spectra or not.
    "max_loop": 300,  # Max number of loops allowed. By default, 100.
    "tol": 1e-20,  # Fitting tolerance, only applicable for "lbfgsb" method.
}

"""

For the fitting method, you can try one among 17 different fitting methods and algorithms of LMFIT,
introduced in `LMFIT method list <https://lmfit.github.io/lmfit-py/fitting.html#choosing-different-fitting-methods>`.

You can see the benchmark result of these algorithms here:
`RADIS Newfitting Algorithm Benchmark <https://github.com/radis/radis-benchmark/blob/master/manual_benchmarks/plot_newfitting_comparison_algorithm.py>`.

"""


# -------------------- Step 3. Run the fitting and retrieve results -------------------- #


# Conduct the fitting process!
s_best, result, log = fit_spectrum(
    s_exp=s_experimental,
    fit_params=fit_parameters,
    bounds=bounding_ranges,
    model=experimental_conditions,
    pipeline=fit_properties,
)


# Now investigate the log

print("\nResidual history: \n")
print(log["residual"])

print("\nFitted values history: \n")
for fit_val in log["fit_vals"]:
    print(fit_val)

print("\nTotal fitting time: ")
print(log["time_fitting"], end=" s\n")
