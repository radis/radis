# -*- coding: utf-8 -*-
"""

================================================================================
Compare performance between old 1T fitting example and new fitting module
================================================================================

This example will try to fit an experimental spectrum using two fitting pipelines:

- Old fitting module of :py:func:`~radis.tools.fitting.fit_spectrum`, which is the
  current RADIS fitting module. You can find the gallery example featuring it at
  `1 temperature fit <https://radis.readthedocs.io/en/latest/auto_examples/plot_1T_fit.html>`.
- New fitting module of :py:func:`~radis.tools.new_fitting.fit_spectrum`, a new
  fitting interface for a more practical and interactive fitting experience. See
  `1 temperature fit using new module <https://radis.readthedocs.io/en/latest/auto_examples/plot_newfitting_Tgas.html>`

With this, you can compare their overall performance, including number of loops,
fitting time, and final residuals between experimental and best-fit spectra, as
an indicator of fit accuracy.

"""

import time

import astropy.units as u

from radis import SpectrumFactory, load_spec, plot_diff
from radis.test.utils import getTestFile, setup_test_line_databases
from radis.tools.fitting import LTEModel
from radis.tools.new_fitting import fit_spectrum

# ------------------------------------ OLD 1-TEMPERATURE FITTING EXAMPLE ------------------------------------ #

setup_test_line_databases()

# fit range
wlmin = 4167
wlmax = 4180

s_exp = (
    load_spec(getTestFile("CO2_measured_spectrum_4-5um.spec"))
    .crop(wlmin, wlmax, "nm")
    .normalize()
    .sort()
    .offset(-0.2, "nm")
)


def LTEModel_withslitnorm(factory, fit_parameters, fixed_parameters):
    s = LTEModel(factory, fit_parameters, fixed_parameters)
    # we could also have added a fittable parameter, such as an offset,
    # or made the slit width a fittable parameter.
    # ... any paramter in model_input will be fitted.
    # s.offset(model_input["offset"], 'nm')
    s.apply_slit(1.4, "nm")
    return s.take("radiance").normalize()


sf = SpectrumFactory(
    wlmin * u.nm,
    wlmax * u.nm,
    wstep=0.001,  # cm-1
    pressure=1 * 1e-3,  # bar
    cutoff=1e-25,
    isotope="1,2",
    path_length=10,  # cm-1
    mole_fraction=1,
    truncation=1,  # cm-1
)
sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
sf.warnings["HighTemperatureWarning"] = "ignore"
sf.load_databank("HITRAN-CO2-TEST")

begin_time_mark = time.time()

s_best, best = sf.fit_spectrum(
    s_exp.take("radiance"),
    model=LTEModel_withslitnorm,
    fit_parameters={
        "Tgas": 300,
        # "offset": 0
    },
    bounds={
        "Tgas": [300, 2000],
        # "offset": [-1, 1],
    },
    plot=True,
    solver_options={
        "maxiter": 15,  # 👈 increase to let the fit converge
        "ftol": 1e-15,
    },
    verbose=2,
)

end_time_mark = time.time()

plot_diff(s_exp, s_best)

# Information to report later in comparison result
oldfitting_residual = best.fun
oldfitting_loops = best.nfev
oldfitting_time = end_time_mark - begin_time_mark

# ------------------------------------------------------------------------------------------------------------ #


# -------------------------------------------- NEW FITTING MODULE -------------------------------------------- #


# ------------------------------------ Step 1. Load experimental spectrum ------------------------------------ #


# Load an experimental spectrum. You can prepare yours, or fetch one of them in the radis/test/files directory.
my_spec = getTestFile("CO2_measured_spectrum_4-5um.spec")
s_experimental = load_spec(my_spec).offset(-0.2, "nm")


# ------------------------------------ Step 2. Fill ground-truths and data ------------------------------------ #


# Experimental conditions which will be used for spectrum modeling. Basically, these are known ground-truths.
experimental_conditions = {
    "molecule": "CO2",  # Molecule ID
    "isotope": "1,2",  # Isotope ID, can have multiple at once
    "wmin": 4167
    * u.nm,  # Starting wavelength/wavenumber to be cropped out from the original experimental spectrum.
    "wmax": 4180 * u.nm,  # Ending wavelength/wavenumber for the cropping range.
    "mole_fraction": 1,  # Species mole fraction, from 0 to 1.
    "pressure": 1
    * 1e-3
    * u.bar,  # Total pressure of gas, in "bar" unit by default, but you can also use Astropy units.
    "path_length": 10
    * u.cm,  # Experimental path length, in "cm" unit by default, but you can also use Astropy units.
    "slit": "1.4 nm",  # Experimental slit, must be a blank space separating slit amount and unit.
    "wstep": 0.001,  # Resolution of wavenumber grid, in cm-1.
    "databank": "hitran",  # Databank used for the spectrum calculation. Must be stated.
}

# List of parameters to be fitted, accompanied by their initial values.
fit_parameters = {
    "Tgas": 1150,  # Gas temperature, in K.
}

# List of bounding ranges applied for those fit parameters above.
# You can skip this step and let it use default bounding ranges, but this is not recommended.
# Bounding range must be at format [<lower bound>, <upper bound>].
bounding_ranges = {
    "Tgas": [
        300,
        2000,
    ],
}

# Fitting pipeline setups.
fit_properties = {
    "method": "lbfgsb",  # Preferred fitting method. By default, "leastsq".
    "fit_var": "radiance",  # Spectral quantity to be extracted for fitting process, such as "radiance", "absorbance", etc.
    "normalize": True,  # Either applying normalization on both spectra or not.
    "max_loop": 150,  # Max number of loops allowed. By default, 200.
    "tol": 1e-15,  # Fitting tolerance, only applicable for "lbfgsb" method.
}

"""

For the fitting method, you can try one among 17 different fitting methods and algorithms of LMFIT,
introduced in `LMFIT method list <https://lmfit.github.io/lmfit-py/fitting.html#choosing-different-fitting-methods>`.

You can see the benchmark result of these algorithms here:
`RADIS Newfitting Algorithm Benchmark <https://github.com/radis/radis-benchmark/blob/master/manual_benchmarks/plot_newfitting_comparison_algorithm.py>`.

"""


# ------------------------------------ Step 3. Run the fitting and retrieve results ------------------------------------ #


# Conduct the fitting process!
s_bestfit, result, log = fit_spectrum(
    s_exp=s_experimental,  # Experimental spectrum.
    fit_params=fit_parameters,  # Fit parameters.
    bounds=bounding_ranges,  # Bounding ranges for those fit parameters.
    model=experimental_conditions,  # Experimental ground-truths conditions.
    pipeline=fit_properties,  # Fitting pipeline references.
    verbose=False,  # If you want a clean result, stay False. If you want to see more about each loop, True.
)

# Information to report later in comparison result
newfitting_residual = log["residual"][-1]
newfitting_loops = result.nfev
newfitting_time = log["time_fitting"]


# ---------------------------------------------------------------------------------------------------------------------- #


# ---------------------------------- PERFORMANCE COMPARISON BETWEEN 2 FITTING METHODS ---------------------------------- #

print(
    "\n\n\n====================  PERFORMANCE COMPARISON BETWEEN 2 FITTING METHODS  ===================="
)

print("\n1. LAST RESIDUAL\n")
print(f"- Old 1T fitting example:   \t{oldfitting_residual}")
print(f"- New fitting module:       \t{newfitting_residual}")

print("\n2. NUMBER OF FITTING LOOPS\n")
print(f"- Old 1T fitting example:   \t{oldfitting_loops} loops")
print(f"- New fitting module:       \t{newfitting_loops} loops")

print("\n3. TOTAL TIME TAKEN (s)\n")
print(f"- Old 1T fitting example:   \t{oldfitting_time} s")
print(f"- New fitting module:       \t{newfitting_time} s")

print(
    "\n==========================================================================================\n"
)


"""

From the comparative result above, which includes last residual, total time fitting and total number of
fitting loops of each fitting method, we can see that under exactly the same ground-truth conditions,
fitting method (L-BFGS-B) and fitting pipeline, we can see that the new fitting module provides a better
accuracy, while requiring less fitting loops and time for execution.

You are free to change the experimental spectrum and its accompanied ground-truth conditions, but please
make sure to keep the same inputs between the two for a transparent comparative result.

"""
