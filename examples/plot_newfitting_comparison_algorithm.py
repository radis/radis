# -*- coding: utf-8 -*-
"""

================================================================================
Fitting benchmark - Performance comparison between different fitting algorithms
================================================================================

Now with the new fitting module implemented as :py:func:`~radis.tools.new_fitting.fit_spectrum`
function on September 6 2022, we can experience an easier-to-use, robust and complete spectrum
fitting.

This example features performance comparison between 17 different fitting methods and algorithms
of LMFIT, introduced in `LMFIT method list <https://lmfit.github.io/lmfit-py/fitting.html#choosing-different-fitting-methods>`.
There are 2 variables to be compared: the last residual, and total number of loops required for
the fitting process. A good fitting method should have small number of loops (robustness) and has
low residual (fitting accuracy).

"""


import matplotlib.pyplot as plt

from radis import load_spec
from radis.test.utils import getTestFile
from radis.tools.new_fitting import fit_spectrum

# ========================================== PREPARE FITTING PROCESS ========================================== #

spec_file = getTestFile("synth-NH3-1-500-2000cm-P10-mf0.01-p1.spec")
s_experimental = load_spec(spec_file)

# Experimental conditions which will be used for spectrum modeling. Basically, these are known ground-truths.
experimental_conditions = {
    "molecule": "NH3",  # Molecule ID
    "isotope": "1",  # Isotope ID, can have multiple at once
    "wmin": 1000,  # Starting wavelength/wavenumber to be cropped out from the original experimental spectrum.
    "wmax": 1050,  # Ending wavelength/wavenumber for the cropping range.
    "wunit": "cm-1",  # Accompanying unit of those 2 wavelengths/wavenumbers above.
    "mole_fraction": 0.01,  # Species mole fraction, from 0 to 1.
    "pressure": 10,  # Total pressure of gas, in "bar" unit by default, but you can also use Astropy units.
    "path_length": 1,  # Experimental path length, in "cm" unit by default, but you can also use Astropy units.
    "slit": "1 nm",  # Experimental slit, must be a blank space separating slit amount and unit.
    "offset": "-0.2 nm",  # Experimental offset, must be a blank space separating offset amount and unit.
    "wstep": 0.003,  # Resolution of wavenumber grid, in cm-1.
    "databank": "hitran",  # Databank used for the spectrum calculation. Must be stated.
}

# List of parameters to be fitted, accompanied by their initial values.
fit_parameters = {
    "Tgas": 700,  # Gas temperature, in K.
}

# List of bounding ranges applied for those fit parameters above.
bounding_ranges = {
    "Tgas": [
        500,
        2000,
    ],
}

# Fitting pipeline setups.
fit_properties = {
    "fit_var": "radiance",  # Spectral quantity to be extracted for fitting process, such as "radiance", "absorbance", etc.
    "normalize": False,  # Either applying normalization on both spectra or not.
    "max_loop": 150,  # Max number of loops allowed. By default, 200.
    "tol": 1e-15,  # Fitting tolerance, only applicable for "lbfgsb" method.
}


# ============================================ BENCHMARKING PROCESS ============================================ #

list_methods = [
    "leastsq",  # Levenberg-Marquardt
    "least_squares",  # Least-Squares minimization, using Trust Region Reflective method
    "differential_evolution",  # Differential evolution
    "brute",  # Brute force
    "basinhopping",  # Basin-hopping
    "ampgo",  # Adaptive Memory Programming for Global Optimization
    "nelder",  # Nelder-Mead method
    "lbfgsb",  # Limited-memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS-B)
    "powell",  # Powell method
    "cg",  # Conjugate-gradient
    "cobyla",  # Cobyla method
    "bfgs",  # Broyden–Fletcher–Goldfarb–Shanno (BFGS)
    "tnc",  # Truncated Newton
    "trust-constr",  # Trust-region for constrained optimization
    "slsqp",  # Sequential Linear Squares Programming
    "shgo",  # Simplicial Homology Global Optimization
    "dual_annealing",  # Dual Annealing optimization
]

# Commence fitting processes of all above fitting methods

list_residual = []  # Store last residual
list_loops = []  # Store number of fitting loops

for method in list_methods:

    fit_properties["method"] = method

    # Commence fitting process for this particular method
    _, result, log = fit_spectrum(
        s_exp=s_experimental,  # Experimental spectrum.
        fit_params=fit_parameters,  # Fit parameters.
        bounds=bounding_ranges,  # Bounding ranges for those fit parameters.
        model=experimental_conditions,  # Experimental ground-truths conditions.
        pipeline=fit_properties,  # Fitting pipeline references.
        verbose=False,  # If you want a clean result, stay False. If you want to see more about each loop, True.
        show_plot=False,  # Since benchmarking involves multiple fitting processes, turn off plot showing.
    )

    # Store the performance result
    list_residual.append(log["residual"][-1])
    list_loops.append(result.nfev)

# Print the result as a table

print("\n======================== BENCHMARKING RESULT ========================\n")
print("||           METHOD          ||          RESIDUAL         || LOOPS ||")
print("||---------------------------||---------------------------||-------||")

for i in range(len(list_methods)):
    print(
        "|| {0: <25} || {1: <25} || {2: <5} ||".format(
            list_methods[i], list_residual[i], list_loops[i]
        )
    )

print("||---------------------------||---------------------------||-------||\n")

# Plot the benchmarking result as 2D scatter plot

fig, ax = plt.subplots()

# Scatter the data points
ax.scatter(list_loops, list_residual)
# Mark method name for each corresponding data point
for i in range(len(list_methods)):
    ax.annotate(list_methods[i], (list_loops[i], list_residual[i]))
# Add axis labels
ax.set_xlabel("Number of loops")
ax.set_ylabel("Residual")

ax.grid(True)
plt.show()


"""

As seen from the data above, we have several methods like basinhopping, ampgo, dual_annealing and
differential_evolution jumping out of the loop limit of 150, making them unstable fitting methods.
Then, in this particular case, if you zoom into clustered section below, you can see several good
methods with excellently low number of loops and residual. They are leastsq, lbfgsb, tnc, bfgs,
nelder and trust-constr.

Further benchmarking attempts with different other spectra have shown that, leastsq and lbfgsb
have consistently good performance. Generally speaking, while leastsq tends to have lower residual
and thus greater accuracy, lbfgsb on the other hand often requires less fitting loops and also
less time in total. As the discrepancy between these 2 methods are not really significant, I decided
to set leastsq as default fitting method, but users are encouraged to switch to lbfgsb as one of
first attempts in case the fitting result is unsatisfied.

"""
