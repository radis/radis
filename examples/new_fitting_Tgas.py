# -*- coding: utf-8 -*-
"""

===================================================================
Fit a spectrum with Tgas as fit parameter using new fitting modules
===================================================================

RADIS has its own fitting feature, as shown in 1 temperature fit example, where you have
to manually create the spectrum model, input the experimental spectrum and other ground-truths
into numerous RADIS native functions, as well as adjust the fitting pipeline yourself.

Now with the new fitting module released, all you have to do is to prepare a .spec file
containing your experimental spectrum, fill some JSON forms describing the ground-truth conditions
just like how you fill your medical checkup paper, call the function fit_spectrum() and let it do
all the work! If you are not satisfied with the result, you can simply adjust the parameters in your
JSON, such as slit and path_length , then recall the function again, until the results are satisfied.

This is way easier and more convenient than dwelling into RADIS documentation to find out and
learn how to use the current example, especially for new RADIS users. This new fitting module
aims to provide an end-to-end fitting experience, with minimum amount of RADIS knowledge needed.

Instruction:

- Step 1: prepare a .spec file. Create a .spec file containing your experimental spectrum. You can
  do it with RADIS by saving a Spectrum object with s.store(). If your current data is not a Spectrum
  object, you can convert it to a Spectrum object from Python arrays or from text files, and then save
  it as .spec file using s.store() as mentioned above.

- Step 2: fill the JSON forms. There are 4 JSON forms you need to fill: `experimental_conditions` with
  ground-truth data about your experimental environment, `fit_parameters` with the parameters you need
  to fit (such as Tgas, mole fraction, etc.), `bounding_ranges` with fitting ranges for parameters you
  listed in `fit_parameters`, and `fit_properties` for some fitting pipeline references.

- Step 3: call `fit_spectrum()` with the experimental spectrum and 4 JSON forms and see the result.

This example features fitting an experimental spectrum with Tgas, using new fitting modules.


"""

from radis.tools.new_fitting import fit_spectrum
from radis.test.utils import getTestFile
from radis import load_spec

# Load an experimental spectrum. You can prepare yours, or fetch one of them in the radis/test/files directory.
my_spec = getTestFile("synth-NH3-1-500-2000cm-P10-mf0.01-p1.spec")
s_experimental = load_spec(my_spec)

# Experimental conditions which will be used for spectrum modeling. Basically, these are known ground-truths.
experimental_conditions = {
    "molecule" : "NH3",         # Molecule ID
    "isotope" : "1",            # Isotope ID, can have multiple at once
    "wmin" : 1000,              # Starting wavelength/wavenumber to be cropped out from the original experimental spectrum.
    "wmax" : 1100,              # Ending wavelength/wavenumber for the cropping range.
    "wunit" : "cm-1",           # Accompanying unit of those 2 wavelengths/wavenumbers above.
    "mole_fraction" : 0.01,     # Species mole fraction, from 0 to 1.
    "pressure" : 10,            # Partial pressure of gas, in "bar" unit.
    "path_length" : 1,          # Experimental path length, in "cm" unit.
    "slit" : "1 nm",            # Experimental slit, must be a blank space separating slit amount and unit.
    "offset" : "-0.2 nm",       # Experimental offset, must be a blank space separating offset amount and unit.
    "databank" : "hitran"       # Databank used for the spectrum calculation. Must be stated.
}

# List of parameters to be fitted, accompanied by their initial values.
fit_parameters = {
    "Tgas" : 700,               # Gas temperature, in K.
}

# List of bounding ranges applied for those fit parameters above.
# You can skip this step and let it use default bounding ranges, but this is not recommended.
bounding_ranges = {
    "Tgas" : [500, 2000],       # Bounding range must be at format [<lower bound>, <upper bound>]
}

# Fitting pipeline setups.
fit_properties = {
    "method" : "lbfgsb",        # Preferred fitting method from the 17 confirmed methods of LMFIT stated in week 4 blog. By default, "leastsq".
    "fit_var" : "radiance",     # Spectral quantity to be extracted for fitting process, such as "radiance", "absorbance", etc.
    "normalize" : False,        # Either applying normalization on both spectra or not.
    "max_loop" : 150,           # Max number of loops allowed. By default, 200.
    "tol" : 1e-15               # Fitting tolerance, only applicable for "lbfgsb" method.
}


# Conduct the fitting process!
s_best, result, log = fit_spectrum(
    s_exp = s_experimental,
    fit_params = fit_parameters,
    bounds = bounding_ranges,
    model = experimental_conditions,
    pipeline = fit_properties
)


# Now investigate the result logs

print("\nResidual history: \n")
print(log["residual"])

print("\nFitted values history: \n")
for fit_val in log["fit_vals"]:
    print(fit_val)

print("\nTotal fitting time: ")
print(log["time_fitting"], end = " s\n")