# -*- coding: utf-8 -*-
"""
.. _example_spectrum_for_multiple_molecules:

====================================================
Calculate and Compare Spectra for Multiple Molecules
====================================================

This example demonstrates how to:

- Compute a spectrum with multiple molecules simultaneously
- Combine precomputed spectra of individual molecules
- Compare both approaches and highlight performance trade-offs

"""

import matplotlib.pyplot as plt

from radis import calc_spectrum, config

# Configure RADIS to handle missing broadening coefficients
config["MISSING_BROAD_COEF"] = "air"

# Approach 1: Direct calculation with multiple molecules
direct_h2o_co2 = calc_spectrum(
    wavelength_min=4165,
    wavelength_max=5000,
    Tgas=300,
    path_length=0.1,
    mole_fraction={"H2O": 0.19, "CO2": 0.095},
    isotope="1,2",
    verbose=False,
)

# Approach 2: Precompute individual spectra and combine
s_h2o = calc_spectrum(
    wavelength_min=4165,
    wavelength_max=5000,
    Tgas=300,
    path_length=0.1,
    mole_fraction={"H2O": 0.1},
    isotope="1,2",
    verbose=False,
)

s_co2 = calc_spectrum(
    wavelength_min=4165,
    wavelength_max=5000,
    Tgas=300,
    path_length=0.1,
    mole_fraction={"CO2": 0.1},
    isotope="1,2",
    verbose=False,
)

# Rescale and combine spectra
pre_computed_then_combined = s_h2o.rescale_mole_fraction(
    0.19
) + s_co2.rescale_mole_fraction(0.095)

# Plot the results
plt.figure(figsize=(8, 5))
direct_h2o_co2.plot(label="Direct Calculation")
pre_computed_then_combined.plot(label="Combined Precomputed Spectra", linestyle="--")
plt.legend()
plt.title("Comparison of Multi-Molecule Spectrum Approaches")
plt.show()
