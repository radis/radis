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

from radis import calc_spectrum, config

# Configure RADIS to handle missing broadening coefficients
config["MISSING_BROAD_COEF"] = "air"

# Approach 1: Direct calculation with multiple molecules
direct_h2o_co2 = calc_spectrum(
    wavelength_min=4800,
    wavelength_max=5000,
    Tgas=300,
    path_length=0.1,
    mole_fraction={"H2O": 0.2, "CO2": 0.8},
    isotope="1,2",
    verbose=False,
)

# Approach 2: Precompute individual spectra and combine
s_h2o = calc_spectrum(
    wavelength_min=4800,
    wavelength_max=5000,
    Tgas=300,
    path_length=0.1,
    mole_fraction={"H2O": 0.1},
    isotope="1,2",
    verbose=False,
)

s_co2 = calc_spectrum(
    wavelength_min=4800,
    wavelength_max=5000,
    Tgas=300,
    path_length=0.1,
    mole_fraction={"CO2": 0.1},
    isotope="1,2",
    verbose=False,
)

# Select abscoeff data from precomputed spectra
s_h2o_selected = s_h2o.take("abscoeff")
s_co2_selected = s_co2.take("abscoeff")

# Rescale and combine spectra
pre_computed_then_combined = s_h2o_selected.rescale_mole_fraction(
    0.2
) + s_co2_selected.rescale_mole_fraction(0.8)

# Plot the results
direct_h2o_co2.plot("abscoeff", label="Direct Calculation")
pre_computed_then_combined.plot(label="Combined Precomputed Spectra", linestyle="--")
