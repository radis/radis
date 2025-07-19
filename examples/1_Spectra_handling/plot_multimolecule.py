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

import time

from radis import calc_spectrum, config, plot_diff

# Configure RADIS to handle missing broadening coefficients
config["MISSING_BROAD_COEF"] = "air"

# Direct Calculation
start_time_direct = time.perf_counter()
direct_h2o_co2 = calc_spectrum(
    wavelength_min=4800,
    wavelength_max=5000,
    Tgas=300,
    path_length=0.1,
    mole_fraction={"H2O": 0.2, "CO2": 0.8},
    isotope="1,2",
    verbose=False,
)
time_direct = time.perf_counter() - start_time_direct

# Precompute and Combine Individual Spectra
start_time_precomputed = time.perf_counter()
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
s_h2o_selected = s_h2o.take("abscoeff")
s_co2_selected = s_co2.take("abscoeff")
combined_spectrum = s_h2o_selected.rescale_mole_fraction(
    0.2
) + s_co2_selected.rescale_mole_fraction(0.8)
time_precomputed = time.perf_counter() - start_time_precomputed

print("Single Iteration Performance Results:")
print(f"Direct Calculation time: {time_direct:.4f} seconds")
print(f"Precomputed then Combined time: {time_precomputed:.4f} seconds")

# Convert combined_spectrum from wavenumber to wavelength in nm
combined_spectrum_nm = combined_spectrum.resample(
    combined_spectrum.get_wavelength(), "nm", inplace=False
)
# Plot the difference between the combined spectrum and the directly calculated spectrum.
plot_diff(combined_spectrum_nm, direct_h2o_co2)
