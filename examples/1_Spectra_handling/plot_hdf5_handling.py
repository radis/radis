# -*- coding: utf-8 -*-
r"""
======================
HDF5 Spectrum Handling
======================

This code example illustrates the process of saving and loading spectrum data in HDF5 format, covering three key scenarios: saving a spectrum to an HDF5 file, loading a full spectrum from a file, and loading a partial spectrum for a specific wavenumber range or specific quantities.
"""

from radis.test.utils import getTestFile
from radis.tools.database import load_spec

# Load an example spectrum
my_file = getTestFile("synth-NH3-1-500-2000cm-P10-mf0.01-p1.spec")
s1 = load_spec(my_file)

#%%
# Load retrieved spectrum

# Save to HDF5 format
s1.to_hdf5("spectrum_data.h5")

s2 = s1.from_hdf5("spectrum_data.h5")

#%%
# Load and plot partial spectrum (specific wavenumber range) from retrieved spectrum
s_partial = s1.from_hdf5("spectrum_data.h5", wmin=1020, wmax=1040)
s_partial.plot("radiance", label="Partial (1020-1040 cm-1)")

#%%
# Load specific quantities and compare
s_specific = s1.from_hdf5("spectrum_data.h5", columns=["radiance"])
s_specific.plot("radiance", label="Specific (radiance)")

#%%
# Print metadata to verify it's preserved
print("\nOriginal Spectrum conditions:")
print(s1.conditions)

print("\nLoaded HDF5 Spectrum conditions:")
print(s2.conditions)

# Print available quantities in original and retrieved spectrum
print("\nAvailable quantities in original spectrum:")
print(s1.get_vars())

print("\nAvailable quantities in retrieved spectrum:")
print(s2.get_vars())
