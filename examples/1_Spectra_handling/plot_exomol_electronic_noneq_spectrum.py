# -*- coding: utf-8 -*-
"""
Example: Calculate and plot the non-equilibrium electronic spectrum of an ExoMol molecule (e.g., OH) using RADIS.

This script demonstrates:
- Calculating non-equilibrium spectra with Telec != Tvib = Trot
- Manually adjusting electronic band intensities (band scaling)
- Accessing and plotting populations for each electronic state

You can adapt this script for any ExoMol molecule by changing the database and .states file
"""

from radis import SpectrumFactory

sf = SpectrumFactory(
    wavelength_min=300,
    wavelength_max=340,
    molecule="OH",
    isotope="1",
    pressure=1,
    path_length=1,
    mole_fraction=1e-3,
    wstep=0.001,
    verbose=True,
    self_absorption=True
)

sf.fetch_databank("exomol", "MoLLIST-OH")
sf.misc.export_rovib_fraction = True

s = sf.non_eq_spectrum(
    Tvib=400,
    Trot=400,
    Telec=10000,
    mole_fraction=1e-3,
    path_length=1,
    pressure=1,
    # band_scaling={"A": 2.0, "X": 1.0}
)
s.apply_slit(0.1)
s.plot('radiance')
 