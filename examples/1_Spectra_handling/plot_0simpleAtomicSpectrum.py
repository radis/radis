# -*- coding: utf-8 -*-
"""
.. _example_atomic_spectra:

========================
Atomic spectrum #0: Calculate an atomic spectrum
========================

In this example, the NIST and Kurucz databases are used to compute an atomic spectrum.
The same partition functions are employed. The Lorentzian broadening is assumed equal for all lines.

"""
#%%
def broad_arbitrary(**kwargs):
    """An arbitrary broadening formula for the Lorentzian component"""
    HWHM = 1 * (296 / kwargs["Tgas"]) ** 0.8
    shift = None
    return HWHM, shift


import radis
from radis import SpectrumFactory, plot_diff

radis.config["ALLOW_OVERWRITE"] = True

sf = SpectrumFactory(
    # wavelength_min=777,
    # wavelength_max=777.8,
    # species="O_I",
    wavelength_min=655,
    wavelength_max=657,
    species="H_I",
    wstep=0.01,
    path_length=1,  # cm
    pressure=1,  # atm
    mole_fraction=0.1,
    verbose=0,  # to keep output quiet
    lbfunc=broad_arbitrary,
    pfsource="nist",
)

# NIST
sf.fetch_databank("nist")
s_NIST = sf.eq_spectrum(Tgas=10000, name="NIST")

# Kurucz
sf.fetch_databank("kurucz")
s_KURUCZ = sf.eq_spectrum(Tgas=10000, name="Kurucz")

plot_diff(s_NIST, s_KURUCZ, "radiance_noslit", wunit="cm-1")
