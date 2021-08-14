# -*- coding: utf-8 -*-
"""
================================
Calculate a spectrum from ExoMol
================================

Auto-download and compute a SiO spectrum from the ExoMol database ([ExoMol-2020]_)

"""

from radis import calc_spectrum

s = calc_spectrum(
    1080,
    1320,  # cm-1
    molecule="SiO",
    isotope="1",
    pressure=1.01325,  # bar
    Tgas=1000,  # K
    mole_fraction=0.1,
    path_length=1,  # cm
    broadening_method="fft",  # @ dev: Doesn't work with 'voigt'
    databank="exomol",  # uses the recommended database. Use ('exomol', "EBJT") for a specific database ("EBJT")
)
s.apply_slit(1, "cm-1")  # simulate an experimental slit
s.plot("radiance")
