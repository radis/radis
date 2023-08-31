# -*- coding: utf-8 -*-
"""
================================
Calculate a spectrum from ExoMol
================================

Auto-download and compute a SiO spectrum from the ExoMol database ([ExoMol-2020]_)

"""

from radis import calc_spectrum

s = calc_spectrum(
    2000,
    5000,  # cm-1
    species="Mg_I",
    isotope="1",
    Tgas=5000,  # K
    databank="kurucz",
    warnings={"AccuracyError": "ignore", "AccuracyWarning": "ignore"}
)

s.plot("radiance_noslit",wunit="nm")
