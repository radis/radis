# -*- coding: utf-8 -*-
"""
================================
Calculate a spectrum from Kurucz
================================

Auto-download and compute an atomic spectra from the Kurucz database

Naming Conventions:
- You can use either the conventional naming or its simplified form.
- For instance, 'Mg_I' can be used for 'Mg' and 'Ca_II' for 'Ca+'.

"""

from radis import calc_spectrum

s = calc_spectrum(
    2000,
    5000,  # cm-1
    species="Mg_I", # Enter species name 
    Tgas=5000,  # K
    databank="kurucz",
    warnings={"AccuracyError": "ignore", "AccuracyWarning": "ignore"}
)

s.plot("radiance_noslit",wunit="nm")
