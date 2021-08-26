# -*- coding: utf-8 -*-
"""
.. _example_calculate_hitemp:

================================
Calculate a spectrum from HITEMP
================================

Auto-download and compute a CO spectrum from HITEMP ([HITEMP-2010]_).

Compare with HITRAN-2016 ([HITRAN-2016]_).

"""

from radis import calc_spectrum

s = calc_spectrum(
    2135,
    2170,  # cm-1
    molecule="CO",
    isotope="1",
    pressure=3,  # bar
    Tgas=2000,  # K
    mole_fraction=0.1,
    path_length=1,  # cm
    databank="hitemp",  # latest version is fetched
    name="HITEMP-2020",
)

from radis.test.utils import setup_test_line_databases

setup_test_line_databases()

s2 = calc_spectrum(
    2135,
    2170,  # cm-1
    molecule="CO",
    isotope="1",
    pressure=3,  # bar
    Tgas=2000,  # K
    mole_fraction=0.1,
    path_length=1,  # cm
    databank="HITRAN-CO-TEST",  # we could also have fetched the latest HITRAN with "databank='hitran'"
    name="HITRAN-2016",
)

# Compare:
from radis import plot_diff

plot_diff(s, s2)
