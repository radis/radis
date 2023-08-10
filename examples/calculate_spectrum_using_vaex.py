# -*- coding: utf-8 -*-
"""
.. _calculate_vaex:

===============================
Calculate spectrum using Vaex
===============================

Calculating a spectrum using pandas can be memory-intensive, as it loads data into memory
for the computation. This can restrict the number of lines that can be processed on a
device, especially for larger datasets. In such cases, Vaex is recommended as the engine.
Vaex is a table management tool similar to PANDAS, but it utilizes memory mapping, a zero-
memory-copy policy, and lazy computations for optimal performance.
"""

# To compute spectrum using vaex add an argument engine= "vaex" in calc_spectrum, similarly add output="vaex" in fetch_hitran method to fetch data in vaex format.

from radis import calc_spectrum

s = calc_spectrum(
    1900,
    2300,  # cm-1
    molecule="CO",
    isotope="1,2,3",
    pressure=1.01325,  # bar
    Tgas=700,  # K
    mole_fraction=0.1,
    path_length=1,  # cm
    databank="hitran",  # or 'hitemp', 'geisa', 'exomol'
    engine="vaex",
)
s.apply_slit(0.5, "nm")  # simulate an experimental slit
s.plot("radiance")

# To compute spectrum using pandas change engine="pandas" (default is also pandas)
from radis import calc_spectrum

s = calc_spectrum(
    1900,
    2300,  # cm-1
    molecule="CO",
    isotope="1,2,3",
    pressure=1.01325,  # bar
    Tgas=700,  # K
    mole_fraction=0.1,
    path_length=1,  # cm
    databank="hitran",  # or 'hitemp', 'geisa', 'exomol'
    engine="pandas",
)
s.apply_slit(0.5, "nm")  # simulate an experimental slit
s.plot("radiance")
