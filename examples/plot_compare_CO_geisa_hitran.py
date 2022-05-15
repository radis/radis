# -*- coding: utf-8 -*-
"""
=========================================================
Compare CO spectrum from the GEISA and HITRAN database
=========================================================

GEISA Database has been newly implemented in RADIS 0.13 release on May 15, 2022. This is among
the very first attemps to compare the spectra generated from the two databases.

Auto-download and calculate CO spectrum from the GEISA database, and the HITRAN database.

Output should be similar, but not exactly! By default these two databases provide different
broadening coefficients. However, the Einstein coefficients & linestrengths should be
approximately the same, therefore the integrals under the lines should be similar.

You can see it by running the code below.

For your interest, GEISA and HITRAN lines can be downloaded and accessed separately using
:py:func:`~radis.io.geisa.fetch_geisa` and :py:func:`~radis.io.hitran.fetch_hitran`

"""
import astropy.units as u

from radis import calc_spectrum, plot_diff

conditions = {
    "wmin": 2002 / u.cm,
    "wmax": 2300 / u.cm,
    "molecule": "CO",
    "pressure": 1.01325,  # bar
    "Tgas": 1000,  # K
    "mole_fraction": 0.1,
    "path_length": 1,  # cm
    "verbose": True,
}

s_geisa = calc_spectrum(**conditions, databank="geisa", name="GEISA's CO")

s_hitran = calc_spectrum(
    **conditions,
    databank="hitran",
    name="HITRAN's CO",
)

"""

In :py:func:`~radis.io.geisa.fetch_geisa`, you can choose to additionally plot the
absolute difference (method='diff') by default, or the ratio (method='ratio'), or both.

"""

plot_diff(s_geisa, s_hitran, method=["diff", "ratio"])
