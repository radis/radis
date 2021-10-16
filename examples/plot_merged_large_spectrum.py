# -*- coding: utf-8 -*-
"""
==========================
Calculate a large spectrum
==========================

Calculate a full-range spectrum on several partial spectral ranges, and combine them.

Uses the :py:func:`~radis.los.slabs.MergeSlabs` function
"""

from radis import MergeSlabs, calc_spectrum

spectra = []
for (wmin, wmax) in [(50, 3000), (3000, 7000), (7000, 10000)]:

    spectra.append(
        calc_spectrum(
            wmin,
            wmax,
            Tgas=1000,
            pressure=10,  # bar
            molecule="OH",
            path_length=1,
            mole_fraction=0.1,
            wstep=0.1,
            databank="hitemp",
            verbose=False,
        )
    )

s = MergeSlabs(*spectra, resample="full", out="transparent")
print(s)
s.plot("transmittance_noslit", wunit="nm")
