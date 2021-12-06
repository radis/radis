# -*- coding: utf-8 -*-
"""
.. _example_large_range_by_part:

==================================
Calculate a large spectrum by part
==================================

In may be faster to calculate a full-range spectrum on several partial spectral ranges,
and combine them.

Uses the :py:func:`~radis.los.slabs.MergeSlabs` function for that.

Starting from 0.11, RADIS introduces a sparse waverange implementation
that should make it possible to directly compute a full range spectrum.
See the :ref:`HITRAN full-range example <example_hitran_full_range>`

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

import matplotlib.pyplot as plt

plt.ylim(0, 1)
