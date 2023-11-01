# -*- coding: utf-8 -*-
"""
.. _example_plot_hitran_xsc:

==========================
Plot HITRAN cross-sections
==========================

Parse HITRAN ([HITRAN-2016]_) cross-sections manually downloaded from https://hitran.org/xsc/
turn them into a Radis Spectrum, and plot the resulting transmittance
rescaled for arbitrary mole fractions or path length

"""

from radis import Spectrum
from radis.phys.units import Unit as u

# Here we use a downloaded cross-section file from the test suite :
from radis.test.utils import getTestFile

datafile = getTestFile("CH3COCH3_233.4_375.2_700.0-1780.0_13.xsc")

# Temperature and pressure are given by the cross-section file used (above: 233.4 K; 375.2 Torr)
# Mole fraction and path length can be recomputed arbitrarily :
conditions = {"mole_fraction": 0.1, "path_length": 20 * u("mm")}

s = Spectrum.from_xsc(datafile, conditions)

s.plot("transmittance_noslit", lw=2)

# Rescale for other path lengths, and plot :
s.rescale_path_length(100 * u("mm")).plot(nfig="same")

import matplotlib.pyplot as plt

plt.title(s.c["molecule"])
