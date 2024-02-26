# -*- coding: utf-8 -*-
"""
=================
Remove a baseline
=================

Remove a baseline from an experimental spectrum :

"""


# Get a spectrum:
from radis import load_spec
from radis.test.utils import getTestFile

s_exp = (
    load_spec(getTestFile(r"CO2_measured_spectrum_4-5um.spec"), binary=True)
    .crop(4300, 5000)
    .sort()
)


sb = s_exp.get_baseline(algorithm="als")

s_exp.plot()
sb.plot(nfig="same", lw=3)

# Plot the spectrum minus the baseline :
s = s_exp - sb
s.plot(nfig="same", lw=2, zorder=0)
