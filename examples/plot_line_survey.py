# -*- coding: utf-8 -*-
"""
===========
Line Survey
===========

Plot details of every single line in a spectrum.

Uses the :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey` function.

"""

from radis import calc_spectrum

s = calc_spectrum(
    wavenum_min=2380,
    wavenum_max=2400,
    mole_fraction=400e-6,
    path_length=100,  # cm
    Tgas=1500,
    molecule="CO2",
    isotope=[1],
    databank="hitran",
    export_lines=True,
)
s.apply_slit(2, "nm")
s.line_survey(overlay="radiance", barwidth=0.01)
