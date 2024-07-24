# -*- coding: utf-8 -*-
"""
===================
Blackbody radiation
===================

Compute Planck's blackbody emission and return a RADIS
:py:class:`~radis.spectrum.spectrum.Spectrum` object
for easier post-processing.

Uses :py:class:`~radis.phys.blackbody.sPlanck`.

"""


from radis.phys.blackbody import sPlanck

sPlanck(wavelength_min=135, wavelength_max=3000, T=4000).plot()
sPlanck(wavelength_min=135, wavelength_max=3000, T=5000).plot(nfig="same")
sPlanck(wavelength_min=135, wavelength_max=3000, T=6000).plot(nfig="same")
sPlanck(wavelength_min=135, wavelength_max=3000, T=7000).plot(nfig="same")
