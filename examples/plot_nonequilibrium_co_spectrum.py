# -*- coding: utf-8 -*-
"""
============================================
Calculate non-LTE spectra of carbon-monoxide
============================================

Compute a CO spectrum with the temperature of the vibrational
mode different from the temperature of the rotational mode.

This example uses the :py:func:`~radis.lbl.calc.calc_spectrum` function,
the [HITRAN-2016]_ line database to derive the line positions
and intensities, and the default RADIS spectroscopic
constants to compute nonequilibrium energies and populations,
but it can be extended to other line databases and other sets
of spectroscopic constants.

"""


from astropy import units as u

from radis import calc_spectrum

s2 = calc_spectrum(
    1900 / u.cm,
    2300 / u.cm,
    molecule="CO",
    isotope="1,2,3",
    pressure=1.01325 * u.bar,
    Tvib=700 * u.K,
    Trot=300 * u.K,
    mole_fraction=0.1,
    path_length=1 * u.cm,
    databank="hitran",  # or use 'hitemp'
)
s2.plot("radiance_noslit")

# Apply a (large) instrumental slit function :
s2.apply_slit(10, "nm")
s2.plot("radiance", nfig="same", lw=2)  # compare with previous
