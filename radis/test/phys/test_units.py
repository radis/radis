# -*- coding: utf-8 -*-
"""
Test conversion functions
"""

import astropy.units as u

from radis.phys.units_astropy import convert_and_strip_units


def test_convert():

    # Dimensionned input: should be converted
    assert convert_and_strip_units(4500 * u.nm, u.um) == 4.5
    assert convert_and_strip_units(4.5 * u.um, u.nm) == 4500  # uses round(10)
    assert convert_and_strip_units(2000 * (1 / u.cm), 1 / u.m) == 200000

    # Convert using the spectral equivalencies (astropy.units.equivalencies)
    assert convert_and_strip_units(5 * u.um, 1 / u.cm) == 2000
    assert convert_and_strip_units(60 * u.THz, 1 / u.cm) == 2001.3845711889

    # Covert using temperature equivalencies
    assert convert_and_strip_units(0 * u.deg_C, u.K) == 273.15
    assert convert_and_strip_units(300 * u.K, u.imperial.deg_F) == 80.33

    # Undimensionned input : shouldnt change
    assert convert_and_strip_units(4500, u.um) == 4500

    # None type : should return None
    assert convert_and_strip_units(None) is None
    assert convert_and_strip_units(None, u.nm) is None


if __name__ == "__main__":
    test_convert()
