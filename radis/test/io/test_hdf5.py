# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 20:36:38 2021

@author: erwan
"""

import pytest

from radis.io.hdf5 import hdf2df
from radis.io.hitemp import fetch_hitemp
from radis.misc.config import getDatabankEntries


@pytest.mark.needs_connection
def test_local_hdf5_lines_loading(*args, **kwargs):
    """
    We use the OH HITEMP line database to test :py:func:`~radis.io.hitemp.fetch_hitemp`
    and :py:func:`~radis.io.hdf5.hdf2df`

    - Partial loading (only specific wavenumbers)
    - Only certain isotopes
    - Only certain columns

    """

    fetch_hitemp("OH")  # to initialize the database

    path = getDatabankEntries("HITEMP-OH")["path"]

    # Initialize the database
    fetch_hitemp("OH")
    path = getDatabankEntries("HITEMP-OH")["path"][0]
    df = hdf2df(path)
    wmin, wmax = df.wav.min(), df.wav.max()
    assert wmin < 2300  # needed for next test to be valid
    assert wmax > 2500  # needed for next test to be valid
    assert len(df.columns) > 5  # many columns loaded by default
    assert len(df.iso.unique()) > 1

    # Test loading only certain columns
    df = hdf2df(path, columns=["wav", "int"])
    assert len(df.columns) == 2 and "wav" in df.columns and "int" in df.columns

    # Test loading only certain isotopes
    df = hdf2df(path, isotope="2")
    assert df.iso.unique() == 2

    # Test partial loading of wavenumbers
    df = hdf2df(path, load_wavenum_min=2300, load_wavenum_max=2500)
    assert df.wav.min() >= 2300
    assert df.wav.max() <= 2500

    # Test with only one
    assert hdf2df(path, load_wavenum_min=2300).wav.min() >= 2300

    # Test with the other
    assert hdf2df(path, load_wavenum_max=2500).wav.max() <= 2500


if __name__ == "__main__":
    test_local_hdf5_lines_loading()
