# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 13:51:40 2021

@author: erwan
"""


import pytest

from radis.io.hitemp import fetch_hitemp
from radis.misc.config import getDatabankList


@pytest.mark.needs_connection
def test_fetch_hitemp(*args, **kwargs):

    fetch_hitemp("OH")

    assert "HITEMP-OH" in getDatabankList()


@pytest.mark.needs_connection
def test_partial_loading(*args, **kwargs):
    """ Assert that using partial loading of the database works """

    wmin, wmax = 2500, 4500

    # First ensures that the database is wider than this (else test is irrelevant) :
    df = fetch_hitemp("OH")
    assert df.wav.min() < wmin
    assert df.wav.max() > wmax

    # Now fetch with partial loading
    df = fetch_hitemp("OH", load_wavenum_min=wmin, load_wavenum_max=wmax)
    assert df.wav.min() >= wmin
    assert df.wav.max() <= wmax


@pytest.mark.needs_connection
def test_calc_hitemp_spectrum(*args, **kwargs):
    """
    Test direct loading of HDF5 files
    """

    from astropy import units as u

    from radis import calc_spectrum

    calc_spectrum(
        wavenum_min=2500 / u.cm,
        wavenum_max=4500 / u.cm,
        molecule="OH",
        Tgas=600,
        databank="hitemp",  # test by fetching directly
        verbose=False,
    )

    calc_spectrum(
        wavenum_min=2500 / u.cm,
        wavenum_max=4500 / u.cm,
        molecule="OH",
        Tgas=600,
        databank="HITEMP-OH",  # test by loading the downloaded database
        verbose=False,
    )

    return


if __name__ == "__main__":
    test_fetch_hitemp()
    test_partial_loading()
    test_calc_hitemp_spectrum()
