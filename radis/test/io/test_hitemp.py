# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 13:51:40 2021

@author: erwan
"""


import pytest

from radis.io.hitemp import HITEMP_SOURCE_FILES, INFO_HITEMP_LINE_COUNT, fetch_hitemp
from radis.misc.config import getDatabankList


@pytest.mark.needs_connection
def test_fetch_hitemp_OH(verbose=True, *args, **kwargs):
    """Test proper download of HITEMP OH database.

    Good to test fetch_hitemp.
    ``13_HITEMP2020.par.bz2`` is only 900 kb so it can be run on online
    tests without burning the planet 🌳

    ⚠️ if using the default `chunksize=100000`, it uncompresses in one pass
    (only 57k lines for OH) and we cannot test that appending to the same HDF5
    works. So here we use a smaller chunksize of 20,000.
    .

    """

    df = fetch_hitemp("OH", cache="regen", chunksize=20000, verbose=3 * verbose)

    assert "HITEMP-OH" in getDatabankList()

    assert len(df) == 57019

    # Load again and make sure it works (ex: metadata properly loaded etc.):
    fetch_hitemp("OH")


@pytest.mark.needs_connection
@pytest.mark.download_large_databases
@pytest.mark.parametrize(
    "molecule", [mol for mol, url in HITEMP_SOURCE_FILES.items() if url]
)
def test_fetch_hitemp_all_molecules(molecule, verbose=False, *args, **kwargs):
    """Test fetch HITEMP for all molecules whose download URL is available.

    ..warning::
        this downloads gigabytes of data. It is unselected by default by Pytest
        (see radis/setup.cfg)

    The bz2 compression factor gives about 17 MB / million lines :

    - OH (57 k lines) is only 900 kb
    - CO (0.7 M lines) is only ~14 Mb
    - CH4 (114 M lines) is 435 MB

    If it fails, check the databases downloaded in ~/.radisdb



    Notes
    -----

    Performance tests of chunksize tested on CO :
        - chunksize=1000000  > 22s  , 1 iteration ~ 22s
        - chunksize=100000 > 18s,  , 1 iteration ~ 4s
        - chunksize=50000 > 19s   ,1 iteration ~ 2s
        - chunksize=1000 --> 90s  ,  1 iteration << 1s
    """

    df = fetch_hitemp(molecule, verbose=verbose)

    assert f"HITEMP-{molecule}" in getDatabankList()

    assert len(df) == INFO_HITEMP_LINE_COUNT[molecule]


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


@pytest.mark.needs_connection
def test_calc_hitemp_CO_noneq(verbose=True, *args, **kwargs):
    """Test proper download of HITEMP CO database.

    Good to test noneq calculations with direct-download of HITEMP.

    ``05_HITEMP2020.par.bz2`` is about 14 Mb. We still download it
    (once per machine) because it allows to compute & test noneq spectra,
    which failed with direct HITEMP download in RADIS 0.9.28.

    Approximate cost for the planet 🌳 :  ~ 14 gr CO2
    (hard to evaluate !).


    """

    from astropy import units as u

    from radis import calc_spectrum

    calc_spectrum(
        wavenum_min=2000 / u.cm,
        wavenum_max=2300 / u.cm,
        molecule="CO",
        isotope="1,2",
        Tvib=1500,
        Trot=300,
        databank="hitemp",  # test by fetching directly
    )

    # Recompute with (now) locally downloaded database [note : this failed on 0.9.28]
    calc_spectrum(
        wavenum_min=2000 / u.cm,
        wavenum_max=2300 / u.cm,
        molecule="CO",
        isotope="1,2",
        Tvib=1500,
        Trot=300,
        databank="HITEMP-CO",  # registered in ~/.radis
    )


if __name__ == "__main__":
    test_fetch_hitemp_OH()
    test_partial_loading()
    test_calc_hitemp_spectrum()
    test_fetch_hitemp_all_molecules("OH")
    test_fetch_hitemp_all_molecules("CO")
    test_fetch_hitemp_all_molecules("N2O", verbose=3)
    test_fetch_hitemp_all_molecules("NO", verbose=3)
    test_calc_hitemp_CO_noneq()
