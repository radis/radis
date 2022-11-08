# -*- coding: utf-8 -*-
"""Test parsers.

Notes
-----

Runs tests for radis/io so that they can be accessed by pytest (and hopefully
the CI test suite)

Examples
--------

Run all tests::

        pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)::

        pytest -m fast

-------------------------------------------------------------------------------
"""

import os
from os.path import getmtime
from warnings import warn

import numpy as np
import pytest

from radis.api.cdsdapi import cdsd2df
from radis.api.hitranapi import hit2df
from radis.misc.warning import IrrelevantFileWarning
from radis.test.utils import getTestFile, setup_test_line_databases


@pytest.mark.fast
def test_hitran_names_match(verbose=True, warnings=True, *args, **kwargs):
    """Compare that HITRAN species defined in :mod:`radis.api.hitranapi` match
    the nomenclature dictionary : :py:data:`radis.api.hitranapi.trans`.

    This should be ensured by developers when adding new species.
    """
    from radis.db.classes import (
        HITRAN_CLASS1,
        HITRAN_CLASS2,
        HITRAN_CLASS3,
        HITRAN_CLASS4,
        HITRAN_CLASS5,
        HITRAN_CLASS6,
        HITRAN_CLASS7,
        HITRAN_CLASS8,
        HITRAN_CLASS9,
        HITRAN_CLASS10,
        HITRAN_MOLECULES,
    )
    from radis.misc.basics import compare_lists

    all_hitran = (
        HITRAN_CLASS1
        + HITRAN_CLASS2
        + HITRAN_CLASS3
        + HITRAN_CLASS4
        + HITRAN_CLASS5
        + HITRAN_CLASS6
        + HITRAN_CLASS7
        + HITRAN_CLASS8
        + HITRAN_CLASS9
        + HITRAN_CLASS10
    )
    all_hitran = list(set(all_hitran))

    # All species in HITRAN groups should be valid HITRAN_MOLECULES names
    for m in all_hitran:
        if not m in HITRAN_MOLECULES:
            raise ValueError(
                "{0} is defined in HITRAN groups but has no HITRAN id".format(m)
            )

    # Species in 'HITRAN_MOLECULES' should be classified in groups, else nonequilibrium
    # calculations are not possible.
    if warnings and all_hitran != HITRAN_MOLECULES:
        warn(
            "Difference between HITRAN groups (left) and HITRAN id "
            + "dictionary (right). Some HITRAN species are not classified in "
            + "groups. Nonequilibrium calculations wont be possible for these!:\n"
            + "{0}".format(
                compare_lists(
                    all_hitran, HITRAN_MOLECULES, verbose=False, return_string=True
                )[1]
            )
        )
    return


@pytest.mark.fast
def test_local_hitran_co(verbose=True, warnings=True, **kwargs):
    """Analyse some default files to make sure everything still works."""

    # 1. Load
    df = hit2df(getTestFile("hitran_CO_fragment.par"), cache="regen")
    if verbose:
        print("Read hitran_CO_fragment.par")
        print("---------------------------")
        print(df.head())

    # 2. Test
    assert list(df.loc[0, ["vu", "vl"]]) == [4, 4]
    assert df.dtypes["vu"] == np.int64
    assert df.dtypes["vl"] == np.int64

    return True


def test_local_hitran_co2(verbose=True, warnings=True, **kwargs):

    # 1. Load
    df = hit2df(getTestFile("hitran_CO2_fragment.par"), cache="regen")
    if verbose:
        print("Read hitran_CO2_fragment.par")
        print("----------------------------")
        print(df.head())

    # 2. Test
    assert list(
        df.loc[0, ["v1u", "v2u", "l2u", "v3u", "v1l", "v2l", "l2l", "v3l"]]
    ) == [4, 0, 0, 0, 0, 0, 0, 1]
    assert df.dtypes["v1l"] == np.int64
    assert df.dtypes["v3u"] == np.int64

    return True


def test_local_hitran_h2o(verbose=True, warnings=True, **kwargs):

    # 1. Load
    df = hit2df(getTestFile("hitran_2016_H2O_2iso_2000_2100cm.par"), cache="regen")
    if verbose:
        print("Read hitran_2016_H2O_2iso_2000_2100cm.par")
        print("-----------------------------------------")
        print(df.head())

    # 2. Test
    assert list(df.loc[0, ["v1u", "v2u", "v3u", "v1l", "v2l", "v3l"]]) == [
        0,
        2,
        0,
        0,
        1,
        0,
    ]

    assert df.loc[26, "ju"] == 5  # in .par : line 27, column 99-100
    assert df.loc[27, "ju"] == 18  # in .par : line 28, column 99-100

    assert df.dtypes["v1l"] == np.int64
    assert df.dtypes["v3u"] == np.int64

    assert df.dtypes["ju"] == np.int64
    assert df.dtypes["Kau"] == np.int64
    assert df.dtypes["Kcu"] == np.int64
    assert df.dtypes["jl"] == np.int64
    assert df.dtypes["Kal"] == np.int64
    assert df.dtypes["Kcl"] == np.int64

    return True


def test_local_hitemp_file(verbose=True, warnings=True, **kwargs):
    """Analyse some default files to make sure everything still works."""

    # 1. Load
    df = cdsd2df(
        getTestFile("cdsd_hitemp_09_header.txt"), cache="regen", drop_non_numeric=True
    )
    if verbose:
        print(df.head())

    # 2. Tests
    assert df.wav[3] == 2250.00096
    # make sure P Q R is correctly replaced by drop_non_numeric:
    assert "branch" in df
    assert df["branch"].iloc[0] == 0  # Q
    assert df["branch"].iloc[1] == 1  # R

    return True


def test_irrelevant_file_loading(*args, **kwargs):
    """check that irrelevant files (irrelevant wavenumber) are not loaded"""

    # For cdsd-hitemp files :

    # Ensures file is not empty
    df = cdsd2df(getTestFile("cdsd_hitemp_09_header.txt"), cache="regen")
    assert len(df) > 0

    # Now test nothing is loaded if asking for outside the wavenumber range
    with pytest.raises(IrrelevantFileWarning):
        df = cdsd2df(getTestFile("cdsd_hitemp_09_header.txt"), load_wavenum_min=100000)
    with pytest.raises(IrrelevantFileWarning):
        df = cdsd2df(getTestFile("cdsd_hitemp_09_header.txt"), load_wavenum_max=0.5)

    # For HITRAN files :

    # Ensures file is not empty
    df = hit2df(getTestFile("hitran_2016_H2O_2iso_2000_2100cm.par"), cache="regen")
    assert len(df) > 0

    # Now test nothing is loaded if asking for outside the wavenumber range
    with pytest.raises(IrrelevantFileWarning):
        df = hit2df(
            getTestFile("hitran_2016_H2O_2iso_2000_2100cm.par"), load_wavenum_min=100000
        )
    with pytest.raises(IrrelevantFileWarning):
        df = hit2df(
            getTestFile("hitran_2016_H2O_2iso_2000_2100cm.par"), load_wavenum_max=0.5
        )


def _run_example(verbose=False):
    from radis import SpectrumFactory

    setup_test_line_databases(
        verbose=verbose
    )  # add HITEMP-CO2-TEST in ~/radis.json if not there
    sf = SpectrumFactory(
        wavelength_min=4165,
        wavelength_max=4200,
        path_length=0.1,
        pressure=20,
        molecule="CO2",
        isotope="1",
        cutoff=1e-25,  # cm/molecule
        truncation=5,  # cm-1
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.load_databank("HITRAN-CO2-TEST")  # this database must be defined in ~/radis.json


def test_cache_regeneration(verbose=True, warnings=True, **kwargs):
    """Checks that if a line database is manually updated (last edited time changes),
    its cache file will be regenerated automatically
    """

    # Initialisation
    # --------------
    # Run case : we generate a first cache file
    _run_example(verbose=verbose)
    # get time when line database file was last modified
    file_last_modification = getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.par")
    )
    # get time when cache file was last modified
    cache_last_modification = getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.h5")
    )

    # To be sure, re run-example and make sure the cache file was not regenerated.
    _run_example(verbose=verbose)
    assert cache_last_modification == getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.h5")
    )

    # Test
    # ----
    # Now we fake the manual editing of the database.
    # change the time when line database was last modified :
    stinfo = os.stat(getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.par"))
    access_time = stinfo.st_atime
    os.utime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.par"),
        (file_last_modification + 1, access_time + 1),
    )

    file_last_modification_again = getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.par")
    )

    assert file_last_modification_again > file_last_modification

    # Run case : this should re-generated the cache file
    _run_example(verbose=verbose)

    cache_last_modification_again = getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.h5")
    )

    assert cache_last_modification_again > cache_last_modification


def _run_testcases(verbose=True, *args, **kwargs):

    test_hitran_names_match(verbose=verbose, *args, **kwargs)
    test_local_hitran_co(verbose=verbose, *args, **kwargs)
    test_local_hitran_co2(verbose=verbose, *args, **kwargs)
    test_local_hitran_h2o(verbose=verbose, *args, **kwargs)
    test_local_hitemp_file(verbose=verbose, *args, **kwargs)
    test_irrelevant_file_loading()
    test_cache_regeneration(verbose=verbose, *args, **kwargs)
    return True


if __name__ == "__main__":
    print("Testing test_hitran_cdsd.py: ", _run_testcases(verbose=True))
    test_cache_regeneration(verbose=3)
