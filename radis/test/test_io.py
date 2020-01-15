# -*- coding: utf-8 -*-
""" Test parsers

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

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.io.hitran import hit2df
from radis.io.cdsd import cdsd2df
from radis.test.utils import getTestFile
import pytest
import numpy as np
from warnings import warn


@pytest.mark.fast
def test_hitran_names_match(verbose=True, warnings=True, *args, **kwargs):
    """ Compare that HITRAN species defined in :mod:`radis.io.hitran` match
    the nomenclature dictionary : :py:data:`radis.io.hitran.trans`.
    
    This should be ensured by developers when adding new species. 
    """
    from radis.io.hitran import (
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
    )
    from radis.io.hitran import HITRAN_MOLECULES
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
def test_hitran_co(verbose=True, warnings=True, **kwargs):
    """ Analyse some default files to make sure everything still works"""

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


def test_hitran_co2(verbose=True, warnings=True, **kwargs):

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


def test_hitran_h2o(verbose=True, warnings=True, **kwargs):

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


def test_hitemp(verbose=True, warnings=True, **kwargs):
    """ Analyse some default files to make sure everything still works"""

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


def _run_testcases(verbose=True, *args, **kwargs):

    test_hitran_names_match(verbose=verbose, *args, **kwargs)
    test_hitran_co(verbose=verbose, *args, **kwargs)
    test_hitran_co2(verbose=verbose, *args, **kwargs)
    test_hitran_h2o(verbose=verbose, *args, **kwargs)
    test_hitemp(verbose=verbose, *args, **kwargs)

    return True


if __name__ == "__main__":
    print("Testing io.py: ", _run_testcases(verbose=True))
