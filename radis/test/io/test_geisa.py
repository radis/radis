# -*- coding: utf-8 -*-
"""Test cases for GEISA parser. Currently able to test parsing of local files.
Notes
-----
Runs tests for radis/io so that they can be accessed by pytest (and hopefully
the CI test suite).
Examples
--------
Run all tests::
        pytest       (in command line, in project folder)

-------------------------------------------------------------------------------
"""

import pytest

from radis.api.geisaapi import gei2df
from radis.test.utils import getTestFile

files = {
    "CO": "geisa_CO_fragment.par",
    "H2O": "geisa_H2O_fragment.par",
    "CO2": "geisa_CO2_fragment.par",
}


@pytest.mark.fast
def test_local_geisa_co(verbose=True, warnings=True, **kwargs):

    # 1. Load local file
    fileName = files["CO"]
    df = gei2df(getTestFile(fileName), cache="regen")
    if verbose:
        print("Read {0}".format(fileName))
        print("-------------------------------------")
        print(df)

    # 2. Check if the parser works correctly
    # Check for accurate parsing of parameters for equilibrium calculations
    assert list(
        df.loc[3, ["wav", "int", "airbrd", "selbrd", "Pshft", "isoG", "iso", "mol"]]
    ) == [
        3.747902,
        1.1440e-27,
        0.0797,
        0.086,
        -0.000265,
        27,
        4,
        5,
    ]

    print("GEISA parsing process of molecule CO works normally.\n")

    return True


@pytest.mark.fast
def test_local_geisa_h2o(verbose=True, warnings=True, **kwargs):

    # 1. Load local file
    fileName = files["H2O"]
    df = gei2df(getTestFile(fileName), cache="regen")
    if verbose:
        print("Read {0}".format(fileName))
        print("-------------------------------------")
        print(df)

    # 2. Check if the parser works correctly
    # Check for accurate parsing of parameters for equilibrium calculations
    assert list(
        df.loc[3, ["wav", "int", "airbrd", "selbrd", "Pshft", "isoG", "iso", "mol"]]
    ) == [
        3.412985,
        2.4320e-32,
        0.0368,
        0.2214,
        0.000000,
        171,
        3,
        1,
    ]

    print("GEISA parsing process of molecule H2O works normally.\n")

    return True


@pytest.mark.fast
def test_local_geisa_co2(verbose=True, warnings=True, **kwargs):

    # 1. Load local file
    fileName = files["CO2"]
    df = gei2df(getTestFile(fileName), cache="regen")
    if verbose:
        print("Read {0}".format(fileName))
        print("-------------------------------------")
        print(df)

    # 2. Check if the parser works correctly
    # Check for accurate parsing of parameters for equilibrium calculations
    assert list(
        df.loc[3, ["wav", "int", "airbrd", "selbrd", "Pshft", "isoG", "iso", "mol"]]
    ) == [
        397.727316,
        4.0280e-29,
        0.0748,
        0.1019,
        -0.000580,
        626,
        1,
        2,
    ]

    print("GEISA parsing process of molecule CO2 works normally.\n")

    return True


# Run all test cases
def _run_testcases(verbose=True, *args, **kwargs):

    test_local_geisa_co(verbose=verbose, *args, **kwargs)
    test_local_geisa_h2o(verbose=verbose, *args, **kwargs)
    test_local_geisa_co2(verbose=verbose, *args, **kwargs)

    return True


if __name__ == "__main__":
    print("Running test_geisa.py")

    _run_testcases(verbose=True)
