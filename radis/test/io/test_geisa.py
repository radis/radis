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

import astropy.units as u
import numpy as np
import pytest

from radis.api.geisaapi import gei2df
from radis.test.utils import getTestFile

files = {
    "CO": "geisa_CO_fragment.par",
    "H2O": "geisa_H2O_fragment.par",
    "CO2": "geisa_CO2_fragment.par",
}
conditions = {
    "wmin": 2002 / u.cm,
    "wmax": 2300 / u.cm,
    "molecule": "CO",
    "isotope": "1",
    "pressure": 1.01325,  # bar
    "mole_fraction": 0.1,
    "path_length": 1,  # cm
    "verbose": True,
}


@pytest.mark.fast
def test_local_geisa_co(verbose=True, warnings=True, **kwargs):

    # 1. Load local file
    fileName = files["CO"]
    df = gei2df(getTestFile(fileName), cache="regen")
    if verbose:
        print(f"Read {fileName}")
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


@pytest.mark.fast
def test_local_geisa_h2o(verbose=True, warnings=True, **kwargs):

    # 1. Load local file
    fileName = files["H2O"]
    df = gei2df(getTestFile(fileName), cache="regen")
    if verbose:
        print(f"Read {fileName}")
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


@pytest.mark.fast
def test_local_geisa_co2(verbose=True, warnings=True, **kwargs):

    # 1. Load local file
    fileName = files["CO2"]
    df = gei2df(getTestFile(fileName), cache="regen")
    if verbose:
        print(f"Read {fileName}")
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


@pytest.mark.fast
@pytest.mark.needs_connection
def test_calc_geisa_spectrum(plot=False):
    """
    Auto-fetch and calculate a CO spectrum from the Geisa database
    """
    from radis import SpectrumFactory

    sf = SpectrumFactory(**conditions)

    # Geisa
    sf.fetch_databank(
        source="geisa",
    )
    s_geisa = sf.eq_spectrum(Tgas=700, path_length=1)

    sf.fetch_databank(
        source="hitran",
    )
    s_hitran = sf.eq_spectrum(Tgas=700, path_length=1)

    # Broadening coefficients are different but areas under the lines should be the same:
    assert np.isclose(
        s_geisa.get_integral("abscoeff"), s_hitran.get_integral("abscoeff"), rtol=0.03
    )
    if plot:
        from radis import plot_diff

        plot_diff(s_geisa, s_hitran, "absorbance")


# Run all test cases
def _run_testcases(verbose=True, *args, **kwargs):

    test_local_geisa_co(verbose=verbose, *args, **kwargs)
    test_local_geisa_h2o(verbose=verbose, *args, **kwargs)
    test_local_geisa_co2(verbose=verbose, *args, **kwargs)
    test_calc_geisa_spectrum(plot=True)


if __name__ == "__main__":
    print("Running test_geisa.py")

    _run_testcases(verbose=True)
