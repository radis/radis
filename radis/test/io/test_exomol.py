# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 22:44:39 2021

@author: erwan
"""

import astropy.units as u
import numpy as np
import pytest

from radis.io.exomol import get_exomol_database_list, get_exomol_full_isotope_name

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


@pytest.mark.needs_connection
def test_exomol_parsing_functions(verbose=True, *args, **kwargs):
    """Test functions used to parse ExoMol website"""

    assert get_exomol_full_isotope_name("H2O", 1) == "1H2-16O"
    assert get_exomol_full_isotope_name("H2O", 4) == "1H-2H-16O"
    assert get_exomol_full_isotope_name("CO2", 3) == "16O-12C-18O"
    assert get_exomol_full_isotope_name("H2O", 1) == "1H2-16O"
    assert get_exomol_full_isotope_name("NH3", 1) == "14N-1H3"
    assert get_exomol_full_isotope_name("H2S", 1) == "1H2-32S"
    assert get_exomol_full_isotope_name("FeH", 1) == "56Fe-1H"
    assert get_exomol_full_isotope_name("CH4", 1) == "12C-1H4"
    assert get_exomol_full_isotope_name("C2H4", 1) == "12C2-1H4"
    assert get_exomol_full_isotope_name("C2H2", 1) == "12C2-1H2"

    #  Note : this may change if new databases are added ; test would have to
    # be updated in that case.
    databases, recommended = get_exomol_database_list("CH4", "12C-1H4")
    assert sorted(databases) == sorted(["xsec-MM", "YT10to10", "YT34to10", "MM"])
    assert recommended == "MM"

    # Test that these databases are found
    KNOWN_EXOMOL_DATABASE_NAMES = {
        "H2O": {"1H2-16O": ["POKAZATEL"]},
        "NH3": {"14N-1H3": ["CoYuTe"]},
        "H2S": {"1H2-32S": ["AYT2"]},
        "FeH": {"56Fe-1H": ["MoLLIST"]},
    }
    for molecule, val in KNOWN_EXOMOL_DATABASE_NAMES.items():
        for isotope_name, lookup_databases in val.items():
            known_databases, _ = get_exomol_database_list(molecule, isotope_name)
            for dat in lookup_databases:
                assert dat in known_databases


@pytest.mark.needs_connection
@pytest.mark.needs_HITRAN_credentials
def test_calc_exomol_vs_hitemp(verbose=True, plot=False, *args, **kwargs):
    """Auto-fetch and calculate a CO spectrum from the ExoMol database
    (HITEMP lineset) and compare to HITEMP (on the HITRAN server)

    https://github.com/radis/radis/pull/320#issuecomment-884508206
    """

    from radis import SpectrumFactory

    sf = SpectrumFactory(**conditions)

    # ExoMol
    sf.fetch_databank(
        source="exomol",
        broadf=False,
        broadf_download=False,  # accelerates the test!
    )
    s_exomol = sf.eq_spectrum(Tgas=1000, path_length=1)

    sf.fetch_databank(
        source="hitemp",
    )
    s_hitemp = sf.eq_spectrum(Tgas=1000, path_length=1, name="HITEMP (Air broadened)")

    if plot:
        s_exomol.plot(
            lw=3,
        )
        s_hitemp.plot(lw=1, nfig="same")
        import matplotlib.pyplot as plt

        plt.legend()

    # Broadening coefficients are different but areas under the lines should be the same:
    assert np.isclose(
        s_exomol.get_integral("abscoeff"), s_hitemp.get_integral("abscoeff"), rtol=0.001
    )


if __name__ == "__main__":
    # test_exomol_parsing_functions()
    # test_calc_exomol_spectrum()
    test_calc_exomol_vs_hitemp()
