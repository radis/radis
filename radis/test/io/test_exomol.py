# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 22:44:39 2021

@author: erwan
"""

from radis.io.exomol import get_exomol_database_list, get_exomol_full_isotope_name


def test_exomol_parsing_functions(verbose=True, *args, **kwargs):
    """Test functions used to parse ExoMol website"""

    assert get_exomol_full_isotope_name("H2O", 4) == "1H-2H-16O"
    assert get_exomol_full_isotope_name("CO2", 3) == "16O-12C-18O"
    assert get_exomol_full_isotope_name("H2O", 1) == "1H2-16O"
    assert get_exomol_full_isotope_name("NH3", 1) == "14N-1H3"
    assert get_exomol_full_isotope_name("H2S", 1) == "1H2-32S"
    assert get_exomol_full_isotope_name("FeH", 1) == "56Fe-1H"

    #  Note : this may change if new databases are added ; test would have to
    # be updated in that case.
    databases, recommended = get_exomol_database_list("CH4", "12C-1H4")
    assert databases == ["xsec-YT10to10", "YT10to10", "YT34to10"]
    assert recommended == "YT34to10"

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


if __name__ == "__main__":
    test_exomol_parsing_functions()
