from radis.api.exomolapi import exact_molname_exomol_to_simple_molname


def test_exact_molname_exomol_to_simple_molname():
    assert exact_molname_exomol_to_simple_molname("12C-1H4") == "CH4"
    assert exact_molname_exomol_to_simple_molname("23Na-16O-1H") == "NaOH"
    assert exact_molname_exomol_to_simple_molname("HeH_p") == "HeH_p"


def test_exact_molname_exomol_to_simple_molname_HDO():
    """test for HDO
    See https://github.com/HajimeKawahara/exojax/issues/528
    """

    assert exact_molname_exomol_to_simple_molname("1H-2H-16O") == "H2O"
