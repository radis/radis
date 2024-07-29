import pytest

from radis.db.classes import is_atom, to_conventional_name


@pytest.mark.fast
def test_is_atom():

    assert is_atom("Fe_I")
    assert is_atom("H_II")
    assert not is_atom("CO2")


@pytest.mark.fast
def test_to_conventional_name():
    assert to_conventional_name("N") == "N_I"
    assert to_conventional_name(801) == "O_II"
    assert to_conventional_name(2700) == "Co_I"
    assert to_conventional_name("H+") == "H_II"
    assert to_conventional_name("Y++++") == "Y_V"
    assert to_conventional_name("26.00") == "Fe_I"
    assert to_conventional_name("06.01") == "C_II"
    assert to_conventional_name("6.01") == "C_II"
    assert to_conventional_name("K_IV") == "K_IV"
    assert to_conventional_name("K IV") == "K_IV"

    # confirm it is friendly to molcules:
    assert to_conventional_name("CO2") == "CO2"
    assert to_conventional_name("OH") == "OH"
    assert to_conventional_name(49) == 49  # ID for COCl2


if __name__ == "__main__":
    test_is_atom()
    test_to_conventional_name()
