import pytest

from radis.db.classes import is_atom, to_conventional_name
from radis.db.conventions import get_convention
from radis.db.degeneracies import gi, gs


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


@pytest.mark.fast
def test_get_convention():
    """Test the get_convention function that determines if coefficients
    use Herzberg or Dunham notation."""
    # Herzberg coefficients
    assert get_convention(["wexe"]) == "herzberg"
    assert get_convention(["wexe1"]) == "herzberg"  # with mode number suffix
    assert get_convention(["Be", "De"]) == "herzberg"
    assert get_convention(["we", "wexe", "Be"]) == "herzberg"

    # Dunham coefficients
    assert get_convention(["Y01", "Y11"]) == "dunham"
    assert get_convention(["Y10", "Y20"]) == "dunham"

    # Mixed conventions should fail
    with pytest.raises(ValueError):
        get_convention(["wexe", "Y01", "Y11"])

    # Unknown coefficient should fail
    with pytest.raises(ValueError):
        get_convention(["unknown_coeff"])


@pytest.mark.fast
def test_degeneracies_gi():
    """Test state-independent degeneracy function."""
    # CO2 isotopes
    assert gi(2, 1) == 1  # 626
    assert gi(2, 2) == 2  # 636
    assert gi(2, 3) == 1  # 628
    assert gi(2, 4) == 6  # 627

    # CO isotopes
    assert gi(5, 1) == 1  # 26
    assert gi(5, 2) == 2  # 36

    # Unknown molecule should fail
    with pytest.raises(NotImplementedError):
        gi(999, 1)

    # Unknown isotope should fail
    with pytest.raises(NotImplementedError):
        gi(2, 999)


@pytest.mark.fast
def test_degeneracies_gs():
    """Test state-dependent degeneracy function."""
    # CO2 isotopes - symmetric isotopes have tuple (1, 0)
    assert gs(2, 1) == (1, 0)  # 626
    assert gs(2, 2) == (1, 0)  # 636
    # Asymmetric CO2 isotopes have scalar 1
    assert gs(2, 3) == 1  # 628
    assert gs(2, 4) == 1  # 627

    # CO isotopes
    assert gs(5, 1) == 1  # 26
    assert gs(5, 2) == 1  # 36

    # Unknown molecule should fail
    with pytest.raises(NotImplementedError):
        gs(999, 1)
