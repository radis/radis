# -*- coding: utf-8 -*-
"""Tests for :py:mod:`radis.db.hitran_isotopes`."""

import pytest


def test_hitran_isotopes_vs_molparam():
    """Ensure HITRAN_ISOTOPES is consistent with isotope_name_dict.

    Every (mol_id, iso_num) pair in ``isotope_name_dict`` must appear
    in ``HITRAN_ISOTOPES`` under the correct molecule name.
    """
    from radis.db.classes import get_molecule
    from radis.db.hitran_isotopes import HITRAN_ISOTOPES
    from radis.db.molparam import isotope_name_dict

    for mol_id, iso_num in isotope_name_dict:
        mol_name = get_molecule(mol_id)
        assert (
            mol_name in HITRAN_ISOTOPES
        ), f"Molecule '{mol_name}' (id={mol_id}) missing from HITRAN_ISOTOPES"
        assert (
            iso_num in HITRAN_ISOTOPES[mol_name]
        ), f"Isotope {iso_num} of '{mol_name}' missing from HITRAN_ISOTOPES"


def test_hitran_isotopes_known_molecules():
    """Spot-check well-known molecules."""
    from radis.db.hitran_isotopes import HITRAN_ISOTOPES

    assert HITRAN_ISOTOPES["CO"] == [1, 2, 3, 4, 5, 6]
    assert HITRAN_ISOTOPES["CO2"] == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    assert HITRAN_ISOTOPES["H2O"] == [1, 2, 3, 4, 5, 6, 7]
    assert HITRAN_ISOTOPES["O2"] == [1, 2, 3]


def test_get_hitran_isotopes():
    """Test the helper function."""
    from radis.db.hitran_isotopes import get_hitran_isotopes

    assert get_hitran_isotopes("CO") == [1, 2, 3, 4, 5, 6]

    with pytest.raises(KeyError):
        get_hitran_isotopes("NOT_A_MOLECULE")


if __name__ == "__main__":
    test_hitran_isotopes_vs_molparam()
    test_hitran_isotopes_known_molecules()
    test_get_hitran_isotopes()
    print("All hitran_isotopes tests passed!")
