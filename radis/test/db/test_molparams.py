# -*- coding: utf-8 -*-
""" """

import numbers

import pytest

from radis.api.exomolapi import get_exomol_full_isotope_name
from radis.db.classes import EXOMOL_MOLECULES, HITRAN_MOLECULES, trans
from radis.db.molparam import MOLPARAMS_EXTRA_PATH, MolParams, get_extra_molparams


def test_molparams(verbose=True, *args, **kwargs):

    molpar = MolParams()

    assert molpar.get(2, 2, "abundance") == 0.0110574
    assert molpar.get("CO2", 1, "molar_mass") == 43.98983

    # Test errors :
    import pytest

    with pytest.raises(NotImplementedError):
        molpar.get("C20H25N3O", 1, "molar_mass")

    # Test reading of extra parameters
    import pytest

    molpar = MolParams(extra_file_json=None)
    with pytest.raises(NotImplementedError):
        molpar.get("trans-P2H2", 1, "molar_mass") == 43.98983

    molpar = MolParams(extra_file_json=MOLPARAMS_EXTRA_PATH)
    assert "abundance" in molpar.extra_molparams
    assert "molar_mass" in molpar.extra_molparams
    assert molpar.get("trans-P2H2", 1, "molar_mass") == 63.9631740604


@pytest.mark.fast
def test_hitran_molecules_have_molparams():
    """Every HITRAN molecule needs a row in molparam.txt, and every row a name in isotope_name_dict"""
    molpar = MolParams()
    for hitran_id, molecule in trans.items():
        row = (int(hitran_id), 1)
        assert (
            row in molpar.df.index
        ), f"{molecule} (HITRAN id {hitran_id}) has no row in molparam.txt"
    unnamed = list(molpar.df.index[molpar.df["isotope_name"].isna()])
    assert not unnamed, f"molparam.txt rows missing from isotope_name_dict: {unnamed}"


def _is_number(value):
    # hjson reads 1.0 as the int 1 (e.g. AlF), so accept ints; reject strings and bools
    return isinstance(value, numbers.Real) and not isinstance(value, bool)


@pytest.mark.fast
@pytest.mark.parametrize("molecule", EXOMOL_MOLECULES)
def test_exomol_molecules_have_molparams(molecule):
    """fetch_exomol needs isotope 1's full name, abundance and molar mass, for every listed molecule"""
    assert get_exomol_full_isotope_name(molecule, 1)
    molpar = MolParams(extra_file_json=MOLPARAMS_EXTRA_PATH)
    abundance = molpar.get(molecule, 1, "abundance")
    molar_mass = molpar.get(molecule, 1, "molar_mass")
    assert _is_number(abundance) and 0 < abundance <= 1
    assert _is_number(molar_mass) and molar_mass > 0


@pytest.mark.fast
def test_molparams_extra_file_is_consistent():
    extra = get_extra_molparams()
    assert list(extra) == ["isotope_full_name", "abundance", "molar_mass"]
    names = extra["isotope_full_name"]
    for section in extra.values():
        assert list(section) == sorted(section), "keep molecules in sorted order"
        assert section.keys() == names.keys()
    for molecule, isotopes in names.items():
        abundances = extra["abundance"][molecule]
        molar_masses = extra["molar_mass"][molecule]
        assert abundances.keys() == isotopes.keys() == molar_masses.keys()
        for iso in isotopes:
            assert _is_number(abundances[iso]) and 0 < abundances[iso] <= 1, (
                molecule,
                iso,
            )
            assert _is_number(molar_masses[iso]) and molar_masses[iso] > 0, (
                molecule,
                iso,
            )
        if molecule not in HITRAN_MOLECULES:
            # (HITRAN molecules only store isotopologues beyond HITRAN's here, e.g. CO 7-9)
            assert list(isotopes) == [str(i) for i in range(1, len(isotopes) + 1)]
            values = [abundances[iso] for iso in isotopes]
            assert values == sorted(
                values, reverse=True
            ), "isotope 1 must be the most abundant"


if __name__ == "__main__":
    test_molparams()
