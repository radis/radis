# -*- coding: utf-8 -*-
""" """

import pytest

from radis.db.classes import trans
from radis.db.molparam import MOLPARAMS_EXTRA_PATH, MolParams


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


if __name__ == "__main__":
    test_molparams()
