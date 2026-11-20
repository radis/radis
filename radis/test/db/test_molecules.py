# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 18:00:54 2018

@author: erwan
"""

import pytest

from radis.db.molecules import Molecules, getMolecule
from radis.db.utils import (
    _get_rovib_coefficients,
    get_default_jsonfile,
    get_herzberg_coefficients,
)


@pytest.mark.fast  # this is a fast test. Run fast tests only with 'pytest -m fast'
def test_herzberg_coefficients_CO2_json(verbose=True, *args, **kwargs):
    """Check that default coefficients are found correctly"""
    # @dev: update this test if we ever are to change the default Herzberg
    # coeffs in CO2/molecules_data.json

    assert set(get_herzberg_coefficients("CO2", 1, "X1SIGu+")) == set(
        ["we1", "we2", "we3", "wexe1", "wexe2", "wexe3", "Be", "De", "He"]
    )

    if verbose:
        print(
            (
                "Herzberg coefficients defined for CO2(X)(iso=1):",
                get_herzberg_coefficients("CO2", 1, "X1SIGu+"),
            )
        )


@pytest.mark.fast  # this is a fast test. Run fast tests only with 'pytest -m fast'
def test_getMolecule(verbose=True, *args, **kwargs):

    # Input
    molec = "CO"
    iso = 1
    state = "X"

    # get predefined
    S = getMolecule(molec, iso, state)

    assert S.get_fullname() == "CO(X1Σ+)(iso1)"


def test_ZPE(verbose=True, *args, **kwargs):
    """Test that Zero-Point-Energy stored in the default RADIS database
    corresponds to the Energy calculated for vi=0, J=0 for all molecules and isotopes
    and electronic states in the RADIS :py:data:`~radis.db.molecules.Molecules`
    list.
    """
    for M, isotopes in Molecules.items():
        if M != "CO2":
            continue
        for iso, elecstates in isotopes.items():
            for state, ElecState in elecstates.items():
                jsonfile = get_default_jsonfile(M)
                coeffs = _get_rovib_coefficients(
                    M,
                    iso,
                    ElecState.get_statename_utf(),
                    jsonfile=jsonfile,
                    remove_trailing_cm1=True,
                )
                if "ZPE" in coeffs:
                    # Calculate ZPE by setting all parameters to 0. Note that
                    # the number of parameters in ElecState.Erovib depends on
                    # the molecule. @dev: we use inspect to get the number of
                    # parameters
                    from radis.misc.utils import getarglist

                    params = [
                        k
                        for k in getarglist(ElecState._Erovib)
                        if k not in ["remove_ZPE", "coeff_dict"]
                    ]
                    # Calculate ZPE:
                    params = [0] * len(params)
                    assert ElecState.Erovib(*params, remove_ZPE=False) == coeffs["ZPE"]
                    if verbose:
                        print(
                            f"ZPE calculated for {ElecState.get_fullname()} matches value in database ({coeffs['ZPE']:.2f} cm-1)"
                        )
                else:
                    if verbose:
                        print(f"No ZPE defined for {ElecState.get_fullname()}")


def _run_testcases(verbose=True, *args, **kwargs):
    #
    test_herzberg_coefficients_CO2_json(verbose=verbose, *args, **kwargs)
    test_getMolecule(verbose=verbose, *args, **kwargs)
    test_ZPE(verbose=verbose, *args, **kwargs)
    #    test_CO_energies_Herzberg_vs_Dunham(verbose=verbose, *args, **kwargs)


if __name__ == "__main__":

    _run_testcases()
