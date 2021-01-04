# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:24:27 2018

@author: erwan
"""

from radis.levels.dunham import EvJ


def test_EvJ(verbose=True, *args, **kwargs):

    # Test hardcoded values with arbitrary Dunham coefficients
    assert EvJ(1, 2, Y10=0.2e4, Y20=-0.1e2, Y30=0.1e-1) == 2977.53375

    # Test **kwargs arguments
    assert EvJ(1, 2, Y10=0.2e4, Y20=-0.1e2, Y30=0.1e-1) == EvJ(
        1, 2, Y20=-0.1e2, Y30=0.1e-1, **{"Y10_cm-1": 0.2e4}
    )


def test_dunham_co(verbose=True, *args, **kwargs):

    molecule = "CO"
    isotope = 1

    # %% Check molecule data JSON
    from radis.db.utils import get_dunham_coefficients

    dunham_coeffs = get_dunham_coefficients(molecule, isotope, "X1SIG+")

    #    # Compare with results from Evj calculation in Herzberg notation
    from radis.db.molecules import Molecules

    CO_X = Molecules[molecule][isotope]["X"]

    # %% Calculate

    v, J = 0, 0
    if verbose:
        print(("Energies for CO(X) v, J=", v, J))

        # ... calculate with Molecular Data from JSON (new format)
        print(("... from JSON: {0:.5f} cm-1".format(EvJ(v, J, **dunham_coeffs))))

        # ... calculate with hardcoded Dunham expansion (legacy)
        print(
            (
                "... from hardcoded Herzberg constants {0:.5f} cm-1".format(
                    CO_X.Erovib(v, J, remove_ZPE=False)
                )
            )
        )

    import numpy as np

    assert np.isclose(EvJ(v, J, **dunham_coeffs), CO_X.Erovib(v, J, remove_ZPE=False))


def _run_all_tests(verbose=True, *args, **kwargs):

    test_EvJ(verbose=verbose, *args, **kwargs)
    test_dunham_co(verbose=verbose, *args, **kwargs)


if __name__ == "__main__":
    _run_all_tests()
