# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 18:00:54 2018

@author: erwan
"""

import numpy as np
from radis.db.molecules import getMolecule
from radis.db.utils import get_herzberg_coefficients
import pytest

@pytest.mark.fast       # this is a fast test. Run fast tests only with 'pytest -m fast'    
def test_herzberg_coefficients_CO2_json(verbose=True, *args, **kwargs):
    ''' Check that default coefficients are found correctly
    '''
    # @dev: update this test if we ever are to change the default Herzberg
    # coeffs in CO2/molecules_data.json
    
    assert set(get_herzberg_coefficients('CO2', 1, 'X1SIGu+')) == set(
            ['we1', 'we2', 'we3', 'wexe1', 'wexe2', 'wexe3', 'Be', 'De', 'He'])

    if verbose:
        print('Herzberg coefficients defined for CO2(X)(iso=1):', get_herzberg_coefficients('CO2', 1, 'X1SIGu+'))


@pytest.mark.fast       # this is a fast test. Run fast tests only with 'pytest -m fast'    
def test_getMolecule(verbose=True, *args, **kwargs):

    # Input
    molec = 'CO'
    iso = 1
    state = 'X'

    # get predefined
    S = getMolecule(molec, iso, state)

    # Print CO Dunham development calculation (can be accessed from console)
    if verbose:
        print('Dunham development')
        print('------------\n')
    S.get_ref('Erovib')

    return True

# Test defined in neq 0.9.24. Deactivated in the move to RADIS 1.0
#@pytest.mark.fast       # this is a fast test. Run fast tests only with 'pytest -m fast'    
#def test_CO_energies_Herzberg_vs_Dunham(verbose=True, *args, **kwargs):
#
#    from radis.db.molecules import CO_X_iso1_Herzberg, CO_X_iso1
#    
#    for (v, J) in [(0, 0), (1, 10), (10, 10)]:
#        if verbose:
#            print('v={0}, J={1}'.format(v,J))
#            print('CO_X_iso1_Herzberg: {0:.3f}cm-1'.format(CO_X_iso1_Herzberg.Erovib(v,J)))
#            print('CO_X_iso1_from_json: {0:.3f}cm-1'.format(CO_X_iso1.Erovib(v,J)))
#        assert np.isclose(CO_X_iso1_Herzberg.Erovib(v,J), CO_X_iso1.Erovib(v,J),
#                          rtol=1e-4)


def _run_testcases(verbose=True, *args, **kwargs):

    test_herzberg_coefficients_CO2_json(verbose=verbose, *args, **kwargs)
    test_getMolecule(verbose=verbose, *args, **kwargs)
#    test_CO_energies_Herzberg_vs_Dunham(verbose=verbose, *args, **kwargs)

    return True



if __name__ == '__main__':

    _run_testcases()
