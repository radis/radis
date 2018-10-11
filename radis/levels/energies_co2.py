# -*- coding: utf-8 -*-
"""
Functions specific to the calculation of CO2 energy levels

Notes
-----

@dev: you can define Python functions here, with hardcoded spectroscopic coefficients. 
Or, you can use the spectroscopic coefficients given in ``radis/db/CO2`` under
JSON format. 

In both cases, you will need to overwrite the :meth:`~radis.db.classes.ElectronicState.Erovib`
method. See ``radis.db.molecules.py``


----------------------------------------------------------------------------------------------

"""


from __future__ import absolute_import
from radis.levels.vibrating_rotor import EvJ_uncoupled_vibrating_rotor, EvJah_uncoupled_vibrating_rotor

def EvJ_co2(v1, v2, l2, v3, J, coeff_dict, remove_ZPE=True):
    ''' Rovibrational energies.
    See :py:func:`~radis.levels.vibrating_rotor.EvJ_uncoupled_vibrating_rotor` 
    '''
    return EvJ_uncoupled_vibrating_rotor(v1, v2, l2, v3, J, coeff_dict=coeff_dict,
                                         gv1=1, gv2=2, gv3=1, remove_ZPE=remove_ZPE)


def EvJah_co2(v1, v2, l2, v3, J, coeff_dict, remove_ZPE=True):
    ''' Harmonic and anharmonic parts of rovibrational energies.
        See :py:func:`~radis.levels.vibrating_rotor.EvJah_uncoupled_vibrating_rotor` 
    '''
    return EvJah_uncoupled_vibrating_rotor(v1, v2, l2, v3, J, coeff_dict=coeff_dict,
                                         remove_ZPE=remove_ZPE)
