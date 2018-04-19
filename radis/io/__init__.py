# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:52:15 2015

Erwan Pannier
EM2C, CentraleSupélec, 2015
CNRS UPR 288

"""

from __future__ import absolute_import, division, print_function, unicode_literals

from .hitran import hit2df, get_molecule, get_molecule_identifier
from .cdsd import cdsd2df


# %% Get list of supported molecules
def _get_supported_molecules_equilibrium():
    ''' Molecules supported in RADIS equilibrium calculations
    Basically, all HITRAN species
    '''
    # Add all HITRAN species
    from radis.io.hitran import (HITRAN_CLASS1, HITRAN_CLASS2, HITRAN_CLASS3,
                                 HITRAN_CLASS4, HITRAN_CLASS5, HITRAN_CLASS6,
                                 HITRAN_CLASS7, HITRAN_CLASS8, HITRAN_CLASS9,
                                 HITRAN_CLASS10)
    
    return (HITRAN_CLASS1+HITRAN_CLASS2+HITRAN_CLASS3+HITRAN_CLASS4+HITRAN_CLASS5+
            HITRAN_CLASS6+HITRAN_CLASS7+HITRAN_CLASS8+HITRAN_CLASS9+HITRAN_CLASS10)
    
def _get_supported_molecules_nonequilibrium():
    ''' Molecules supported in RADIS non equilibrium calculations without need
    for extra databases. 
    Basically, molecules whose spectroscopic constants are built-in RADIS database
    (see radis.db)
    '''
    
    # Hardcoded for the moment. 
    # TODO Look up radis.db once it's merged here
    
    return ['CO', 'CO2']
    
    
__supported_molecules_equilibrium__ = _get_supported_molecules_equilibrium()
__supported_molecules_nonequilibrium__ = _get_supported_molecules_nonequilibrium()
