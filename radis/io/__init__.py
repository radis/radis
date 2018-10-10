# -*- coding: utf-8 -*-
"""

-------------------------------------------------------------------------------

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

    return (HITRAN_CLASS1+HITRAN_CLASS2+HITRAN_CLASS3+HITRAN_CLASS4+HITRAN_CLASS5 +
            HITRAN_CLASS6+HITRAN_CLASS7+HITRAN_CLASS8+HITRAN_CLASS9+HITRAN_CLASS10)


def _get_supported_molecules_nonequilibrium():
    ''' Molecules supported in RADIS non equilibrium calculations without need
    for extra databases. 
    Basically, molecules whose spectroscopic constants are built-in RADIS database
    (see radis.db)
    '''

    from radis.db.molecules import Molecules

    return list(Molecules.keys())


MOLECULES_LIST_EQUILIBRIUM = __supported_molecules_equilibrium__ = _get_supported_molecules_equilibrium()
''' list: molecules that can be calculated in RADIS at equilibrium. 
Basically, all HITRAN species are available. '''
MOLECULES_LIST_NONEQUILIBRIUM = __supported_molecules_nonequilibrium__ = _get_supported_molecules_nonequilibrium()
''' list: molecules that can be calculated in RADIS at nonequilibrium. Built-in
spectroscopic constants to calculate energy levels are needed 

See Also
--------

:py:data:`~radis.db.molecules.Molecules`
'''
