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
    Basically, all [HITRAN-2016]_ species
    '''
    # Add all HITRAN species
    from radis.io.hitran import HITRAN_MOLECULES

    return HITRAN_MOLECULES


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
All [HITRAN-2016]_ species are available:

- 1 	``'H2O'`` : 	Water 	

.. image:: https://github.com/radis/radis-examples/blob/master/hitran_spectra/out/0%20-%20H2O%20infrared%20spectrum.png
    :width: 600
    :target: https://github.com/radis/radis-examples/blob/master/hitran_spectra/out/0%20-%20H2O%20infrared%20spectrum.png
    :alt: Water H2O infrared spectrum

- 2 	``'CO2'`` : 	Carbon Dioxide 	
- 3 	``'O3'`` : 	Ozone 	
- 4 	``'N2O'`` : 	Nitrogen oxide 	
- 5 	``'CO'`` : 	Carbon Monoxide 	
- 6 	``'CH4'`` : 	Methane 	
- 7 	``'O2'`` : 	Oxygen 	
- 8 	``'NO'`` : 	Nitric Oxide 	
- 9 	``'SO2'`` : 	Sulfur Dioxide 	
- 10 	``'NO2'`` : 	Nitrogen Dioxide 	
- 11 	``'NH3'`` : 	Ammonia 	
- 12 	``'HNO3'`` : 	Nitric Acid 	
- 13 	``'OH'`` : 	Hydroxyl 	
- 14 	``'HF'`` : 	Hydrogen Fluoride 	
- 15 	``'HCl'`` : 	Hydrogen Chloride 	
- 16 	``'HBr'`` : 	Hydrogen Bromide 	
- 17 	``'HI'`` : 	Hydrogen Iodide 	
- 18 	``'ClO'`` : 	Chlorine Monoxide 	
- 19 	``'OCS'`` : 	Carbonyl Sulfide 	
- 20 	``'H2CO'`` : 	Formaldehyde 	
- 21 	``'HOCl'`` : 	Hypochlorous Acid 	
- 22 	``'N2'`` : 	Nitrogen 	
- 23 	``'HCN'`` : 	Hydrogen Cyanide 	
- 24 	``'CH3Cl'`` : 	Methyl Chloride 	
- 25 	``'H2O2'`` : 	Hydrogen Peroxide 	
- 26 	``'C2H2'`` : 	Acetylene 	
- 27 	``'C2H6'`` : 	Ethane 	
- 28 	``'PH3'`` : 	Phosphine 	
- 29 	``'COF2'`` : 	Carbonyl Fluoride 	
- 30 	``'SF6'`` : 	Sulfur Hexafluoride 	
- 31 	``'H2S'`` : 	Hydrogen Sulfide 	
- 32 	``'HCOOH'`` : 	Formic Acid 	
- 33 	``'HO2'`` : 	Hydroperoxyl 	
- 34 	``'O'`` : 	Oxygen Atom 	
- 35 	``'ClONO2'`` : 	Chlorine Nitrate 	
- 36 	``'NO+'`` : 	Nitric Oxide Cation 	
- 37 	``'HOBr'`` : 	Hypobromous Acid 	
- 38 	``'C2H4'`` : 	Ethylene 	
- 39 	``'CH3OH'`` : 	Methanol 	
- 40 	``'CH3Br'`` : 	Methyl Bromide 	
- 41 	``'CH3CN'`` : 	Acetonitrile 	
- 42 	``'CF4'`` : 	CFC-14 	
- 43 	``'C4H2'`` : 	Diacetylene 	
- 44 	``'HC3N'`` : 	Cyanoacetylene 	
- 45 	``'H2'`` : 	Hydrogen 	
- 46 	``'CS'`` : 	Carbon Monosulfide 	
- 47 	``'SO3'`` : 	Sulfur trioxide 	
- 48 	``'C2N2'`` : 	Cyanogen 	
- 49 	``'COCl2'`` : 	Phosgene 	

See Also
--------

:py:data:`~radis.db.molecules.Molecules`,
:py:func:`~radis.db.molecules.getMolecule`

'''

MOLECULES_LIST_NONEQUILIBRIUM = __supported_molecules_nonequilibrium__ = _get_supported_molecules_nonequilibrium()
''' list: molecules that can be calculated in RADIS at nonequilibrium. 
Spectroscopic constants to calculate energy levels are needed. 

RADIS features some built-in :ref:`spectroscopic constants <label_db_spectroscopic_constants>` 
for the following species ([HITRAN-2016]_ nomenclature):

- 2 	``'CO2'`` : 	Carbon Dioxide 	
- 5 	``'CO'`` : 	Carbon Monoxide 	


See Also
--------

:py:data:`~radis.db.molecules.Molecules`,
:py:func:`~radis.db.molecules.getMolecule`
'''
