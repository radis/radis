# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:47:36 2015

@author: Erwan
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os
import json
from collections import OrderedDict


def getFile(*relpath):
    ''' Converts the relative path of a database file in a the full path.
    Used by processing script not to worry about where the database is stored

    Examples of use::
        
        radis.db.getFile('CN','CN_Violet_vib.dat')
        radis.db.getFile('CN\CN_Violet_vib.dat')    
        radis.db.getFile('CN\\CN_Violet_vib.dat')

    '''
    
#    return os.path.join(os.path.dirname(__file__), *relpath)
    from radis.misc.utils import getProjectRoot
    return os.path.join(getProjectRoot(), 'db', *relpath)



def check_molecule_data_structure(fname, verbose=True):
    ''' Check that ``fname`` has a valid JSON structure for molecular data
    
    Parameters
    ----------
    
    fname: str
        molecular data JSON file 
        
        
    Notes
    -----
    
    Order in the json doesnt matter, however, the order of the 
    electronic_level_names matters: the index of all levels should
    match the `index` key of these states. 
    
    '''

    with open(fname) as f:
        db = json.load(f) #, object_pairs_hook=OrderedDict)
            
    for molecule, db_molec in db.items():
        # ... Check number of isotopes is correct
        isotope_names = db_molec['isotopes_names']
        isotopes = db_molec['isotopes']
        if len(isotopes) != len(isotope_names):
            raise ValueError(
                    'In molecule {0}: isotope names '.format(molecule)+\
                    '({0}) dont match the number of isotopes ({1})'.format(
                    isotope_names, list(isotopes.keys())))
            
        # ... Check number of electronic states is correct
        for isotope, db_iso in db_molec['isotopes'].items():
            elec_states_names = db_iso['electronic_levels_names']
            elec_states = db_iso['electronic_level']
            if len(elec_states_names) != len(elec_states):
                raise ValueError(
                    'In molecule {0}, isotope {1}: electronic '.format(molecule, isotope)+\
                    'levels names ({0}) dont match the number of levels ({1})'.format(
                    elec_states_names, list(elec_states.keys())))
        
        # ... Check they are properly ordered
            for state, db_state in elec_states.items():
                if elec_states_names.index(db_state['name'])+1 != db_state['index']:  
                    # ... + 1 because Python index start at 0 and FORTRAN at 1 unless allocated with 0:N
                    raise ValueError(
                        'In molecule {0}, isotope {1}: index of electronic '.format(
                                molecule, isotope, db_state['index'])+\
                        'state {0} ({1}): {2} does not match the list of states: {3}. '.format(
                        state, db_state['name'], db_state['index'], elec_states_names))
        
    if verbose: print('Structure of {0} looks correct'.format(fname))



def get_dunham_coefficients(molecule, isotope, electronic_state, jsonfile='default'):
    ''' Returns Dunham coefficients ``Yij`` for ``molecule``, ``isotope``, ``electronic_state`` 
    by parsing a JSON file of molecule data. 
    
    Dunham coefficients are identified as starting with ``Y`` (i.e. ``alpha_e``
    won't be recognized)
    
    Parameters
    ----------
    
    molecule: str
        molecule name
        
    isotope: int
        isotope number
        
    electronic_state: str
        electronic state name
        
    jsonfile: str, or ``default``
        path to json file. If ``default``, the ``molecules_data`` JSON file 
        in the ``radis.db`` database is used::
        
            radis\db\[molecule]\molecules_data.json
    
    '''

    if jsonfile == 'default':
        jsonfile = getFile('{0}/molecules_data.json'.format(molecule))
        
    check_molecule_data_structure(jsonfile, verbose=False)
    
    with open(jsonfile) as f:
        db = json.load(f, object_pairs_hook=OrderedDict)

    # Get Dunham coefficients in 001 state (X)
    elec_state_names = db[molecule]['isotopes'][str(isotope)]['electronic_levels_names']
    
    try:
        elec_state_index = elec_state_names.index(electronic_state)
    except ValueError:
        raise ValueError("{0} not in the electronic state list for {1}(iso={2}): {3}".format(
                electronic_state, molecule, isotope, elec_state_names))
    else:
        elec_state_index = '{:03d}'.format(elec_state_index+1)  # 1 based index
    
    dunham_coeffs = db[molecule]['isotopes'][str(isotope)]['electronic_level'][elec_state_index]
    # Only get Dunham coeffs, i.e, these that start with Y
    dunham_coeffs = {k:v for (k, v) in dunham_coeffs.items() if k.startswith('Y')}

    if len(dunham_coeffs) == 0:
        raise ValueError('No Dunham coefficients found for {0}{1}(iso={2})'.format(
                molecule, isotope, electronic_state))

    return dunham_coeffs

def get_herzberg_coefficients(molecule, isotope, electronic_state, jsonfile='default'):
    ''' Returns spectroscopic coefficients with Herzberg conventions for 
    ``molecule``, ``isotope``, ``electronic_state``  by parsing a JSON file of molecule data. 
    
    Herzberg coefficients are the usual:
        
    .. math::
        
        \\omega_e, \\alpha_e, B_e, D_e, etc.
    
    The extensive list of parameters considered as Herzberg coefficients is found in 
    :py:data:`~radis.db.conventions.herzberg_coefficients`
    
    Parameters
    ----------
    
    molecule: str
        molecule name
        
    isotope: int
        isotope number
        
    electronic_state: str
        electronic state name
        
    jsonfile: str, or ``default``
        path to json file. If ``default``, the ``molecules_data`` JSON file 
        in the ``radis.db`` database is used::
        
            radis\db\[molecule]\molecules_data.json
    
    See Also
    --------
    
    :py:data:`~radis.db.conventions.herzberg_coefficients`
    
    '''
    
    from radis.db.conventions import herzberg_coefficients

    if jsonfile == 'default':
        jsonfile = getFile('{0}/molecules_data.json'.format(molecule))
        
    check_molecule_data_structure(jsonfile, verbose=False)
    
    with open(jsonfile) as f:
        db = json.load(f, object_pairs_hook=OrderedDict)

    # Get Dunham coefficients in 001 state (X)
    elec_state_names = db[molecule]['isotopes'][str(isotope)]['electronic_levels_names']
    
    try:
        elec_state_index = elec_state_names.index(electronic_state)
    except ValueError:
        raise ValueError("{0} not in the electronic state list for {1}(iso={2}): {3}".format(
                electronic_state, molecule, isotope, elec_state_names))
    else:
        elec_state_index = '{:03d}'.format(elec_state_index+1)  # 1 based index
    
    def remove_cm1(coef):
        ''' Remove the trailing '_cm-1' in the spectroscopic coefficient name, 
        if defined '''
        if coef.endswith('_cm-1'):
            coef = coef[:-5]
        return coef
        
    herzberg_coeffs = {remove_cm1(k):v for (k,v) in 
                     db[molecule]['isotopes'][str(isotope)]['electronic_level'][elec_state_index].items()}
#    # Only get Herzberg coeffs, i.e, these that are defined in :py:data:`~radis.db.conventions.herzberg_coefficients`
    
    def ignore_trailing_number(coef):
        ''' Used so that ``wexe1`` matches ``wexe`` as a well defined 
        Herzberg coefficient '''
        if str.isdigit(coef[-1]):
            coef = coef[:-1]
        return coef
    
    herzberg_coeffs = {k:v for (k, v) in herzberg_coeffs.items() if ignore_trailing_number(k) in herzberg_coefficients}

    if len(herzberg_coeffs) == 0:
        raise ValueError('No Herzberg spectroscopic coefficients found for {0}{1}(iso={2})'.format(
                molecule, isotope, electronic_state))

    return herzberg_coeffs
