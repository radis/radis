#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Wrapper to fetch line database from HITRAN with Astroquery [1]_ 

Notes
-----

Astroquery [1]_ is itself based on HAPI [2]_ 
    
Created on Sun Feb 18 18:29:58 2018

@author: erwan
 

References
----------

.. [1] `Astroquery <https://astroquery.readthedocs.io>`_ 

.. [2] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_

"""

from radis.io.hitran import get_molecule, get_molecule_identifier
from radis.misc import is_float
from astroquery.hitran import read_hitran_file, cache_location, download_hitran
from os.path import join
import numpy as np

def fetch_astroquery(molecule, isotope, wmin, wmax):
    ''' Wrapper to Astroquery [1]_ fetch function to download a line database 
    
    Notes
    -----
    
    Astroquery [1]_ is itself based on HAPI [2]_ 
    
    Parameters
    ----------
    
    molecule: str, or int
        molecule name or identifier
        
    isotope: int
        isotope number 
        
    wmin, wmax: float  (cm-1)
        wavenumber min and max 
        
    References
    ----------
    
    .. [1] `Astroquery <https://astroquery.readthedocs.io>`_ 
    
    .. [2] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_
    
    
    '''
    
    # Check input
    if not is_float(molecule):
        mol_id = get_molecule_identifier(molecule)
    else:
        mol_id = molecule
        molecule = get_molecule(mol_id)
    assert is_float(isotope)

    # Download using the astroquery library
    download_hitran(mol_id, isotope, wmin, wmax)
    tbl = read_hitran_file(join(cache_location, molecule+'.data'))
    
    df = tbl.to_pandas()

    # Rename columns from Astroquery to RADIS format 
    rename_columns = {'molec_id':'id',
                      'local_iso_id':'iso',
                      'nu':'wav',
                      'sw':'int',
                      'a':'A',
                      'gamma_air':'airbrd',
                      'gamma_self':'selbrd',
                      'elower':'El',
                      'n_air':'Tdpair',
                      'delta_air':'Pshft',
                      'global_upper_quanta':'globu',
                      'global_lower_quanta':'globl',
                      'local_upper_quanta':'locu',
                      'local_lower_quanta':'locl',
                      'line_mixing_flag':'lmix',
                      'gp':'gp',
                      'gpp':'gpp',
                      }
    df = df.rename(columns=rename_columns)
    
    # Cast type to float64
    cast_type = {'wav':np.float64,
                 'int':np.float64,
                 'A':np.float64,
                 'airbrd':np.float64,
                 'selbrd':np.float64,
                 'El':np.float64,
                 'Tdpair':np.float64,
                 'Pshft':np.float64,
                 }
    for c, typ in cast_type.items():
        df[c] = df[c].astype(typ)

    return df

