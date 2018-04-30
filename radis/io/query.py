#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Wrapper to fetch line database from HITRAN with Astroquery [1]_ 

Notes
-----

Astroquery [R1]_ is itself based on HAPI [R2]_ 
    
Created on Sun Feb 18 18:29:58 2018

@author: erwan
 

References
----------

.. [R1] `Astroquery <https://astroquery.readthedocs.io>`_ 

.. [R2] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_

"""

from __future__ import absolute_import
from __future__ import print_function
from radis.io.hitran import get_molecule, get_molecule_identifier
from radis.misc import is_float
from astroquery.hitran import read_hitran_file, cache_location, download_hitran
from os.path import join
import numpy as np
import pandas as pd

def fetch_astroquery(molecule, isotope, wmin, wmax, verbose=True):
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
        
    Other Parameters
    ----------------
    
    verbose: boolean
        Default True
        
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

    empty_range = False

    # Download using the astroquery library
    try:
        download_hitran(mol_id, isotope, wmin, wmax)
    except:
        # Maybe there are just no lines for this species in this range
        # In that case we usually end up with errors like:
        
            # (<class 'Exception'>, Exception('Query failed: 404 Client Error: 
            # Not Found for url: http://hitran.org/lbl/api?numax=25000&numin=19000&iso_ids_list=69\n',), 
            # <traceback object at 0x7f0967c91708>)
            
        import sys
        _err_class, _err_details, _err_obj = sys.exc_info()
        
        if 'Not Found for url:' in str(_err_details):
            # Let's bet it's just that there are no lines in this range
            empty_range = True
            if verbose: 
                print(('Not lines for {0} (id={1}), iso={2} between {3:.2f}-{4:.2f}cm-1. '.format(
                        molecule, mol_id, isotope, wmin, wmax)))
        else:
            raise ValueError('An error occured during the download of HITRAN files '+\
                             'for {0} (id={1}), iso={2} between {3:.2f}-{4:.2f}cm-1. '.format(
                                     molecule, mol_id, isotope, wmin, wmax)+\
                             'See details of the error below: {0}'.format(_err_details))
                      
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
    
    if not empty_range:
        tbl = read_hitran_file(join(cache_location, molecule+'.data'))
        df = tbl.to_pandas()
        df = df.rename(columns=rename_columns)
    else:
        df = pd.DataFrame(columns=list(rename_columns.values()))
    
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

if __name__ == '__main__':
    from radis.test.io.test_query import _run_testcases
    _run_testcases(verbose=True)