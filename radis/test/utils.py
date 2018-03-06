# -*- coding: utf-8 -*-
""" Tools to test RADIS library

Notes
-----

Created on Thu May 28 14:47:36 2015

@author: Erwan

Tools to test RADIS library

Todo
---- 

add all files in a separetaly downloadble folder. getTestFile should check
the file exist else suggest to download the folder

Examples
--------

Run all tests:
    
    cd radis/test
    pytest
    
Run only "fast" tests (tests that have "fast" in their name, and should be 
a few seconds only):
    
    cd radis/test
    pytest -k fast
    
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os
from warnings import warn
from radis.misc.config import (getDatabankList, getDatabankEntries, addDatabankEntries,
                             diffDatabankEntries)
from radis.misc.utils import FileNotFoundError
from radis.misc.basics import compare_dict, compare_paths
from os.path import join, dirname

TEST_FOLDER_PATH = join(dirname(dirname(__file__)), 'test')
IGNORE_MISSING_DATABASES = True        # TODO: move in ~/.radis

def getTestFile(file):
    ''' Return the full path of a test file. Used by test functions not to
    worry about the project architecture'''

    return join(TEST_FOLDER_PATH, 'files', file)

def getValidationCase(file):
    ''' Return the full path of a validation case file. Used by test functions not to
    worry about the project architecture'''

    return join(TEST_FOLDER_PATH, 'validation', file)

try: # Python 3.6 only 
    getTestFile.__annotations__['file'] = os.listdir(join(TEST_FOLDER_PATH, 'files'))
    getValidationCase.__annotations__['file'] = os.listdir(join(TEST_FOLDER_PATH, 'validation'))
except: 
    pass

# %% Comparison functions

def testEqual(a, b, info=''):
    if a != b:
        print('Mismatch', info, ':', a, '!=', b)
    return a == b

# %% Test Databases

TEST_DATABASES = {
        'HITRAN-CO2-TEST':{
            'info':'HITRAN 2016 database, CO2, 1 main isotope (CO2-626), bandhead: '+\
                   '2380-2398 cm-1 (4165-4200 nm)',
            'path':[getTestFile(r'hitran_co2_626_bandhead_4165_4200nm.par')],
            'format':'hitran',
            'parfuncfmt':'hapi',
            'levelsfmt':'neq', # TODO: replace with 'radis'
            }, 
        'HITRAN-CO-TEST':{
            'info':'HITRAN 2016 database, CO, 3 main isotopes (CO-26, 36, 28), '+\
                   '2000-2300 cm-1',
            'path':[getTestFile(r'hitran_co_3iso_2000_2300cm.par')],
            'format':'hitran',
            'parfuncfmt':'hapi',
            'levelsfmt':'neq', # TODO: replace with 'radis'
            },
        }

# %% Utils to test spec module

def build_test_databases(verbose=True):
    ''' Build test databases and add them in ~/.radis. Generate the file if it 
    doesnt exist
    
    In particular:
    
    - HITRAN-CO2-TEST: CO2, HITRAN 2016, 4165-4200 nm 
    - HITRAN-CO-TEST: CO, HITRAN 2016, 2000-2300 cm-1
    
    These test databases are used to run the different test routines. They can
    obviously be used by Users to run simulations, but we suggest Users to download
    their own line databases files and add them to ~/.radis so they have more control
    on it
    
    '''
    
    # Get list of databases
    try:
        dbnames = getDatabankList()
    except FileNotFoundError:
        dbnames = []
        
    # %% Add test databases
    
    def add_to_parser(config, name, dic):
        for k, v in dic.items():
            config[name][k] = v
        if verbose: print("Adding '{0}' database in ~/.radis".format(name))
        
    for dbname, dbentries in TEST_DATABASES.items():
        
        if dbname in dbnames:  # Check entries are correct
#            for k 
            diff = diffDatabankEntries(getDatabankEntries(dbname), dbentries,
                                       verbose=False)
            if diff is not None:
                raise ValueError('{0}'.format(diff)+\
                                 '\nIn ~/.radis\n----------\n{0}'.format(getDatabankEntries(dbname))+\
                                 '\n\nExpected\n---------\n{0}\n\n'.format(dbentries)+\
                                 'Test Database {0} doesnt match expected '.format(dbname)+\
                                 'entries for key `{0}`. See comparison above'.format(diff))

        else:  #  add them (create ~/.radis file if doesnt exist yet)
            addDatabankEntries(dbname, dbentries)
    
    return


# %% Deal with missing databases
def _failsafe_if_no_db(testcase, *args, **kwargs):
    '''finally not implemented?'''
    from radis.misc.utils import DatabankNotFound
    try:
        testcase(*args, **kwargs)
    except DatabankNotFound:
        import sys
        print((sys.exc_info()))
        print(('Testing {0}: Database not defined. \n'.format(testcase.__name__)+\
                       'Ignoring the test'))
        return True
    
def IgnoreMissingDatabase(err, file='', warnings=True):
    if IGNORE_MISSING_DATABASES:
        if warnings:
            import sys
            print(sys.exc_info())
            print('In {0}: Database not defined: {1}'.format(file, err.filename)+\
                  '\n Ignoring the test')
        return True
    else:
        raise err

if __name__ == '__main__':
#    run_tests()
    build_test_databases()