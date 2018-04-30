# -*- coding: utf-8 -*-
""" Functions to parse ~/.radis file

Notes
-----

@author: Erwan

Functions to parse ~/.radis file

Create a ~/.radis file in your HOME that contains all machine-specific information
(e.g: path to databanks)

"""

from __future__ import print_function, absolute_import, division, unicode_literals

from radis.misc.utils import FileNotFoundError, DatabankNotFound, configparser
from os.path import expanduser, join, exists, dirname
from six import string_types
from radis.misc.basics import compare_lists, compare_dict, stdpath

# %% Functions to parse ~/.radis file

DBFORMAT = ("""
--------------------------

[CDSD]                           #  your databank name
info = CDSD-HITEMP databank      #  whatever you want
path =                           #  no "", multipath allowed
       D:\Databases\CDSD-HITEMP\cdsd_hitemp_07
       D:\Databases\CDSD-HITEMP\cdsd_hitemp_08
       D:\Databases\CDSD-HITEMP\cdsd_hitemp_09
format = cdsd                    #  'hitran' or 'cdsd' (no ")
                                 # databank text file format. More info in
                                 # SpectrumFactory.load_databank function.
parfuncfmt:                      #  'cdsd', 'hapi', etc.
                                 # format to read tabulated partition function 
                                 # file. If `hapi`, then HAPI (HITRAN Python 
                                 # interface) is used to retrieve them (valid if
                                 # your databank is HITRAN data). HAPI is embedded 
                                 # into RADIS. Check the version.            
# Optional
# ----------
parfunc:                         #  path to tabulated partition function to use.
                                 # If `parfuncfmt` is `hapi` then `parfunc` 
                                 # should be the link to the hapi.py file. If 
                                 # not given, then the hapi.py embedded in RADIS 
                                 # is used (check version)
levels_iso1                      #  path to energy levels (needed for non-eq 
                                 # calculations). Default None
levels_iso2                      # etc
levels_iso4                      # etc
levelsfmt:                       #  'cdsd', etc. 
                                 # how to read the previous file. Default None.
levelsZPE:                       #  zero-point-energy (cm-1): offset for all level 
                                 # energies. Default 0 (if not given)

--------------------------""")

CONFIG_PATH = join(expanduser("~"),".radis")
    

def getConfig():
    ''' Read config file and returns it

    Config file name is harcoded: `~/.radis`
    '''

    config = configparser.ConfigParser()
    configpath = CONFIG_PATH

    # Test ~/.radis exists
    if not exists(configpath):

        raise FileNotFoundError("Create a `.radis` file in {0} to store links to ".format(
                                        dirname(configpath))+\
                                "your local databanks. Format must be:\n {0}".format(
                                        DBFORMAT)+\
                                "\n(it can be empty too)")
    config.read(configpath)

    return config
#

def getDatabankEntries(dbname):
    ''' Read ~/.radis config file and returns a dictionary of entries.


    Notes
    -----
    
    Databank format:
    
        [CDSD]                           # your databank name
        info = CDSD-HITEMP databank      # whatever you want
        path =                           # no "", multipath allowed
               D:\Databases\CDSD-HITEMP\cdsd_hitemp_07
               D:\Databases\CDSD-HITEMP\cdsd_hitemp_08
               D:\Databases\CDSD-HITEMP\cdsd_hitemp_09
        format = cdsd                    # 'hitran' or 'cdsd' (no ")
                                         # Databank text file format. More info in
                                         # SpectrumFactory.load_databank function.
    
        # Optional:
    
        parfunc                          # path or 'USE_HAPI'
                                         # path to tabulated partition functions. If
                                         # `USE_HAPI`, then HAPI (HITRAN Python
                                         interface) [1]_ is used to retrieve them (valid
                                         if your databank is HITRAN data). HAPI
                                         is embedded into RADIS. Check the version.
    
        parfuncfmt:                      # 'cdsd'
                                         # format to read tabulated partition function
                                         # file. If `USE_HAPI` is given as `parfunc`
                                         # parameter then this line should not be used.
    
        levels_iso1                      # path to energy levels (needed for non-eq)
                                         # calculations.
        levels_iso2                      # etc
        levels_iso4                      # etc
    
        levelsfmt                        # 'cdsd'
                                         # how to read the previous file.

    References
    ----------

    .. [1] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_


    '''

    config = getConfig()

    # Make sure it looks like a databank (path and format are given)
    try:
        config.get(dbname, 'path')
        config.get(dbname, 'format')
    except configparser.NoSectionError:
        msg = ("{1}\nDBFORMAT\n{0}\n".format(DBFORMAT, dbname)+\
               "No databank named {0} in `{1}`. ".format(dbname, CONFIG_PATH) +
               "Available databanks: {0}. ".format(getDatabankList())+\
               "See databank format above")
        raise DatabankNotFound(msg)

    entries = dict(config.items(dbname))
    
    # Parse paths correctly
    entries['path'] = entries['path'].strip('\n').split('\n')
    
    # Merge all isotope-dependant levels into one dict
    iso_list = [k for k in entries if k.startswith('levels_iso')]
    if len(iso_list) > 0:
        levels = {}
        for k in iso_list:
            iso = float(k[10:])
            levels[iso] = entries.pop(k)
        entries['levels'] = levels
    else:
        if 'levels' in entries:
            raise SyntaxError('in {0}: `levels` replaced with '.format(CONFIG_PATH)+\
                              '`levels_iso#` where # is the isotope number')

    return entries

def getDatabankList():
    ''' Get all databanks available in ~/.radis'''

    config = getConfig()

    # Get databank path and format
    validdb = []
    for dbname in config.sections():
        try:
            config.get(dbname, 'path')
            config.get(dbname, 'format')
        except configparser.NoSectionError:
            # not a db
            continue
        except configparser.NoOptionError:
            # not a db
            continue            
        # looks like a db. add to db
        validdb.append(dbname)
    return validdb

def addDatabankEntries(dbname, dict_entries, verbose=True):
    ''' Add database dbname with entries from dict_entries. If database 
    already exists in ~/.radis, raises an error
    '''
        
    # Get ~/.radis if exists, else create it
    try:
        dbnames = getDatabankList()
    except FileNotFoundError:
        # generate ~/.radis:
        dbnames = []
        open(CONFIG_PATH, 'a').close()
        if verbose: print('Created ~/.radis in {0}'.format(dirname(CONFIG_PATH)))
        
    # Check database doesnt exist
    if dbname in dbnames:
        raise ValueError('Database already exists: {0}'.format(dbname)+\
                         '. Cant add it')
        
    # Add entries to parser
    config = configparser.ConfigParser()
    config[dbname] = {}
    
    if 'info' in dict_entries: # optional
        config[dbname]['info'] = dict_entries.pop('info')
    
    # ... Parse paths correctly
    if dict_entries['path'] in string_types:
        config[dbname]['path'] = dict_entries.pop('path')
    else:  # list 
        config[dbname]['path'] = '\n       '.join(dict_entries.pop('path'))
    
    config[dbname]['format'] = dict_entries.pop('format')
    config[dbname]['parfuncfmt'] = dict_entries.pop('parfuncfmt')
    
    # Optional:
    # ... Split all isotopes in separate keys
    levels_dict = dict_entries.pop('levels', {})
    for iso, levels_iso in levels_dict.items():
        dict_entries['levels_iso{0}'.format(iso)] = levels_iso
        
    if 'levelsfmt' in dict_entries:
        config[dbname]['levelsfmt'] = dict_entries.pop('levelsfmt')
    
    # Check nothing is left
    if dict_entries != {}:
        raise ValueError('Unexpected keys: {0}'.format(list(dict_entries.keys())))
        
    # Write to ~/.radis
    # ... Note: what if there is a PermissionError here? Try/except pass? 
    with open(CONFIG_PATH, 'a') as configfile:
        configfile.write('\n')
        config.write(configfile)
    if verbose: print("Added {0} database in {1}".format(dbname, CONFIG_PATH))

    return

def diffDatabankEntries(dict_entries1, dict_entries2, verbose=True):
    ''' Compare two Databank entries under dict format (i.e: output of 
    getDatabankEntries)
    
    Returns None if no differences are found, or the first different key 
    '''
    
    k = None
    
    try:
        assert len(dict_entries1) == len(dict_entries2)
        assert compare_lists(list(dict_entries1.keys()), list(dict_entries2.keys()),
                             verbose=verbose) == 1
        for k in dict_entries1.keys():
            v1 = dict_entries1[k]
            v2 = dict_entries2[k]
            if k in ['info', 'format', 'parfuncfmt', 'levelsfmt']:
                assert v1 == v2
            elif k in ['path']:
                assert compare_lists([stdpath(path1) for path1 in v1],
                                     [stdpath(path2) for path2 in v2],
                                     verbose=verbose) == 1
            elif k in ['levels']:
                assert compare_dict(v1, v2, compare_as_paths=list(v1.keys()),
                                    verbose=verbose) == 1
            else:
                raise ValueError('Unexpected key:', k)
                
        return None
    
    except AssertionError:
        if verbose: print('Key doesnt match:', k)
        return k
    
def printDatabankEntries(dbname, crop=200):
    ''' Print databank info
    
    
    Parameters    
    ----------
    
    dbname: str
        database name in ~/.radis
        
    crop: int
        if > 0, cutoff entries larger than that
        
    '''
    entries = getDatabankEntries(dbname)
    print(dbname,'\n-------')
    for k, v in entries.items():
        # Add extra arguments
        args = []
        if k == 'levelszpe':
            args.append('cm-1')
        v = '{0}'.format(v)
        if len(v) > crop and crop > 0:
            v = v[:crop] + '...'
        # Print item
        print(k,':',v, *args)

def printDatabankList():
    ''' Print all databanks available in ~/.radis '''
    try:
        print('Databanks in {0}: '.format(CONFIG_PATH), ','.join(getDatabankList()))
        for dbname in getDatabankList():
            print('\n')
            printDatabankEntries(dbname)
    except FileNotFoundError:
        print('No config file {0}'.format(CONFIG_PATH))
        # it's okay

# %% Test

def _test(*args, **kwargs):

    printDatabankList()

    return True

if __name__=='__main__':
    _test()
