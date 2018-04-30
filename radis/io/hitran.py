# -*- coding: utf-8 -*-
'''
HITRAN database parser 
'''

from __future__ import print_function, absolute_import, division, unicode_literals

import sys
import pandas as pd
import numpy as np
from collections import OrderedDict
import os
from os.path import exists, splitext
from six.moves import range
from six.moves import zip


# %% Hitran groups and classes
# As defined in Rothman et al, "The HITRAN 2004 molecular spectroscopic database"
# Tables 3 and 4

# Groups (define local quanta)

HITRAN_GROUP1 = ['H2O', 'O3', 'SO2', 'NO2', 'HNO3', 'H2CO', 'HOCl', 'H2O2', 'COF2',
               'H2S', 'HO2', 'HCOOH', 'ClONO2', 'HOBr', 'C2H4'] # asymmetric rotors

HITRAN_GROUP2 = ['CO2', 'N2O', 'CO', 'HF', 'HCl', 'HBr', 'HI', 'OCS', 'N2', 'HCN',
                 'C2H2', 'NO+']   # diatomic and linear molecules

HITRAN_GROUP3 = ['SF6', 'CH4']   # Spherical rotors

HITRAN_GROUP4 = ['CH3D', 'CH3Cl', 'C2H6', 'NH3', 'PH3', 'CH3OH']  # symmetric rotors

HITRAN_GROUP5 = ['O2']  # Triplet-Sigma ground electronic states

HITRAN_GROUP6 = ['NO', 'OH', 'ClO']  # Doublet-Pi ground electronic states 

# Classes (define global quanta)

HITRAN_CLASS1 = ['CO', 'HF', 'HCl', 'HBr', 'HI', 'N2', 'NO+']

HITRAN_CLASS2 = ['O2']  # Diatomic molecules with different electronic levels

HITRAN_CLASS3 = ['NO', 'OH', 'ClO']  # Diatomic molecules with doublet-Pi electronic state
    
HITRAN_CLASS4 = ['N2O', 'OCS', 'HCN']  # Linear triatomic

HITRAN_CLASS5 = ['CO2']  # Linear triatomic with large Fermi resonance

HITRAN_CLASS6 = ['H2O', 'O3', 'SO2', 'NO2', 'HOCl', 'H2S', 'HO2', 'HOBr'] # Non-linear triatomic

HITRAN_CLASS7 = ['C2H2']  # Linear tetratomic

HITRAN_CLASS8 = ['NH3', 'PH3']  # Pyramidal tetratomic

HITRAN_CLASS9 = ['H2CO', 'H2O2', 'COF2']  # Non-linear tetratomic

HITRAN_CLASS10 = ['CH4', 'CH3D', 'CH3Cl', 'C2H6', 'HNO3', 'SF6', 'HCOOH', 
                  'ClONO2', 'C2H4', 'CH3OH']  # Pentatomic or greater polyatomic

# %% Parsing functions

# General case
columns_2004 = OrderedDict([ (
               # name    # format # type  # description                                 # unit 
               'id',     ('a2',   int,   'Molecular number'                               ,''                      )),(
               'iso',    ('a1',   int,   'isotope number'                                 ,''                      )),(
               'wav',    ('a12',  float, 'vacuum wavenumber'                              ,'cm-1'                  )),(
               'int',    ('a10',  float, 'intensity at 296K'                              ,'cm-1/(molecule/cm-2)', )),(
               'A',      ('a10',  float, 'Einstein A coefficient'                         ,'s-1'                   )),(
               'airbrd', ('a5',   float, 'air-broadened half-width at 296K'               ,'cm-1.atm-1'            )),(
               'selbrd', ('a5',   float, 'self-broadened half-width at 296K'              ,'cm-1.atm-1'            )),(
               'El',     ('a10',  float, 'lower-state energy'                             ,'cm-1'                  )),(
               'Tdpair', ('a4',   float, 'temperature-dependance exponent for Gamma air'  ,''                      )),(
               'Pshft',  ('a8',   float, 'air pressure-induced line shift at 296K'        ,'cm-1.atm-1'            )),( 
               'globu',  ('a15',  str,   'electronic and vibrational global upper quanta' ,''                      )),(
               'globl',  ('a15',  str,   'electronic and vibrational global lower quanta' ,''                      )),(
               'locu',   ('a15',  str,   'electronic and vibrational local upper quanta'  ,''                      )),(
               'locl',   ('a15',  str,   'electronic and vibrational local lower quanta'  ,''                      )),(
               'ierr',   ('a6',   str,   'ordered list of indices corresponding to uncertainty estimates of transition parameters'    ,''   )),(
               'iref',   ('a12',  str,   'ordered list of reference identifiers for transition parameters'                    ,''           )),(
               'lmix',   ('a1',   str,   'flag indicating the presence of additional data and code relating to line-mixing'   ,''           )),(
               'gp',     ('a7',   float, 'upper state degeneracy'                         ,''                      )),(
               'gpp',    ('a7',   float, 'lower state degeneracy'                         ,''                      ))
               ])


# quick fix # TODO: proper implementation of HITRAN classes, and groups
# Update: classes are now implemented properly, groups remain to be done 


def _generate_cache_file(fname, df):
    ''' Generate cache file for pandas Dataframe df under name fname, which 
    must be a HDF5 
    
    Parameters
    ----------
    
    fname: str
        HDF5 filename
        
    df: pandas Dataframe
        line database to store
        
    Note
    ----
    
    Performance of different compression librairies / methods:
        
    - storing as np.float32 can save up to 50% space, but we can loose some 
    information such as low linestrength numbers (1e-80, 1e-100) that becomes 0
    
    - storing with compression (complevel=9, complib='blosc') saves up to 50%
    space, but reading speed is doubled 
    
    - storing with low compression 'complevel=1', complib='blosc' can still
    save up to 50%, but reading speed is only increased by about 20%
    
    - format='fixed' is always much faster than format='table' (about a factor 5!)
    
    Eventually we use ``format='fixed', mode='w', complevel=1, complib='blosc')`` 

    
    '''
    
    assert fname.endswith('h5')
    
    df.to_hdf(fname, 'df', format='fixed', mode='w', complevel=1, complib='blosc')
#    df.to_hdf(fname, 'df', format='fixed', mode='w')  # would be 10-20% faster, but take 2x more space
    


# %% Hitran global quanta classes

def _parse_HITRAN_class1(df):
    ''' Diatomic molecules: CO, HF, HCl, HBr, HI, N2, NO+
    
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax:
    
    >>>      v1
    >>>  13x I2
    
    '''
    dgu = df['globu'].str.extract(
            '[ ]{13}(?P<v1u>[\d ]{2})',
            expand=True)
    dgl = df['globl'].str.extract(
            '[ ]{13}(?P<v1l>[\d ]{2})',
            expand=True)
    dgu = dgu.apply(pd.to_numeric)
    dgl = dgl.apply(pd.to_numeric)
    return pd.concat([df, dgu, dgl], axis=1)

def _parse_HITRAN_class2(df):
    ''' Diatomic molecules with different electronic levels: O2
    
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax:
    
    >>>      X  v1
    >>>  12x A1 I2
    
    '''
    print('parse_global_quanta not implemented for molecules of this HITRAN class')
    return df

def _parse_HITRAN_class3(df):
    ''' Diatomic molecules with doublet-Pi electronic state: NO, OH, ClO
    
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -------
    
    HITRAN syntax:
    
    >>>      X i     v1
    >>>  7x A1 A3 2x I2
    '''
    print('parse_global_quanta not implemented for molecules of this HITRAN class')
    return df

def _parse_HITRAN_class4(df):
    ''' Parse linear triatomic class in HITRAN [1]_: N2O, OCS, HCN
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    Notes
    -----
    
    HITRAN syntax:
    
    >>>     v1 v2 l2 v3
    >>>  7x I2 I2 I2 I2
    
    Note: I2 in regexp: [\d ]{2} 
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    '''
    dgu = df['globu'].str.extract(
            '[ ]{7}(?P<v1u>[\d ]{2})(?P<v2u>[\d ]{2})(?P<l2u>[\d ]{2})(?P<v3u>[\d ]{2})',
            expand=True)
    dgl = df['globl'].str.extract(
            '[ ]{7}(?P<v1l>[\d ]{2})(?P<v2l>[\d ]{2})(?P<l2l>[\d ]{2})(?P<v3l>[\d ]{2})',
            expand=True)
    dgu = dgu.apply(pd.to_numeric)
    dgl = dgl.apply(pd.to_numeric)
    return pd.concat([df, dgu, dgl], axis=1)

def _parse_HITRAN_class5(df):
    ''' Parse linear triatomic with large Fermi resonance in HITRAN [1]_: CO2
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    Notes
    -----
    
    HITRAN syntax:
    
    >>>     v1 v2 l2 v3 r
    >>>  6x I2 I2 I2 I2 I1
    
    Note: I2 in regexp: [\d ]{2} 
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    '''
    dgu = df['globu'].str.extract(
            '[ ]{6}(?P<v1u>[\d ]{2})(?P<v2u>[\d ]{2})(?P<l2u>[\d ]{2})(?P<v3u>[\d ]{2})(?P<ru>\d)',
            expand=True)
    dgl = df['globl'].str.extract(
            '[ ]{6}(?P<v1l>[\d ]{2})(?P<v2l>[\d ]{2})(?P<l2l>[\d ]{2})(?P<v3l>[\d ]{2})(?P<rl>\d)',
            expand=True)
    dgu = dgu.apply(pd.to_numeric)
    dgl = dgl.apply(pd.to_numeric)
    return pd.concat([df, dgu, dgl], axis=1)

def _parse_HITRAN_class6(df):
    ''' Parse non-linear triatomic in HITRAN [1]_: H2O, O3, SO2, NO2, HOCl, H2S, HO2, HOBr
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    Notes
    -----
    
    HITRAN syntax:
    
    >>>     v1 v2 v3
    >>>  9x I2 I2 I2
    
    Note: I2 in regexp: [\d ]{2} 
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    '''
    dgu = df['globu'].str.extract(
            '[ ]{9}(?P<v1u>[\d ]{2})(?P<v2u>[\d ]{2})(?P<v3u>[\d ]{2})',
            expand=True)
    dgl = df['globl'].str.extract(
            '[ ]{9}(?P<v1l>[\d ]{2})(?P<v2l>[\d ]{2})(?P<v3l>[\d ]{2})',
            expand=True)
    dgu = dgu.apply(pd.to_numeric)
    dgl = dgl.apply(pd.to_numeric)
    return pd.concat([df, dgu, dgl], axis=1)

def _parse_HITRAN_class7(df):
    ''' Parse linear tetratomic in HITRAN [1]_: C2H2
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    Notes
    -----
    
    HITRAN syntax:
    
    >>>
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    '''
    print('parse_global_quanta not implemented for molecules of this HITRAN class')
    return df

def _parse_HITRAN_class8(df):
    ''' Pyramidal tetratomic in HITRAN [1]_: NH3, PH3
    
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax:
    
    >>>
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    '''
    print('parse_global_quanta not implemented for molecules of this HITRAN class')
    return df

def _parse_HITRAN_class9(df):
    ''' Non-linear tetratomic in HITRAN [1]_: H2CO, H2O2, COF2
    
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax:
    
    >>>
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    '''
    print('parse_global_quanta not implemented for molecules of this HITRAN class')
    return df

def _parse_HITRAN_class10(df):
    ''' Pentatomic or greater polyatomic in HITRAN [1]_
    
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax:
    
    >>>
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    '''
    print('parse_global_quanta not implemented for molecules of this HITRAN class')
    return df
    
# %% HITRAN Local quanta
    
    
def _parse_HITRAN_group1(df):
    ''' 
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax: [1]_
    
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    
    '''
    print('parse_local_quanta not implemented for molecules of this HITRAN group')
    return df

    
    
def _parse_HITRAN_group2(df):
    ''' 
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax: [1]
    
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    
    '''
    
    dgu = df['locu'].str.extract(
            '[ ]{10}(?P<Fu>.{5})',
            expand=True)
    dgl = df['locl'].str.extract(
            '[ ]{5}(?P<branch>[\S]{1})(?P<jl>[\d ]{3})(?P<sym>.)(?P<Fl>.{5})',
            expand=True)
    
    #dgl['jl'] = dgl.jl.apply(pd.to_numeric)
    dgl['jl'] = pd.to_numeric(dgl.jl)
    
    del df['locu']
    del df['locl']
    
    return pd.concat([df, dgu, dgl], axis=1)

    
#                # 10X included in Fu
#               'Fu',     ('a15',  str,  'upper state total angular momentum including nuclear spin'  ,''         )),(
#               #'locl',   ('a15',  str,   'electronic and vibrational local lower quanta'  ,''                      )),(
#               'branch', ('a6',   str,     'O, P, Q, R, S branch symbol'                  ,''                      )),(
#               'jl',     ('a3',   int,    'lower state rotational quantum number'         ,''                      )),(
#               'sym',    ('a1',   str,     'symmetry'                                     ,''                      )),(
#               'Fl',     ('a5',   str,   'lower state total angular momentum including nuclear spin', ''         )),(

    
def _parse_HITRAN_group3(df):
    ''' 
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax [1]_:
        
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    
    '''
    print('parse_local_quanta not implemented for molecules of this HITRAN group')
    return df

    
def _parse_HITRAN_group4(df):
    ''' 
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax [1]_:
        
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    
    '''
    print('parse_local_quanta not implemented for molecules of this HITRAN group')
    return df

    
def _parse_HITRAN_group5(df):
    ''' 
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax [1]_:
        
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    
    '''
    print('parse_local_quanta not implemented for molecules of this HITRAN group')
    return df

    
def _parse_HITRAN_group6(df):
    ''' 
    
    Parameters
    ----------
    
    df: pandas Dataframe
        lines read from a HITRAN-like database
        
    
    Notes
    -----
    
    HITRAN syntax [1]_:
    
        
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    
    '''
    print('parse_local_quanta not implemented for molecules of this HITRAN group')
    return df


# %% Reading function


def _format_dtype(dtype):
    ''' Format dtype from specific columns. Crash with hopefully helping error message '''
    
    try:
        dt = np.dtype([(str(k), c) for k, c in dtype])
        # Note: dtype names cannot be `unicode` in Python2. Hence the str()
    except TypeError:
        # Cant read database. Try to be more explicit for user
        print('Data type')
        print('-'*30)
        for (k, c) in dtype:
            print(str(k), '\t', c)
        print('-'*30)
        raise
    return dt

def _cast_to_dtype(data, dtype):
    ''' Cast array to certain type, crash with hopefull helping error message.
    Return casted data
    
    
    Parameters    
    ----------
    
    data: array to cast
    
    dtype: (ordered) list of (param, type)
    
    '''
    
    dt = _format_dtype(dtype)
    
    try:
        data = np.array(data, dtype=dt)
    except ValueError:
        try:
            # Cant read database. Try to be more explicit for user
            print('Cant cast data to specific dtype. Trying column by column:')
            print('-'*30)
            for i in range(len(data[0])):
                print(dtype[i], '\t', np.array(data[0][i], dtype=dt[i]))
            print('-'*30)
        except ValueError:
            print('>>> Next param:', dtype[i], '. Value:', data[0][i], '\n')
            raise ValueError('Cant cast data to specific dtype. Tried column by column. See results above')

    return data

def hit2df(fname, count=-1, cache=False, verbose=True):
    ''' Convert a HITRAN/HITEMP [1]_ file to a Pandas dataframe 
    
    
    Parameters    
    ----------
    
    fname: str
        HITRAN-HITEMP file name 
        
    count: int
        number of items to read (-1 means all file)
        
    cache: boolean
        if True, a pandas-readable HDF5 file is generated on first access, 
        and later used. This saves on the datatype cast and conversion and
        improves performances a lot (but changes in the database are not 
        taken into account). If False, no database is used. If 'regen', temp
        file are reconstructed. Default False. 
    
    
    Returns
    -------
    
    df: pandas Dataframe
        dataframe containing all lines and parameters
        
    
    
    References
    ----------

    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    
    
    Notes
    -----
    
    Performances: see CDSD-HITEMP parser
    
    '''
    
    columns = columns_2004

    if cache: # lookup if cached file exist. 
#        fcache = fname+'.cached'
        fcache = splitext(fname)[0]+'.h5'
        if exists(fcache):
            if cache == 'regen':
                os.remove(fcache)
                if verbose: print('Deleted h5 cache file : {0}'.format(fcache))
            else:
                if verbose: print('Using h5 file: {0}'.format(fcache))
    #            return pd.read_csv(fcache)
                return pd.read_hdf(fcache, 'df')
        
    # Detect the molecule by reading the start of the file
    with open(fname) as f:
        mol = get_molecule(int(f.read(2)))
        
    # %% Start reading the full file
    
    # To be faster, we read file totally in bytes mode with fromfiles. But that
    # requires to properly decode the line return character:
    
    # problem arise when file was written in an OS and read in another OS (for instance,
    # line return characters are not converted when read from .egg files). Here
    # we read the first line and infer the line return character for it 
    
    # ... Create a dtype with the binary data format and the desired column names
    dtype = [(k, c[0]) for (k, c) in columns.items()]+[('_linereturn','a2')]   
    # ... _linereturn is to capture the line return symbol. We delete it afterwards
    dt = _format_dtype(dtype)
    data = np.fromfile(fname, dtype=dt, count=1)   # just read the first line 
    
    # get format of line return
    from radis.misc.basics import to_str
    linereturn = to_str(data[0][-1])
    if to_str('\n\r') in linereturn:
        linereturnformat = 'a2'
    elif to_str('\n') in linereturn or to_str('\r') in linereturn:
        linereturnformat = 'a1'
    else:
        raise ValueError('Line return format unknown: {0}. Please update RADIS'.format(linereturn))
        
    # Now re-read with correct line return character

    # ... Create a dtype with the binary data format and the desired column names
    dtype = [(k, c[0]) for (k, c) in columns.items()]+[('_linereturn',linereturnformat)]   
    # ... _linereturn is to capture the line return symbol. We delete it afterwards
    dt = _format_dtype(dtype)
    data = np.fromfile(fname, dtype=dt, count=count)
    
    # ... Cast to new type
    # This requires to recast all the data already read, but is still the fastest
    # method I found to read a file directly (for performance benchmark see 
    # CDSD-HITEMP parser)
    newtype = [c[0] if (c[1]==str) else c[1] for c in columns.values()]
    dtype = list(zip(list(columns.keys()), newtype))+[('_linereturn',linereturnformat)]
    data = _cast_to_dtype(data, dtype)
    
    # %% Create dataframe    
    df = pd.DataFrame(data.tolist(), columns=list(columns.keys())+['_linereturn'])

    # assert one molecule per database only. Else the groupbase data reading 
    # above doesnt make sense
    nmol = len(set(df['id']))
    if nmol == 0:
        raise ValueError('Databank looks empty')        
    elif nmol!=1:
        # Crash, give explicity error messages
        try:
            secondline = df.iloc[1]
        except IndexError:
            secondline = ''
        raise ValueError('Multiple molecules in database ({0}). Current '.format(nmol)+\
                         'spectral code only computes 1 species at the time. Use MergeSlabs. '+\
                         'Verify the parsing was correct by looking at the first row below: '+\
                         '\n{0}'.format(df.iloc[0])+'\n----------------\nand the second row '+\
                         'below: \n{0}'.format(secondline))
    
    for k, c in columns.items():
        if c[1] == str:
            df[k] = df[k].str.decode("utf-8")
    
    # %% Add local quanta attributes, based on the HITRAN group
    df = parse_local_quanta(df, mol)
    
    # %% Add global quanta attributes, based on the HITRAN class
    df = parse_global_quanta(df, mol)

    # Strip whitespaces around PQR columns (due to 2 columns jumped)
    if 'branch' in df:
        df['branch'] = df.branch.str.strip()
    
    # Delete dummy column than handled the line return character
    del df['_linereturn']
    
    if cache: # cached file mode but cached file doesn't exist yet (else we had returned)
        if verbose: print('Generating cached file: {0}'.format(fcache))
        try:
#            df.to_csv(fcache)
            _generate_cache_file(fcache, df)
        except:
            if verbose:
                print(sys.exc_info())
                print('An error occured in cache file generation. Lookup access rights')
            pass
        
    return df 

def parse_local_quanta(df, mol):
    '''
    Parameters
    ----------
    
    df: pandas Dataframe
    
    mol: str
        molecule name
    '''
    
    # had some problems with bytes types
    df['locu'] = df.locu.astype(str)
    df['locl'] = df.locl.astype(str)
    
    if mol in HITRAN_GROUP1:
        df = _parse_HITRAN_group1(df)
    elif mol in HITRAN_GROUP2:
        df = _parse_HITRAN_group2(df)
    elif mol in HITRAN_GROUP3:
        df = _parse_HITRAN_group3(df)
    elif mol in HITRAN_GROUP4:
        df = _parse_HITRAN_group4(df)
    elif mol in HITRAN_GROUP5:
        df = _parse_HITRAN_group5(df)
    elif mol in HITRAN_GROUP6:
        df = _parse_HITRAN_group6(df)
    else:
        raise ValueError('Unknown group for molecule {0}. Cant parse local quanta'.format(
                mol))
        
    return df
            
def parse_global_quanta(df, mol):
    ''' 
    
    Parameters
    ----------
    
    df: pandas Dataframe
    
    mol: str
        molecule name
    '''
    
    # had some problems with bytes types
    df['globu'] = df.globu.astype(str)
    df['globl'] = df.globl.astype(str)
    
    if mol in HITRAN_CLASS1:
        df = _parse_HITRAN_class1(df)
    elif mol in HITRAN_CLASS2:
        df = _parse_HITRAN_class2(df)
    elif mol in HITRAN_CLASS3:
        df = _parse_HITRAN_class3(df)
    elif mol in HITRAN_CLASS4:
        df = _parse_HITRAN_class4(df)
    elif mol in HITRAN_CLASS5:
        df = _parse_HITRAN_class5(df)
    elif mol in HITRAN_CLASS6:
        df = _parse_HITRAN_class6(df)
    elif mol in HITRAN_CLASS7:
        df = _parse_HITRAN_class7(df)
    elif mol in HITRAN_CLASS8:
        df = _parse_HITRAN_class8(df)
    elif mol in HITRAN_CLASS9:
        df = _parse_HITRAN_class9(df)
    elif mol in HITRAN_CLASS10:
        df = _parse_HITRAN_class10(df)
    else:
        raise ValueError('Unknown class for molecule {0}. Cant parse global quanta'.format(
                mol))
        
    return df
            

def get_molecule_identifier(molecule_name):
    '''
    For a given input molecular formula, return the corresponding HITRAN molecule 
    identifier number [1]_.
    
    
    Parameters
    ----------
    molecular_formula : str
        The string describing the molecule.
        
        
    Returns
    -------
    M : int
        The HITRAN molecular identified number.
        
        
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    Function is from from https://github.com/nzhagen/hitran/blob/master/hitran.py
    
    '''

    trans = { '1':'H2O',    '2':'CO2',   '3':'O3',      '4':'N2O',   
              '5':'CO',    '6':'CH4',   '7':'O2',     '8':'NO',
              '9':'SO2',   '10':'NO2',  '11':'NH3',    '12':'HNO3', 
              '13':'OH',   '14':'HF',   '15':'HCl',   '16':'HBr',
             '17':'HI',    '18':'ClO',  '19':'OCS',    '20':'H2CO', 
             '21':'HOCl', '22':'N2',   '23':'HCN',   '24':'CH3Cl',
             '25':'H2O2',  '26':'C2H2', '27':'C2H6',   '28':'PH3',  
             '29':'COF2', '30':'SF6',  '31':'H2S',   '32':'HCOOH',
             '33':'HO2',   '34':'O',    '35':'ClONO2', '36':'NO+',  
             '37':'HOBr', '38':'C2H4', '39':'CH3OH', '40':'CH3Br',
             '41':'CH3CN', '42':'CF4',  '43':'C4H2',   '44':'HC3N', 
             '45':'H2',   '46':'CS',   '47':'SO3'}
    ## Invert the dictionary.
    trans = {v:k for k,v in trans.items()}

    return int(trans[molecule_name])

def get_molecule(molecule_id):
    '''
    For a given input molecular identifier, return the corresponding HITRAN 
    molecule name [1]_.
    
    
    Parameters    
    ----------
    
    molecular_id : str
        Hitran identifier of the molecule.
        
    
    References
    ----------
    
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__
    
    '''
    
    # assert str
    id = '{0}'.format(molecule_id)

    trans = { '1':'H2O',    '2':'CO2',   '3':'O3',      '4':'N2O',   
              '5':'CO',    '6':'CH4',   '7':'O2',     '8':'NO',
              '9':'SO2',   '10':'NO2',  '11':'NH3',    '12':'HNO3', 
              '13':'OH',   '14':'HF',   '15':'HCl',   '16':'HBr',
             '17':'HI',    '18':'ClO',  '19':'OCS',    '20':'H2CO', 
             '21':'HOCl', '22':'N2',   '23':'HCN',   '24':'CH3Cl',
             '25':'H2O2',  '26':'C2H2', '27':'C2H6',   '28':'PH3',  
             '29':'COF2', '30':'SF6',  '31':'H2S',   '32':'HCOOH',
             '33':'HO2',   '34':'O',    '35':'ClONO2', '36':'NO+',  
             '37':'HOBr', '38':'C2H4', '39':'CH3OH', '40':'CH3Br',
             '41':'CH3CN', '42':'CF4',  '43':'C4H2',   '44':'HC3N', 
             '45':'H2',   '46':'CS',   '47':'SO3'}

    return trans[id]

## ======================================================
# %% Test

if __name__ == '__main__':
    from radis.test.test_io import test_hitran_parser
    print('Testing HITRAN parsing: ', test_hitran_parser())
        