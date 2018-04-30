# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 12:59:37 2014

@author: Erwan

Small functions used in other procedures
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import os
from os.path import join, abspath, normcase, normpath
import numpy as np
import sys
import fileinput
from six.moves import filter
from six.moves import range
from six import binary_type
if sys.version_info[0] == 2:
    from itertools import ifilterfalse as filterfalse
else:
    from itertools import filterfalse
from itertools import tee
import pandas as pd

verbose = True

# %%
#==============================================================================
# Basic Functions
#==============================================================================

def all_in(keys, L):
    ''' Returns whether all items in keys are in list L '''
    return all([k in L for k in keys])
        
def any_in(keys, L):
    ''' Returns whether any of the items in keys are in list L '''
    return any([k in L for k in keys])
        
def key_max_val(d):
    '''Return the dictionary key with max value'''
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]


def exec_file(afile, globalz=None, localz=None):
    with open(afile, "r") as fh:
        exec(fh.read(), globalz, localz)

def remove_duplicates(l):
    ''' Remove duplicates from a list, without changing the order. Note that 
    if the order doesn't matter you could just do set(l) '''

    l1 = []
    for e in l:
        if not e in l1:
            l1.append(e)
    return l1

def partition(pred, iterable):
    ''' Use a predicate to partition entries into false entries and true entries
    

    Returns
    -------
    
    Returns two lists: positive, and negative
    
    
    Example
    -------
    
     >>> partition(is_odd, range(10)) 
     --> [0 2 4 6 8], [1 3 5 7 9]
    
    '''
    t1, t2 = tee(iterable)
    return list(filter(pred, t2)), list(filterfalse(pred, t1))

# Just for performance comparison. Below: naive implementation of partition
# %timeit partition(lambda x:x>0.5, np.random.rand(10))
# ... partition         2.67 µs per loop
# ... partition_naive   28.5 µs per loop
#def partition_naive(pred, iterable):
#    trues = []
#    falses = []
#    for item in iterable:
#        if pred(item):
#            trues.append(item)
#        else:
#            falses.append(item)
#    return trues, falses
    
# %% Compare / merge tools

def compare_dict(d1, d2, verbose=True, compare_as_paths=[]):
    ''' Returns ratio of equal keys [0-1]
    If verbose, also print all keys and values on 2 columns
    
    
    Parameters    
    ----------
    
    d1, d2: dict
    
    
    Other Parameters
    ----------------
    
    compare_as_paths: list of keys
        compare the values corresponding to given keys as path (irrespective of
        forward / backward slashes, or case )
    
    verbose: boolean, or 'if_different'
        'if_different' means results will be shown only if there is a difference. 
        function is called twice
        

    Returns
    -------
    
    out: float [0-1]
        ratio of matching keys 
        
    '''
    
    if verbose == 'if_different':
        restart = True
        verbose = False
    else:
        restart = False
    
    if verbose: print('{0:15}{1}\t{2}'.format('Key', 'Left', 'Right'))
    if verbose: print('-'*40)
    all_keys = set(list(d1.keys())+list(d2.keys()))
    s = 0       # counter of all matching keys
    for k in all_keys:
        if k in d1 and k in d2:   # key in both dicts. Let's compare values
            # Deal with path case
            if k in compare_as_paths:
                if not compare_paths(d1[k], d2[k]):
                    print('{0:15}{1}\t{2}'.format(k, d1[k], d2[k]))
                else:
                    s += 1 
            # Other cases
            else:
                if d1[k] != d2[k]:
                    if verbose: print('{0:15}{1}\t{2}'.format(k, d1[k], d2[k]))
                else:
                    s += 1 
        elif k in d1: 
            if verbose: print('{0:15}{1}\tN/A'.format(k, d1[k]))
        else:
            if verbose: print('{0:15}N/A\t{1}'.format(k, d2[k]))
    if verbose: print('-'*40)
    
    if len(all_keys) == 0:
        out = 1 
    else:
        out = s/len(all_keys)

    # Exit
    if restart and out != 1: 
        return compare_dict(d1, d2, verbose=True, compare_as_paths=compare_as_paths)
    else:
        return out

def compare_lists(l1, l2, verbose=True):
    ''' Compare 2 lists of elements that may not be of the same length, irrespective
    of order. Returns the ratio of elements [0-1] present in both lists. If verbose, 
    prints the differences 
    
    
    Parameters    
    ----------
    
    l1, l2: list-like 
    
    verbose: boolean, or 'if_different'
        'if_different' means results will be shown only if there is a difference. 
        function is called twice
        

    Returns
    -------
    
    out: float [0-1]
        ratio of matching keys 
        
    '''
    
    if verbose == 'if_different':
        restart = True
        verbose = False
    else:
        restart = False
    
    if verbose: print('{0:20}\t{1}'.format('Left', 'Right'))
    if verbose: print('-'*40)
    all_keys = set(list(l1)+list(l2))
    s = 0       # counter of all matching keys
    for k in all_keys:
        if k in l1 and k in l2:   # key in both lists
            s += 1 
        elif k in l1: 
            if verbose: print('{0:20}\tN/A'.format('{0} ({1})'.format(k, type(k))))
        else:
            if verbose: print('{0:20}\t{1} ({2})'.format('N/A', k, type(k)))
    if verbose: print('-'*40)
    
    if len(all_keys) == 0:
        out = 1 
    else:
        out = s/len(all_keys)
        
    # Exit
    if restart and out != 1: 
        return compare_lists(l1, l2, verbose=True)
    else:
        return out

def stdpath(p):
    ''' Convert path p in standard path (irrespective of slash / backslash,
    or case)
    '''

    return normpath(normcase(abspath(p)))

def compare_paths(p1, p2):
    ''' Compare 2 paths p1 and p2 '''
    return stdpath(p1) == stdpath(p2)
    
def merge_lists(lists):
    ''' Merge a list of lists and return a list with unique elements '''
    return list(set(sum([l for l in lists], [])))

#==============================================================================
# %% Pandas specific
#==============================================================================
    
def merge_rename_columns(df, columns1, columns2, merged_names):
    ''' Merge all columns under easier names. Only keep the useful ones
    Returns a new dataframe
    
    
    Parameters    
    ----------
    
    df: pandas Dataframe
    
    columns1: list
        list of columns names
        
    columns2: list
        list of columns names, whose index match columns 1
        
    merged_names: list
        new names
        
       
    Example
    -------
    
    df = merge_rename_columns(df1, ['lvl_u', 'ju', 'Eu', 'nu', 'gu', 'grotu'],
                                   ['lvl_l', 'jl', 'El', 'nl', 'gl', 'grotl'],
                                   ['lvl',   'j',  'E',  'n',  'g',  'grot']
                                   )
    '''
    
    assert all_in(columns1, list(df.keys()))
    assert all_in(columns2, list(df.keys()))
 
    df1 = df.loc[:,columns1]
    df2 = df.loc[:,columns2]
    df1.rename(columns={columns1[i]: merged_names[i] for i in range(len(merged_names))}, inplace=True)
    df2.rename(columns={columns2[i]: merged_names[i] for i in range(len(merged_names))}, inplace=True)
    df = pd.concat((df1, df2), ignore_index=True)

    return df.drop_duplicates()

def print_series(a):
    ''' Print a pandas series `a` , explicitely showing all rows'''
    
    for i, k in enumerate(a.keys()):
        print(k, '\t', a.values[0][i])

# %%
#==============================================================================
# Types 
#==============================================================================

def list_if_float(a):
    if type(a) is list:
        return a
    else:
        return [a]

def is_list(a):
    ''' Returns True if a has list-like type: list, np.array, tuple, set, etc.)'''
    return type(a) in [list, np.ndarray, tuple, set]

def is_float(a):
    ''' Returns True if a has float-like type: float, np.float64, np.int64, etc.)'''
    return type(a) in [float, np.float64, np.int32, np.float32, int, np.int64]

def to_str(a):
    if isinstance(a, binary_type):
        return a.decode('utf-8')
    else:
        return a
    
    