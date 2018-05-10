# -*- coding: utf-8 -*-
""" Test parsers

Notes
-----

Created on Tue Aug 15 14:56:25 2017

@author: erwan

Runs tests for radis/io so that they can be accessed by pytest (and hopefully 
the CI test suite)

Examples
--------

Run all tests::

	pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have 'fast' in their name)::

	pytest -k fast

"""

from __future__ import print_function, absolute_import, division, unicode_literals

from radis.io.hitran import hit2df
from neq.io.hitran.hiparser import hit2dfTAB
from radis.test.utils import getTestFile
from time import time
import pytest

@pytest.mark.fast
def test_hitran_parser(verbose=True, warnings=True, **kwargs):
    ''' Analyse some default files to make sure everything still works'''
    
    t0 = time()
    df = hit2df(getTestFile('hitran_CO_fragment.par'))
    if verbose: print('File loaded in {0:.0f}s'.format(time()-t0))
    if verbose: print(df.head())
    assert (list(df.loc[0, ['v1u', 'v1l']]) == [4, 4])

    t0 = time()
    df = hit2df(getTestFile('hitran_CO2_fragment.par'))
    if verbose: print('File loaded in {0:.0f}s'.format(time()-t0))
    if verbose: print(df.head())
    assert (list(df.loc[0, ['v1u', 'v2u', 'l2u', 'v3u', 'v1l', 'v2l', 'l2l', 'v3l']]) ==
              [4, 0, 0, 0, 0, 0, 0, 1])
    
    return True

def _run_testcases(verbose=True, *args, **kwargs):

    test_hitran_parser(verbose=verbose,*args, **kwargs)
    
    return True

def test_hitran_tab_parser(verbose=True):
    '''Opens a test file and tries to read some key parameters in it'''
    
    t0 = time()
    df = hit2dfTAB(getTestFile('CO_in_air_tabtest.out'))
    if verbose: print('File loaded in {0:.0f}s'.format(time()-t0))
    assert (df.statep[0] == 'ElecStateLabel=X;v=2;J=26')
    #TODO: change it when you wrote the extraction of J and v
    
    return True

if __name__ == '__main__':
    print('Testing io.py: ', _run_testcases(verbose=True))
    
