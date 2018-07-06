# -*- coding: utf-8 -*-
""" Test parsers

Notes
-----

Runs tests for radis/io so that they can be accessed by pytest (and hopefully 
the CI test suite)

Examples
--------

Run all tests::

	pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)::

	pytest -m fast

-------------------------------------------------------------------------------

"""

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.io.hitran import hit2df
from radis.io.cdsd import cdsd2df
from radis.test.utils import getTestFile
import pytest
import numpy as np


@pytest.mark.fast
def test_hitran_co(verbose=True, warnings=True, **kwargs):
    ''' Analyse some default files to make sure everything still works'''

    # 1. Load
    df = hit2df(getTestFile('hitran_CO_fragment.par'), cache='regen')
    if verbose:
        print('Read hitran_CO_fragment.par')
        print('---------------------------')
        print(df.head())
        
    # 2. Test
    assert (list(df.loc[0, ['vu', 'vl']]) == [4, 4])
    assert df.dtypes['vu'] == np.int64
    assert df.dtypes['vl'] == np.int64
    
    return True

def test_hitran_co2(verbose=True, warnings=True, **kwargs):
    
    # 1. Load
    df = hit2df(getTestFile('hitran_CO2_fragment.par'), cache='regen')
    if verbose:
        print('Read hitran_CO2_fragment.par')
        print('----------------------------')
        print(df.head())
        
    # 2. Test
    assert (list(df.loc[0, ['v1u', 'v2u', 'l2u', 'v3u', 'v1l', 'v2l', 'l2l', 'v3l']]) ==
            [4, 0, 0, 0, 0, 0, 0, 1])
    assert df.dtypes['v1l'] == np.int64
    assert df.dtypes['v3u'] == np.int64

    return True

def test_hitran_h2o(verbose=True, warnings=True, **kwargs):
    
    # 1. Load
    df = hit2df(getTestFile('hitran_2016_H2O_2iso_2000_2100cm.par'), cache='regen')
    if verbose:
        print('Read hitran_2016_H2O_2iso_2000_2100cm.par')
        print('-----------------------------------------')
        print(df.head())
        
    # 2. Test
    assert (list(df.loc[0, ['v1u', 'v2u', 'v3u', 'v1l', 'v2l', 'v3l']]) ==
            [0, 2, 0, 0, 1, 0])
    assert df.dtypes['v1l'] == np.int64
    assert df.dtypes['v3u'] == np.int64
    
    return True


def test_hitemp(verbose=True, warnings=True, **kwargs):
    ''' Analyse some default files to make sure everything still works'''

    # 1. Load
    df = cdsd2df(getTestFile('cdsd_hitemp_09.txt'), cache='regen',
                 drop_non_numeric=True)
    if verbose:
        print(df.head())
    
    # 2. Tests
    assert df.wav[3] == 2250.00096
    # make sure P Q R is correctly replaced by drop_non_numeric: 
    assert 'branch' in df    
    assert df['branch'].iloc[0] in ['Q', 0]
    assert df['branch'].iloc[1] in ['R', 1]
    
    return True

def _run_testcases(verbose=True, *args, **kwargs):

    test_hitran_co(verbose=verbose, *args, **kwargs)
    test_hitran_co2(verbose=verbose, *args, **kwargs)
    test_hitran_h2o(verbose=verbose, *args, **kwargs)
    test_hitemp(verbose=verbose, *args, **kwargs)

    return True


if __name__ == '__main__':
    print('Testing io.py: ', _run_testcases(verbose=True))
