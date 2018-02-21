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
from radis.test.utils import getTestFile
from time import time


def test_hitran_parser__fast(verbose=True, warnings=True, **kwargs):
    ''' Analyse some default files to make sure everything still works'''
    
    b = True
    
    t0 = time()
    df = hit2df(getTestFile('hitran_CO_fragment.par'))
    print('File loaded in {0:.0f}s'.format(time()-t0))
    if verbose: print(df.head())
    b *= (list(df.ix[0, ['v1u', 'v1l']]) == [4, 4])

    t0 = time()
    df = hit2df(getTestFile('hitran_CO2_fragment.par'))
    print('File loaded in {0:.0f}s'.format(time()-t0))
    if verbose: print(df.head())
    b *= (list(df.ix[0, ['v1u', 'v2u', 'l2u', 'v3u', 'v1l', 'v2l', 'l2l', 'v3l']]) ==
              [4, 0, 0, 0, 0, 0, 0, 1])
    
    return (bool(b))

def _run_testcases(verbose=True, *args, **kwargs):

    b1 = test_hitran_parser__fast(verbose=verbose,*args, **kwargs)
    
    if verbose:
        print('>>> test_hitran_parser__fast:', b1)
        
    b = bool(b1)
    
    return bool(b)

if __name__ == '__main__':
    print('Testing io.py: ', _run_testcases(verbose=True))
    
