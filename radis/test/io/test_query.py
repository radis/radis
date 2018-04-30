# -*- coding: utf-8 -*-
"""
Test query functions
"""

import pytest
from radis.io.query import fetch_astroquery

@pytest.mark.needs_connection   # ignored by pytest with argument -m "not needs_connection" 
def test_fetch_astroquery(verbose=True, *args, **kwargs):
    ''' Test astroquery '''
    fetch_astroquery('CO2', 1, 2200, 2400, verbose=verbose)
    
    if verbose:
        print('test_fetch_astroquery: no test defined')
    
    return

@pytest.mark.needs_connection   # ignored by pytest with argument -m "not needs_connection" 
def test_fetch_astroquery_empty(verbose=True, *args, **kwargs):
    ''' Test astroquery: get a spectral range where there are no lines'''
    df = fetch_astroquery(2, 1, 25000, 50000, verbose=verbose)  # 200-400 nm
    
    assert len(df) == 0
    
    return

def _run_testcases(verbose=True, *args, **kwargs):

    test_fetch_astroquery(verbose=verbose,*args, **kwargs)
    test_fetch_astroquery_empty(verbose=verbose,*args, **kwargs)
    
    return True

if __name__ == '__main__':
    print('test_query.py: ', _run_testcases(verbose=True))
    