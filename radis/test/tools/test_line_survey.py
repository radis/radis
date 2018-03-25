# -*- coding: utf-8 -*-
"""
Test that line survey works
"""

from radis.misc.utils import DatabankNotFound
from radis.test.utils import getTestFile
from radis.tools.database import load_spec
import os
from os.path import exists

def test_line_survey__fast(verbose=True, plot=False, warnings=True, *args, **kwargs):
    ''' Test line survey '''

    _temp_file = 'radis_test_line_survey.html'
    if exists(_temp_file):
        os.remove(_temp_file)
    
    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'), binary=True)
    s.line_survey(overlay='abscoeff', display=plot, filename=_temp_file)
    
    assert exists(_temp_file)
    
    if verbose:
        print('test_line_survey: html file was correctly generated')
    
    if not plot:
        # clean after use
        os.remove(_temp_file)
        
    return True

def _run_testcases(plot=True, verbose=True, *args, **kwargs):

    # Show media line_shift
    test_line_survey__fast(plot=plot, verbose=verbose, *args, **kwargs)

    return True

if __name__ == '__main__':

    print('Testing line survey functions:', _run_testcases(plot=True))