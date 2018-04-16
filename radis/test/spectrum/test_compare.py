# -*- coding: utf-8 -*-
"""

Examples
--------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have 'fast' in their name)

>>> pytest -k fast

"""

from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from radis.test.utils import getTestFile
from radis.tools.database import load_spec
from radis.spectrum.compare import get_distance, plot_diff

# Test routines

def test_compare_methods(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    ''' Just run all Spectrum compare methods to check they work'''

    if plot and close_plots:
        import matplotlib.pyplot as plt
        plt.close('all')
    
    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'))
    s.resample(np.linspace(2193, 2193.8, 100))   # limits to a single line, because get_distance() 
                                                   # is very computationaly heavy 
    s.update('radiance_noslit')
    s_noabsorption = s.copy()
    s.name = 'solve RTE'
    s_noabsorption.name = 'optically thin'
    
    # rescale, normal 
#    s.rescale_mole_fraction(10)
    s.rescale_path_length(10)
    
    # rescale, optically thin mode
    s_noabsorption.conditions['self_absorption'] = False
#    s_noabsorption.rescale_mole_fraction(10)
    s_noabsorption.rescale_path_length(10)
    
    # Compare 
    get_distance(s, s_noabsorption, 'radiance_noslit')  # should be added in an example with experimental fit of bandhead
    title = 'CO x={0:.2f}, L={1:.2f}cm'.format(s.conditions['mole_fraction'], s.conditions['path_length'])
    if plot:
        plot_diff(s, s_noabsorption, method='diff', title=title)
        plot_diff(s, s_noabsorption, method='ratio', normalize=True, title=title)
    
    
    
# %%

def _run_testcases(plot=True, verbose=True, warnings=True, *args, **kwargs):
    ''' Test procedures
    '''
    
    # Test all Spectrum compare methods 
    # ----------------------------------
    test_compare_methods(verbose=verbose, plot=plot, *args, **kwargs)
    return True


if __name__ == '__main__':
    print('Test_compare.py: ', _run_testcases())
