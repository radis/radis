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
from radis.los.slabs import MergeSlabs, SerialSlabs

def test_merge_slabs(verbose=True, plot=True, close_plots=True, warnings=True, debug=False, 
                      *args, **kwargs):
    ''' Merge 10 slabs with 1/10 of concentration, and compare the results. 
    Ensure error is < 0.1%
    
    Note that we won't have exactly the same results as the broadening is not
    updated in MergeSlabs 
    '''
    
    from radis.tools.database import load_spec
    import matplotlib.pyplot as plt
    from radis.test.utils import getTestFile
    if plot:
        plt.ion()   # dont get stuck with Matplotlib if executing through pytest
        if close_plots: plt.close('all')

    for optically_thin in [True, False]:
    
        # Get Some spectra
        s1 = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'))
        s2 = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.5.spec'))
        s1.update('all')
        s2.update('all')
        
        # Merge 50 times s1 in the same slab
        s50 = [s1]*50
        s1N = MergeSlabs(*s50, optically_thin=optically_thin, debug=debug)
        
        if plot:
            for k in ['radiance_noslit']: #, 'transmittance_noslit']:
                s2.plot(k, lw=3,label='1x[CO=0.5]')
                s1N.plot(k, nfig='same', label='50x[CO=0.01]')
                plt.legend()
                plt.title('Optically thin: {0}'.format(optically_thin))
                plt.tight_layout()
            
        if verbose:
            print('test_merge_slabs')
            print('... Compare 50x[CO=0.01] vs 1x[CON=0.5] (optically thin: {0})'.format(
                    optically_thin))
            print('... Difference: {0:.2f}%'.format(abs(s1N.get_power()/s2.get_power()-1)*100))
        assert np.isclose(s2.get_power(), s1N.get_power(), 1.5e-2)
    
    return True
#        
        
def _run_testcases(verbose=True, plot=True, close_plots=True, debug=False, warnings=True, *args, **kwargs):
    
    test_merge_slabs(verbose=verbose, plot=plot, close_plots=close_plots, debug=debug, warnings=warnings,
                           *args, **kwargs)

    # Todo: add a test with _serial_slabs ... 
    if verbose: print('WARNING: no test for SerialSlabs defined yet')
        
    return True

if __name__ == '__main__':
    print('test_slabs: ', _run_testcases(verbose=True))
