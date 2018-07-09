# -*- coding: utf-8 -*-
"""

Examples
--------

Run all tests::

    pytest      # (in command line, in project folder)

Run only fast tests (i.e: tests that a  'fast' label)::

    pytest -m fast

-------------------------------------------------------------------------------

"""

from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from radis.los.slabs import MergeSlabs, SerialSlabs
import pytest

@pytest.mark.fast
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
        if close_plots:
            plt.close('all')

    for optically_thin in [True, False]:

        # Get Some spectra
        s1 = load_spec(getTestFile(
            'CO_Tgas1500K_mole_fraction0.01.spec'), binary=True)
        s2 = load_spec(getTestFile(
            'CO_Tgas1500K_mole_fraction0.5.spec'), binary=True)
        s1.conditions['thermal_equilibrium'] = True   # fix: condition not given in radis < 0.2.2
        s2.conditions['thermal_equilibrium'] = True   # fix: condition not given in radis < 0.2.2
        s1.update('all')
        s2.update('all')

        # Merge 50 times s1 in the same slab
        s50 = [s1]*50
        s1N = MergeSlabs(*s50, optically_thin=optically_thin, debug=debug)

        if plot:
            for k in ['radiance_noslit']:  # , 'transmittance_noslit']:
                s2.plot(k, lw=3, label='1x[CO=0.5]')
                s1N.plot(k, nfig='same', label='50x[CO=0.01]')
                plt.legend()
                plt.title('Optically thin: {0}'.format(optically_thin))
                plt.tight_layout()

        if verbose:
            print('test_merge_slabs')
            print('... Compare 50x[CO=0.01] vs 1x[CON=0.5] (optically thin: {0})'.format(
                optically_thin))
            print('... Difference: {0:.2f}%'.format(
                abs(s1N.get_power()/s2.get_power()-1)*100))
        assert np.isclose(s2.get_power(), s1N.get_power(), 1.5e-2)

    return True
#


@pytest.mark.fast
def test_serial_slabs(verbose=True, plot=True, warnings=True, debug=False,
                     *args, **kwargs):
    ''' Add some slabs in serie, ensure that overall transmittance decreases 
    Also check that some quantities are propagated (like path length?)

    '''

    from radis.tools.database import load_spec
    import matplotlib.pyplot as plt
    from radis.test.utils import getTestFile
    if plot:
        plt.ion()   # dont get stuck with Matplotlib if executing through pytest

    # Get Some spectra
    s = load_spec(getTestFile(
        'CO_Tgas1500K_mole_fraction0.5.spec'), binary=True)
#    s.update('all')
    assert 'transmittance_noslit' not in s.get_vars()
    
    s_los = SerialSlabs(s, s, s)
    
    w, T = s_los.get('transmittance_noslit')  
    # note that it was calculated despite transmittance not in the initial
    # spectrum
    
    s.update('transmittance_noslit')   # get from abscoeff
    w1, T1 = s.get('transmittance_noslit')
    
    # Check it looks alright
    assert (T <= T1).all()
    assert (T != T1).any()
    assert (w == w1).all()
    
    # because they were the same slabs, intensive conditions should have been propagated
    assert(s.conditions['Tgas'] == s_los.conditions['Tgas'])
    # but extensive shouldnt:
    assert(s.conditions['path_length']*3 == s_los.conditions['path_length'])
    
    return True
#

from radis.misc.basics import all_in

def test_serial_slabs_updates(verbose=True, warnings=True, *args, **kwargs):
    ''' Test how quantities are updated in a SerialSlab operation

    '''

    from radis.tools.database import load_spec
    from radis.test.utils import getTestFile

    # Get Some spectra
    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.5.spec'), binary=True)
    s.conditions['thermal_equilibrium'] = True   # fix: condition not given in radis < 0.2.2
    
    # Ensure Update works
    assert s.get_vars() == ['abscoeff']
    s.update()
    assert all_in(['radiance_noslit', 'transmittance_noslit', 'emisscoeff', 'absorbance'],
                  s.get_vars())
    
    # Ensure Update() doesnt modify existing quantities
    # ... Hack: kill radiance (see how it's used to generate a Pure Transmittance
    # ... spectrum in validation/compare_torch
    s.conditions['thermal_equilibrium'] = False  # ensure radiance not recomputed from abscoeff + Kirchoff
    s._q['radiance_noslit'] *= 0
    s.update('radiance_noslit')
    assert (s._q['radiance_noslit'] == 0).all()
    


def _run_testcases(verbose=True, plot=True, close_plots=True, debug=False, warnings=True, *args, **kwargs):

    test_merge_slabs(verbose=verbose, plot=plot, close_plots=close_plots, debug=debug, warnings=warnings,
                     *args, **kwargs)

    test_serial_slabs(verbose=verbose, plot=plot, debug=debug, warnings=warnings,
                     *args, **kwargs)

    test_serial_slabs_updates(verbose=verbose, warnings=warnings,
                     *args, **kwargs)
    
    return True


if __name__ == '__main__':
    print('test_slabs: ', _run_testcases(verbose=True))
