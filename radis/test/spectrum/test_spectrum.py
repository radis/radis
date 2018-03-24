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
from radis.spectrum import Spectrum, calculated_spectrum
from radis.phys.convert import nm2cm
#from radis.misc.utils import DatabankNotFound
#from radis.test.utils import IgnoreMissingDatabase, build_test_databases
import numpy as np
from numpy import allclose, linspace
from os.path import basename

fig_prefix = basename(__file__)+': '

# %% Test routines

def test_spectrum_get_methods(verbose=True, *args, **kwargs):
    ''' Test all spectrum methods on a Spectrum generated in Specair '''
    
    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec
    from radis.tools.slit import get_FWHM

    s = load_spec(getTestFile('N2C_specair_380nm.spec'))
    
    if verbose: print(s)
    assert s.get_name() == 'N2C_specair_380nm'
    assert all(s.get_radiance_noslit(Iunit='W/m2/sr/nm') == 
               s.get('radiance_noslit', Iunit='W/m2/sr/nm')[1])
    assert all(nm2cm(s.get_wavelength(medium='vacuum')) == 
               s.get_wavenumber())
    assert s.get_power(unit='W/cm2/sr') == 2631.6288408588148
    assert s.get_waveunit() == 'nm'
    assert s.get_power(unit='W/cm2/sr') == s.get_integral('radiance_noslit', Iunit='W/cm2/sr/nm')
    assert s.get_conditions()['Tgas'] == 1500
    assert len(s.get_vars()) == 2
    assert s.is_at_equilibrium() == False
    assert s.is_optically_thin() == False
    
    # Check applied slit has the correct width
    s.apply_slit(0.5)
    wslit, Islit = s.get_slit()
    wstep = np.diff(wslit)[0]
    assert np.isclose(get_FWHM(*s.get_slit()), 0.5, atol=1.1*wstep)
    
    if verbose:
        print('Tested Spectrum methods:')
        print('...print(Spectrum)')
        print('.. get_name()')
        print('.. get_radiance_noslit() vs get()')
        print('.. get_wavelength() vs get_wavenumber')
        print('.. get_power()')
        print('.. get_waveunit()')
        print('.. get_power() vs get_integral()')
        print('.. get_conditions()')
        print('.. get_vars()')
        print('.. is_at_equilibrium()')
        print('.. is_optically_thin()')
        print('.. get_slit()')

def test_copy(verbose=True, *args, **kwargs):
    ''' Test that a Spectrum is correctly copied 
    
    We compare a Spectrum that has:
    - all available spectral quantities
    - a slit 
    - many calculation conditions
    - no populations
    - no lines 
    '''
    
    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec
    
    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'))
    s.update()
    s.apply_slit(1.5)
    
    s2 = s.copy()
    
    assert s == s2
    assert s is not s2
    
    if verbose:
        print('Tested that s2 == s after Spectrum copy')
        print('Warning: current test doesnt include populations nor lines')
    

def test_intensity_conversion__fast(verbose=True, *args, **kwargs):
    ''' Test conversion of intensity cm-1 works'''

    from radis import planck, planck_wn

    w_nm = linspace(300, 3000)
    I_nm = planck(w_nm, T=6000, unit='mW/sr/cm2/nm')

    s = calculated_spectrum(w_nm, I_nm, wunit='nm', Iunit='mW/sr/cm2/nm', 
                            conditions={'medium':'vacuum'})
    w, I = s.get('radiance_noslit', Iunit='mW/sr/cm2/cm_1')

    w_cm = nm2cm(w_nm)
    I_cm = planck_wn(w_cm, T=6000, unit='mW/sr/cm2/cm_1')

    assert allclose(I_cm, I, rtol=1e-3)

def test_rescaling_function__fast(verbose=True, *args, **kwargs):
    ''' Test rescaling functions '''

    from radis.test.utils import getTestFile

    s = Spectrum.from_txt(getTestFile('calc_N2C_spectrum_Trot1200_Tvib3000.txt'),
                          quantity='radiance_noslit', 
                          waveunit='nm',
                          unit='mW/cm2/sr/µm',   # Specair units: mW/cm2/sr/µm
                          conditions={'Tvib':3000, 'Trot':1200,
                                        'path_length':1,  # arbitrary
                                        'medium':'air',
                                        },
                          populations={'molecules':{'N2C':1e13}}, # arbitrary
                                                            # (just an example)
                          )
    s.update(optically_thin=True)
    s.rescale_path_length(10)

    assert np.isclose(s.get_radiance_noslit(Iunit='mW/cm2/sr/nm')[0], 352.57305783248)


def _run_testcases(plot=False, verbose=True, debug=False, warnings=True, *args, **kwargs):
    ''' Test procedures

    Input
    ------

    debug: boolean
        swamps the console namespace with local variables. Default False

    '''
    
    # Test all Spectrum get methods 
    # -----------------------------
    test_spectrum_get_methods(debug=debug, verbose=verbose, *args, **kwargs)
    
    
    # Test populations
    # ----------
#    test_populations__fast(verbose=verbose, *args, **kwargs)

    # Test conversion of intensity cm-1 works
    # -------------
    test_intensity_conversion__fast(debug=debug, verbose=verbose, *args, **kwargs)

    # Test updating / rescaling functions (no self absorption)
    # ---------
    test_rescaling_function__fast(debug=debug, *args, **kwargs)
#    test_rescaling_path_length__fast(plot=plot, verbose=verbose, debug=debug,
#                                     warnings=warnings, *args, **kwargs)
#    test_rescaling_mole_fraction__fast(plot=plot, verbose=verbose, debug=debug,
#                                     warnings=warnings, *args, **kwargs)
    
    # Test propagating medium
#    test_medium__fast(plot=plot, verbose=verbose, debug=debug,
#                           warnings=warnings, *args, **kwargs)

    return True


if __name__ == '__main__':
    print(('Test spectrum: ', _run_testcases(debug=False)))
