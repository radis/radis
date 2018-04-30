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
import os
from os.path import basename, exists
import pytest

fig_prefix = basename(__file__)+': '

# %% Test routines

def test_spectrum_get_methods(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    ''' Test all spectrum methods on a Spectrum generated in Specair '''
    
    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec
    from radis.tools.slit import get_FWHM


    if plot and close_plots:
        import matplotlib.pyplot as plt
        plt.close('all')

    s = load_spec(getTestFile('N2C_specair_380nm.spec'))
    
    # general methods
    if verbose: print(s)
    dir(s)
    
    # access properties
    assert s.get_name() == 'N2C_specair_380nm'
    assert all(s.get_radiance_noslit(Iunit='W/m2/sr/nm') == 
               s.get('radiance_noslit', Iunit='W/m2/sr/nm')[1])
    assert all(nm2cm(s.get_wavelength(medium='vacuum')) == 
               s.get_wavenumber())
    assert s.get_power(unit='W/cm2/sr') == 2631.6288408588148
    assert s.get_waveunit() == 'nm'
    assert s.get_power(unit='W/cm2/sr') == s.get_integral('radiance_noslit', Iunit='W/cm2/sr/nm')
    assert s.get_conditions()['Tgas'] == 1500
    assert s.get_medium() == 'air'
    assert len(s.get_vars()) == 2
    assert s.is_at_equilibrium() == False
    assert s.is_optically_thin() == False
    
    # Check applied slit has the correct width
    s.apply_slit(0.5)
    wslit, Islit = s.get_slit()
    wstep = np.diff(wslit)[0]
    assert np.isclose(get_FWHM(*s.get_slit()), 0.5, atol=1.1*wstep)
    
    if plot:
        s.plot_slit()
    
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
    
def test_populations(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    ''' Test that populations in a Spectrum are correctly read 
    '''
    
    if plot:
        import matplotlib.pyplot as plt
        plt.ion()   # dont get stuck with Matplotlib if executing through pytest
    

    if plot and close_plots:
        import matplotlib.pyplot as plt
        plt.close('all')

    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec
    import pytest
    
    # get a spectrum
    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'))
    
    # Get all populations
    pops = s.get_populations(molecule='CO', isotope=1)
    assert np.isclose(pops['Ia'], 0.986544)
    
    # Get vib levels
    df_vib = s.get_vib_levels()
    assert len(df_vib) == 49
    
    # Get rovib levels (dont exist: should fail!)
    with pytest.raises(KeyError):  # expected behavior
        s.get_rovib_levels()
        
    if plot:
        s.plot_populations()
    
    if verbose:
        print('test_populations: OK')
    
def test_store_functions(verbose=True, *args, **kwargs):
    ''' Test some store / retrieve functions '''
    
    from radis.spectrum.spectrum import transmittance_spectrum
    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec
    
    temp_file = 'test_radis_tempfile_transmittance.txt'
    assert not exists(temp_file)
         
    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'), binary=True)
    s.update()
    try:
        s.savetxt(temp_file, 'transmittance_noslit', wunit='nm', medium='vacuum')
        w, T = np.loadtxt(temp_file).T
    finally:
        os.remove(temp_file)

    s2 = transmittance_spectrum(w, T, wunit='nm', conditions={'medium':'vacuum'})       
    assert s.compare_with(s2, spectra_only='transmittance_noslit', plot=False)
    
    return True


@pytest.mark.fast
def test_intensity_conversion(verbose=True, *args, **kwargs):
    ''' Test conversion of intensity cm-1 works:
        
    - conversion of mW/sr/cm2/nm -> mW/sr/cm2/cm-1
    
    '''

    from radis import planck, planck_wn

    w_nm = linspace(300, 3000)
    w_cm = nm2cm(w_nm)
    I_nm = planck(w_nm, T=6000, unit='mW/sr/cm2/nm')

    s = calculated_spectrum(w_nm, I_nm, wunit='nm', Iunit='mW/sr/cm2/nm', 
                            conditions={'medium':'vacuum'})
    
    # mW/sr/cm2/nm -> mW/sr/cm2/cm-1
    w, I = s.get('radiance_noslit', Iunit='mW/sr/cm2/cm_1')
    I_cm = planck_wn(w_cm, T=6000, unit='mW/sr/cm2/cm_1')
    assert allclose(I_cm, I, rtol=1e-3)
    
#def test_emissivity_conversion__fast(verbose=True, *args, **kwargs):
#    ''' Test conversion of intensity cm-1 works:
#        
#    - conversion of mW/sr/cm2/nm -> mW/sr/cm2/cm-1
#    
#    '''
#    
#    from radis.phys.units import convert_emi2nm
#    from radis.test.utils import getTestFile
#    from radis.tools.database import load_spec
#    
#    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'), binary=True)
#    s.update()
#
#    # mW/sr/cm3/nm -> mW/sr/cm3/cm-1
#    
#    w, I = s.get('emissivity_noslit', Iunit='mW/sr/cm3/cm_1')
#    I_cm = convert_emi2nm(I, 'mW/sr/cm3/cm_1', 'mW/sr/m3/µm')
#    
    
# TODO: finish implementing emissivity_conversino above


def test_rescaling_function(verbose=True, *args, **kwargs):
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

def test_resampling_function(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    ''' Test resampling functions 
    
    Get a Spectrum calculated in cm-1, then resample on a smaller range in cm-1, 
    and in approximately the same range (but in nm). Check that all 3 overlap 
    '''
# %%
    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec
    from radis.spectrum import get_residual

    if plot and close_plots:
        import matplotlib.pyplot as plt
        plt.close('all')

    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'), binary=True)
    s.name = 'original'
    s2 = s.copy()
    s2b = s.copy()
    s3 = s.copy()
    s2.resample(np.linspace(4500, 4700, 10000), unit='nm', medium='vacuum')
    s2b.resample(np.linspace(4500, 4700, 10000), unit='nm', medium='air')
    s3.resample(np.linspace(2127.2, 2227.7, 10000), unit='cm-1')
    s2.name = 'resampled in nm (vacuum)'
    s2b.name = 'resampled in nm (air)'
    s3.name = 'resampled in cm-1'
    
    s.compare_with(s2, plot=plot, title='Residual: {0:.2g}'.format(
            get_residual(s, s2, 'abscoeff', ignore_nan=True)))
    s.compare_with(s2b, plot=plot, title='Residual: {0:.2g}'.format(
            get_residual(s, s2b, 'abscoeff', ignore_nan=True)))
    s.compare_with(s3, plot=plot, title='Residual: {0:.2g}'.format(
            get_residual(s, s3, 'abscoeff', ignore_nan=True)))
    
    assert get_residual(s, s2, 'abscoeff', ignore_nan=True) < 1e-4
    assert get_residual(s, s2b, 'abscoeff', ignore_nan=True) < 1e-3
    assert get_residual(s, s3, 'abscoeff', ignore_nan=True) < 1e-5
    
# %%

def _run_testcases(plot=True, close_plots=False, verbose=True, debug=False, warnings=True, *args, **kwargs):
    ''' Test procedures

    Input
    ------

    debug: boolean
        swamps the console namespace with local variables. Default False

    '''
    
    # Test all Spectrum methods 
    # -------------------------
    test_spectrum_get_methods(debug=debug, verbose=verbose, plot=plot, close_plots=close_plots, *args, **kwargs)
    test_copy(verbose=verbose, *args, **kwargs)
    test_populations(verbose=verbose, plot=plot, close_plots=close_plots, *args, **kwargs)
    test_store_functions(verbose=verbose, *args, **kwargs)
    
    
    # Test populations
    # ----------
#    test_populations__fast(verbose=verbose, *args, **kwargs)

    # Test conversion of intensity cm-1 works
    # -------------
    test_intensity_conversion(debug=debug, verbose=verbose, *args, **kwargs)

    # Test updating / rescaling functions (no self absorption)
    # ---------
    test_rescaling_function(debug=debug, *args, **kwargs)
    test_resampling_function(debug=debug, plot=plot, close_plots=close_plots, *args, **kwargs)
#    test_rescaling_path_length__fast(plot=plot, verbose=verbose, debug=debug,
#                                     warnings=warnings, *args, **kwargs)
#    test_rescaling_mole_fraction__fast(plot=plot, verbose=verbose, debug=debug,
#                                     warnings=warnings, *args, **kwargs)
    
    # Test propagating medium
#    test_medium__fast(plot=plot, verbose=verbose, debug=debug,
#                           warnings=warnings, *args, **kwargs)

    return True


if __name__ == '__main__':
    print(('Test spectrum: ', _run_testcases(debug=False, close_plots=False)))
