# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 17:14:01 2018

@author: erwan

Test broadening against HAPI and tabulated data

We're looking at CO(0->1) line 'R1' at 2150.86 cm-1

"""

from __future__ import unicode_literals, print_function, absolute_import, division
from radis.lbl.factory import SpectrumFactory
from radis.spectrum.spectrum import Spectrum
from radis import plot_diff, get_residual_integral, get_residual
from radis.misc.utils import DatabankNotFound
from radis.test.utils import IgnoreMissingDatabase, setup_test_line_databases
from radis.misc.printer import printm
from os.path import join, dirname
import matplotlib.pyplot as plt
import pytest


@pytest.mark.fast
def test_broadening(rtol=1e-2, verbose=True, plot=False, *args, **kwargs):
    '''
    Test broadening against HAPI and tabulated data

    We're looking at CO(0->1) line 'R1' at 2150.86 cm-1
    '''
    from radis.io.hapi import db_begin, fetch, tableList, absorptionCoefficient_Voigt

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/.radis if not there

    # Conditions
    T = 3000
    p = 0.0001
    wstep = 0.001
    wmin = 2150  # cm-1
    wmax = 2152  # cm-1
    broadening_max_width = 10  # cm-1

    # %% HITRAN calculation
    # -----------

    # Generate HAPI database locally
    db_begin(join(dirname(__file__), __file__.replace('.py', '_HAPIdata')))
    if not 'CO' in tableList():   # only if data not downloaded already
        fetch('CO', 5, 1, wmin-broadening_max_width /
              2, wmax+broadening_max_width/2)
        # HAPI doesnt correct for side effects

    # Calculate with HAPI
    nu, coef = absorptionCoefficient_Voigt(SourceTables='CO',
                                           Environment={'T': T,  # K
                                                        'p': p/1.01325,  # atm
                                                        },
                                           WavenumberStep=wstep,
                                           HITRAN_units=False)

    s_hapi = Spectrum.from_array(nu, coef, 'abscoeff', 'cm-1', 'cm_1',
                                 conditions={'Tgas': T},
                                 name='HAPI')

    # %% Calculate with RADIS
    # ----------
    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=1,
        path_length=1,   # doesnt change anything
        wstep=wstep,
        pressure=p,
        broadening_max_width=broadening_max_width,
        isotope=[1],
        warnings={'MissingSelfBroadeningWarning':'ignore',
                  'NegativeEnergiesWarning':'ignore',
                  'HighTemperatureWarning':'ignore',
                  'GaussianBroadeningWarning':'ignore'}
        )  # 0.2)
    sf.load_databank('HITRAN-CO-TEST')
#    s = pl.non_eq_spectrum(Tvib=T, Trot=T, Ttrans=T)
    s = sf.eq_spectrum(Tgas=T, name='RADIS')
    
    
    if plot:  # plot broadening of line of largest linestrength
        sf.plot_broadening(i=sf.df1.S.idxmax())

    # Plot and compare
    res = abs(get_residual_integral(s, s_hapi, 'abscoeff'))
    if plot:
        plot_diff(s, s_hapi, var='abscoeff',
                  title='{0} bar, {1} K (residual {2:.2g}%)'.format(p, T, res*100), show_points=False)
        plt.xlim((wmin, wmax))
    if verbose:
        printm('residual:', res)
    assert res < rtol


@pytest.mark.fast
def test_voigt_broadening_methods(verbose=True, plot=False, *args, **kwargs):
    '''
    Test direct Voigt broadening vs convolution of Gaussian x Lorentzian
    '''

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/.radis if not there

    # Conditions
    T = 3000
    p = 1
    wstep = 0.1
    wmin = 2150  # cm-1
    wmax = 2152  # cm-1
    broadening_max_width = 10  # cm-1
    
    for i, wstep in enumerate([0.01, 0.1, 0.5]):

        # %% Calculate with RADIS
        # ----------
        sf = SpectrumFactory(
            wavenum_min=wmin,
            wavenum_max=wmax,
            mole_fraction=1,
            path_length=1,   # doesnt change anything
            wstep=wstep,
            pressure=p,
            broadening_max_width=broadening_max_width,
            isotope='1',
            verbose=False,
            warnings={'MissingSelfBroadeningWarning':'ignore',
                      'NegativeEnergiesWarning':'ignore',
                      'HighTemperatureWarning':'ignore',
                      'GaussianBroadeningWarning':'ignore'}
            )  # 0.2)
        sf.load_databank('HITRAN-CO-TEST')
    #    s = pl.non_eq_spectrum(Tvib=T, Trot=T, Ttrans=T)
        sf._broadening_method = 'voigt'
        s_voigt = sf.eq_spectrum(Tgas=T, name='direct')
        
        sf._broadening_method = 'convolve'
        s_convolve = sf.eq_spectrum(Tgas=T, name='convolve')
    
        res = get_residual(s_voigt, s_convolve, 'abscoeff')
        
        if verbose:
            print('Residual:', res)

        # plot the last one    
        if plot:
            plot_diff(s_voigt, s_convolve, 'abscoeff', nfig='test_voigt_broadening_methods'+str(i),
                      title='P {0} bar, T {1} K, wstep {2} cm-1'.format(p, T, wstep))
        
        assert res < 2e-4
        

#@pytest.mark.needs_config_file
#@pytest.mark.needs_db_HITEMP_CO2_DUNHAM
@pytest.mark.needs_connection
def test_abscoeff_continuum(plot=False, verbose=2, warnings=True, *args, **kwargs):
    ''' 
    Test calculation with pseudo-continuum

    Assert results on abscoeff dont change
    
    
    Notes
    -----
    
    Uses HITRAN so it can deployed and tested on [Travis]_, but we should switch
    to HITEMP if some HITEMP files can be downloaded automatically at the 
    execution time. 

    '''

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    try:
        if verbose:
            printm('>>> test_abscoeff_continuum')

        sf = SpectrumFactory(wavelength_min=4200, wavelength_max=4500,
                             parallel=False, bplot=False, cutoff=1e-23,
                             molecule='CO2',
                             isotope='1,2', 
                             db_use_cached=True,
                             broadening_max_width=10,
                             path_length=0.1,
                             mole_fraction=1e-3,
                             medium='vacuum',
                             verbose=verbose)
        sf.warnings.update({
            'MissingSelfBroadeningWarning': 'ignore',
            'NegativeEnergiesWarning': 'ignore',
            'LinestrengthCutoffWarning': 'ignore',
            'HighTemperatureWarning':'ignore'})
        sf.fetch_databank()   # uses HITRAN: not really valid at this temperature, but runs on all machines without install
#        sf.load_databank('HITEMP-CO2-DUNHAM')       # to take a real advantage of abscoeff continuum, should calculate with HITEMP
        sf._export_continuum = True   # activate it 

        # Calculate one without pseudo-continuum
        sf.params.pseudo_continuum_threshold = 0
        s1 = sf.eq_spectrum(Tgas=2000)
        s1.name = 'All lines resolved ({0}) ({1:.1f}s)'.format(s1.conditions['lines_calculated'],
                                                               s1.conditions['calculation_time'])
        assert s1.conditions['pseudo_continuum_threshold'] == 0

        # Calculate one with pseudo-continuum
        sf.params.pseudo_continuum_threshold = 0.05
        s2 = sf.eq_spectrum(Tgas=2000)
        s2.name = 'Semi-continuum + {0} lines ({1:.1f}s)'.format(s2.conditions['lines_calculated'],
                                                                 s2.conditions['calculation_time'])
        assert s2.conditions['pseudo_continuum_threshold'] == 0.05
        assert 'abscoeff_continuum' in s2.get_vars()

        # Plot
        if plot:
            plot_diff(s1, s2, 'radiance_noslit', Iunit='µW/cm2/sr/nm',
                      nfig='test_abscoeff_continuum: diff')

            s2.plot('abscoeff', label=s2.name,
                    nfig='test_abscoeff_continuum: show continuum')
            s2.plot('abscoeff_continuum', nfig='same', label='Pseudo-continuum (aggreg. {0:g} lines)'.format(
                    s2.conditions['lines_in_continuum']),
                    force=True)

        # Compare
        res = get_residual(s1, s2, 'abscoeff')
        if verbose:
            printm('residual:', res)

        assert res < 1e-6

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)

#@pytest.mark.needs_config_file
#@pytest.mark.needs_db_HITEMP_CO2_DUNHAM
@pytest.mark.needs_connection
def test_noneq_continuum(plot=False, verbose=2, warnings=True, *args, **kwargs):
    ''' 
    Test calculation with pseudo-continuum under nonequilibrium

    Assert results on emisscoeff dont change

    
    Notes
    -----
    
    Uses HITRAN so it can deployed and tested on [Travis]_, but we should switch
    to HITEMP if some HITEMP files can be downloaded automatically at the 
    execution time. 

    '''
   
    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    try:
        if verbose:
            printm('>>> test_noneq_continuum')

        sf = SpectrumFactory(wavelength_min=4200, wavelength_max=4500,
                             parallel=False, bplot=False, cutoff=1e-23,
                             molecule='CO2',
                             isotope='1,2', db_use_cached=True,
                             broadening_max_width=10,
                             path_length=0.1,
                             mole_fraction=1e-3,
                             medium='vacuum',
                             verbose=verbose)
        sf.warnings.update({
            'MissingSelfBroadeningWarning': 'ignore',
            'NegativeEnergiesWarning': 'ignore',
            'LinestrengthCutoffWarning': 'ignore',
            'HighTemperatureWarning':'ignore'})
        sf.fetch_databank()   # uses HITRAN: not really valid at this temperature, but runs on all machines without install
#        sf.load_databank('HITEMP-CO2-DUNHAM')       # to take a real advantage of abscoeff continuum, should calculate with HITEMP
        sf._export_continuum = True   # activate it 

        # Calculate one without pseudo-continuum
        sf.params.pseudo_continuum_threshold = 0
        s1 = sf.non_eq_spectrum(Tvib=2000, Trot=1000)
        s1.name = 'All lines resolved ({0}) ({1:.1f}s)'.format(s1.conditions['lines_calculated'],
                                                               s1.conditions['calculation_time'])
        assert s1.conditions['pseudo_continuum_threshold'] == 0

        # Calculate one with pseudo-continuum
        sf.params.pseudo_continuum_threshold = 0.05
        s2 = sf.non_eq_spectrum(Tvib=2000, Trot=1000)
        s2.name = 'Semi-continuum + {0} lines ({1:.1f}s)'.format(s2.conditions['lines_calculated'],
                                                                 s2.conditions['calculation_time'])
        assert s2.conditions['pseudo_continuum_threshold'] == 0.05
        assert 'abscoeff_continuum' in s2.get_vars()
        assert 'emisscoeff_continuum' in s2.get_vars()

        # Plot
        if plot:
            plot_diff(s1, s2, 'radiance_noslit', Iunit='µW/cm2/sr/nm',
                      nfig='test_noneq_continuum: diff')

            s2.plot('emisscoeff', label=s2.name,
                    nfig='test_noneq_continuum: show continuum')
            s2.plot('emisscoeff_continuum', nfig='same', label='Pseudo-continuum (aggreg. {0:g} lines)'.format(
                    s2.conditions['lines_in_continuum']), force=True)

        # Compare
        res = get_residual(s1, s2, 'abscoeff') + get_residual(s1, s2, 'emisscoeff')
        if verbose:
            printm('residual:', res)

        assert res < 5e-6

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)



def _run_testcases(plot=False, verbose=True, *args, **kwargs):

    # Test broadening
    test_broadening(plot=plot, verbose=verbose, *args, **kwargs)
#    test_voigt_broadening_methods(plot=plot, verbose=verbose, *args, **kwargs)
#
#    # Test pseudo-continuum
#    test_abscoeff_continuum(plot=plot, verbose=verbose, *args, **kwargs)
#    test_noneq_continuum(plot=plot, verbose=verbose, *args, **kwargs)

    return True


if __name__ == '__main__':
    printm('test_broadening: ', _run_testcases(
        plot=True, verbose=True, debug=False))
