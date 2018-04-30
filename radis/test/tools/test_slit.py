# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 10:18:48 2017

@author: erwan

Examples
--------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have 'fast' in their name)

>>> pytest -k fast

"""

from __future__ import print_function, absolute_import, division, unicode_literals
#from neq.spec import SpectrumFactory
from radis.spectrum.spectrum import calculated_spectrum, transmittance_spectrum
from radis.tools.database import load_spec
from radis.los.slabs import SerialSlabs
from radis.tools.slit import (gaussian_slit, triangular_slit, trapezoidal_slit,
                           import_experimental_slit, convolve_with_slit, 
                           get_FWHM, get_effective_FWHM)
from radis.phys.units import is_homogeneous
from radis.phys.convert import dcm2dnm, dnm2dcm
#from radis.misc.utils import DatabankNotFound
#from radis.test.utils import IgnoreMissingDatabase, build_test_databases
import matplotlib.pyplot as plt 
import numpy as np
from numpy import sqrt, linspace, abs, trapz
from os.path import basename

fig_prefix = basename(__file__)+': '

# %% --------------------------------------------------------------------------
#  Test cases
# -----------------------------------------------------------------------------

def test_all_slit_shapes__fast(FWHM=0.4, verbose=True, plot=True, close_plots=True, *args, **kwargs):
    ''' Test all slit generation functions and make sure we get the expected FWHM'''

    if plot:
        plt.ion()   # dont get stuck with Matplotlib if executing through pytest
        if close_plots: plt.close('all')
    
    # get spectrum
    from radis.test.utils import getTestFile
    from radis.spectrum.spectrum import Spectrum
    s = Spectrum.from_txt(getTestFile('calc_N2C_spectrum_Trot1200_Tvib3000.txt'),
                          quantity='radiance_noslit', waveunit='nm', unit='mW/cm2/sr/µm',
                          conditions={'medium':'air'})
    wstep = np.diff(s.get_wavelength())[0]
    
    # Plot all slits
    # ... gaussian
    s.apply_slit(FWHM, unit='nm', shape='gaussian', plot_slit=plot)
    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=2*wstep)
    
    # ... triangular
    s.apply_slit(FWHM, unit='nm', shape='triangular', plot_slit=plot)
    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=2*wstep)
    
    # ... trapezoidal
    s.apply_slit((FWHM*0.9, FWHM*1.1), unit='nm', shape='trapezoidal', plot_slit=plot)
    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=2*wstep)
    
#    # ... trapezoidal
#    s.apply_slit(FWHM, unit='nm', shape='trapezoidal', plot_slit=plot, norm='max')
#    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=1.1*wstep)
    
    # ... experimental
    s.apply_slit(getTestFile('slitfunction.txt'), unit='nm', plot_slit=plot)
    assert np.isclose(get_effective_FWHM(*s.get_slit()), FWHM, atol=0.01)
    # note that we're applying a slit function measured at 632.5 nm to a Spectrum 
    # at 4.7 µm. It's just good for testing the functions
    
#    # ... experimental, convolve with max
#    s.apply_slit(getTestFile('slitfunction.txt'), unit='nm', norm_by='max', plot_slit=plot)
#    assert np.isclose(get_FWHM(*s.get_slit()), FWHM, atol=1.1*wstep)
    
    if verbose: print('\n>>> _test_all_slits yield correct FWHM (+- wstep) : OK\n')
    
    return True # nothing defined yet

def test_slit_unit_conversions_spectrum_in_cm(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    ''' Test that slit is consistently applied for different units
    
    Assert that:
    
    - calculated FWHM is the one that was applied
    
    '''
    
    from radis.test.utils import getTestFile
    from radis.tools.database import load_spec
    
    if plot:    # dont get stuck with Matplotlib if executing through pytest
        plt.ion()
        if close_plots: plt.close('all')
    
    # %% Get a Spectrum (stored in cm-1)
    s_cm = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'), binary=True)
    s_cm.rescale_mole_fraction(1)   # just because it makes better units
    s_cm.update()
    wstep = s_cm.conditions['wstep']
    
    assert s_cm.get_waveunit() == 'cm-1'       # ensures it's stored in cm-1
    
    for shape in ['gaussian', 'triangular']:
        
        # Apply slit in cm-1 
        slit_cm = 2    
        s_cm.name = 'Spec in cm-1, slit {0:.2f} cm-1'.format(slit_cm)
        s_cm.apply_slit(slit_cm, unit='cm-1', shape=shape, mode='same')  
        # ... mode=same to keep same output length. It helps compare both Spectra afterwards
        fwhm = get_FWHM(*s_cm.get_slit())   # in cm-1 as that's s.get_waveunit()
        assert np.isclose(slit_cm, fwhm, atol=2*wstep)
            
        # Apply slit in nm this time
        s_nm = s_cm.copy()
        w_cm = s_nm.get_wavenumber(which='non_convoluted')
        slit_nm = dcm2dnm(slit_cm, w_cm[len(w_cm)//2])
        s_nm.name = 'Spec in cm-1, slit {0:.2f} nm'.format(slit_nm)
        s_nm.apply_slit(slit_nm, unit='nm', shape=shape, mode='same')
        
        plotargs = {}
        if plot: plotargs['title'] = 'test_slit_unit_conversions: {0} ({1} cm-1)'.format(shape, slit_cm)
        s_cm.compare_with(s_nm, spectra_only='radiance', rtol=1e-3, 
                          verbose=verbose, plot=plot, **plotargs)
    
def test_slit_unit_conversions_spectrum_in_nm(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    ''' Test that slit is consistently applied for different units
    
    Assert that:
    
    - calculated FWHM is the one that was applied
    
    '''
    
    from radis.test.utils import getTestFile
    from radis.spectrum.spectrum import Spectrum
 
    
    if plot:    # dont get stuck with Matplotlib if executing through pytest
        plt.ion()
        if close_plots: plt.close('all')
   
    # %% Get a Spectrum (stored in nm)
    s_nm = Spectrum.from_txt(getTestFile('calc_N2C_spectrum_Trot1200_Tvib3000.txt'),
                          quantity='radiance_noslit', waveunit='nm', unit='mW/cm2/sr/µm',
                          conditions={'medium':'air',
                                      'self_absorption':False})
    s_nm.rescale_path_length(1, 0.001)   # just because it makes better units
    wstep = np.diff(s_nm.get_wavelength())[0]
    
    assert s_nm.get_waveunit() == 'nm'       # ensures it's stored in cm-1
    
    for shape in ['gaussian', 'triangular']:
        
        # Apply slit in nm
        slit_nm = 0.5    
        s_nm.name = 'Spec in nm, slit {0:.2f} nm'.format(slit_nm)
        s_nm.apply_slit(slit_nm, unit='nm', shape=shape, mode='same')  
        # ... mode=same to keep same output length. It helps compare both Spectra afterwards
        fwhm = get_FWHM(*s_nm.get_slit())   # in cm-1 as that's s.get_waveunit()
        assert np.isclose(slit_nm, fwhm, atol=2*wstep)
    
        # Apply slit in nm this time
        s_cm = s_nm.copy()
        w_nm = s_nm.get_wavelength(which='non_convoluted')
        slit_cm = dnm2dcm(slit_nm, w_nm[len(w_nm)//2])
        s_cm.name = 'Spec in nm, slit {0:.2f} cm-1'.format(slit_cm)
        s_cm.apply_slit(slit_cm, unit='cm-1', shape=shape, mode='same')
        
        plotargs = {}
        if plot: plotargs['title'] = 'test_slit_unit_conversions: {0} ({1} nm)'.format(shape, slit_nm)
        s_nm.compare_with(s_cm, spectra_only='radiance', rtol=1e-3, 
                          verbose=verbose, plot=plot, **plotargs)
    
    
    # %% 
    
def test_against_specair_convolution__fast(plot=True, close_plots=True, verbose=True, debug=False, *args, **kwargs):

    if plot:
        plt.ion()   # dont get stuck with Matplotlib if executing through pytest
        if close_plots: plt.close('all')
    
    # Test
    from radis.test.utils import getTestFile
    
    # Plot calculated vs convolved with slit
    w, I = np.loadtxt(getTestFile('calc_N2C_spectrum_Trot1200_Tvib3000.txt')).T   # Specair units: mW/cm2/sr/µm
    s = calculated_spectrum(w, I, conditions={'Tvib':3000, 'Trot':1200, 
                                              'medium':'air'},
                             Iunit='mW/cm2/sr/µm')

    if plot:
        fig = plt.figure(fig_prefix+'SPECAIR convoluted vs SPECAIR non convoluted')
        s.plot('radiance_noslit', nfig=fig.number, label='calc')
    slit_nm = 0.1
    s.apply_slit(slit_nm, shape='triangular', norm_by='area')
    if plot:
        s.plot('radiance', nfig=fig.number, label='slit {0}nm'.format(slit_nm),
               color='r', lw=2)
        plt.legend()

    # Compare with Specair slit function
    s.apply_slit(slit_nm, norm_by='max')
    # Note unit conversion from mW/cm2/sr/um*nm is properly done!
    if plot:
        fig = plt.figure(fig_prefix+'convoluted RADIS vs convoluted SPECAIR')
        s.plot('radiance', nfig=fig.number, label='slit {0}nm, norm with max'.format(slit_nm),
               color='r', lw=2, Iunit='mW/cm2/sr')
    # ... Get initial spectrum convoluted with Specair
    ws, Is = np.loadtxt(getTestFile('calc_N2C_spectrum_Trot1200_Tvib3000_slit0.1.txt')).T  # Specair units: mW/cm2/sr
    if plot:
        plt.plot(ws, Is, 'k', label='normalized in Specair (sides not cropped)')

    # ... Test output is the same
    wu, Iu = s.get('radiance', Iunit='mW/cm2/sr')

    As = np.trapz(Is, ws)  # Todo one day: replace with get_power() function
    Au = np.trapz(Iu, wu)
    if verbose:
        print(('Integrals should match: {0:.2f} vs {1:.2f} ({2:.2f}% error)'.format(
                As, Au, 100*abs(As-Au)/As)))
    assert np.isclose(As, Au, rtol=1e-2)

    # Test resampling
    s.resample(linspace(376,380.6,3000), unit='nm', if_conflict_drop='convoluted')
    s.apply_slit(slit_nm, norm_by='max')
    if plot:
        s.plot('radiance', Iunit='mW/cm2/sr', lw=3, nfig=fig.number, color='b',
               zorder=-3, label='Resampled')
        plt.legend()

    if verbose: print('\n>>>Testing spectrum slit matches Specair: OK')
    
    return True

def test_normalisation_mode__fast(plot=True, close_plots=True, verbose=True, *args, **kwargs):
    ''' Test norm_by = 'area' vs norm_by = 'max' '''

    from radis.test.utils import getTestFile
    
    if plot:
        plt.ion()   # dont get stuck with Matplotlib if executing through pytest
        if close_plots: plt.close('all')
    
        from publib import set_style
        set_style('origin')
    
    # %% Compare spectra convolved with area=1 and max=1
    # Slit in nm
    # Spectrum in nm
    w, I = np.loadtxt(getTestFile('calc_N2C_spectrum_Trot1200_Tvib3000.txt')).T   # Specair units: mW/cm2/sr/µm
    s = calculated_spectrum(w, I, conditions={'Tvib':3000, 'Trot':1200},
                             Iunit='mW/cm2/sr/µm')
    
    FWHM = 2
    
    s.apply_slit(FWHM, norm_by='area')   # spectrum convolved with area=1
    w_area, I_area = s.get('radiance')
    if plot:
        fig = plt.figure(fig_prefix+'Spectrum in nm + slit in nm')
        fig.clear()
        ax = fig.gca()
        s.plot(nfig=fig.number, wunit='nm', label='norm_by: area', lw=3)
    s.apply_slit(FWHM, norm_by='max')    # spectrum convolved with max=1     
    w_max, I_max = s.get('radiance', wunit='nm')
    if plot:
        ax.plot(w_max, I_max/FWHM, 'r', label='(norm_by:max)/FWHM')
        ax.legend(loc='best')
    assert np.allclose(I_area, I_max/FWHM)
    if verbose: 
        print("equivalence of normalisation mode for spectrum in 'nm': OK") 
    
    # %% Compare spectra convolved with area=1 and max=1
    # Slit in nm
    # Spectrum in cm-1
    
    s = load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'))
    s.update()
    s.apply_slit(FWHM, norm_by='area', plot_slit=plot)   # spectrum convolved with area=1
    w_area, I_area = s.get('radiance')
    if plot:
        fig = plt.figure(fig_prefix+'Spectrum in cm-1 + slit in nm')
        fig.clear()
        ax = fig.gca()
        s.plot(nfig=fig.number, wunit='nm', label='norm_by: area', lw=3)
    s.apply_slit(FWHM, norm_by='max', plot_slit=plot)    # spectrum convolved with max=1
    w_max, I_max = s.get('radiance', wunit='nm')
    if plot:
        ax.plot(w_max, I_max/FWHM, 'r', label='(norm_by:max)/FWHM')
        ax.legend(loc='best')
    assert np.allclose(I_area, I_max/FWHM)
    if verbose: 
        print("equivalence of normalisation mode for spectrum in 'cm-1': {0}: OK")
    assert is_homogeneous(s.units['radiance'], 'mW/cm2/sr')
    if verbose: 
        print(("radiance unit ({0}) is homogeneous to 'mW/cm2/sr': OK".format(s.units['radiance'])))
    
    return True


def test_slit_energy_conservation__fast(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    ''' Convoluted and non convoluted quantities should have the same area
    (difference arises from side effects if the initial spectrum is not 0 on 
    the sides '''
    
    from radis.test.utils import getTestFile
    
    if plot:
        import matplotlib.pyplot as plt
        plt.ion()   # dont get stuck with Matplotlib if executing through pytest
    if close_plots: plt.close('all')
    
    if verbose: print('\n>>> _test_slit_energy_conservation\n')

    s = calculated_spectrum(*np.loadtxt(getTestFile('calc_N2C_spectrum_Trot1200_Tvib3000_slit0.1.txt')).T, 
                            wunit='nm', Iunit='mW/cm2/sr/nm')   # arbitrary)
    
    P = s.get_power(unit='mW/cm2/sr')
    s.apply_slit(1.5, norm_by='area')
    w, I = s.get('radiance', wunit='nm', Iunit='mW/cm2/sr/nm')
    Pc = abs(np.trapz(I, x=w))  # mW/cm2/sr
    
    b = np.isclose(P,Pc, 3e-2)   

    if plot:
        fig = plt.figure(fig_prefix+'energy conservation during resampling')
        s.plot(nfig=fig.number, label='{0:.1f} mW/cm2/sr'.format(P))
        s.plot('radiance_noslit', nfig=fig.number, label='{0:.1f} mW/cm2/sr'.format(Pc))
        plt.title('Energy conservation: {0}'.format(b)) 
        plt.legend()
        plt.tight_layout()
        
    assert np.isclose(P,Pc, 3e-2)   
        
    return True
    

def test_slit_function_effect__fast(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    ''' A test case to show the effect of wavelength dispersion (cf spectrometer
    reciprocal function) on the slit function '''

    from radis.test.utils import getTestFile
    from publib import fix_style
    from numpy import pi, tan, cos

    if plot:
        plt.ion()   # dont get stuck with Matplotlib if executing through pytest
        if close_plots: plt.close('all')
    
    # Effect of Slit dispersion

    def dirac(w0, width=20, wstep=0.009):
        ''' Return a Dirac in w0. Plus some space on the side '''

        w = np.arange(w0-width, w0+width+wstep, wstep)
        I = np.zeros_like(w)
        I[len(I)//2] = 1/wstep

        return w, I

    def linear_dispersion(w, f=750, phi=-6, m=1, gr=300):
        ''' dlambda / dx
        Default values correspond to Acton 750i

        Input
        -------
        f: focal length (mm)
             default 750 (SpectraPro 2750i)

        phi: angle in degrees (°)
            default 9

        m: order of dispersion
            default 1

        gr: grooves spacing (gr/mm)
            default 300
        '''
        # correct units:
        phi *= 2*pi/360
        d = 1e-3/gr
        disp = w/(2*f)*(tan(phi)+sqrt((2*d/m/(w*1e-9)*cos(phi))**2-1))
        return disp  # to nm/mm

    w_slit, I_slit = import_experimental_slit(getTestFile('slitfunction.txt'))

    if plot:
        plt.figure(fig_prefix+'Linear dispersion effect')
        plt.plot(w_slit, I_slit, '--k', label='Exp: FWHM @{0}nm: {1:.3f} nm'.format(632.8,
                 get_effective_FWHM(w_slit, I_slit)))
        
    # Test how slit function FWHM scales with linear_dispersion
    for w0, FWHM in zip([380, 1000, 4200, 5500],
                        [0.396, 0.388, 0.282, 0.188]):
        w, I = dirac(w0)
        wc, Ic = convolve_with_slit(w, I, w_slit, I_slit, norm_by='area', 
                                    slit_dispersion=linear_dispersion)
        assert np.isclose(FWHM, get_effective_FWHM(wc, Ic), atol=0.001)

        if plot: 
            plt.plot(wc, Ic, label='FWHM @{0:.2f} nm: {1:.3f} nm'.format(w0, 
                     get_effective_FWHM(wc, Ic)))
    
    if plot:
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Dirac $x$ slit function')
        plt.legend(loc='best', prop={'size':15})
        fix_style('article')

    return True # nothing defined yet

def _run_testcases(plot=True, close_plots=False, verbose=True, *args, **kwargs):

    test_against_specair_convolution__fast(plot=plot, close_plots=close_plots, verbose=verbose,
                                           *args, **kwargs)
    test_normalisation_mode__fast(plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs)
    test_slit_energy_conservation__fast(plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs)
    test_slit_function_effect__fast(plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs)
#    test_constant_source(plot=plot, verbose=verbose, *args, **kwargs)
    test_all_slit_shapes__fast(plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs)
    test_slit_unit_conversions_spectrum_in_cm(verbose=verbose, plot=plot, close_plots=close_plots, *args, **kwargs)
    test_slit_unit_conversions_spectrum_in_nm(verbose=verbose, plot=plot, close_plots=close_plots, *args, **kwargs)
#    test_resampling__fast(plot=plot, verbose=verbose, *args, **kwargs)

    return True


if __name__ == '__main__':
    print(('Testing slit.py: ', _run_testcases(plot=True)))
    
