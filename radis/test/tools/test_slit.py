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

from neq.spec import SpectrumFactory
from radis.spectrum.spectrum import calculated_spectrum, transmittance_spectrum
from radis.tools.database import load_spec
from radis.los.slabs import SerialSlabs
from radis.tools.slit import (gaussian_slit, triangular_slit, trapezoidal_slit,
                           import_experimental_slit, convolve_with_slit, get_FWHM)
from radis.phys.units import is_homogeneous
from radis.misc.utils import DatabankNotFound
from neq.test.utils import IgnoreMissingDatabase, build_test_databases
import matplotlib.pyplot as plt 
import numpy as np
from numpy import sqrt, linspace, abs, trapz
from os.path import basename

fig_prefix = basename(__file__)+': '

# %% --------------------------------------------------------------------------
#  Test cases
# -----------------------------------------------------------------------------

def test_all_slits__fast(FWHM=2, wstep=0.01, verbose=True, plot=False, *args, **kwargs):
    ''' Test all slit generation functions and make sure we get the expected FWHM'''

    if plot:
        plt.figure(fig_prefix+'all_slits')
        
    w, I = gaussian_slit(FWHM, wstep, calc_range=4)
    if plot:
        plt.plot(w, I, label='Gaussian slit (FWHM={0})'.format(get_FWHM(w, I)))
    assert np.isclose(get_FWHM(w, I), FWHM, atol=1.1*wstep)

    w, I = triangular_slit(FWHM, wstep)
    if plot:
        plt.plot(w, I, label='Triangular slit (FWHM={0})'.format(get_FWHM(w, I)))
    assert np.isclose(get_FWHM(w, I), FWHM, atol=1.1*wstep)

    w, I = trapezoidal_slit(FWHM*0.9, FWHM*1.1, wstep)
    if plot:
        plt.plot(w, I, label='Trapezoidal slit (FWHM={0})'.format(get_FWHM(w, I)))
        plt.legend()
    assert np.isclose(get_FWHM(w, I), FWHM, atol=1.1*wstep)

    if verbose: print('\n>>> _test_all_slits yield correct FWHM (+- wstep) : OK\n')
    
    return True # nothing defined yet

def test_against_specair_convolution__fast(plot=False, verbose=True, debug=False, *args, **kwargs):

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
        print('Integrals should match: {0:.2f} vs {1:.2f} ({2:.2f}% error)'.format(
                As, Au, 100*abs(As-Au)/As))
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

def test_normalisation_mode__fast(plot=False, verbose=True, *args, **kwargs):
    ''' Test norm_by = 'area' vs norm_by = 'max' '''

    from radis.test.utils import getTestFile
    
    # Test
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
        print("radiance unit ({0}) is homogeneous to 'mW/cm2/sr': OK".format(s.units['radiance']))
    
    return True


def test_slit_energy_conservation__fast(verbose=True, plot=False, *args, **kwargs):
    ''' Convoluted and non convoluted quantities should have the same area
    (difference arises from side effects if the initial spectrum is not 0 on 
    the sides '''
    
    from radis.test.utils import getTestFile
    
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
    

def test_slit_function_effect__fast(verbose=True, plot=False, *args, **kwargs):
    ''' A test case to show the effect of wavelength dispersion (cf spectrometer
    reciprocal function) on the slit function '''

    from neq.test.utils import getTestFile
    from publib import fix_style
    from numpy import pi, tan, cos

    # Effect of Slit dispersion

    def dirac(w0, width=20, wstep=0.009):
        ''' Return a Dirac in w0. Plus some space on the side '''

        w = np.arange(w0-width, w0+width+wstep, wstep)
        I = np.zeros_like(w)
        I[len(I)//2] = 1/wstep

        return w, I

    def FWHM(w, I):
        ''' Returns effective FWHM, which is the area when peak is normalized
        to 1 '''
        return np.trapz(I/I.max(), x=w)

    def linear_dispersion(w, f=750, phi=6, m=1, gr=300):
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


    if plot:
        plt.figure(fig_prefix+'Linear dispersion effect')

        w_slit, I_slit = import_experimental_slit(getTestFile('slitfunction.txt'))

        plt.plot(w_slit, I_slit, '--k', label='Exp: FWHM @{0}nm: {1:.3f} nm'.format(632.8,
                 FWHM(w_slit, I_slit)))
        for w0 in [380, 1000, 4200, 5500]:
            w, I = dirac(w0)
            wc, Ic = convolve_with_slit(w, I, w_slit, I_slit, slit_dispersion=linear_dispersion)
            plt.plot(wc, Ic, label='FWHM @{0:.2f} nm: {1:.3f} nm'.format(w0, FWHM(wc, Ic)))
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Dirac $x$ slit function')
        plt.legend(loc='best', prop={'size':15})
        fix_style(str('article'))

    if verbose: print('\n>>> _test_slit_function_effect : NOT DEFINED\n')
    return True # nothing defined yet

def _run_testcases(plot=False, verbose=True, *args, **kwargs):

    test_against_specair_convolution__fast(plot=plot, verbose=verbose,
                                           *args, **kwargs)
    test_normalisation_mode__fast(plot=plot, verbose=verbose, *args, **kwargs)
    test_slit_energy_conservation__fast(plot=plot, verbose=verbose, *args, **kwargs)
    test_slit_function_effect__fast(plot=plot, verbose=verbose, *args, **kwargs)
#    test_constant_source(plot=plot, verbose=verbose, *args, **kwargs)
    test_all_slits__fast(plot=plot, verbose=verbose, *args, **kwargs)
#    test_resampling__fast(plot=plot, verbose=verbose, *args, **kwargs)

    return True


if __name__ == '__main__':
    print('Testing slit.py: ', _run_testcases(plot=True))
    