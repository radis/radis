#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Handy functions for manipulation of :class:`~radis.spectrum.spectrum.Spectrum` objects

Routine Listing
---------------

- :func:`~radis.spectrum.operations.substract`

-------------------------------------------------------------------------------


"""

from __future__ import print_function, absolute_import, division, unicode_literals
#from radis.misc.curve import curve_substract, curve_add
from radis.spectrum import Spectrum

# %% Filter Spectra


def Transmittance(s):
    ''' Makes a new Spectrum with only the transmittance part of Spectrum ``s``
    
    Parameters
    ----------
    
    s: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object
        
    Returns
    -------
    
    s_tr: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object, with only the ``transmittance``,
        ``absorbance`` and/or ``abscoeff`` part of ``s``, where ``radiance_noslit`` ,
        ``emisscoeff`` and ``emissivity_noslit`` (if they exist) have been set to 0
        
    Examples
    --------
    
    '''
    
    s_tr = s.copy()
    for k in ['radiance_noslit', 'emisscoeff', 'emissivity_noslit']:
        if k in s_tr._q:
            s_tr._q[k] *= 0
            
    for k in ['radiance', 'emissivity']:
        if k in s_tr._q_conv:
            s_tr._q_conv[k] *= 0

    # Deactivate equilibrium conditions (so Kirchoff's law cannot be used anymore)
    s_tr.conditions['thermal_equilibrium'] = False
    
    s_tr.name = 'Transmittance({0})'.format(s.get_name())
    
    return s_tr
    

# %% Algebric operations on Spectra
def _multiply(s, coef, var='radiance', wunit='nm', name='None'):
    '''Multiply s[var] by the float 'coef'

    Parameters    
    ----------
    s: Spectrum objects
        The spectra to multiply.
    coef: Float
        Coefficient of the multiplication.
    var: str
        'radiance', 'transmittance', ...
    wunit: str
        'nm'or 'cm-1'
    name: str
        name of output spectrum
    Returns    
    ----------
    mult : Spectrum object where intensity of s['var'] is multiplied by coef


    '''
    w, I = s.get(var, wunit=wunit)
    mult = Spectrum.from_array(w, coef*I, var,
                              waveunit=wunit,
                              unit=s.units[var],
                              conditions={
                                  'medium': s.conditions['medium'], 'waveunit': wunit},
                              name=name)
    return mult

def _add_constant(s, cst, var='radiance', wunit='nm', name='None'):
    '''Return a new spectrum with a constant added to s[var] 

    Parameters    
    ----------
    s: Spectrum objects
        Spectrum you want to modify
    cst: Float
        Constant to add.
    var: str
        'radiance', 'transmittance', ...
    wunit: str
        'nm'or 'cm-1'
    name: str
        name of output spectrum
    Returns    
    ----------
    add : Spectrum object where cst is added to intensity of s['var']

    Note    
    ----------
    Use only for rough work. If you want to work properly with spectrum 
    objects, see MergeSlabs.
    '''
    w, I = s.get(var, wunit=wunit)
    add = Spectrum.from_array(w, cst+I, var,
                              waveunit=wunit,
                              unit=s.units[var],
                              conditions={
                                  'medium': s.conditions['medium'], 'waveunit': wunit},
                              name=name)
    return add

def _sub_baseline(s, left, right, var='radiance', wunit='nm', name='None'):
    '''Return a new spectrum with a baseline substracted to s[var] 
    
     Parameters    
    ----------
    s: Spectrum objects
        Spectrum you want to modify
    left: Float
        Constant to substract on the left of the spectrum.
    right: Float
        Constant to substract on the right of the spectrum.
    var: str
        'radiance', 'transmittance', ...
    wunit: str
        'nm'or 'cm-1'
    name: str
        name of output spectrum
    Returns    
    ----------
    output : Spectrum object where the baseline was substracted to intensity of s['var']

    Note    
    ----------
    Use only for rough work. 
    '''
    
    w, I = s.get(var, wunit=wunit)
    I_final = I - np.linspace(left, right, num=np.size(I))
    output = Spectrum.from_array(w, I_final, var,
                              waveunit=wunit,
                              unit=s.units[var],
                              conditions={
                                  'medium': s.conditions['medium'], 'waveunit': wunit},
                              name=name)
    return output
    
def _offset(s, offset, var='radiance', wunit='nm', name='None'):
    w, I = s.get(var, wunit=wunit)
    w_final = w + offset
    output = Spectrum.from_array(w_final, I, var,
                              waveunit=wunit,
                              unit=s.units[var],
                              conditions={
                                  'medium': s.conditions['medium'], 'waveunit': wunit},
                              name=name)
    return output

#def _offset(s, offset, var='radiance', wunit='nm', name='None'):
#     '''Return a new spectrum with a baseline substracted to s[var] 
#    
#     Parameters    
#    ----------
#    s: Spectrum objects
#        Spectrum you want to modify
#    offset: Float
#        Constant to add to the wavelength.
#    right: Float
#        Constant to substract on the right of the spectrum.
#    var: str
#        'radiance', 'transmittance', ...
#    wunit: str
#        'nm'or 'cm-1'
#    name: str
#        name of output spectrum
#    Returns    
#    ----------
#    output : Spectrum object shifted in wavelength
#
#    Note    
#    ----------
#    Use only for rough work. 
#    '''
#    
#    w, I = s.get(var, wunit=wunit)
#    w_final = w + offset
#    output = Spectrum.from_array(w_final, I, var,
#                              waveunit=wunit,
#                              unit=s.units[var],
#                              conditions={
#                                  'medium': s.conditions['medium'], 'waveunit': wunit},
#                              name=name)
#    return output


def _test_multiplyAndAddition(s):
    s_bis = _add_constant(s, 1)
    diff = get_diff(s_bis, s, 'radiance') 
    test = diff[1]-1
    assert np.all(test<1e-10)
    
    s_ter = _multiply(_multiply(s, 50), 1/50)
#    plot_diff(s_ter, s_5)
    diff = get_diff(s_ter, s, 'radiance') 
    ratio = abs(np.trapz(diff[1], x=diff[0])/s.get_integral('radiance'))
    assert ratio<1e-10
    
def _test_visualTestBaseline(s):
    from radis import plot_diff
    s2 = _sub_baseline(s, 2e-4, -2e-4, name = 'sub_arb_baseline')
    plot_diff(s, s2)

def _test_visualTestOffset(s):
    from radis import plot_diff
    s2 = _offset(s, 10, wunit='nm', name = 'sub_arb_baseline')
    plot_diff(s, s2)
    
    
if __name__ == '__main__':
    from radis import load_spec, get_diff, get_residual
    from radis.test.utils import getTestFile
    import numpy as np
#    s_5 = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.5.spec"))
    s_01 = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"))
    s_01.conditions['self_absorption']=False
    s_01.update()
    s_01.apply_slit(0.1)
    _test_multiplyAndAddition(s_01)
#    _test_visualTestBaseline(s_01)
    _test_visualTestOffset(s_01)
    
    
    
    # Added by @erwan
    
    # An implement of Spectrum Algebra
    # Reload:
    s=load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'))
    s.update()
    
    # Test addition
    s.plot(lw=2, nfig='Addition (Merge): s+s')
    (s+s).plot(nfig='same')
    
    # TODO: add test
    # @EP: the test fails at the moment because multiply only works with radiance,
    # and MergeSlabs only works with non convoluted quantities
    # Do we want that? Up for discussion...
    
#    # This should be very small if the spectrum is optically thin (which it is)
#    assert get_residual(2*s, s+s, 'radiance_noslit') < 1e-3

    s.apply_slit(0.1, 'nm')
    # TODO: make 2*s  (multiply)    crash if there are more than 1 spectral quantity?
    
    # Test multiplication with float
    s.plot(lw=2, nfig='Multiplication (Serial): s*s', wunit='nm')
#    (s*s).plot(nfig='same')
    (2*s).plot(nfig='same', wunit='nm')
