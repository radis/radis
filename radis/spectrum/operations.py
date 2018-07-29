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
from radis.misc.curve import curve_substract, curve_add
from radis.spectrum.compare import get_diff, get_sum
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
def checkImput(s1, s2, Iunit, var):
    if s1.get_conditions()['medium'] != s2.get_conditions()['medium']:
        raise ValueError('Not the same medium for s1 and s2')
    # Check inputs, get defaults
    # ----
    if Iunit == 'default':
        try:
            Iunit = s1.units[var]
        except KeyError:  # unit not defined in dictionary
            raise KeyError('Iunit not defined in spectrum. Cant plot')
    return Iunit
    
def substract(s1, s2, var='radiance', wunit='nm', Iunit='default',
              resample=True, name='default'):
    '''Substract s2 to s1

    Parameters    
    ----------
    s1, s2: Spectrum objects
        The spectra to substract.
    var: str, optional
        Spectral quantity (ex: 'radiance', 'transmittance'...).
    wunit: str, optional 
        Wave unit ('nm', 'cm-1', ...) to use for output. If default, s1 wunit is taken.        
    Iunit: str
        If 'default' s1 unit is used for variable var.
    name: str, optional
        If not given will be s1 and s2 names separated by '-'.
    Returns    
    ----------
    sub : Spectrum object where intensities of s1 and s2 are substracted

    Note    
    ----------
    Wavelength of s1 is taken as reference.
    '''
    Iunit = checkImput(s1, s2, Iunit, var)
    w1, I3 = get_diff(s1, s2, var=var, wunit=wunit, Iunit=Iunit,
                      resample=resample)
    if name != 'default':
        name = s1.get_name()+'-'+s2.get_name()
    sub = Spectrum.from_array(w1, I3, var,
                              waveunit=wunit,
                              unit=Iunit,
                              conditions={
                                  'medium': s1.conditions['medium'], 'waveunit': wunit},
                              name=name)
    return sub

def add(s1, s2, var='radiance', wunit='nm', Iunit='default',
              resample=True, name='default'):
    '''Add s2 to s1

    Parameters    
    ----------
    s1, s2: Spectrum objects
        The spectra to add.
    var: str, optional
        Spectral quantity (ex: 'radiance', 'transmittance'...).
    wunit: str, optional 
        Wave unit ('nm', 'cm-1', ...) to use for output. If default, s1 wunit is taken.        
    Iunit: str
        If 'default' s1 unit is used for variable var.
    name: str, optional
        If not given will be s1 and s2 names separated by '+'.
    Returns    
    ----------
    add : Spectrum object where intensities of s1 and s2 are added

    Note    
    ----------
    Wavelength of s1 is taken as reference.
    '''
    Iunit = checkImput(s1, s2, Iunit, var)
    w1, I3 = get_sum(s1, s2, var=var, wunit=wunit, Iunit=Iunit,
                      resample=resample)
    if name != 'default':
        name = s1.get_name()+'+'+s2.get_name()
    add = Spectrum.from_array(w1, I3, var,
                              waveunit=wunit,
                              unit=Iunit,
                              conditions={
                                  'medium': s1.conditions['medium'], 'waveunit': wunit},
                              name=name)
    return add


def multiply(s, coef, var='radiance', wunit='nm', name='None'):
    '''Multiply the spectrum by the float 'coef'

    Parameters    
    ----------
    s: Spectrum objects
        The spectra to multiply.
    coef: Float
        Coefficient of the multiplication.
    Returns    
    ----------
    sub : Spectrum object where intensity of s is multiplied by coef

    Note    
    ----------
    Good for fittings without absolute calibration. No unit in output !
    '''
    w, I = s.get(var, wunit=wunit)
    mult = Spectrum.from_array(w, coef*I, var,
                              waveunit=wunit,
                              unit='None',
                              conditions={
                                  'medium': s.conditions['medium'], 'waveunit': wunit},
                              name=name)
    return mult

def Nionized(qty):
    s = one_slab(wavelength_min=400,   # nm
            wavelength_max=600,
            mole_fractions={'N II':qty},
            Tgas=30000,            # K
            path_length=0.1,      # cm
            pressure=1,           # bar
            ) 
    return s

if __name__ == '__main__':
    from pyspecair import one_slab
    from radis import plot_diff, get_residual_integral
    s_3 = Nionized(0.3)
    s_6 = Nionized(0.6)

    
    s_3.apply_slit(1)
    s_6.apply_slit(1)
    
    s_3_bis = substract(s_6, s_3)
    s_6_bis = add(s_3, s_3)
    plot_diff(s_6, s_6_bis)
    plot_diff(s_3, s_3_bis)
    
    A = get_residual_integral(s_6, s_3, 'radiance')
    B = get_residual_integral(s_6, s_3_bis, 'radiance')
    
    C = get_residual_integral(s_6_bis, s_3, 'radiance')
    
    print('ratio of A/B ='+str(A/B))
    print('ratio of A/C ='+str(A/C))