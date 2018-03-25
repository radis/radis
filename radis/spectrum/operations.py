#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 23:44:47 2018
Handy functions for manipulation of spectra
@author: minou
"""

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.misc.curve import curve_substract, curve_add
from radis import get_diff, Spectrum

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
        If not given will be s1 and s2 names separated by '+'.
    Returns    
    ----------
    sub : Spectrum object where intensities of s1 and s2 are substracted
    
    Note    
    ----------
    Wavelength of s1 is taken as reference.
    '''
    if s1.get_conditions()['medium'] != s2.get_conditions()['medium']:
        raise ValueError('Not the same medium for s1 and s2')
        # Check inputs, get defaults
    # ----
    if Iunit == 'default':
        try:
            Iunit = s1.units[var]
        except KeyError:  # unit not defined in dictionary
            raise KeyError('Iunit not defined in spectrum. Cant plot') 
            
    w1, I3 = get_diff(s1, s2, var=var, wunit=wunit, Iunit=Iunit, 
                                  resample=resample)
    if name != 'default':
        name = s1.get_name()+'+'+s2.get_name()
    sub = Spectrum.from_array(w1, I3, var, 
                               waveunit=wunit, 
                               unit=Iunit,
                               conditions={'medium' : s1.conditions['medium'], 'waveunit':wunit}, 
                               name=name)
    return sub
