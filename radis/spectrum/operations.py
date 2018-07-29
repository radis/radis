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


if __name__ == '__main__':
    from radis import plot_diff, get_residual_integral, load_spec
    from radis.test.utils import getTestFile
    s_5 = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.5.spec"))
    s_01 = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"))
    
    