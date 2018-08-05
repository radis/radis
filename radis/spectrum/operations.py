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
from radis.spectrum.compare import get_diff
from radis.spectrum import Spectrum
from radis.phys.convert import cm2nm, nm2cm, cm2nm_air, nm_air2cm, air2vacuum, vacuum2air

# %% Filter Spectra


def Transmittance(s:Spectrum): # -> Spectrum:
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
    
    ::
        
        from radis import load_spec, Transmittance
        s = load_spec('file.spec')
        s_tr = Transmittance(s)
        
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
    
    

# %% Change wavelength
    

def crop(s:Spectrum, wmin:float, wmax:float, wunit:str, medium=None,
         inplace=False): # -> Spectrum:
    ''' Crop spectrum to ``wmin-wmax`` range in ``wunit``
    
    Parameters
    ----------
    
    s: Spectrum object
        object to crop
        
    wmin, wmax: float
        boundaries of waverange
        
    wunit: 'nm', 'cm-1'
        which wavespace to use for ``wmin, wmax``
    
    medium: 'air', vacuum'
        necessary if cropping in 'nm'
        
    Other Parameters
    ----------------
    
    inplace: bool
        if ``True``, modifiy ``s`` directly. Else, returns a copy.
    
    Returns
    -------
    
    s_crop: Spectrum
        a cropped Spectrum.
        if using ``inplace``, then ``s_crop`` and ``s`` are still the same object
    
    Examples
    --------
    
    ::
        
        crop(s, 420, 480, 'nm', 'air')
        
    Or in ``cm-1``::
        
        crop(s, 2000, 2300, 'cm-1')
    
    '''
    
    # Check inputs
    assert wunit in ['nm', 'cm-1']
    if wunit == 'nm' and medium is None:
        raise ValueError("Precise give wavelength medium with medium='air' or "+\
                         "medium='vacuum'")
    if wmin >= wmax:
        raise ValueError('wmin should be > wmax')
    
    if len(s._q)>0 and len(s._q_conv)>0:
        raise NotImplementedError('Cant crop this Spectrum as there are both convoluted '+\
                                  'and not convoluted quantities stored')
        # Could bring unexpected errors... For instance, if cropping both 
        # with slit and without slit  quantities to the same waverange, 
        # reapplying the slit in 'valid' mode would reduce the wavelength range 
        # of the convoluted quantities 
        # Implementation: better ask User to drop some of the quantities themselves
    
    if not inplace:
        s = s.copy()
    
    # Convert wmin, wmax to Spectrum wavespace    
    # (deal with cases where wavelength are given in 'air' or 'vacuum')
    waveunit = s.wavespace()
    if wunit == 'nm' and waveunit == 'cm-1':
        if medium == 'air':
            wmin, wmax = nm_air2cm(wmax), nm_air2cm(wmin)   # reverted
        else:
            wmin, wmax = nm2cm(wmax), nm2cm(wmin)   # reverted
    elif wunit == 'cm-1' and waveunit == 'nm':
        if s.get_medium() == 'air':
            wmin, wmax = cm2nm_air(wmax), cm2nm_air(wmin)   # nm in air
        else:
            wmin, wmax = cm2nm(wmax), cm2nm(wmin)   # get nm in vacuum
    elif wunit == 'nm' and waveunit == 'nm':
        if s.get_medium() == 'air' and medium == 'vacuum':
            # convert from given medium ('vacuum') to spectrum medium ('air')
            wmin, wmax = vacuum2air(wmin), vacuum2air(wmax)
        elif s.get_medium() == 'vacuum' and medium == 'air':
            # the other way around
            wmin, wmax = air2vacuum(wmin), air2vacuum(wmax)
    else:
        assert wunit == waveunit           # correct wmin, wmax
    
    # Crop non convoluted
    if len(s._q)>0:
        b = (wmin <= s._q['wavespace']) & (s._q['wavespace'] <= wmax)
        for k, v in s._q.items():
            s._q[k] = v[b]
  
    # Crop convoluted
    if len(s._q_conv)>0:
        b = (wmin <= s._q_conv['wavespace']) & (s._q_conv['wavespace'] <= wmax)
        for k, v in s._q_conv.items():
            s._q_conv[k] = v[b]

    return s
    



# %% Algebric operations on Spectra

def substract(s1:Spectrum, s2:Spectrum, var='radiance', wunit='nm', Iunit='default',
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
    -------
    sub : Spectrum 
        Spectrum object where intensities of s1 and s2 are substracted

    Notes
    -----
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
                              conditions={
                                  'medium': s1.conditions['medium'], 'waveunit': wunit},
                              name=name)
    return sub



def multiply(s:Spectrum, coef:float, var='radiance', wunit='nm', name='None'):
    '''Multiply the spectrum by the float 'coef'

    Parameters    
    ----------
    
    s: Spectrum objects
        The spectra to multiply.
        
    coef: Float
        Coefficient of the multiplication.
        
    Returns    
    -------
    
    sub : Spectrum object where intensity of s is multiplied by coef

    Notes
    -----
    
    Godd for fittings without absolute calibration. No unit in output !
    '''
    w, I = s.get(var, wunit=wunit)
    mult = Spectrum.from_array(w, coef*I, var,
                              waveunit=wunit,
                              unit='None',
                              conditions={
                                  'medium': s.conditions['medium'], 'waveunit': wunit},
                              name=name)
    return mult