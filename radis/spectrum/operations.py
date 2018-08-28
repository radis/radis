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
from radis.phys.convert import (cm2nm, nm2cm, cm2nm_air, nm_air2cm, air2vacuum, vacuum2air,
                                dcm2dnm, dnm2dcm)

# %% Filter Spectra

def Transmittance(s):
    # type: (Spectrum) -> Spectrum
    ''' Returns a new Spectrum with only the ``transmittance`` component of ``s`` 
    
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
 
    '''

    return s.copy(copy_lines=True, quantity='transmittance')

def Transmittance_noslit(s):
    ''' Returns a new Spectrum with only the ``transmittance_noslit`` component of ``s`` 
    
    Parameters
    ----------
    
    s: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object
        
    Returns
    -------
    
    s_tr: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object, with only 
        ``transmittance_noslit`` defined
 
    '''

    return s.copy(copy_lines=True, quantity='transmittance_noslit')

def Radiance(s):
    ''' Returns a new Spectrum with only the ``radiance`` component of ``s`` 
    
    Parameters
    ----------
    
    s: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object
        
    Returns
    -------
    
    s_tr: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object, with only ``radiance``
        defined
 
    '''

    return s.copy(copy_lines=True, quantity='radiance')

def Radiance_noslit(s):
    ''' Returns a new Spectrum with only the ``radiance_noslit`` component of ``s`` 
    
    Parameters
    ----------
    
    s: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object
        
    Returns
    -------
    
    s_tr: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object, with only ``radiance_noslit``
        defined
 
    '''

    return s.copy(copy_lines=True, quantity='radiance_noslit')

# Useful for determining line-of-sight contributions:

def PerfectAbsorber(s):
    ''' Makes a new Spectrum with the same transmittance/absorbance as Spectrum
    ``s``, but with radiance set to 0. 
    Useful to get contribution of different slabs in line-of-sight 
    calculations)
    
    .. note:
        
        formerly named "Transmittance", but "Transmittance(s)" wouldnt 
        return the Transmittance exactly
    
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
    
    s_tr.name = 'PureAbsorber({0})'.format(s.get_name())
    
    return s_tr
    

# %% Change wavelength
    

def crop(s, wmin, wmax, wunit, medium=None,
         inplace=False):
    # type: (Spectrum, float, float, str, str, bool) -> Spectrum
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
    
def _get_unique_var(s):
    ''' Returns the unique spectral quantity in the Spectrum ``s``. If there are more, 
    raises an error.
    '''
    # TODO: make sure units are okay
    quantities = s.get_vars()
    if len(quantities)>1:
        raise KeyError('There is an ambiguity with the Spectrum algebraic operation. '+\
                        'There should be only one var in Spectrum {0}. Got {1}\n'.format(
                                s.get_name(), s.get_vars())+\
                        "Think about using 'Transmittance(s)', 'Radiance(s)', etc.")
    elif len(quantities)==0:
        raise KeyError('No spectral quantity defined in Spectrum {0}'.format(s.get_name()))
    else:
        return quantities[0]


def _multiply(s, coef, var=None, name='None', inplace=False):
    '''Multiply s[var] by the float 'coef'

    Parameters    
    ----------
    s: Spectrum objects
        The spectra to multiply.
    coef: Float
        Coefficient of the multiplication.
    var: str, or ``None``
        'radiance', 'transmittance', ... If ``None``, get the unique spectral
        quantity of ``s`` or raises an error if there is any ambiguity
    name: str
        name of output spectrum
    inplace: bool
        if ``True``, modifies ``s`` directly. Else, returns a copy. 
        Default ``False``
        
    Returns    
    -------
    s : Spectrum
        Spectrum object where intensity of s['var'] is multiplied by coef
        If ``inplace=True``, ``s`` has been modified directly.
        
    '''
    if var is None:
        var = _get_unique_var(s)
        
    if not inplace:
        s = s.copy(quantity=var)
    
    # Multiply inplace       ( @dev: we have copied already if needed )
    w, I = s.get(var, wunit=s.get_waveunit(), copy=False)  
    I *= coef         # @dev: updates the Spectrum directly because of copy=False

    return s

def _add_constant(s, cst, unit=None, var=None, wunit='nm', name='None', inplace=False):
    '''Return a new spectrum with a constant added to s[var] 

    Parameters    
    ----------
    s: Spectrum objects
        Spectrum you want to modify
    cst: Float
        Constant to add.
    unit: str
        unit for ``cst``. If ``None``, uses the default unit in ``s`` for 
        variable ``var``.
    var: str, or ``None``
        'radiance', 'transmittance', ... If ``None``, get the unique spectral
        quantity of ``s`` or raises an error if there is any ambiguity
    wunit: str
        'nm'or 'cm-1'
    name: str
        name of output spectrum
    inplace: bool
        if ``True``, modifies ``s`` directly. Else, returns a copy. 
        Default ``False``

    Returns    
    -------
    s : Spectrum
        Spectrum object where cst is added to intensity of s['var']
        If ``inplace=True``, ``s`` has been modified directly.

    Notes   
    -----
    Use only for rough work. If you want to work properly with spectrum 
    objects, see MergeSlabs.
    '''
    
    if var is None:
        var = _get_unique_var(s)
        
    # Convert to Spectrum unit
    if unit is not None:
        Iunit = s.units[var]
        if unit != Iunit:
            from radis.phys.convert import conv2
            cst = conv2(cst, unit, Iunit)
        
    if not inplace:
        s = s.copy(quantity=var)
    
    # Add inplace       ( @dev: we have copied already if needed )
    w, I = s.get(var, wunit=s.get_waveunit(), copy=False)  
    I += cst
    # @dev: updates the Spectrum directly because of copy=False

    return s

def _sub_baseline(s, left, right, unit=None, var=None, wunit='nm', name='None', 
                  inplace=False):
    '''Return a new spectrum with a baseline substracted to s[var] 
    
    Parameters    
    ----------
    
    s: Spectrum objects
        Spectrum you want to modify
    left: Float
        Constant to substract on the left of the spectrum.
    right: Float
        Constant to substract on the right of the spectrum.
    unit: str
        unit for ``cst``. If ``None``, uses the default unit in ``s`` for 
        variable ``var``.
    var: str
        'radiance', 'transmittance', ...  If ``None``, get the unique spectral
        quantity of ``s`` or raises an error if there is any ambiguity
    wunit: str
        'nm'or 'cm-1'
    name: str
        name of output spectrum
    inplace: bool
        if ``True``, modifies ``s`` directly. Else, returns a copy. 
        Default ``False``
        
    Returns    
    -------
    
    s: Spectrum
        Spectrum object where the baseline was substracted to intensity of s['var']
        If ``inplace=True``, ``s`` has been modified directly.

    Notes    
    -----
    Use only for rough work. 
    '''
    
    if var is None:
        var = _get_unique_var(s)
        
    # Convert to Spectrum unit
    if unit is not None:
        Iunit = s.units[var]
        if unit != Iunit:
            from radis.phys.convert import conv2
            left = conv2(left, unit, Iunit)
            right = conv2(right, unit, Iunit)

    if not inplace:
        s = s.copy(quantity=var)
       
    # @EP: 

    # Substract inplace       ( @dev: we have copied already if needed )
    w, I = s.get(var, wunit=s.get_waveunit(), copy=False)  
    I -= np.linspace(left, right, num=np.size(I))
    # @dev: updates the Spectrum directly because of copy=False

    return s
    
from warnings import warn

def _substract_spectra(s1, s2, var=None, wunit='nm', name=None):
    '''Return a new spectrum with s2 substracted from s1
    
    Parameters    
    ----------
    
    s1, s2: Spectrum objects
        Spectrum you want to substract
    var: str
        quantity to manipulate: 'radiance', 'transmittance', ... If ``None``, 
        get the unique spectral quantity of ``s1``, or the unique spectral
        quantity of ``s2``, or raises an error if there is any ambiguity
    wunit: str
        'nm'or 'cm-1'
    name: str
        name of output spectrum
        
    Returns    
    -------
    
    s: Spectrum
        Spectrum object with the same units and waveunits as ``s1``
    '''

    from radis import get_diff
    
    # Get variable
    if var is None:
        try:
            var = _get_unique_var(s2)    # unique variable of 2nd spectrum  
        except KeyError:
            var = _get_unique_var(s1)    # if doesnt exist, unique variable of 1st spectrum
            # if it fails, let it fail
    # Make sure it is in both Spectra
    if var not in s1.get_vars():
        raise KeyError('Variable {0} not in Spectrum {1}'.format(var, s1.get_name()))
    if var not in s2.get_vars():
        raise KeyError('Variable {0} not in Spectrum {1}'.format(var, s1.get_name()))
    
    # Use same units 
    Iunit = s1.units[var]            
    wunit = s1.get_waveunit()
    w1, Idiff = get_diff(s1, s2, var=var, wunit=wunit, Iunit=Iunit, resample=True)
    
    name = s1.get_name()+'-'+s2.get_name()
    sub = Spectrum.from_array(w1, Idiff, var, 
                               waveunit=wunit, 
                               unit=Iunit,
                               conditions={'medium' : s1.conditions['medium']}, 
                               name=name)
    warn("Conditions of the left spectrum were copied in the substraction.", Warning)
    return sub

def offset(s, offset, unit, name=None, inplace=False):
    # type: (Spectrum, float, str, str, bool) -> Spectrum
    '''Offset the spectrum by a wavelength or wavenumber 

    Parameters    
    ----------
    s: Spectrum
        Spectrum you want to modify
    offset: float
        Constant to add to all quantities in the Spectrum.
    unit: 'nm' or 'cm-1'
        unit for ``offset``.
    name: str
        name of output spectrum
    inplace: bool
        if ``True``, modifies ``s`` directly. Else, returns a copy. 
        Default ``False``

    Returns    
    -------
    s : Spectrum
        Spectrum object where cst is added to intensity of s['var']
        If ``inplace=True``, ``s`` has been modified directly.
    '''
    
    # Convert to correct unit:
    if unit == 'nm' and s.get_waveunit() == 'cm-1':
        # Note @EP: technically we should check the medium is air or vacuum...  TODO
        # Note @EP: here we're offsetting by a constant value in 'cm-1', which is
        # not a constant value in 'nm'. We end of with offset as an array 
        offset_q = - dnm2dcm(offset, s.get_wavelength(which='non_convoluted'))  # this is an array
        offset_qconv = - dnm2dcm(offset, s.get_wavelength(which='convoluted'))  # this is an array
    elif unit == 'cm-1' and s.get_waveunit() == 'nm':
        offset_q = - dcm2dnm(offset, s.get_wavenumber(which='non_convoluted'))  # this is an array
        offset_qconv = - dcm2dnm(offset, s.get_wavenumber(which='convoluted'))  # this is an array
    else:
        assert unit == s.get_waveunit()
        offset_q = offset
        offset_qconv = offset
        
    if not inplace:
        s = s.copy()
        
    # Update all variables
    if 'wavespace' in s._q:
        s._q['wavespace'] += offset_q
        # @dev: updates the Spectrum directly because of copy=False
    if 'wavespace' in s._q_conv:
        s._q_conv['wavespace'] += offset_qconv
        
    if name:
        s.name = name
        
    return s


# %% Tests

def test_multiplyAndAddition(*args, **kwargs):

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.update('radiance_noslit', verbose=False)
    s.apply_slit(0.1)
    s = Radiance(s)
    assert s.units['radiance'] == 'mW/cm2/sr/nm'

    s_bis = _add_constant(s, 1, 'mW/cm2/sr/nm')
    w, Idiff = get_diff(s_bis, s, 'radiance') 
    test = Idiff[1]-1
    assert np.all(test<1e-10)

    s_ter = _multiply(_multiply(s, 50), 1/50)
#    plot_diff(s_ter, s_5)
    diff = get_diff(s_ter, s, 'radiance') 
    ratio = abs(np.trapz(diff[1], x=diff[0])/s.get_integral('radiance'))
    assert ratio<1e-10

def test_visualTestBaseline(plot=True, *args, **kwargs):

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.update('radiance_noslit', verbose=False)
    s.apply_slit(0.1)
    s = Radiance_noslit(s)
    assert s.units['radiance'] == 'mW/cm2/sr/nm'
    
    s2 = _sub_baseline(s, 2e-4, -2e-4, name = 'sub_arb_baseline')
    if plot:
        plot_diff(s, s2)
    # TODO: add Test on intensity on both sides?

def test_offset(plot=True):
    
    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.update('radiance_noslit', verbose=False)
    s.apply_slit(0.1)
    
    s2 = offset(s, 10, 'nm', name = 'offset_10nm')
    if plot:
        plot_diff(s, s2)
    assert np.allclose(s2.get_wavelength(which='convoluted'), 
                       s.get_wavelength(which='convoluted')+10)
    assert np.allclose(s2.get_wavelength(which='non_convoluted'), 
                       s.get_wavelength(which='non_convoluted')+10)
    
    # Test inplace version
    s.offset(10, 'nm')
    assert np.allclose(s2.get_wavelength(which='convoluted'), 
                       s.get_wavelength(which='convoluted'))
    
def test_invariants():
    ''' Ensures adding 0 or multiplying by 1 does not change the spectra '''
    from radis import load_spec
    from radis.test.utils import getTestFile
    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"))
    s.update()
    s = Radiance_noslit(s)

    assert s.compare_with(_add_constant(s, 0, 'W/cm2/sr/nm'), spectra_only='radiance_noslit')
    assert s.compare_with(_multiply(s, 1), spectra_only='radiance_noslit')
    

if __name__ == '__main__':
    
    from radis import load_spec, get_diff, plot_diff
    from radis.test.utils import getTestFile
    import numpy as np
    import pytest
   
    
    test_invariants()
    test_multiplyAndAddition()
    test_offset()
    
    
    
    # Added by @erwan
    
    # An implement of Spectrum Algebra
    # Reload:
    s=load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'))
    s.update()
        
    # Test addition of Spectra
    s.plot(lw=2, nfig='Merge: s//s')
    (s//s).plot(nfig='same')
    
    # Test substraction of Spectra
    s_tr = Transmittance_noslit(s)
    assert (s_tr-1.0*s_tr).get_integral('transmittance_noslit') == 0

    # TODO: add test
    # @EP: the test fails at the moment because multiply only works with radiance,
    # and MergeSlabs only works with non convoluted quantities
    # Do we want that? Up for discussion...
    
    # There should be an error if algebraic operations are used when 
    # multiple quantities are defined:
    with pytest.raises(KeyError):
        2*s
        
    s.apply_slit(0.1, 'nm')
    s_rad = Radiance(s)
    
    # Test multiplication with float
    s.plot(lw=2, nfig='Multiplication (by scalar): 2*s', wunit='nm')
#    (s*s).plot(nfig='same')
    (2*s_rad).plot(nfig='same', wunit='nm')
    
    # Test Serial:
    s.rescale_path_length(20)
    s.plot('radiance_noslit', lw=2, nfig='Line of sight (SerialSlabs): s > s > s')
    (s>s>s).plot('radiance_noslit', nfig='same')

