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
from radis.phys.convert import cm2nm, nm2cm, cm2nm_air, nm_air2cm, air2vacuum, vacuum2air

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
    
def test_invariants():
    from radis import load_spec
    from radis.test.utils import getTestFile
    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"))
    s.update()
    s.apply_slit(0.1, 'nm')
    
    
    # Note @EP: this test doesnt work for the moment as constant & multiply
    # seem to return a Spectrum in a difference medium (vacuum vs air) -> 
    # creates a ~1 nm offset for CO IR 
    
    assert s == _add_constant(s, 0)
    assert s == _multiply(s, 1)
    
    
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
    
    test_invariants()   #not working

    
    # Added by @erwan
    
    # An implement of Spectrum Algebra
    # Reload:
    s=load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'))
    s.update()
    
    # Test addition of Spectra
    s.plot(lw=2, nfig='Addition (Merge): s+s')
    (s//s).plot(nfig='same')
    
    # Test substraction of Spectra
    # TODO : Does not work because of Transmittance @Erwan
    s2 = Transmittance(s)
    s_test = s2-s2
    assert s_test.get_integral('abscoeff') == 0
    
    
    # TODO: add test
    # @EP: the test fails at the moment because multiply only works with radiance,
    # and MergeSlabs only works with non convoluted quantities
    # Do we want that? Up for discussion...
    
#    # This should be very small if the spectrum is optically thin (which it is)
#    assert get_residual(2*s, s+s, 'radiance_noslit') < 1e-3

    s.apply_slit(0.1, 'nm')
    # TODO: make 2*s  (multiply)    crash if there is more than 1 spectral quantity?
    
    # Test multiplication with float
    s.plot(lw=2, nfig='Multiplication (by scalar): 2*s', wunit='nm')
#    (s*s).plot(nfig='same')
    (2*s).plot(nfig='same', wunit='nm')
    
    
    # Test Serial:
    s.rescale_path_length(20)
    s.plot('radiance_noslit', lw=2, nfig='Line of sight (SerialSlabs): s > s > s')
    (s>s>s).plot('radiance_noslit', nfig='same')

