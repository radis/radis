#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 15:59:06 2017

This module contains:
- a function to convolve with a slit function, that can resample and correct
for slit dilation
- predefined slit functions generators (triangular, gaussian, etc... see SLIT_SHAPES)
and experimental slit function importer
- functions to manipulate slit functions (plot, get_FWHM, get_effective_FWHM,
recenter_slit, crop_slit)


"""

from __future__ import print_function, absolute_import, division, unicode_literals

import matplotlib.pyplot as plt
import numpy as np
from numpy import exp, sqrt, trapz
from numpy import log as ln
from scipy.interpolate import splrep, splev
from warnings import warn
from radis.misc.arrays import evenly_distributed
from radis.misc.basics import is_float
from radis.misc.signal import resample_even
from radis.misc.debug import printdbg
from radis.phys.convert import cm2nm, nm2cm, dnm2dcm, dcm2dnm
from radis.spectrum.spectrum import cast_waveunit
from six import string_types

SLIT_SHAPES = ['triangular', 'trapezoidal', 'gaussian']

# %% Get slit function

def get_slit_function(slit_function, unit='nm', norm_by='area', shape='triangular',
                      center_wavespace=None, return_unit='same', wstep=None,
                       plot=False, resfactor=2,
                       *args, **kwargs):
    ''' Import or generate slit function in correct wavespace
    Give a file path to import, or a float / tuple to generate arbitrary shapes

    Warning with units: read about unit and return_unit parameters.

    See :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit` and 
    :func:`~radis.tools.slit.convolve_with_slit` for more info

    
    Parameters    
    ----------

    slit_function: tuple, or str
        If float:
            generate slit function with FWHM of `slit_function` (in `unit`)
        If .txt file:
            import experimental slit function (in `unit`): format must be 2-columns
            with wavelengths and intensity (doesn't have to be normalized)

    unit: 'nm' or 'cm-1'
        unit of slit_function FWHM, or unit of imported file

    norm_by: 'area', 'max', or None
        how to normalize. `area` conserves energy. With `max` the slit is normalized
        so that its maximum is one (that is what is done in Specair: it changes
        the outptut spectrum unit, e.g. from 'mW/cm2/sr/µm' to 'mW/cm2/sr')
        None doesnt normalize. Default 'area'

    shape: 'triangular', 'trapezoidal', 'gaussian'
        which shape to use when generating a slit. Default 'triangular'

    center_wavespace: float, or None
        center of slit when generated (in unit). Not used if slit is imported.

    return_unit: 'nm', 'cm-1', or 'same'
        if not 'same', slit is converted to the given wavespace.

    wstep: float
        which discretization step to use (in return_unit) when generating a slit
        function. Not used if importing


    Other Parameters
    ----------------

    resfactor: int
        resolution increase when resampling from nm to cm-1, or the other way
        round. Default 2.

    energy_threshold: float
         tolerance fraction. Only used when importing experimental slit as the
         theoretical slit functions are directly generated in spectrum wavespace
         Default 1e-3 (0.1%)
         

    Returns
    -------
    
    wslit, Islit: array
        wslit is in `return_unit` . Islit is normalized according to  `norm_by`
         
        
    Examples
    --------
    
    >>> wslit, Islit = get_slit_function(1, 'nm', shape='triangular', 
    center_wavespace=600, wstep=0.01)     
    
    Returns a triangular slit function of FWHM = 1 nm, centered on 600 nm, with
    a step of 0.01 nm
    
    >>> wslit, Islit = get_slit_function(1, 'nm', shape='triangular', 
    center_wavespace=600, return_unit='cm-1', wstep=0.01)   
    
    Returns a triangular slit function expressed in cm-1, with a FWHM = 1 nm 
    (converted in equivalent width in cm-1 at 600 nm), centered on 600 nm, with 
    a step of 0.01 cm-1 (!)   
    

    Notes
    -----

    In norm_by 'max' mode, slit is normalized by slit max. In RADIS, this is done
    in the spectrum wavespace (to avoid errors that would be caused by interpolating
    the spectrum). 
    
    A problem arise if the spectrum wavespace is different from the slit wavespace:
    typically, slit is in 'nm' but a spectrum calculated by RADIS is stored in 
    'cm-1': in that case, the convoluted spectrum is multiplied by /int(Islit*dν)
    instead of /int(Islit*dλ). The output unit is then [radiance]*[spec_unit]
    instead of [radiance]*[slit_unit], i.e, typically, [mW/cm2/sr/nm]*[cm-1]
    instead of [mW/cm2/sr/nm]*[nm]=[mW/cm2/sr]
    
    While this remains true if the units are taken into account, this is not the
    expected behaviour. For instance, Specair users are used to having a FWHM(nm) 
    factor between spectra convolved with slit normalized by max and slit normalized
    by area.
    
    The norm_by='max' behaviour adds a correction factor
    `/int(Islit*dλ)/int(Islit*dν)` to maintain an output spectrum in [radiance]*[slit_unit]
    
    See Also
    --------
    
    :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`, 
    :func:`~radis.tools.slit.convolve_with_slit`
    
    '''

    if 'waveunit' in kwargs:
        assert return_unit == 'same'  # default
        return_unit = kwargs.pop('waveunit')
        warn(DeprecationWarning('waveunit renamed return_unit'))
    if 'slit_unit' in kwargs:
        assert unit == 'nm' # default
        unit = kwargs.pop('slit_unit')
        warn(DeprecationWarning('slit_unit renamed unit'))

    energy_threshold = kwargs.pop('energy_threshold', 1e-3)   # type: float
                    # tolerance fraction
                    # when resampling (only used in experimental slit as the)
                    # theoretical slit functions are directly generated in
                    # spectrum wavespace

    def check_input_gen():
        if center_wavespace is None:
            raise ValueError('center_wavespace has to be given when generating '+\
                             'slit function')
        if wstep is None:
            raise ValueError('wstep has to be given when generating '+\
                             'slit function')

    # Cast units
    if return_unit == 'same':
        return_unit = unit
    unit = cast_waveunit(unit)
    return_unit = cast_waveunit(return_unit)
    scale_slit = 1   # in norm_by=max mode, used to keep units in [Iunit]*return_unit in [Iunit]*unit
                     # not used in norm_by=area mode

    # First get the slit in return_unit space
    if is_float(slit_function):  # Generate slit function (directly in return_unit space)

        check_input_gen()

        # ... first get FWHM in return_unit  (it is in `unit` right now)
        FWHM = slit_function
        if return_unit == 'cm-1' and unit == 'nm':
            # center_wavespace ~ nm, FWHM ~ nm
            FWHM = dnm2dcm(FWHM, center_wavespace)   # wavelength > wavenumber
            center_wavespace = nm2cm(center_wavespace)
            if norm_by=='max': scale_slit = slit_function/FWHM         # [unit/return_unit]
        elif return_unit == 'nm' and unit == 'cm-1':
            # center_wavespace ~ cm-1, FWHM ~ cm-1
            FWHM = dcm2dnm(FWHM, center_wavespace)  # wavenumber > wavelength
            center_wavespace = cm2nm(center_wavespace)
            if norm_by=='max': scale_slit = slit_function/FWHM         # [unit/return_unit]
        else:
            pass # correct unit already
        # Now FWHM is in 'return_unit'

        # ... now, build it (in our wavespace)
        if __debug__: printdbg('get_slit_function: {0} FWHM {1:.2f}{2}, center {3:.2f}{2}, norm_by {4}'.format(
                shape, FWHM, return_unit, center_wavespace, norm_by))

        if shape == 'triangular':
            wslit, Islit = triangular_slit(FWHM, wstep,
                                           center=center_wavespace, bplot=plot,
                                           norm_by = norm_by,
                                           waveunit=return_unit,
                                           scale=scale_slit,
                                           *args, **kwargs)

        # Insert other slit shapes here
        # ...

        elif shape == 'gaussian':
            wslit, Islit = gaussian_slit(FWHM, wstep,
                                         center=center_wavespace, bplot=plot,
                                         norm_by = norm_by,
                                         waveunit=return_unit,
                                         scale=scale_slit,
                                         *args, **kwargs)

        elif shape == 'trapezoidal':
            raise TypeError('A (top, base) tuple must be given with a trapezoidal slit')

        else:
            raise TypeError('Slit function ({0}) not in known slit shapes: {1}'.format(
                    shape, SLIT_SHAPES))

    elif isinstance(slit_function, tuple):

        check_input_gen()

        try:
            top, base = slit_function
        except:
            raise TypeError('Wrong format for slit function: {0}'.format(slit_function))
        if shape == 'trapezoidal':
            pass
        elif shape == 'triangular': # it's the default
            warn('Triangular slit given with a tuple: we used trapezoidal slit instead')
            shape = 'trapezoidal'
        else:
            raise TypeError('A (top, base) tuple must be used with a trapezoidal slit')

        # ... first get FWHM in our wavespace unit
        if return_unit == 'cm-1' and unit == 'nm':
            # center_wavespace ~ nm, FWHM ~ nm
            top = dnm2dcm(top, center_wavespace)   # wavelength > wavenumber
            base = dnm2dcm(base, center_wavespace)   # wavelength > wavenumber
            center_wavespace = nm2cm(center_wavespace)
            if norm_by=='max': scale_slit = sum(slit_function)/(top+base)  # [unit/return_unit]
        elif return_unit == 'nm' and unit == 'cm-1':
            # center_wavespace ~ cm-1, FWHM ~ cm-1
            top = dcm2dnm(top, center_wavespace)  # wavenumber > wavelength
            base = dcm2dnm(base, center_wavespace)  # wavenumber > wavelength
            center_wavespace = cm2nm(center_wavespace)
            if norm_by=='max': scale_slit = sum(slit_function)/(top+base)  # [unit/return_unit]
        else:
            pass # correct unit already
            
        FWHM = (top+base)/2

        # ... now, build it (in our wavespace)
        if __debug__: printdbg('get_slit_function: {0}, FWHM {1:.2f}{2}, center {3:.2f}{2}, norm_by {4}'.format(
                shape, FWHM, return_unit, center_wavespace, norm_by))

        wslit, Islit = trapezoidal_slit(top, base, wstep,
                                     center=center_wavespace, bplot=plot,
                                     norm_by = norm_by,
                                     waveunit=return_unit,
                                     scale=scale_slit,
                                     *args, **kwargs)

    elif isinstance(slit_function, string_types): # import it
        if __debug__: printdbg('get_slit_function: {0} in {1}, norm_by {2}, return in {3}'.format(
                slit_function, unit, norm_by, return_unit))
        wslit, Islit = import_experimental_slit(slit_function,norm_by=norm_by, # norm is done later anyway
                                                waveunit=unit,
                                                bplot=False, # we will plot after resampling
                                                *args, **kwargs)
        # ... get unit 
        # Normalize
        if norm_by == 'area': # normalize by the area
    #        I_slit /= np.trapz(I_slit, x=w_slit)
            Iunit = '1/{0}'.format(unit)
        elif norm_by == 'max': # set maximum to 1
            Iunit = '1'
        elif norm_by is None:
            Iunit=None
        else:
            raise ValueError('Unknown normalization type: `norm_by` = {0}'.format(norm_by))

        # ... check it looks correct
        unq, counts = np.unique(wslit, return_counts=True)
        dup = counts > 1
        if dup.sum() > 0:
            raise ValueError('Not all wavespace points are unique: slit function '+
                 'format may be wrong. Duplicates for w={0}'.format(unq[dup]))

        # ... resample if needed
        if return_unit == 'cm-1' and unit == 'nm': # wavelength > wavenumber
            wold, Iold = wslit, Islit
            wslit, Islit = resample_even(nm2cm(wslit), Islit, resfactor=resfactor,
                                    energy_threshold=energy_threshold,
                                    print_conservation=True)
            scale_slit=trapz(Iold, wold)/trapz(Islit, wslit)    # [unit/return_unit]
            renormalize = True
        elif return_unit == 'nm' and unit == 'cm-1': # wavenumber > wavelength
            wold, Iold = wslit, Islit
            wslit, Islit = resample_even(cm2nm(wslit), Islit, resfactor=resfactor,
                                    energy_threshold=energy_threshold,
                                    print_conservation=True)
            scale_slit=trapz(Iold, wold)/trapz(Islit, wslit)    # [unit/return_unit]
            renormalize = True
        else:  # return_unit == unit
            renormalize = False
        # Note: if wstep dont match with quantity it's alright as it gets
        # interpolated in the `convolve_with_slit` function

        # re-Normalize if needed (after changing units)
        if renormalize:
            if __debug__: printdbg('get_slit_function: renormalize')
            if norm_by == 'area': # normalize by the area
                Islit /= abs(np.trapz(Islit, x=wslit))
                Iunit = '1/{0}'.format(return_unit)
            elif norm_by == 'max': # set maximum to 1
                Islit /= abs(np.max(Islit))
                Islit *= scale_slit
                Iunit = '1'
                if scale_slit != 1:
                    Iunit += 'x{0}'.format(scale_slit)
#            elif norm_by == 'max2': # set maximum to 1    # removed this mode for simplification
#                Islit /= abs(np.max(Islit))
            elif norm_by is None:
                Iunit=None
            else:
                raise ValueError('Unknown normalization type: `norm_by` = {0}'.format(norm_by))

        if plot: # (plot after resampling / renormalizing)
            # Plot slit
            plot_slit(wslit, Islit, waveunit=return_unit, Iunit=Iunit)

    else:
        raise TypeError('Unexpected type for slit function: {0}'.format(type(slit_function)))

    return wslit, Islit

try:  # Python >3.6 only
    get_slit_function.__annotations__['shape'] = SLIT_SHAPES
except AttributeError:
    pass  # old Python version

# %% Convolve with slit function

def convolve_with_slit(w, I, w_slit, I_slit, norm_by = 'area', 
                       mode='valid', slit_dispersion=None,
                       k=1, bplot=False, conv_type=None, verbose=True):
    ''' Convolves spectrum (w,I) with instrumental slit function (w_slit, I_slit)
    Returns a convolved spectrum on a valid range.

    
    Parameters    
    ----------

    w, I: array
        theoretical spectrum (wavespace unit (nm / cm-1))

    w_slit, I_slit: array
        instrumental slit function  (wavespace unit (nm / cm-1))
        Both wavespaces have to be the same!

    norm_by: 'area', 'max', or None
        how to normalize. `area` conserves energy. `max` is what is done in
        Specair and changes spectrum units, e.g. from 'mW/cm2/sr/µm' to 'mW/cm2/sr'
        None doesnt normalize. Default 'area'

    mode: 'valid', 'same'
        'same' returns output of same length as initial spectra, 
        but boundary effects are still visible. 'valid' returns 
        output of length len(spectra) - len(slit) + 1, for 
        which lines outside of the calculated range have
        no impact. Default 'valid'. 

    slit_dispersion: func of (lambda), or None
        spectrometer reciprocal function : dλ/dx(λ)
        If not None, then the slit_dispersion function is used to correct the
        slit function for the whole range. Can be important if slit function
        was measured far from the measured spectrum  (e.g: a slit function
        measured at 632.8 nm will look broader at 350 nm because the spectrometer
        dispersion is higher at 350 nm. Therefore it should be corrected)
        Default None

        Warning: slit dispersion is not unit aware: if your spectrum is stored
        in cm-1 the slit function is converted in cm-1 but the slit dispersion
        is not changed, so that may result in errors
        # TODO. If slit dispersion first force slit function to be given in nm ?
        # Else it's not relevant

        a Python implementation:

        >>> def f(lbd):
        >>>    return  w/(2*f)*(tan(Φ)+sqrt((2*d/m/(w*1e-9)*cos(Φ))^2-1))

        Theoretical / References:

        >>> dλ/dx ~ d/mf    # at first order
        >>> dλ/dx = w/(2*f)*(tan(Φ)+sqrt((2*d/m/(w)*cos(Φ))^2-1))  # cf

        with:

        - Φ: spectrometer angle (°)
        - f: focal length (mm)
        - m: order of dispersion
        - d: grooves spacing (mm)   = 1/gr  with gr in (gr/mm)

        See Laux 1999 "Experimental study and modeling of infrared air plasma
        radiation" for more information

        slit_dispersion is assumed to be constant on the whole measured range,
        and corrected for the center wavelength. If there is an error >1% on
        the whole range a warning is raised.

    k: int
        order of spline interpolation. 3: cubic, 1: linear. Default 1.

    bplot: boolean
        if True, plots the convolve slit function (for debugging)

    conv_type: (deprecated)
        former norm_by name

    verbose: boolean
        blabla


    Returns
    -------

    w_conv, I_conv: arrays
        convoluted quantity. w_conv is defined on a valid range 
        (sides are cut) if mode='valid'

    Notes
    -----
    
    Implementation is done in 7 steps:
    - Check input
    - Correct for Slit dispersion
    - Interpolate the slit function on the spectrum grid, resample it if not
      evenly spaced
    - Normalize
    - Check slit aspect, plot slit if asked for
    - Convolve!
    - Remove boundary effects

    See Also
    --------
    
    :func:`~radis.tools.slit.get_slit_function`, 
    :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`

    '''

    # 1. Check input
    # --------------

    if conv_type is not None:
        raise ValueError('`conv_type` deprecated. Use `norm_by`')

    # Assert slit function is thin enough
    try:
        assert(abs(w[-1]-w[0])>abs(w_slit[-1]-w_slit[0]))
    except AssertionError:
        raise AssertionError('Slit function is broader ({0:.1f}nm) than spectrum'.format(
                            abs(w_slit[-1]-w_slit[0])) + \
                            ' ({0:.1f}nm). No valid range.'.format(abs(w[-1]-w[0])))


    # 2. Correct for Slit dispersion
    # --------------

    if slit_dispersion is not None:
        w0 = w[len(w)//2]
        wslit0 = w_slit[len(w_slit)//2]

        # Check that slit dispersion is about constant (<1% change) on the calculated range
        threshold = 0.01
        if not 1-threshold < slit_dispersion(w.max())/slit_dispersion(w.min()) < 1+threshold:
            warn('Slit dispersion changes slightly ({2:.2f}%) between {0:.3f} and {1:.3f}nm'.format(
                    w.min(), w.max(), abs(slit_dispersion(w.max())/slit_dispersion(w.min())-1
                          )*100)+'. Consider splitting your spectrum')

        # Offset slit and correct for dispersion
        w_slit = w0 + slit_dispersion(w0)/slit_dispersion(wslit0)*(w_slit-wslit0)


    # 3. Interpolate the slit function on the spectrum grid, resample it if not
    #    evenly spaced
    # --------------

    # ... Resample if not evenly spaced
    # TODO: add a criteria based based on FWHM rather than absolute?
    wstep = abs(np.diff(w)).min()   # spectrum wavelength spacing
    if not evenly_distributed(w, tolerance=wstep*1e-3):
        # TODO: automatically find a resampling factor?
        warn('Spectrum not evenly spaced. Resampling')
        w, I = resample_even(w, I, resfactor=2, print_conservation=True)
        wstep = abs(np.diff(w)).min()   # new spectrum wavelength spacing

    # ... Check reversed (interpolation requires objects are sorted)
    reverse = (w_slit[-1] < w_slit[0])
    if reverse:
        w_slit = w_slit[::-1]
        I_slit = I_slit[::-1]

    if not np.allclose(np.diff(w_slit),wstep):
        if verbose: print('interpolating slit function over spectrum grid') # numerical errors can be produced
        try:
            tck = splrep(w_slit, I_slit, k=k)
        except ValueError:
            # Probably error on input data. Print it before crashing.
            print('\nValueError - Input data below:')
            print('-'*5)
            print(w_slit)
            print(I_slit)
            print('Check figure')
            plot_slit(w_slit, I_slit, waveunit='')
            raise

        w_slit_int = np.arange(w_slit[0],w_slit[-1]+wstep,wstep)  # can be a bug here if wstep has the wrong sign.
                                                              # more test needed.
        I_slit_int = splev(w_slit_int, tck)
    else:
        w_slit_int = w_slit
        I_slit_int = I_slit

    # 4. Normalize
    # --------------

    if norm_by == 'area': # normalize by the area
        I_slit_int /= abs(np.trapz(I_slit_int, x=w_slit_int))
    elif norm_by == 'max': # set maximum to 1
        I_slit_int /= abs(np.max(I_slit_int))
    elif norm_by is None:
        pass
    else:
        raise ValueError('Unknown normalization type: `norm_by` = {0}'.format(norm_by))

    # 5. Check aspect
    # -------------

    # Check no nan
    if np.isnan(I_slit).sum() > 0:
        raise ValueError('Slit has nan value')

    # check slit is positive
    if not (I_slit>=0).all():
        plot_slit(w_slit, I_slit, waveunit='')
        raise ValueError('Slit is partially negative. Check Figure')

    # Plot slit if asked for
    if bplot:
        plot_slit(w_slit, I_slit, waveunit='')

    # 7. Convolve!
    # --------------

    I_conv = np.convolve(I, I_slit_int, mode='same')*wstep

    # 6. Remove boundary effects
    # --------------

    # Remove boundary effects with the x-axis changed accordingly
    if mode == 'valid':
        la = min(len(I), len(I_slit_int))
        a = int((la-1)/2)
        b = int((la)/2)
        I_conv = I_conv[a:-b]
        w_conv = w[a:-b]
    elif mode == 'same':
        I_conv = I_conv
        w_conv = w
    else:
        raise ValueError('Unexpected mode: {0}'.format(mode))

    # reverse back if needed
    # Todo: add test case for that
    if reverse:
        w_slit = w_slit[::-1]
        I_slit = I_slit[::-1]

    return w_conv, I_conv

#%% Slit function methods

def get_FWHM(w, I, return_index=False):
    ''' Calculate full width half maximum by comparing amplitudes

    
    Parameters    
    ----------

    w, I: arrays

    return_index: boolean
        if True, returns indexes for half width boundaries. Default False


    Returns
    -------

    FWHM: float

    [xmin, xmax]: int


    Todo
    ----

    Linearly interpolate at the boundary? insignificant for large number of points

    '''

    upper = np.argwhere(I>=I.max()/2)

    xmin = upper.min()
    xmax = upper.max()

    if return_index:
        return abs(w[xmax] - w[xmin]), xmin, xmax
    else:
        return abs(w[xmax] - w[xmin])

def get_effective_FWHM(w, I):
    ''' Calculate FWHM of a triangular slit of same area and height 1

    
    Parameters    
    ----------

    w, I: arrays
    '''

    Imax = I.max()

    area = abs(np.trapz(I, w))

    return area/Imax

def plot_slit(w, I=None, waveunit='', plot_unit='same', Iunit=None, warnings=True):
    ''' Plot slit, calculate and display FWHM, and calculate effective FWHM.
    FWHM is calculated from the limits of the range above the half width,
    while FWHM is the equivalent width of a triangular slit with the same area

    
    Parameters    
    ----------

    w, I: arrays    or   (str, None)
        if str, open file directly

    waveunit: 'nm', 'cm-1' or ''
        unit of input w

    plot_unit: 'nm, 'cm-1' or 'same'
        change plot unit (and FWHM units)
        
    Iunit: str, or None
        give Iunit

    warnings: boolean
        if True, test if slit is correctly centered and output a warning if it
        is not. Default True

    '''

    try:
        from neq.plot.toolbar import add_tools     # TODO: move in publib
        add_tools()       # includes a Ruler to measure slit
    except:
        pass

    # Check input
    if isinstance(w, string_types) and I is None:
        w, I = np.loadtxt(w).T
    assert len(w)==len(I)
    if np.isnan(I).sum()>0:
        warn('Slit function has nans')
        w = w[~np.isnan(I)]
        I = I[~np.isnan(I)]
    assert len(I)>0

    # cast units
    waveunit = cast_waveunit(waveunit, force_match=False)
    plot_unit = cast_waveunit(plot_unit, force_match=False)
    if plot_unit == 'same':
        plot_unit = waveunit

    # Convert wavespace unit if needed
    elif waveunit == 'cm-1' and plot_unit == 'nm': # wavelength > wavenumber
        w = cm2nm(w)
        waveunit = 'nm'
    elif waveunit == 'nm' and plot_unit == 'cm-1': # wavenumber > wavelength
        w = nm2cm(w)
        waveunit = 'cm-1'
    else: # same units
        pass
#        raise ValueError('Unknown plot unit: {0}'.format(plot_unit))

    # Recalculate FWHM
    FWHM, xmin, xmax = get_FWHM(w, I, return_index=True)
    FWHM_eff = get_effective_FWHM(w, I)

    # Get labels
    if plot_unit == 'nm':
        xlabel = 'Wavelength (nm)'
    elif plot_unit == 'cm-1':
        xlabel = 'Wavenumber (cm-1)'
    elif plot_unit == '':
        xlabel = 'Wavespace'
    else:
        raise ValueError('Unknown unit for plot_unit: {0}'.format(plot_unit))
    ylabel = 'Slit function'
    if Iunit is not None:
        ylabel += ' ({0})'.format(Iunit)

    fig, ax = plt.subplots()
    ax.plot(w, I, 'o', color='lightgrey')
    ax.plot(w, I, '-k', label='FWHM: {0:.3f} {1}'.format(FWHM, plot_unit)+\
                             '\nEff. FWHM: {0:.3f} {1}'.format(FWHM_eff, plot_unit)+\
                             '\nArea: {0:.3f}'.format(abs(np.trapz(I, x=w))),
                             )

    # Vertical lines on center, and FWHM
    plt.axvline(w[len(w)//2], ls='-', lw=2, color='lightgrey')  # center
    plt.axvline(w[(xmin+xmax)//2], ls='--', color='k', lw=0.5)   # maximum (should be center)
    plt.axvline(w[xmin], ls='--', color='k', lw=0.5)      # FWHM min
    plt.axvline(w[xmax], ls='--', color='k', lw=0.5)      # FWHM max
    plt.axhline(I.max()/2, ls='--', color='k', lw=0.5)      # half maximum

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.legend(loc='best', prop={'size':16})

    # extend axis:
    fig.tight_layout()
    xlmin, xlmax= ax.get_xlim()
    ax.set_xlim((xlmin-0.5, xlmax+0.5))

    if warnings:
        if w[(xmin+xmax)//2] != w[len(w)//2]:
            warn('Slit function doesnt seem centered (center measured with FWHM)'+\
                  ' is not the array center. This can induce offsets!')

        if I[0] != 0 or I[-1] != 0:
            warn('Slit function should have zeros on both sides')

    return fig, ax

def recenter_slit(w_slit, I_slit, verbose=True):
    ''' Recenters the slit on the maximum calculated from the two
    FWHM limits. To recenter, zeros are added on the shorter side (zero padding)
    '''

    _, xmin, xmax = get_FWHM(w_slit, I_slit, return_index=True)
    xcenter = (xmax + xmin)/2
    offset = len(w_slit)/2 - xcenter
    add_zeros = int(2*offset)
    if add_zeros > 0:
        # add zeros on the left (assume regular step)
        wstep = w_slit[1] - w_slit[0]
        w_left = np.linspace(w_slit[0] - wstep*(add_zeros+1), w_slit[0]-wstep,
                             add_zeros)
        w_slit = np.hstack((w_left, w_slit))
        I_slit = np.hstack((np.zeros(add_zeros), I_slit))
        if verbose: print('... added zeros to the slit to center it')
    elif add_zeros < 0:
        add_zeros = -add_zeros
        # add zeros on the right (assume regular step)
        wstep = w_slit[-1] - w_slit[-2]
        w_right = np.linspace(w_slit[-1]+wstep, w_slit[-1] + wstep*(add_zeros+1),
                              add_zeros)
        w_slit = np.hstack((w_slit, w_right))
        I_slit = np.hstack((I_slit, np.zeros(add_zeros)))
        if verbose: print('... added zeros to the slit to center it')
    else:
        pass

    return w_slit, I_slit

def crop_slit(w_slit, I_slit, verbose=True):
    ''' Removes unnecessary zeros on the side for a faster convolution.
    (remove as many zeros on the left as on the right).
    '''

    nzeros_index = np.argwhere(I_slit!=0)
    zeros_left = nzeros_index.min()
    zeros_right = len(I_slit) - nzeros_index.max() - 1
    remove = 0
    if zeros_left > 1 and zeros_right > zeros_left:
        remove = zeros_left-1
    elif zeros_right > 1 and zeros_left > zeros_right:
        remove = zeros_right-1
    if remove > 0:
        w_slit = w_slit[remove:-remove]
        I_slit = I_slit[remove:-remove]
        if verbose: print('... removed {0} zeros to the slit on each side'.format(remove))

    return w_slit, I_slit


#%% Slit function models


def import_experimental_slit(fname, norm_by='area', bplot=False,
                             waveunit='nm', auto_recenter=True, auto_crop=True,
                             center=None, scale=1, verbose=True):
    ''' Import instrumental slit function and normalize it

    
    Parameters    
    ----------

    fname: str
        slit function spectrum

    norm_by: None, 'area', 'max'
        normalisation type. 'area' conserves energy. 'max' is what is
        done in Specair and changes units. None doesnt normalize
        experimental slit. Default 'area'

    bplot: boolean
        plot normalized slit function (for debugging). Default False

    waveunit: 'nm', 'cm-1'
        used for plot only. Slit function is generated assuming you use the
        correct wavespace. No conversions are made here. Default 'nm'

    auto_recenter: boolean
        if True, recenters the slit on the maximum calculated from the two
        FWHM limits. To recenter, zeros are added on the shorter side (zero padding)
        Default True

    auto_crop: boolean
        If True, remove unnecessary zeros on the side for a faster convolution.
        (remove as many zeros on the left as on the right). Default True

    center: float, None
        normally the slit instrumental slit function is centered on where it
        was measured. If not None, center overwrites the center and
        offsets the slit to the given value. Default None

    scale: float
        multiply slit by an arbitrary factor. Default 1.

    verbose: boolean
        Display messages

    '''

    w_slit, I_slit = np.loadtxt(fname).T

    # check sides are zeros
    if (not I_slit[0] == 0 and I_slit[-1] == 0):
        raise ValueError('Slit function must be null on each side. Fix it')

    # recenter if asked for
    if auto_recenter:
        w_slit, I_slit = recenter_slit(w_slit, I_slit, verbose=verbose)
         # note: if auto_crop is true we may be adding zeros just to remove them
         # right way. Fix that. Someday.

    # remove unecessary zeros
    if auto_crop:   # (note that we know I_slit has zeros on each side already)
        w_slit, I_slit = crop_slit(w_slit, I_slit, verbose=verbose)

    # Offset slit if needed
    if center is not None:
        w_slit += center - w_slit[len(w_slit)//2]

    # Normalize
    if norm_by == 'area': # normalize by the area
#        I_slit /= np.trapz(I_slit, x=w_slit)
        I_slit /= abs(np.trapz(I_slit, x=w_slit))
        Iunit = '1/{0}'.format(waveunit)
    elif norm_by == 'max': # set maximum to 1
        I_slit /= abs(np.max(I_slit))
        Iunit = '1'
    elif norm_by is None:
        Iunit=None
    else:
        raise ValueError('Unknown normalization type: `norm_by` = {0}'.format(norm_by))

    # scale
    I_slit *= scale
#    if Iunit is not None and scale != 1:
#        Iunit += 'x{0:.2f}'.format(scale)

    # Plot slit
    if bplot:
        plot_slit(w_slit, I_slit, waveunit=waveunit, Iunit=Iunit)

    return w_slit, I_slit


def triangular_slit(FWHM, wstep, center=0, norm_by='area', bplot=False,
                    waveunit='', scale=1, footerspacing=0):
    r''' Generate normalized slit function

    
    Parameters    
    ----------

    FWHM: (nm)
        full-width at half maximum

    wstep: (nm)
        wavelength step

    center: (nm)
        center wavelength for the wavelength axs of the slit function

    norm_by: 'area', 'max'
        normalisation type. 'area' conserves energy. 'max' is what is
        done in Specair and changes units. Default 'area'

    bplot: boolean
        plot normalized slit function (for debugging). Default False

    wavespace: '', 'nm', 'cm-1'
        used for plot only. Slit function is generated assuming you use the
        correct wavespace. No conversions are made here.

    scale: float
        multiply slit by an arbitrary factor. Default 1.

    footerspacing: int
        spacing (footer) on left and right. Default 10.


    Returns
    -------

    Slit function of shape
    
               ^
              / \
         _ _ /   \ _ _
        :   :  :  :   :
       -a-f -a 0  a  a+f
       

    with FWHM = (a+1/2)*wstep, `f` spacing on left & right


    Notes
    -----
    
    slit is generated with an odd number of elements and centered on its
    maximum. This may result in sightly different FWHM if wstep is too large.
    However, slit is corrected so that effective FWHM is preserved.

    '''

    # Build first half
    slope = 1/(FWHM-wstep)
    a = int(FWHM/wstep)
    I = 1-np.arange(0, a)*slope*wstep    # slope

    # Zeros
    if FWHM % wstep:  # add one extra zero when not a divider
        footerspacing += 1
    f = int(footerspacing)   # add zeros (footer)
    I = np.hstack((I, np.zeros(f)))
    w = np.linspace(0, len(I)*wstep, len(I))

    # Mirror to get second half
    I = np.hstack((I[1:][::-1], I))
    w = np.hstack((-w[1:][::-1], w)) + center

    # Normalize
    if norm_by == 'area': # normalize by the area
        I /= np.trapz(I, x=w)
        Iunit = '1/{0}'.format(waveunit)
    elif norm_by == 'max': # set maximum to 1
        I /= np.max(I)
        Iunit = '1'
    elif norm_by is None:
        Iunit=None
    else:
        raise ValueError('Unknown normalization type: `norm_by` = {0}'.format(norm_by))

    # Scale
    I *= scale
#    if Iunit is not None and scale != 1:
#        Iunit += 'x{0:.2f}'.format(scale)

    # Plot slit
    if bplot:
        plot_slit(w, I, waveunit=waveunit, Iunit=Iunit)

    return w, I

# trapezoidal instrumental broadening function of base base nm and top top nm
def trapezoidal_slit(top, base, wstep, center=0, norm_by='area', bplot=False,
                    waveunit='', scale=1, footerspacing=0):
    r""" Build a trapezoidal slit. Remember that FWHM = (top + base) / 2

    
    Parameters    
    ----------

    top: (nm)
        top of the trapeze

    base: (nm)
        base of the trapeze

    wstep: (nm)
        wavelength step

    center: (nm)
        center wavelength for the wavelength axs of the slit function

    norm_by: 'area', 'max'
        normalisation type. 'area' conserves energy. 'max' is what is
        done in Specair and changes units. Default 'area'

    bplot: boolean
        plot normalized slit function (for debugging). Default False

    wavespace: '', 'nm', 'cm-1'
        used for plot only. Slit function is generated assuming you use the
        correct wavespace. No conversions are made here.

    scale: float
        multiply slit by an arbitrary factor. Default 1.

    footerspacing: int
        spacing (footer) on left and right. Default 10.


    Returns
    -------

    Slit function of shape
                 _______
                /       \
               / :     : \
         _ _ _/  :     :  \____
         0    :  :     :  :   0
         :   -b  :     :  b   :
       -b-f     -t     t     b+f
       

    with `t`and `b` the half-top and half-base in number of elements, `f` spacing
    on left & right. FWHM = (t+b+1)*wstep


    Notes
    -----
    
    slit is generated with an odd number of elements and centered on its
    maximum. This may result in sightly different FWHM if wstep is too large.
    However, slit is corrected so that effective FWHM is preserved.

    """


    if top > base:
        top, base = base, top

    FWHM = (base+top)/2
    b = 2*int(top/wstep//2)+1              # number of points on top (even)

    # Build first half
    slope = 1/(FWHM-b*wstep)
    a = int(FWHM/wstep)-b                  # number of points in slope
    I = 1-np.arange(0, a+1)*slope*wstep    # slope

    if len(I) == 0:
        I = np.ones(1)

    if I[-1] < 0:
        I[-1] = 0
    elif I[-1] > 0:
        # add one extra zero
        footerspacing += 1

    # Zeros
#    if abs((base-top) % wstep) < wstep:  # add one extra zero when a divider
#        footerspacing += 1

    f = int(footerspacing)   # add zeros (footer)
    I = np.hstack((I, np.zeros(f)))
#    w = np.linspace(0, len(I)*wstep, len(I))

    # Mirror to get second half, add top
    I = np.hstack((I[1:][::-1], np.ones(b), I[1:]))
    w = wstep*np.linspace(-len(I)/2, len(I)/2, len(I)) + center

    # Normalize
    if norm_by == 'area': # normalize by the area
        I /= np.trapz(I, x=w)
        Iunit = '1/{0}'.format(waveunit)
    elif norm_by == 'max': # set maximum to 1
        I /= np.max(I)
        Iunit = '1'
    elif norm_by is None:
        Iunit = None
    else:
        raise ValueError('Unknown normalization type: `norm_by` = {0}'.format(norm_by))

    # scale
    I *= scale
#    if Iunit is not None and scale != 1:
#        Iunit += 'x{0:.2f}'.format(scale)

    # Plot slit
    if bplot:
        plot_slit(w, I, waveunit=waveunit, Iunit=Iunit)

    return w, I


def gaussian_slit(FWHM, wstep, center=0, norm_by='area', bplot=False,
                  waveunit='', calc_range=4, scale=1, footerspacing=0):
    r''' Generate normalized slit function


    Parameters    
    ----------

    FWHM: (nm)
        full-width at half maximum

    wstep: (nm)
        wavelength step

    center: (nm)
        center wavelength for the wavelength axs of the slit function

    norm_by: 'area', 'max'
        normalisation type. 'area' conserves energy. 'max' is what is
        done in Specair and changes units. Default 'area'

    bplot: boolean
        plot normalized slit function (for debugging). Default False

    wavespace: '', 'nm', 'cm-1'
        used for plot only. Slit function is generated assuming you use the
        correct wavespace. No conversions are made here.

    calc_range: float   (number of sigma)
        function broadening range expressed in number of standard deviation

    scale: float
        multiply slit by an arbitrary factor. Default 1.

    footerspacing: int
        spacing (footer) on left and right. Default 10


    Returns
    -------

    Slit function of shape
    
                 .-.
                /   \
               /     \
         _ _.-'       `-._ _
         0    [gaussian]   0
         :  :     :     :  :
       -c-f -c    0     c  c+f          (c : calc_range)
       

   with `c` cutoff in number of elements, `f` spacing on left & right


    Notes
    -----
    
    slit is generated with an odd number of elements and centered on its
    maximum. This may result in sightly different FWHM if wstep is too large.

    
    '''

    f = int(footerspacing)   # spacing (footer) on left and right
    sigma = FWHM / 2 / sqrt(2 * ln(2))

    a = 2*int(calc_range*sigma//wstep/2)    # half-base in number of elements (even number)

    # 2 sigma: gaussian non calculated residual: 5%
    # 3 sigma: gaussian non calculated residual: 1%

    w0 = wstep*np.linspace(-a, a, 2*a+1) # centered
    Igauss = exp(-w0**2 / (2 * sigma**2))

    I = np.hstack((np.zeros(f), Igauss, np.zeros(f)))
    w = wstep*np.linspace(-(a+f), (a+f), 2*a+2*f+1) + center

    # Normalize
    if norm_by == 'area': # normalize by the area
        I /= np.trapz(I, x=w)
        Iunit = '1/{0}'.format(waveunit)
    elif norm_by == 'max': # set maximum to 1
        I /= np.max(I)
        Iunit = '1'
    elif norm_by is None:
        Iunit=None
    else:
        raise ValueError('Unknown normalization type: `norm_by` = {0}'.format(norm_by))

    # scale
    I *= scale
#    if Iunit is not None and scale != 1:
#        Iunit += 'x{0:.2f}'.format(scale)

    # Plot slit
    if bplot:
        plot_slit(w, I, waveunit=waveunit, Iunit=Iunit)

    return w, I


# %% Test
if __name__ == '__main__':
    from radis.test.tools.test_slit import _run_testcases
    print('Testing slit.py: ', _run_testcases())
