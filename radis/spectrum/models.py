# -*- coding: utf-8 -*-
"""
Summary
-------

Models built around the :class:`~radis.spectrum.spectrum.Spectrum` class


Routine Listing
---------------

- :func:`~radis.spectrum.models.Transmittance`
- :func:`~radis.spectrum.models.Radiance`

- :func:`~radis.spectrum.models.calculated_spectrum`,
- :func:`~radis.spectrum.models.experimental_spectrum`,
- :func:`~radis.spectrum.models.transmittance_spectrum`,

-------------------------------------------------------------------------------


"""

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.spectrum.spectrum import Spectrum
import numpy as np


# %% Array-to-Spectrum functions


def calculated_spectrum(w, I, wunit='nm', Iunit='mW/cm2/sr/nm',
                        conditions=None, cond_units=None, populations=None,
                        name=None): # -> Spectrum:
    ''' Convert ``(w, I)`` into a :py:class:`~radis.spectrum.spectrum.Spectrum`  
    object that has unit conversion, plotting and slit convolution capabilities


    Parameters    
    ----------

    w, I: np.array
        wavelength and intensity

    wunit: 'nm', 'cm-1'
        wavespace unit

    Iunit: str
        intensity unit (can be 'counts', 'mW/cm2/sr/nm', etc...). Default
        'mW/cm2/sr/nm' (note that non-convoluted Specair spectra are in 'mW/cm2/sr/µm')


    Other Parameters
    ----------------

    conditions: dict
        (optional) calculation conditions to be stored with Spectrum. Default ``None``

    cond_units: dict
        (optional) calculation conditions units. Default ``None``

    populations: dict
        populations to be stored in Spectrum. Default ``None``

    name: str
        (optional) give a name


    See Also
    --------

    :func:`~radis.spectrum.spectrum.transmittance_spectrum`, 
    :func:`~radis.spectrum.spectrum.experimental_spectrum`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_array`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`,
    :func:`~radis.tools.database.load_spec`

    '''

    return Spectrum.from_array(np.array(w), np.array(I), 'radiance_noslit',
                               waveunit=wunit, unit=Iunit,
                               conditions=conditions, cond_units=cond_units,
                               populations=populations,
                               name=name)


def transmittance_spectrum(w, T, wunit='nm', Tunit='I/I0',
                           conditions=None, cond_units=None,
                           name=None): # -> Spectrum:
    ''' Convert ``(w, I)`` into a :py:class:`~radis.spectrum.spectrum.Spectrum`  
    object that has unit conversion, plotting and slit convolution capabilities


    Parameters    
    ----------

    w, I: np.array
        wavelength and transmittance (no slit)

    wunit: ``'nm'``, ``'cm-1'``
        wavespace unit

    Iunit: str
        intensity unit. Default ``'I/I0'``


    Other Parameters
    ----------------

    conditions: dict
        (optional) calculation conditions to be stored with Spectrum

    cond_units: dict
        (optional) calculation conditions units

    name: str
        (optional) give a name


    See Also
    --------

    :func:`~radis.spectrum.models.calculated_spectrum`, 
    :func:`~radis.spectrum.models.experimental_spectrum`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_array`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`,
    :func:`~radis.tools.database.load_spec`

    '''

    return Spectrum.from_array(np.array(w), np.array(T), 'transmittance_noslit',
                               waveunit=wunit, unit=Tunit,
                               conditions=conditions, cond_units=cond_units,
                               name=name)


def experimental_spectrum(w, I, wunit='nm', Iunit='counts', medium='air',
                          conditions={}, cond_units=None, name=None): # -> Spectrum:
    ''' Convert ``(w, I)`` into a :py:class:`~radis.spectrum.spectrum.Spectrum` 
    object that has unit conversion and plotting
    capabilities. Convolution is not available as the spectrum is assumed to
    have be measured experimentally (hence it is already convolved with the slit function)


    Parameters    
    ----------

    w, I: np.array
        wavelength and intensity

    wunit: 'nm', 'cm-1'
        wavespace unit

    Iunit: str
        intensity unit (can be 'counts', 'mW/cm2/sr/nm', etc...). Default
        'counts' (default Winspec output)

    Other Parameters
    ----------------

    medium: 'air', 'vacuum'
        which medium. Default 'air'

    conditions: dict
        (optional) calculation conditions to be stored with Spectrum

    cond_units: dict
        (optional) calculation conditions units

    name: str
        (optional) give a name

    Examples
    --------
    
    Load and plot an experimental spectrum::
    
        from numpy import loadtxt
        from radis import experimental_spectrum
        w, I = loadtxt('my_file.txt').T    # assuming 2 columns 
        s = experimental_spectrum(w, I, Iunit='mW/cm2/sr/nm')
        s.plot()
        

    See Also
    --------

    :func:`~radis.spectrum.spectrum.calculated_spectrum`, 
    :func:`~radis.spectrum.spectrum.transmittance_spectrum`, 
    :meth:`~radis.spectrum.spectrum.Spectrum.from_array`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`,
    :func:`~radis.tools.database.load_spec`

    '''
    
    if 'medium' in conditions:
        assert conditions['medium'] == medium
    else:
        conditions.update({'medium':medium})
        
    return Spectrum.from_array(np.array(w), np.array(I), 'radiance',
                               waveunit=wunit, unit=Iunit,
                               conditions=conditions, cond_units=cond_units,
                               name=name)