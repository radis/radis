# -*- coding: utf-8 -*-
"""
Summary
-------

Models built around the :class:`~radis.spectrum.spectrum.Spectrum` class


Routine Listing
---------------

- :func:`~radis.spectrum.models.Transmittance`
- :func:`~radis.spectrum.models.Radiance`


-------------------------------------------------------------------------------


"""

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.spectrum.spectrum import Spectrum


# %% Array-to-Spectrum functions


def calculated_spectrum(w, I, wunit='nm', Iunit='mW/cm2/sr/nm',
                        conditions=None, cond_units=None, populations=None,
                        name=None): # -> Spectrum:
    ''' Convert (w, I) into a Spectrum object that has unit conversion, plotting
    and slit convolution capabilities


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

    return Spectrum.from_array(w, I, 'radiance_noslit',
                               waveunit=wunit, unit=Iunit,
                               conditions=conditions, cond_units=cond_units,
                               populations=populations,
                               name=name)


def transmittance_spectrum(w, T, wunit='nm', Tunit='I/I0',
                           conditions=None, cond_units=None,
                           name=None): # -> Spectrum:
    ''' Convert (w, I) into a Spectrum object that has unit conversion, plotting
    and slit convolution capabilities


    Parameters    
    ----------

    w, I: np.array
        wavelength and transmittance (no slit)

    wunit: 'nm', 'cm-1'
        wavespace unit

    Iunit: str
        intensity unit. Default 'I/I0'


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

    :func:`~radis.spectrum.spectrum.calculated_spectrum`, 
    :func:`~radis.spectrum.spectrum.experimental_spectrum`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_array`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`,
    :func:`~radis.tools.database.load_spec`

    '''

    return Spectrum.from_array(w, T, 'transmittance_noslit',
                               waveunit=wunit, unit=Tunit,
                               conditions=conditions, cond_units=cond_units,
                               name=name)


def experimental_spectrum(w, I, wunit='nm', Iunit='counts',
                          conditions=None, cond_units=None, name=None): # -> Spectrum:
    ''' Convert (w, I) into a Spectrum object that has unit conversion and plotting
    capabilities. Convolution is not available as the spectrum is assumed to
    be measured experimentally (hence deconvolution of the slit function would
    be required)


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

    conditions: dict
        (optional) calculation conditions to be stored with Spectrum

    cond_units: dict
        (optional) calculation conditions units

    name: str
        (optional) give a name


    See Also
    --------

    :func:`~radis.spectrum.spectrum.calculated_spectrum`, 
    :func:`~radis.spectrum.spectrum.transmittance_spectrum`, 
    :meth:`~radis.spectrum.spectrum.Spectrum.from_array`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`,
    :func:`~radis.tools.database.load_spec`

    '''

    return Spectrum.from_array(w, I, 'radiance',
                               waveunit=wunit, unit=Iunit,
                               conditions=conditions, cond_units=cond_units,
                               name=name)