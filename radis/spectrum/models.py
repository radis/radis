# -*- coding: utf-8 -*-
"""
Summary
-------

Models built around the :class:`~radis.spectrum.spectrum.Spectrum` class


Routine Listing
---------------

- :func:`~radis.spectrum.models.calculated_spectrum`,
- :func:`~radis.spectrum.models.experimental_spectrum`,
- :func:`~radis.spectrum.models.transmittance_spectrum`,

See Also
--------

To extract some spectral quantities from a Spectrum, and create a new Spectrum,
see the functions in :py:mod:`radis.spectrum.operations`:

- :func:`~radis.spectrum.operations.Radiance`
- :func:`~radis.spectrum.operations.Radiance_noslit`
- :func:`~radis.spectrum.operations.Transmittance`
- :func:`~radis.spectrum.operations.Transmittance_noslit`


-------------------------------------------------------------------------------


"""

import numpy as np

from radis.spectrum.spectrum import Spectrum

# %% Array-to-Spectrum functions


def calculated_spectrum(
    w,
    I,
    wunit="nm",
    Iunit="mW/cm2/sr/nm",
    conditions=None,
    cond_units=None,
    populations=None,
    name=None,
) -> Spectrum:
    """Convert ``(w, I)`` into a :py:class:`~radis.spectrum.spectrum.Spectrum`
    object that has unit conversion, plotting and slit convolution
    capabilities.

    Parameters
    ----------
    w: np.array
        wavelength, or wavenumber
    I: np.array
        intensity (no slit)
    wunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
        wavespace unit: wavelength in air (``'nm'``), wavenumber
        (``'cm-1'``), or wavelength in vacuum (``'nm_vac'``). Default ``'nm'``.
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


    Examples
    --------

    ::

        # w, I are numpy arrays for wavelength and radiance
        from radis import calculated_spectrum
        s = calculated_spectrum(w, I, wunit='nm', Iunit='W/cm2/sr/nm')     # creates 'radiance_noslit'



    See Also
    --------

    :func:`~radis.spectrum.models.transmittance_spectrum`,
    :func:`~radis.spectrum.models.experimental_spectrum`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_array`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`,
    :func:`~radis.tools.database.load_spec`
    """

    return Spectrum.from_array(
        np.array(w),
        np.array(I),
        "radiance_noslit",
        wunit=wunit,
        unit=Iunit,
        conditions=conditions,
        cond_units=cond_units,
        populations=populations,
        name=name,
    )


def transmittance_spectrum(
    w, T, wunit="nm", Tunit="", conditions=None, cond_units=None, name=None
) -> Spectrum:
    """Convert ``(w, I)`` into a :py:class:`~radis.spectrum.spectrum.Spectrum`
    object that has unit conversion, plotting and slit convolution
    capabilities.

    Parameters
    ----------
    w: np.array
        wavelength, or wavenumber
    T: np.array
        transmittance (no slit)
    wunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
        wavespace unit: wavelength in air (``'nm'``), wavenumber
        (``'cm-1'``), or wavelength in vacuum (``'nm_vac'``). Default ``'nm'``.
    Iunit: str
        intensity unit. Default ``""`` (adimensionned)


    Other Parameters
    ----------------
    conditions: dict
        (optional) calculation conditions to be stored with Spectrum
    cond_units: dict
        (optional) calculation conditions units
    name: str
        (optional) give a name


    Examples
    --------
    ::

        # w, T are numpy arrays for wavelength and transmittance
        from radis import transmittance_spectrum
        s2 = transmittance_spectrum(w, T, wunit='nm')                       # creates 'transmittance_noslit'


    See Also
    --------

    :func:`~radis.spectrum.models.calculated_spectrum`,
    :func:`~radis.spectrum.models.experimental_spectrum`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_array`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`,
    :func:`~radis.tools.database.load_spec`
    """

    return Spectrum.from_array(
        np.array(w),
        np.array(T),
        "transmittance_noslit",
        wunit=wunit,
        unit=Tunit,
        conditions=conditions,
        cond_units=cond_units,
        name=name,
    )


def experimental_spectrum(
    w, I, wunit="nm", Iunit="count", conditions={}, cond_units=None, name=None
) -> Spectrum:
    """Convert ``(w, I)`` into a :py:class:`~radis.spectrum.spectrum.Spectrum`
    object that has unit conversion and plotting capabilities. Convolution is
    not available as the spectrum is assumed to have be measured experimentally
    (hence it is already convolved with the slit function)

    Parameters
    ----------
    w: np.array
        wavelength, or wavenumber
    I: np.array
        intensity
    wunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
        wavespace unit: wavelength in air (``'nm'``), wavenumber
        (``'cm-1'``), or wavelength in vacuum (``'nm_vac'``). Default ``'nm'``.
    Iunit: str
        intensity unit (can be ``'count'``, 'mW/cm2/sr/nm', etc...). Default
        ``'count'`` (i.e., non calibrated output)

    Other Parameters
    ----------------
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
        w, I = loadtxt('my_file.txt').T    # transpose is often useful, depending on your data.
        s = experimental_spectrum(w, I, Iunit='mW/cm2/sr/nm')             # creates 'radiance'
        s.plot()


    See Also
    --------
    :func:`~radis.spectrum.models.calculated_spectrum`,
    :func:`~radis.spectrum.models.transmittance_spectrum`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_array`,
    :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`,
    :func:`~radis.tools.database.load_spec`
    """

    if np.shape(w) != np.shape(I):
        raise ValueError(
            "Wavelength {0} and intensity {1} do not have the same shape".format(
                np.shape(w), np.shape(I)
            )
        )
    return Spectrum.from_array(
        np.array(w),
        np.array(I),
        "radiance",
        wunit=wunit,
        unit=Iunit,
        conditions=conditions,
        cond_units=cond_units,
        name=name,
    )
