# -*- coding: utf-8 -*-
"""Constants and correlations for air. Wavelength correction due to medium
dispersion.

Routine Listing
---------------

- :func:`~radis.phys.air.air_index_dispersion`
- :func:`~radis.phys.air.vacuum2air`
- :func:`~radis.phys.air.air2vacuum`

-------------------------------------------------------------------------------
"""


def air_index_dispersion(lbd):
    """Return air index dispersion as a function of wavelength with the
    relation of Ciddor [1]_

    Parameters
    ----------

    lbd: array-like (µm)
        wavelength


    References
    ----------

    .. [1] `P. E. Ciddor. "Refractive index of air: new equations for the visible and near infrared", Appl. Optics 35, 1566-1573 (1996) <https://refractiveindex.info/?shelf=other&book=air&page=Ciddor>`__

    Standard air: dry air at 15 °C, 101.325 kPa and with 450 ppm CO2 content
    """

    n = 1 + 0.05792105 / (238.0185 - lbd ** -2) + 0.00167917 / (57.362 - lbd ** -2)

    return n


# Wavelength medium conversion functions


def vacuum2air(wavelength):
    """Converts wavelength as seen in vacuum to wavelength as seen in air.

    Parameters
    ----------

    wavelength: array-like (nm)
        wavelength in vacuum


    Returns
    -------

    wavelength: array-like (nm)
        wavelength in air


    See Also
    --------

    :func:`~radis.phys.air.air2vacuum`
    """

    air_index = air_index_dispersion(wavelength * 1e-3)  # nm > µm
    return wavelength / air_index


def air2vacuum(wavelength):
    """Converts wavelength as seen in air to wavelength as seen in vacuum.

    Parameters
    ----------

    wavelength: array-like (nm)
        wavelength in air


    Returns
    -------

    wavelength: array-like (nm)
        wavelength in vacuum


    Note
    ----

    Not exactly true, as air_index_dispersion is defined for vacuum wavelength
    However, air_index_dispersion doesnt vary much on 1-2 cm-1 (which is typical
    of air index dispersion effect on air/vacuum wavelength in the mid-IR)

    Estimate the error you make with::

        vacuum2air(air2vacuum(w))-w


    See Also
    --------

    :func:`~radis.phys.air.vacuum2air`
    """

    air_index = air_index_dispersion(wavelength * 1e-3)  # nm > µm
    return air_index * wavelength
