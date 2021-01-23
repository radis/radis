# -*- coding: utf-8 -*-
"""Conversion formulas between different units (works with numpy arrays)

Checks the order of magnitudes to help detect conversion errors (comparing with
standard order of magnitudes in plasma physics)

Examples
--------

Get equivalent width in nm of a 10cm-1 width at 380 nm

    >>> from radis import *
    >>> dcm2dnm(10, nm2cm(380))


-------------------------------------------------------------------------------
"""

import numpy as np

from radis.phys.air import air2vacuum, vacuum2air
from radis.phys.constants import c, eV, h, hc_k, k_b

# Make conv2 accessible from .convert:
from radis.phys.units import conv2

assert conv2  # removes flake8 warnings

# %% Energy units


def J2eV(E):
    """J to eV."""
    return E / eV


def J2K(E):
    """J to Kelvin."""
    return E / k_b


def J2cm(E):
    """J to cm-1."""
    return E / (h * c) / 100


def eV2J(E):
    """eV to J."""
    _asserteV(E)
    return E * eV


def eV2nm(E):
    """nm to eV."""
    _asserteV(E)
    return 1 / E * 1e9 * (h * c) / eV


def eV2K(E):
    """eV to K."""
    _asserteV(E)
    return E * eV / k_b


def eV2cm(E):
    """eV to cm."""
    _asserteV(E)
    return E * eV / (h * c) / 100


def K2eV(E):
    """eV to K."""
    _assertK(E)
    return E * k_b / eV


def K2J(E):
    """Kelvin to J."""
    _assertK(E)
    return E * k_b


def K2cm(E):
    """K to cm-1."""
    _assertK(E)
    return E / hc_k


def cm2J(E):
    """cm-1 to J."""
    _assertcm(E)
    return (E * 100) * (h * c)


def cm2K(E):
    """cm-1 to K.

    That's the classical 1.44 K/cm-1. See
    :data:`~radis.phys.constants.hc_k`
    """
    _assertcm(E)
    return E * hc_k


def cm2eV(E):
    """cm-1 to eV."""
    _assertcm(E)
    return (E * 100) * (h * c) / eV


# %% Wavelength, wavenumbers and frequencies


def cm2nm(wl_cm1):
    """cm-1 to (vacuum) nm."""
    return 1 / wl_cm1 * 1e9 / 100


def cm2hz(wl_cm1):
    """wavenumber to frequency wl_cm1, output in hz."""
    return wl_cm1 * c * 100


def nm2cm(wl_nm):
    """(vacuum) nm to cm-1."""
    return 1 / wl_nm * 1e9 / 100


def cm2nm_air(wl_cm1):
    """cm-1 to (air) nm.

    References
    ----------

    :func:`~radis.phys.air.vacuum2air'
    """
    return vacuum2air(cm2nm(wl_cm1))


def nm_air2cm(wl_nm_air):
    """(air) nm to cm-1.

    References
    ----------

    :func:`~radis.phys.air.air2vacuum'
    """
    return nm2cm(air2vacuum(wl_nm_air))


def nm2eV(wl_nm):
    """nm to eV."""
    return 1 / wl_nm * 1e9 * (h * c) / eV


def hz2nm(f_Hz):
    """frequency to wavelength f in Hz, output in nm."""
    return c * 1e9 / f_Hz


def hz2cm(f_Hz):
    """frequency to wavenumber f in Hz, output in cm-1."""
    return f_Hz / c / 100


def nm2hz(lbd_nm):
    """wavelength to frequency lbd in Hz, output in nm."""
    return c * 1e9 / lbd_nm


# Convert Broadenings


def dcm2dnm(delta_nu, nu_0):
    """Converts (ex: FWHM) from Δcm to Δnm.

    Parameters
    ----------

    delta_nu: float (cm-1)
        wavenumber broadening

    nu_0: float (cm-1)
        center wavenumber

    Returns
    -------

    delta_nm: float (nm)
        broadening in wavelength (vacuum)
    """
    return cm2nm(nu_0 - delta_nu / 2) - cm2nm(nu_0 + delta_nu / 2)


def dnm2dcm(delta_lbd, lbd_0):
    """Converts (ex: FWHM) from Δnm to Δcm.

    Parameters
    ----------

    delta_lbd: float (nm)
        wavelength broadening (vacuum)

    lbd_0: float (nm)
        center wavelength (vacuum)

    Returns
    -------

    delta_cm: float (cm-1)
        broadening in wavenumber
    """
    return nm2cm(lbd_0 - delta_lbd / 2) - nm2cm(lbd_0 + delta_lbd / 2)


def dcm2dnm_air(delta_nu, nu_0):
    """Converts (ex: FWHM) from Δcm to Δnm.

    Parameters
    ----------

    delta_nu: float (cm-1)
        wavenumber broadening

    nu_0: float (cm-1)
        center wavenumber

    Returns
    -------

    delta_nm: float (nm)
        broadening in wavelength (air)
    """
    return cm2nm_air(nu_0 - delta_nu / 2) - cm2nm_air(nu_0 + delta_nu / 2)


def dnm_air2dcm(delta_lbd, lbd_0):
    """Converts (ex: FWHM) from Δnm to Δcm.

    Parameters
    ----------

    delta_lbd: float (nm)
        wavelength broadening (air)

    lbd_0: float (nm)
        center wavelength (air)

    Returns
    -------

    delta_cm: float (cm-1)
        broadening in wavenumber
    """
    return nm_air2cm(lbd_0 - delta_lbd / 2) - nm_air2cm(lbd_0 + delta_lbd / 2)


def dhz2dnm(deltaf_hz, f_0):
    """Converts (ex: FWHM) from ΔHz to Δnm.

    Parameters
    ----------

    deltaf_hz: Hz
        frequency broadening

    nu_0: Hz
        center frequency
    """
    return hz2nm(f_0 - deltaf_hz / 2) - hz2nm(f_0 + deltaf_hz / 2)


def dnm2dhz(delta_lbd, lbd_0):
    """Converts (ex: FWHM) from Δnm to Δhz.

    Parameters
    ----------

    delta_lbd: nm
        wavelength broadening

    lbd_0: nm
        center wavelength
    """
    return nm2hz(lbd_0 - delta_lbd / 2) - nm2hz(lbd_0 + delta_lbd / 2)


# %% Pressure units


def torr2bar(p_torr):
    """Torr to bar."""
    return p_torr * 1.01325 / 760


def torr2atm(p_atm):
    return p_atm / 760


def bar2torr(p_bar):
    return p_bar / 1.01325 * 760


def bar2atm(p_bar):
    return p_bar / 1.01325


def atm2torr(p_atm):
    return p_atm * 760


def atm2bar(p_atm):
    return p_atm * 1.01325


# %% Assert functions


def _magn(x):
    return np.round((np.log10(np.abs(x))))


def _assertK(E):
    if np.sum(np.abs(E)) != 0:  # check E != 0 for both floats and arrays
        try:
            m = _magn(E)
            assert ((0 <= m) & (m <= 6)).all()
        except AssertionError:
            print(("Warning. Input values may not be in Kelvin", E, "K?"))


def _assertcm(E):
    if np.sum(np.abs(E)) != 0:  # check E != 0 for both floats and arrays
        try:
            m = _magn(E)
            assert ((1 <= m) & (m <= 5)).all()
        except AssertionError:
            print(("Warning. Input values may not be in cm-1", E, "cm-1?"))


def _asserteV(E):
    if np.sum(np.abs(E)) != 0:  # check E != 0 for both floats and arrays
        try:
            m = _magn(E)
            assert ((0 <= m) & (m <= 2)).all()
        except AssertionError:
            print(("Warning. Input values may not be in eV", E, "eV?"))


# %% Test
if __name__ == "__main__":
    from radis.test.test_phys import test_convert

    print(("Test :", test_convert()))
