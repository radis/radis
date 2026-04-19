# -*- coding: utf-8 -*-
"""
Summary
-------

Runs tests for radis/phys so that they can be accessed by pytest (and hopefully
the CI test suite)

Examples
--------

Run all tests::

    pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)::

    pytest -m fast

-------------------------------------------------------------------------------

"""

import numpy as np
import pytest
from numpy import isclose

from radis.phys.convert import (
    J2K,
    K2J,
    J2cm,
    J2eV,
    K2cm,
    K2eV,
    atm2bar,
    atm2torr,
    bar2atm,
    bar2torr,
    cm2eV,
    cm2hz,
    cm2J,
    cm2K,
    cm2nm,
    cm2nm_air,
    dcm2dnm,
    dcm2dnm_air,
    dhz2dnm,
    div_safe,
    dnm2dcm,
    dnm2dhz,
    dnm_air2dcm,
    eV2cm,
    eV2J,
    eV2K,
    eV2nm,
    hz2cm,
    hz2nm,
    nm2cm,
    nm2eV,
    nm2hz,
    nm_air2cm,
    torr2atm,
    torr2bar,
    zero2nan,
)
from radis.phys.units import conv2, is_homogeneous


@pytest.mark.fast
def test_convert(verbose=True, *args, **kwargs):
    """Test conversions."""

    E = np.linspace(1, 5, 5)  # eV
    assert (J2eV(K2J(J2K(eV2J(E)))) == E).all()

    E = 2150  # cm-1
    assert J2cm(cm2J(E)) == E
    assert K2cm(cm2K(E)) == E
    if verbose:
        print(f"K -> cm: ~ {cm2K(1):.2f} K/cm")

    E = 1  # eV
    assert eV2cm(E) == J2cm(eV2J(1))
    assert round(eV2K(E), 0) == 11605
    assert K2eV(eV2K(E)) == E

    E = 250  # nm
    assert isclose(nm2eV(E), cm2eV(nm2cm(E)))
    assert isclose(eV2nm(nm2eV(E)), E)
    assert hz2nm(nm2hz(E)) == E

    fwhm = 1.5  # nm
    lbd_0 = 632.8  # nm
    assert isclose(fwhm, dcm2dnm(dnm2dcm(fwhm, lbd_0), nm2cm(lbd_0)))

    fwhm = 2e-3  # 0.002 nm
    lbd_0 = 632.8
    fwhm_hz = dnm2dhz(fwhm, lbd_0)  # ~ 1.5 GHz
    if verbose:
        print((f"{fwhm:.2g} nm broadening at {lbd_0} nm = {fwhm_hz * 1e-9:.2g} Ghz"))
    assert isclose(fwhm_hz * 1e-9, 1.4973307983125002)
    assert isclose(dhz2dnm(fwhm_hz, nm2hz(lbd_0)), fwhm)

    assert atm2bar(1) == 1.01325
    assert isclose(torr2bar(atm2torr(bar2atm(1))), 1)
    assert isclose(torr2atm(bar2torr(atm2bar(1))), 1)

    # Hz
    assert isclose(1 / hz2cm(1e9), 30, atol=0.1)  # 1 Ghz is about 30 cm
    assert hz2cm(cm2hz(600)) == 600


@pytest.mark.fast
def test_air_conversions(verbose=True, *args, **kwargs):
    """Test air wavelength conversion functions."""
    wl_cm1 = 15000  # cm-1
    wl_nm_air = cm2nm_air(wl_cm1)
    # Round-trip consistency
    assert isclose(nm_air2cm(wl_nm_air), wl_cm1, rtol=1e-6)
    assert wl_nm_air < cm2nm(wl_cm1)  # air < vacuum
    # Hard-coded reference values (Ciddor 1996 dispersion formula)
    assert isclose(wl_nm_air, 666.483, rtol=1e-4)
    assert isclose(cm2nm(wl_cm1), 666.667, rtol=1e-4)
    # Delta conversions
    fwhm_cm, nu_0 = 10, 15000
    fwhm_nm_air = dcm2dnm_air(fwhm_cm, nu_0)
    assert isclose(dnm_air2dcm(fwhm_nm_air, cm2nm_air(nu_0)), fwhm_cm, rtol=1e-6)
    assert isclose(fwhm_nm_air, 0.4443, rtol=1e-3)


@pytest.mark.fast
def test_safe_division_functions(verbose=True, *args, **kwargs):
    """Test zero2nan and div_safe helper functions."""
    assert np.isnan(zero2nan(0)) and zero2nan(5) == 5
    arr = np.array([0.0, 1.0, 2.0, 0.0, 3.0])
    result = zero2nan(arr)
    assert np.isnan(result[0]) and np.isnan(result[3]) and result[2] == 2.0
    safe_reciprocal = div_safe(lambda x: 1 / x)
    result = safe_reciprocal(np.array([1.0, 2.0, 0.0, 4.0]))
    assert result[0] == 1.0 and result[1] == 0.5 and np.isnan(result[2])


@pytest.mark.fast
def test_units(verbose=True, *args, **kwargs):

    # Test unit-ware arrays
    # RADIS pint-aware array
    #    b = (print(a.to('V * cm^2'))==print(res))

    # Test conversion
    convtable = [
        (500, "nm", 0.5, "µm"),
        (1, "erg/s", 1e-7, "W"),
        (1, "m**2", 10000, "cm**2"),
    ]
    for a, f, r, t in convtable:
        cr = conv2(a, f, t)
        if verbose:
            print(f"{a} {f} = {cr} {t}")
        assert isclose(cr, r)

    # Ensures that an error is raised if units with angles are converted
    # (mathematically correct because angles are dimensionless, but prone
    # to user errors)
    assert not is_homogeneous("mW/cm2/sr/nm", "mW/cm2/nm")
    with pytest.raises(TypeError):
        conv2(1, "mW/cm**2/sr/nm", "mW/cm**2/nm")


def _run_testcases(*args, **kwargs):

    test_convert(*args, **kwargs)
    test_air_conversions(*args, **kwargs)
    test_safe_division_functions(*args, **kwargs)
    test_units(*args, **kwargs)
    return True


if __name__ == "__main__":
    print("testing phys.py:", _run_testcases(verbose=True))
