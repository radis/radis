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
    dcm2dnm,
    dhz2dnm,
    dnm2dcm,
    dnm2dhz,
    eV2cm,
    eV2J,
    eV2K,
    eV2nm,
    hz2cm,
    hz2nm,
    nm2cm,
    nm2eV,
    nm2hz,
    torr2atm,
    torr2bar,
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
        print("K -> cm: ~ {0:.2f} K/cm".format(cm2K(1)))

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
        print(
            (
                "{0:.2g} nm broadening at {1} nm = {2:.2g} Ghz".format(
                    fwhm, lbd_0, fwhm_hz * 1e-9
                )
            )
        )
    assert isclose(fwhm_hz * 1e-9, 1.4973307983125002)
    assert isclose(dhz2dnm(fwhm_hz, nm2hz(lbd_0)), fwhm)

    assert atm2bar(1) == 1.01325
    assert isclose(torr2bar(atm2torr(bar2atm(1))), 1)
    assert isclose(torr2atm(bar2torr(atm2bar(1))), 1)

    # Hz
    assert isclose(1 / hz2cm(1e9), 30, atol=0.1)  # 1 Ghz is about 30 cm
    assert hz2cm(cm2hz(600)) == 600

    return True


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
            print(("{0} {1} = {2} {3}".format(a, f, cr, t)))
        assert isclose(cr, r)

    # Ensures that an error is raised if units with angles are converted
    # (mathematically correct because angles are dimensionless, but prone
    # to user errors)
    assert not is_homogeneous("mW/cm2/sr/nm", "mW/cm2/nm")
    with pytest.raises(TypeError):
        conv2(1, "mW/cm**2/sr/nm", "mW/cm**2/nm")

    return True


def _run_testcases(*args, **kwargs):

    assert test_convert(*args, **kwargs)
    assert test_units(*args, **kwargs)

    return True


if __name__ == "__main__":
    print("testing phys.py:", _run_testcases(verbose=True))
