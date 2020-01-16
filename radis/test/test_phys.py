# -*- coding: utf-8 -*-
"""
Summary
-------

Runs tests for neq/phys so that they can be accessed by pytest (and hopefully 
the CI test suite)

Examples
--------

Run all tests::

    pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)::

    pytest -m fast

-------------------------------------------------------------------------------

"""


from __future__ import absolute_import, unicode_literals, print_function
import numpy as np
from numpy import isclose
from radis.phys.convert import (
    J2eV,
    J2cm,
    cm2J,
    eV2cm,
    eV2K,
    eV2nm,
    nm2eV,
    cm2K,
    K2cm,
    dnm2dcm,
    nm2cm,
    dnm2dhz,
    dhz2dnm,
    cm2eV,
    eV2J,
    K2J,
    J2K,
    dcm2dnm,
    K2eV,
    hz2nm,
    nm2hz,
    hz2cm,
    cm2hz,
    torr2bar,
    torr2atm,
    bar2torr,
    bar2atm,
    atm2torr,
    atm2bar,
)
from radis.phys.units import uarray, Q_, conv2
import pytest


@pytest.mark.fast
def test_convert(verbose=True, *args, **kwargs):
    """ Test conversions  """

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
    assert eV2nm(nm2eV(E)) == E
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
    a = uarray(np.linspace(10, 100, 10), "Td")
    res = Q_(np.linspace(1e-16, 1e-15, 10), "V * cm^2")  # pint definition

    assert (
        np.round(np.array(a.to("V * cm^2")) - np.array(res), 5) == np.zeros_like(res)
    ).all()
    #    b = (print(a.to('V * cm^2'))==print(res))

    # Test conversion
    convtable = [
        (500, "nm", 0.5, "µm"),
        (1, "erg/s", 1e-7, "W"),
        (1, "m2", 10000, "cm2"),
    ]
    for a, f, r, t in convtable:
        cr = conv2(a, f, t)
        if verbose:
            print(("{0} {1} = {2} {3}".format(a, f, cr, t)))
        assert isclose(cr, r)

    # Ensures that an error is raised if units with angles are converted
    # (mathematically correct because angles are dimensionless, but prone
    # to user errors)
    from radis.phys.units import DimensionalityError

    with pytest.raises(DimensionalityError):
        conv2(1, "mW/cm2/sr/nm", "mW/cm2/nm")

    return True


def _run_testcases(*args, **kwargs):

    assert test_convert(*args, **kwargs)
    assert test_units(*args, **kwargs)

    return True


if __name__ == "__main__":
    print("testing phys.py:", _run_testcases(verbose=True))
