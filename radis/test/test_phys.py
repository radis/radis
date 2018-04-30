# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 14:56:25 2017

@author: erwan

Runs tests for neq/phys so that they can be accessed by pytest (and hopefully 
the CI test suite)

Examples
--------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have 'fast' in their name)

>>> pytest -k fast

"""



from __future__ import absolute_import, unicode_literals, print_function

import numpy as np
from numpy import isclose
from radis.phys.convert import (J2eV, J2cm, cm2J, eV2cm, eV2K, eV2nm, nm2eV,
                                dnm2dcm, nm2cm, dnm2dhz, dhz2dnm, cm2eV, eV2J, K2J,
                                J2K, dcm2dnm, K2eV, hz2nm, nm2hz,
                                torr2bar, torr2atm, bar2torr, bar2atm, atm2torr,
                                atm2bar)
from radis.phys.units import uarray, Q_, conv2
import pytest

@pytest.mark.fast
def test_convert(verbose=True, *args, **kwargs):
    ''' Test conversions  '''
    
    E = np.linspace(1, 5, 5)  # eV
    assert (J2eV(K2J(J2K(eV2J(E)))) == E).all()

    E = 2150  # cm-1
    assert J2cm(cm2J(E))

    E = 1  # eV
    assert (eV2cm(E) == J2cm(eV2J(1)))
    assert (round(eV2K(E), 0) == 11605)
    assert K2eV(eV2K(E)) == E
    
    E = 250 # nm
    assert isclose(nm2eV(E),cm2eV(nm2cm(E)))
    assert (eV2nm(nm2eV(E))==E)
    assert hz2nm(nm2hz(E)) == E

    fwhm = 1.5 # nm
    lbd_0 = 632.8 # nm
    assert (isclose(fwhm, dcm2dnm(dnm2dcm(fwhm, lbd_0), nm2cm(lbd_0))))

    fwhm = 2e-3  # 0.002 nm
    lbd_0 = 632.8
    fwhm_hz = dnm2dhz(fwhm, lbd_0)         #  ~ 1.5 GHz
    if verbose: print(('{0:.2g} nm broadening at {1} nm = {2:.2g} Ghz'.format(fwhm, lbd_0, fwhm_hz*1e-9)))
    assert isclose(fwhm_hz*1e-9, 1.4973307983125002)
    assert isclose(dhz2dnm(fwhm_hz, nm2hz(lbd_0)), fwhm)
    
    assert atm2bar(1) == 1.01325
    assert isclose(torr2bar(atm2torr(bar2atm(1))), 1)
    assert isclose(torr2atm(bar2torr(atm2bar(1))), 1)


    return True

def test_units(verbose=True, *args, **kwargs):

    # Test unit-ware arrays
    a = uarray(np.linspace(10, 100, 10), 'Td')          # RADIS pint-aware array
    res = Q_(np.linspace(1e-16, 1e-15, 10), 'V * cm^2')  # pint definition

    assert (np.round(np.array(a.to('V * cm^2')) - np.array(res), 5)
         == np.zeros_like(res)).all()
#    b = (print(a.to('V * cm^2'))==print(res))


    # Test conversion 
    convtable = [
            (500, 'nm', 0.5, 'µm'),
            (1, 'erg/s', 1e-7, 'W'),
            (1, 'm2',10000, 'cm2')
            ]
    for a, f, r, t in convtable:
        cr = conv2(a, f, t)
        if verbose: print(('{0} {1} = {2} {3}'.format(a, f, cr, t)))
        assert isclose(cr,r)

    return True

def _run_testcases(*args, **kwargs):
    
    assert test_convert(*args, **kwargs)
    assert test_units(*args, **kwargs)
    
    return True

if __name__== '__main__':
    print('testing phys.py:', _run_testcases(verbose=True))
