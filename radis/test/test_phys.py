# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 14:56:25 2017

@author: erwan

Runs tests for neq/phys so that they can be accessed by pytest (and hopefully 
the CI test suite)

Use
------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have 'fast' in their name)

>>> pytest -k fast

"""



from __future__ import absolute_import

import numpy as np
from radis.phys.convert import (J2eV, J2cm, cm2J, eV2cm, eV2K, eV2nm, nm2eV,
                                dnm2dcm, nm2cm, dnm2dhz, cm2eV, eV2J, K2J,
                                J2K, dcm2dnm)
from radis.phys.units import uarray, Q_, conv2

def test_convert__fast(*args, **kwargs):
    btest = True

    E = np.linspace(1, 5, 5)  # eV
    btest *= (J2eV(K2J(J2K(eV2J(E)))) == E).all()

    E = 2150  # cm-1
    btest *= J2cm(cm2J(E))

    E = 1  # eV
    btest *= (eV2cm(E) == J2cm(eV2J(1)))
    btest *= (round(eV2K(E), 0) == 11605)
    
    E = 250 # nm
    btest *= np.isclose(nm2eV(E),cm2eV(nm2cm(E)))
    btest *= (eV2nm(nm2eV(E))==E)

    fwhm = 1.5 # nm
    lbd_0 = 632.8 # nm
    btest *= (np.isclose(fwhm, dcm2dnm(dnm2dcm(fwhm, lbd_0), nm2cm(lbd_0))))

    fwhm = 2e-3  # nm
    lbd_0 = 632.8
    fwhm_hz = dnm2dhz(fwhm, lbd_0)
    print(('{0:.2g} nm broadening at {1} nm = {2:.2g} Ghz'.format(fwhm, lbd_0, fwhm_hz*1e-9)))
    btest *= np.isclose(fwhm_hz*1e-9, 1.4973307983125002)

    return bool(btest)

def test_units__fast(verbose=True, *args, **kwargs):

    b = True

    # Test unit-ware arrays
    a = uarray(np.linspace(10, 100, 10), 'Td')          # RADIS pint-aware array
    res = Q_(np.linspace(1e-16, 1e-15, 10), 'V * cm^2')  # pint definition

    b *= (np.round(np.array(a.to('V * cm^2')) - np.array(res), 5)
         == np.zeros_like(res)).all()
#    b = (print(a.to('V * cm^2'))==print(res))


    # Test conversion 
    convtable = [
            (500, 'nm', 0.5, 'µm'),
            (1, 'erg/s', 1e-7, 'W'),
            (1, 'm2',10000,'cm2')
            ]
    for a, f, r, t in convtable:
        cr = conv2(a, f, t)
        print(('{0} {1} = {2} {3}'.format(a, f, cr, t)))
        b *= np.isclose(cr,r)

    return bool(b)

def _run_testcases(*args, **kwargs):
    
    b = True
    
    b *= test_convert__fast(*args, **kwargs)
    b *= test_units__fast(*args, **kwargs)
    
    return bool(b)

if __name__== '__main__':
    print('testing phys.py:', _run_testcases(verbose=True))
