# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 15:27:36 2017

Planck functions:
- planck: planck radiation with wavelength input
- planck_wn: planck radiation with wavenumber input
- see also neq.spec.sPlanck: a NeQ Spectrum blackbody object

"""

from __future__ import absolute_import
from numpy import exp
from radis.phys.constants import k_b, c, h
from radis.phys.constants import k_b_CGS, c_CGS, h_CGS
from radis.phys.units import conv2, Q_

def planck(lmbda, T, eps=1, unit='mW/sr/cm2/nm'):
    ''' Planck function for blackbody radiation

    Input
    ---------

    λ: np.array   (nm)
       wavelength

    T: float    (K)
        equilibrium temperature

    eps: grey-body emissivity
        default 1

    unit: output unit
        default 'mW/sr/cm2/nm'

    Output
    ---------

    planck: np.array   (mW.sr-1.cm-2/nm)
        equilibrium radiance

    '''

    k = k_b
    lbd = lmbda * 1e-9
    iplanck = eps*(2*h*c**2/lbd**5) * 1/(exp(h*c/(lbd*k*T)) - 1)   # S.I  (W.sr-1.m-3)
    iplanck *= 1e-10  # W.sr-1.m-3 >>> mW.sr-1.cm-2.nm-1

    if Q_(unit) != Q_('mW/sr/cm2/nm'):
        iplanck = conv2(iplanck, 'mW/sr/cm2/nm', unit)

    return iplanck

def planck_wn(wavenum, T, eps=1, unit='mW/sr/cm2/cm_1'):
    ''' Planck function for blackbody radiation, wavenumber version

    Input
    ---------

    wavenum: np.array   (cm-1)
       wavenumber

    T: float    (K)
        equilibrium temperature

    eps: grey-body emissivity
        default 1

    unit: str
        output unit. Default 'mW/sr/cm2/cm_1'

    Output
    ---------

    planck: np.array   default (mW/sr/cm2/cm_1)
        equilibrium radiance

    '''

    k = k_b_CGS
    h = h_CGS
    c = c_CGS

    iplanck = eps*(2*h*c**2*wavenum**3) * 1/(exp(h*c*wavenum/(k*T)) - 1) 
    # iplanck in erg/s/sr/cm2/cm-1
    iplanck *= 1e-4     # erg/s/sr/cm2/cm-1 > mW/sr/cm^2/cm-1

    if Q_(unit) != Q_('mW/sr/cm2/cm_1'):
        iplanck = conv2(iplanck, 'mW/sr/cm2/cm_1', unit)

    return iplanck


if __name__ == '__main__':
    pass