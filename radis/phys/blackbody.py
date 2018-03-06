# -*- coding: utf-8 -*-
"""

Notes
-----

Planck functions:
- planck: planck radiation with wavelength input
- planck_wn: planck radiation with wavenumber input
- sPlanck: a RADIS Spectrum blackbody object

Example
-------

Generate Earth blackbody::
    
    s = sPlanck(wavelength_min=3000, wavelength_max=50000,
                T=288, eps=1)
    s.plot()

"""

from __future__ import absolute_import
from numpy import exp, arange, ones_like, zeros_like, inf
from radis.phys.constants import k_b, c, h
from radis.phys.constants import k_b_CGS, c_CGS, h_CGS
from radis.phys.units import conv2, Q_

def planck(lmbda, T, eps=1, unit='mW/sr/cm2/nm'):
    ''' Planck function for blackbody radiation

    
    Parameters    
    ----------

    λ: np.array   (nm)
       wavelength

    T: float    (K)
        equilibrium temperature

    eps: grey-body emissivity
        default 1

    unit: output unit
        default 'mW/sr/cm2/nm'


    Returns
    -------

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

    
    Parameters    
    ----------

    wavenum: np.array   (cm-1)
       wavenumber

    T: float    (K)
        equilibrium temperature

    eps: grey-body emissivity
        default 1

    unit: str
        output unit. Default 'mW/sr/cm2/cm_1'


    Returns
    -------

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

# %% Predefined Spectra objects
        
def sPlanck(wavenum_min=None, wavenum_max=None,
            wavelength_min=None, wavelength_max=None,
            T=None, eps=1, 
            wstep=0.01, medium='air',
            **kwargs):
    ''' Return a RADIS Spectrum object with blackbody radiation. 
    
    It's easier to plug in a MergeSlabs / SerialSlabs config than the Planck
    radiance calculated by iPlanck. And you don't need to worry about units as
    they are handled internally.
    
    See neq.spec.Spectrum documentation for more information
    
    Parameters
    ----------
    
    wavenum_min / wavenum_max: (cm-1)
        minimum / maximum wavenumber to be processed in cm^-1. 

    wavelength_min / wavelength_max: (nm)
        minimum / maximum wavelength to be processed in nm

    T: float (K)
        blackbody temperature
        
    eps: float [0-1]
        blackbody emissivity. Default 1
    
    Other Parameters
    ----------------
    
    wstep: float (cm-1 or nm)
        wavespace step for calculation
        
    **kwargs: other keyword inputs
        all are forwarded to spectrum conditions. For instance you can add
        a 'path_length=1' after all the other arguments

    Example
    -------
    
    Generate Earth blackbody::
        
        s = sPlanck(wavelength_min=3000, wavelength_max=50000,
                    T=288, eps=1)
        s.plot()

    '''
    
    from radis.spectrum.spectrum import Spectrum

    # Check inputs    
    if ((wavelength_min is not None or wavelength_max is not None) and
        (wavenum_min is not None or wavenum_max is not None)):
        raise ValueError('Wavenumber and Wavelength both given... you twart')

    if (wavenum_min is not None and wavenum_max is not None):
        assert(wavenum_min<wavenum_max)
        waveunit = 'cm-1'
    else:
        assert(wavelength_min<wavelength_max)
        waveunit = 'nm'
    
    if T is None:
        raise ValueError('T must be defined')
        
    if not (eps>=0 and eps<= 1):
        raise ValueError('Emissivity must be in [0-1]')
    
    # Test range is correct:
    if waveunit == 'cm-1':
        w = arange(wavenum_min,wavenum_max+wstep,wstep) #generate the vector of wavenumbers (shape M)
        Iunit = 'mW/sr/cm2/cm_1'
        I = planck_wn(w, T, eps=eps, unit=Iunit)
    else:
        w = arange(wavelength_min,wavelength_max+wstep,wstep) #generate the vector of wavenumbers (shape M)
        Iunit = 'mW/sr/cm2/nm'
        I = planck(w, T, eps=eps, unit=Iunit)
        
    conditions = {'wstep':wstep,
                  'medium':medium}
    conditions.update(**kwargs)  # add all extra parameters in conditions (ex: path_length)
    
    return Spectrum(quantities={'radiance_noslit':(w, I),
                              'transmittance_noslit':(w, zeros_like(w)),
                              'absorbance':(w, ones_like(w)*inf)},
                    conditions=conditions,
                    units={'radiance_noslit':Iunit,
                           'transmittance_noslit':'I/I0',
                           'absorbance':'-ln(I/I0)'},
                    cond_units={'wstep':waveunit},
                    waveunit=waveunit)


if __name__ == '__main__':
    from radis.test.phys.test_blackbody import test_planck
    test_planck()