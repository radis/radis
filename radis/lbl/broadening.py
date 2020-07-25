# -*- coding: utf-8 -*-
"""

Summary
-------

A class to handle all broadening related functions (and unload factory.py)

BroadenFactory is inherited by SpectrumFactory eventually

Routine Listing
---------------

Most methods are written in inherited class with the following inheritance scheme:
    
:py:class:`~radis.lbl.loader.DatabankLoader` > :py:class:`~radis.lbl.base.BaseFactory` > 
:py:class:`~radis.lbl.broadening.BroadenFactory` > :py:class:`~radis.lbl.bands.BandFactory` > 
:py:class:`~radis.lbl.factory.SpectrumFactory` > :py:class:`~radis.lbl.parallel.ParallelFactory`

.. inheritance-diagram:: radis.lbl.parallel.ParallelFactory
   :parts: 1

PUBLIC FUNCTIONS - BROADENING

- :py:func:`radis.lbl.broadening.doppler_broadening_HWHM`
- :py:func:`radis.lbl.broadening.gaussian_lineshape`
- :py:func:`radis.lbl.broadening.pressure_broadening_HWHM`
- :py:func:`radis.lbl.broadening.lorentzian_lineshape`
- :py:func:`radis.lbl.broadening.voigt_broadening_HWHM`
- :py:func:`radis.lbl.broadening.voigt_lineshape`

PRIVATE METHODS - BROADENING
(all computational-heavy functions: calculates all lines broadening,
convolve them, apply them on all calculated range)

- :py:func:`radis.lbl.broadening._whiting`
- :py:func:`radis.lbl.broadening._whiting_jit` : precompiled version
- :py:meth:`radis.lbl.broadening.BroadenFactory._calc_broadening_HWHM`
- :py:meth:`radis.lbl.broadening.BroadenFactory._add_voigt_broadening_HWHM`
- :py:meth:`radis.lbl.broadening.BroadenFactory._voigt_broadening`
- :py:meth:`radis.lbl.broadening.BroadenFactory._calc_lineshape`
- :py:meth:`radis.lbl.broadening.BroadenFactory._apply_lineshape`
- :py:meth:`radis.lbl.broadening.BroadenFactory._calc_broadening`
- :py:meth:`radis.lbl.broadening.BroadenFactory._calc_broadening_noneq`
- :py:meth:`radis.lbl.broadening.BroadenFactory._find_weak_lines`
- :py:meth:`radis.lbl.broadening.BroadenFactory._calculate_pseudo_continuum`
- :py:meth:`radis.lbl.broadening.BroadenFactory._add_pseudo_continuum`


Notes
-----

Formula in docstrings generated with :py:func:`~pytexit.pytexit.py2tex` ::
    
    from pytexit import py2tex
    py2tex('...')

----------

"""
from __future__ import print_function, absolute_import, division, unicode_literals
from radis.lbl.base import BaseFactory
from radis.phys.constants import Na
from radis.phys.constants import k_b_CGS, c_CGS
from radis.misc.printer import printg
from radis.misc.basics import is_float
from radis.misc.arrays import is_sorted
from numpy import exp, arange, zeros_like, trapz, pi, sqrt, sin
from numpy import log as ln
from multiprocessing import Pool, cpu_count
from time import time
from warnings import warn
from radis.misc.progress_bar import ProgressBar
from radis.misc.warning import reset_warnings
import numpy as np
import matplotlib.pyplot as plt
import sys
from six.moves import zip
from numba import jit, float64
from radis.misc.debug import printdbg
from six.moves import range

# %% Broadening functions


def doppler_broadening_HWHM(wav, molar_mass, Tgas):
    """ Computes Gaussian (Doppler) broadening HWHM over all lines with [1]_, [2]_

    .. math::
        
        \\frac{w}{c} \\sqrt{\\frac{2N_a k_b T_{gas} \\ln 2}{M}}
        
    with ``k`` and ``c`` in CGS

    Parameters
    ----------

    wav: array like (nm / cm-1)
        transition waverange  [length N = number of lines]

    molar_mass: array like (g/mol)
        molar mass for isotope of given transition, in ``g/mol``   [length N]

    Tgas: float (K)
        (translational) gas temperature

    Returns
    -------

    hwhm_gauss: numpy array    (nm / cm-1)    [shape N]
        gaussian HWHM for all lines

    References
    ----------

    .. math::

        f_{G} = \\frac{1}{\\alpha} \\sqrt{\\frac{\\ln 2}{\\pi}} exp\\left(- \\ln 2 \\left(\\frac{w - w_0}{\\alpha}\\right)^2 \\right)

    with α the full-width half maximum (HWHM) calculated in ``cm-1``, ``nm`` 
    or ``Hz``:

    .. [1] `HITRAN.org <https://hitran.org/docs/definitions-and-units/>`_
           Eqn. (5)  [in cm-1]. ``c`` in CGS, HWHM \\alpha in ``nm``, M in ``g/mol``:
         
        .. math::
            
            \\alpha w_G= \\frac{w}{c_{CGS}} \\sqrt{\\frac{2N_a k_b T \\ln 2}{M}}
               
    .. [2] `Laux et al, 2003, "Optical diagnostics of atmospheric pressure air plasmas" <http://iopscience.iop.org/article/10.1088/0963-0252/12/2/301/meta>`_
           Eqn. (6)  [in nm]. FWHM \\Delta in ``nm``, M in ``g/mol``:
               
        .. math::
            
            \\alpha \\lambda_G=3.58\\times 10^{-07}  \\lambda \\sqrt{\\frac{T}{M}}

    .. [3] `Wikipedia <https://en.wikipedia.org/wiki/Doppler_broadening>`_ [in Hz]

        .. math::

            \\alpha f_G=\\frac{1}{2} \\sqrt{\\frac{8 k T \\ln 2}{m c^2}} f_0

    Here we work in ``cm-1`` with the CGS nomenclature (molar mass M in ``g/mol``)

    Notes
    -----
          
    *equation generated from the Python formula with* :py:func:`~pytexit.pytexit.py2tex`

    See Also
    --------
    
    :py:func:`~radis.lbl.broadening.gaussian_lineshape`

    """

    # Broadening parameter:
    # ... Doppler broadened half-width:
    # ... Reference: either Wikipedia (in cm-1), or former Sean's code (HITRAN based)
    # ... in CGS, or Laux 2003 (in nm). Both give the same results. Here we
    # ... use the CGS nomenclature (reminder: dg.molar_mass in g/mol)
    hwhm_gauss = (wav / c_CGS) * sqrt(
        (2 * Na * k_b_CGS * Tgas * ln(2)) / molar_mass
    )  # HWHM (cm-1)

    return hwhm_gauss  # HWHM (in wav unit)


def gaussian_lineshape(w_centered, hwhm):
    """ Computes Doppler (Gaussian) lineshape over all lines with [1]_, [2]_

    .. math::
        
        \\frac{1}{\\alpha_g} \\sqrt{\\frac{\\ln(2)}{\\pi}} 
        \\operatorname{exp}\\left(-\\ln 2 {\\left(\\frac{w_{centered}}{\\alpha_g}\\right)}^2\\right)

    Parameters
    ----------

    w_centered: 2D array       [one per line: shape W x N]
        waverange (nm / cm-1) (centered on 0, size W = broadening width size)
        
    hwhm:  array   [shape N = number of lines]
        Half-width at half-maximum (HWHM) of Gaussian

    Tgas: K
        (translational) gas temperature

    Returns
    -------

    lineshape: array          [shape N x W]
        line profile  

    References
    ----------

    $$ I_{gaussian} = exp(- ln2 (Δw/HWHM)^2 )) $$

    with HWHM full-width half maximum

    .. [1] `Wikipedia <https://en.wikipedia.org/wiki/Doppler_broadening>`_ (in cm-1)

    .. [2] `Laux et al, 2003, "Optical diagnostics of atmospheric pressure air plasmas" <http://iopscience.iop.org/article/10.1088/0963-0252/12/2/301/meta>`_
           (in nm)

    Both give the same results. 
    
    Notes
    -----
    
    *formula generated from the Python equation with* :py:func:`~pytexit.pytexit.py2tex`

    See Also
    --------
    
    :py:func:`~radis.lbl.broadening.doppler_broadening_HWHM`, 
    :py:func:`~radis.lbl.broadening.lorentzian_lineshape`, 
    :py:func:`~radis.lbl.broadening.voigt_lineshape`, 

    """

    # Calculate broadening
    # ------
    lineshape = 1 / hwhm * sqrt(ln(2) / pi) * exp(-ln(2) * (w_centered / hwhm) ** 2)

    return lineshape


def gaussian_FT(w_centered, hwhm):
    """ Fourier Transform of a Gaussian lineshape 
    
    Parameters
    ----------

    w_centered: 2D array       [one per line: shape W x N]
        waverange (nm / cm-1) (centered on 0)
        
    hwhm:  array   [shape N = number of lines]
        Half-width at half-maximum (HWHM) of Gaussian

    See Also
    --------
    
    :py:func:`~radis.lbl.broadneing.gaussian_lineshape`
    """

    # FT of w*np.sqrt(np.log(2)/np.pi)*np.exp(-np.log(2)*((v-v0)/w)**2)*dv
    n = len(w_centered)
    I = np.zeros(n)
    I[: n // 2 + n % 2] = np.exp(
        -((np.pi * w_centered[: n // 2 + n % 2] * hwhm) ** 2) / (np.log(2))
    )
    I[-(n // 2) :] = I[n // 2 : 0 : -1]
    return I


def pressure_broadening_HWHM(
    airbrd, selbrd, Tdpair, Tdpsel, pressure_atm, mole_fraction, Tgas, Tref
):
    """ Calculates collisional broadening HWHM over all lines by scaling 
    tabulated HWHM for new pressure and mole fractions conditions [1]_

    Note that collisional broadening is computed with a different coefficient
    for air broadening and self broadening. In the case of non-air mixtures,
    then the results should be used with caution, or better, the air broadening
    coefficient be replaced with the gas broadening coefficient

    Parameters
    ----------

    airbrd: array like    [length N]
        half-width half max coefficient (HWHM ) for collisional broadening with air

    selbrd: array like    [length N]
        half-width half max coefficient (HWHM ) for resonant (self) broadening with air

    Tdpair: array like     [length N]
        temperature dependance coefficient for collisional broadening with air

    Tdpsel: array like     [length N]
        temperature dependance coefficient for resonant (self) broadening

    pressure_atm: float  [atm]
        pressure in atmosphere (warning, not bar!)

    mole_fraction: float    [0-1]
        mole fraction

    Tgas: float [K]
        (translational) gas temperature

    Tref: float [K]
        reference temperature at which tabulated HWHM pressure 
        broadening coefficients were tabulated

    Returns
    -------

    lineshape: pandas Series        [shape N x W]
        Lorentzian half-width at half-maximum (FWHM) for each line profile  

    References
    ----------

    .. [1] `Rothman 1998 (HITRAN 1996) eq (A.14) <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`_

    See Also
    --------
    
    :py:func:`~radis.lbl.broadening.lorentzian_lineshape`
    
    """

    # Prepare coefficients, vectorize
    # --------

    # Temperature and pressure dependent half width
    # ... Reference: Rothman 1998 (HITRAN 1996) eq (A.12)
    # ... Hypothesis: we only consider self broadening, and air broadening,
    # ... weighted with mole fractions

    # Note: if not Tdpsel in dg Tdpair is used. Lookup parent function
    # | dev note: in that case we simplify the expression by calculation the
    # | power function once only.

    if Tdpsel is None:
        gamma_lb = ((Tref / Tgas) ** Tdpair) * (
            (airbrd * pressure_atm * (1 - mole_fraction))
            + (selbrd * pressure_atm * mole_fraction)
        )
    else:
        gamma_lb = ((Tref / Tgas) ** Tdpair) * (
            airbrd * pressure_atm * (1 - mole_fraction)
        ) + ((Tref / Tgas) ** Tdpsel) * (selbrd * pressure_atm * mole_fraction)

    return gamma_lb


def lorentzian_lineshape(w_centered, gamma_lb):
    """ Computes collisional broadening over all lines [1]_

    .. math:: 
        
        \\frac{1}{\\pi} \\frac{\\gamma_{lb}}{\\gamma_{lb}^2+w_{centered}^2}

    Parameters
    ----------

    w_centered: 2D array       [one per line: shape W x N]
        waverange (nm / cm-1) (centered on 0)
        
    gamma_lb: array   (cm-1)        [length N]
        half-width half maximum coefficient (HWHM) for pressure broadening 
        calculation

    Returns
    -------

    lineshape: array        [shape N x W]
        line profile  

    References
    ----------

    .. [1] `Rothman 1998 (HITRAN 1996) eq (A.14) <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`_

    Notes
    -----
    
    *formula generated from the Python equation with* :py:func:`~pytexit.pytexit.py2tex`

    See Also
    --------
    
    :py:func:`~radis.lbl.broadening.pressure_broadening_HWHM`, 
    :py:func:`~radis.lbl.broadening.gaussian_lineshape`, 
    :py:func:`~radis.lbl.broadening.voigt_lineshape`

    """

    # Calculate broadening
    # -------
    lineshape = 1 / pi * gamma_lb / ((gamma_lb ** 2) + (w_centered ** 2))

    return lineshape


def lorentzian_FT(w_centered, gamma_lb):
    """ Fourier Transform of a Lorentzian lineshape 
    
    Parameters
    ----------

    w_centered: 2D array       [one per line: shape W x N]
        waverange (nm / cm-1) (centered on 0)
        
    gamma_lb: array   (cm-1)        [length N]
        half-width half maximum coefficient (HWHM) for pressure broadening 
        calculation

    See Also
    --------
    
    :py:func:`~radis.lbl.broadneing.lorentzian_lineshape`
    """

    # FT of (1/np.pi) * w / ((v-v0)**2 + w**2)*dv
    n = len(w_centered)
    I = np.zeros(n)
    I[: n // 2 + n % 2] = np.exp(-np.pi * w_centered[: n // 2 + n % 2] * 2 * gamma_lb)
    I[-(n // 2) :] = I[n // 2 : 0 : -1]
    return I


def voigt_broadening_HWHM(
    airbrd,
    selbrd,
    Tdpair,
    Tdpsel,
    wav,
    molar_mass,
    pressure_atm,
    mole_fraction,
    Tgas,
    Tref,
):
    """ Calculate Voigt profile half-width at half-maximum (HWHM) from the Gaussian and 
    Collisional broadening with the empirical formula of Olivero [1]_

    Gaussian broadening is calculated with [2]_, Collisional broadening with [3]_

    Exact for a pure Gaussian and pure Lorentzian

    Parameters
    ----------

    Line parameters   [length N]
        
    airbrd: np.array  [length N]  (cm-1/atm)
        air broadening half-width half maximum (HWHM)

    selbrd:   np.array  [length N]  (cm-1/atm)
        self broadening half-width half maximum (HWHM)

    Tdpair:   np.array  [length N]  
        temperature dependance of collisional broadening 
        by air

    Tdpsel:   np.array  [length N]  
        temperature dependance of collisional self-broadening 

    wav:   np.array  [length N]    (cm-1)
        transition wavenumber

    molar_mass:   np.array  [length N]   (g/mol)
        molar mass for isotope of given transition 

    Environment parameters

    pressure_atm: float  [atm]
        pressure

    mole_fraction: float [0-1]
        mole fraction

    Tgas: K
        (translational) gas temperature

    Tref: K
        reference temperature at which tabulated HWHM pressure 
        broadening coefficients were tabulated

    Returns
    -------

    gamma_voigt, gamma_lb, gamma_db: numpy array
        Voigt, Lorentz, and Gaussian HWHM    

    References
    ----------
    
    Notes: Olivero [1] uses FWHM. 

    .. [1] `Olivero 1977 "Empirical fits to the Voigt line width: A brief review" <https://www.sciencedirect.com/science/article/pii/0022407377901613>`_
           Also found in NEQAIR96 Manual.

    .. [2] `Wikipedia <https://en.wikipedia.org/wiki/Doppler_broadening>`_ (in cm-1)

    .. [3] `Rothman 1998 (HITRAN 1996) eq (A.12) <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`_

    See Also
    --------
    
    :py:func:`~radis.lbl.broadening.olivero_1977`,
    :py:func:`~radis.lbl.broadening.pressure_broadening_HWHM`, 
    :py:func:`~radis.lbl.broadening.doppler_broadening_HWHM`, 
    :py:func:`~radis.lbl.broadening.voigt_lineshape`

    """

    # Collisional broadening HWHM
    gamma_lb = pressure_broadening_HWHM(
        airbrd, selbrd, Tdpair, Tdpsel, pressure_atm, mole_fraction, Tgas, Tref
    )

    # Doppler Broadening HWHM:
    gamma_db = doppler_broadening_HWHM(wav, molar_mass, Tgas)

    # Calculate broadening
    # -------
    # ... Reference Olivero 1977, also in NEQAIR 96 manual Eqn. (D2)
    # ... note that formula is given in wavelength (nm) [doesnt change anything]
    # ... and uses full width half maximum (FWHM)
    wg = 2 * gamma_db  # HWHM > FWHM
    wl = 2 * gamma_lb  # HWHM > FWHM

    return olivero_1977(wg, wl) / 2, gamma_lb, gamma_db  # FWHM > HWHM


def olivero_1977(wg, wl):
    """ Calculate approximate Voigt FWHM with Olivero 77, also in NEQAIR 96 manual Eqn. (D2)
    Note that formula is given in wavelength (nm) [doesnt change anything]
    and uses full width half maximum (FWHM) of Gaussian and Lorentzian profiles.
    
    For use with the Whiting formula of :py:func:`~radis.lbl.broadening.voigt_lineshape`. 
    
    Parameters
    ----------
    
    wg: numpy array
        Gaussian profile FWHM
        
    wl: numpy array
        Lorentzian profile FWHM
    
    Returns
    -------

    gamma_voigt: numpy array
        Voigt FWHM
        
    See Also
    --------
    
    :py:func:`~radis.lbl.broadening.voigt_broadening_HWHM`, 
    :py:func:`~radis.lbl.broadening.voigt_lineshape`
    """
    #    wv = wl/2 + sqrt((1/4*wl**2+wg**2))
    sd = (wl - wg) / (wl + wg)
    wv = (
        1
        - 0.18121 * (1 - sd ** 2)
        - (0.023665 * exp(0.6 * sd) + 0.00418 * exp(-1.9 * sd)) * sin(pi * sd)
    ) * (wl + wg)
    return wv


def voigt_lineshape(w_centered, hwhm_lorentz, hwhm_voigt, jit=True):
    """ Calculates Voigt lineshape using the approximation of the Voigt profile of
    Whiting [1]_, [2]_ that maintains a good accuracy in the far wings. Exact for a pure 
    Gaussian and pure Lorentzian

    Parameters
    ----------

    w_centered: 2D array       [one per line: shape W x N]
        waverange (nm / cm-1) (centered on 0)
        
    hwhm_lorentz: array   (cm-1)        [length N]
        half-width half maximum coefficient (HWHM) for Lorentzian broadening

    hwhm_voigt: array   (cm-1)        [length N]
        half-width half maximum coefficient (HWHM) for Voigt broadening,
        calculated by :py:func:`~radis.lbl.broadening.voigt_broadening_HWHM`

    Other Parameters
    ----------------
    
    jit: boolean
        if ``True``, use just in time compiler. Usually faster when > 10k lines.
        Default ``True``.

    Returns
    -------

    lineshape: pandas Series        [shape N x W]
        line profile  

    References
    ----------
    
    .. [1] `NEQAIR 1996 User Manual, Appendix D <https://ntrs.nasa.gov/search.jsp?R=19970004690>`_

    .. [2] `Whiting 1968 "An empirical approximation to the Voigt profile", JQSRT <https://www.sciencedirect.com/science/article/pii/0022407368900812>`_

    See Also
    --------
    
     :py:func:`~radis.lbl.broadening.voigt_broadening_HWHM`

    """

    # Note: Whiting and Olivero use FWHM. Here we keep HWHM in all public function
    # arguments for consistency.
    wl = 2 * hwhm_lorentz  # HWHM > FWHM
    wv = 2 * hwhm_voigt  # HWHM > FWHM

    if jit:
        lineshape = _whiting_jit(w_centered, wl, wv)
    else:
        lineshape = _whiting(w_centered, wl, wv)

    # Normalization
    #    integral = wv*(1.065+0.447*(wl/wv)+0.058*(wl/wv)**2)
    # ... approximation used by Whiting, equation (7)
    # ... performance: ~ 6µs vs ~84µs for np.trapz(lineshape, w_centered) ):
    # ... But not used because:
    # ... - it may yield wrong results when the broadening range is not refined enough
    # ... - it is defined for wavelengths only. Here we may have wavenumbers as well
    integral = np.trapz(lineshape, w_centered, axis=0)

    # Normalize
    lineshape /= integral

    assert not np.isnan(lineshape).any()

    return lineshape


# Pseudo-voigts approximations:


def _whiting(w_centered, wl, wv):
    """ 
    Parameters
    ----------
    
    wl: array
        Lorentzian FWHM
        
    wv: array
        Voigt FWHM
        
    w_centered: 2D array
        broadening spectral range for all lines
    
    Notes
    -----
    
    Performances:
        
    using @jit yield a performance increase from 8.9s down to 5.1s 
    on a 50k lines, 250k wavegrid case (performances.py) 
    """
    # Calculate some temporary arrays
    # ... fasten up the calculation by 25% (ex: test on 20 cm-1, ~6000 lines:
    # ... 20.5.s > 16.5s) on the total eq_spectrum calculation
    # ... w_wv is typically a (10.001, 1997) array
    w_wv = w_centered / wv  # w_centered can be ~500 Mb
    w_wv_2 = w_wv ** 2
    wl_wv = wl / wv
    w_wv_225 = np.abs(w_wv) ** 2.25

    # Calculate!  (>>> this is the performance bottleneck <<< : ~ 2/3 of the time spent
    #              on lineshape equation below + temp array calculation above
    #              In particular exp(...) and ()**2.25 are very expensive <<< )
    # ... Voigt 1st order approximation
    lineshape = (
        (1 - wl_wv) * exp(-2.772 * w_wv_2)
        + wl_wv * 1 / (1 + 4 * w_wv_2)
        # ... 2nd order correction
        + 0.016 * (1 - wl_wv) * wl_wv * (exp(-0.4 * w_wv_225) - 10 / (10 + w_wv_225))
    )
    return lineshape


@jit(
    float64[:, :](float64[:, :], float64[:, :], float64[:, :]),
    nopython=True,
    cache=True,
)  # , parallel=True)
def _whiting_jit(w_centered, wl, wv):
    """ 
    Parameters
    ----------
    
    wl: array
        Lorentzian FWHM
        
    wv: array
        Voigt FWHM
        
    w_centered: 2D array
        broadening spectral range for all lines
    
    Notes
    -----
    
    Performances:
        
    using @jit yield a performance increase from 8.9s down to 5.1s 
    on a 50k lines, 250k wavegrid case (performances.py) 
    """
    # Calculate some temporary arrays
    # ... fasten up the calculation by 25% (ex: test on 20 cm-1, ~6000 lines:
    # ... 20.5.s > 16.5s) on the total eq_spectrum calculation
    # ... w_wv is typically a (10.001, 1997) array
    w_wv = w_centered / wv  # w_centered can be ~500 Mb
    w_wv_2 = w_wv ** 2
    wl_wv = wl / wv
    w_wv_225 = np.abs(w_wv) ** 2.25

    # Calculate!  (>>> this is the performance bottleneck <<< : ~ 2/3 of the time spent
    #              on lineshape equation below + temp array calculation above
    #              In particular exp(...) and ()**2.25 are very expensive <<< )
    # ... Voigt 1st order approximation
    lineshape = (
        (1 - wl_wv) * exp(-2.772 * w_wv_2)
        + wl_wv * 1 / (1 + 4 * w_wv_2)
        # ... 2nd order correction
        + 0.016 * (1 - wl_wv) * wl_wv * (exp(-0.4 * w_wv_225) - 10 / (10 + w_wv_225))
    )
    return lineshape


# %% Tools


class BroadenFactory(BaseFactory):
    """ A class that holds all broadening methods, inherited by 
    :class:`~radis.lbl.factory.SpectrumFactory` eventually 
    
    .. inheritance-diagram:: radis.lbl.factory.SpectrumFactory
       :parts: 1
    
    See Also
    --------

    :class:`~radis.lbl.factory.SpectrumFactory`
    """

    def __init__(self):

        super(BroadenFactory, self).__init__()

        # Name variables (initialized later in SpectrumFactory)
        self.wbroad_centered = None

        self.wavenumber = None
        self.wavenumber_calc = None
        self.woutrange = None

        self._broadening_method = "voigt"
        """str: 'voigt', 'convolve', 'fft'
        
        Calculates broadening with a direct voigt approximation ('voigt') or
        by convoluting independantly calculated Doppler and collisional
        broadening ('convolve'). First is much faster, 2nd can be used to
        compare results. Not a user available parameter for the moment, but
        you can edit the SpectrumFactory manually::
            
            sf = SpectrumFactory(...)
            sf._broadening_method = 'voigt'
            
        Fast fourier transform ``'fft'`` is only available if using the DLM lineshape  
        calculation optimisation. Because the DLM convolves all lines at the same time,  
        and thus operates on large arrays, ``'fft'`` becomes more appropriate than
        convolutions in real space (``'voit'``, ``'convolve'`` )
        
        ``'fft'`` is automatically selected if DLM is used. 

        """

        # Predict broadening times (helps trigger warnings for optimization)
        self._broadening_time_ruleofthumb = 1e-7  # s / lines / point
        # Ex: broadening width of 10 cm-1  with resolution 0.01 cm-1:  1000 points
        # Ex: 250k lines: 1e-7 (s/lines/point) * 250k (lines) * 1000 (points) = 25s

    # %% ======================================================================
    # PRIVATE METHODS - BROADENING
    # (all computational-heavy functions: calculates all lines broadening,
    # convolve them, apply them on all calculated range)
    # ---------------------------------
    # _collisional_lineshape
    # _gaussian_lineshape
    # _calc_lineshape
    # plot_broadening  (public, debugging function)
    # _apply_lineshape
    # _broaden_lines
    # _broaden_lines_noneq
    # _calc_broadening
    # _calc_broadening_noneq
    #
    # XXX =====================================================================

    # %% Factory Broadenings methods:
    # Only collisional & gaussian at the moment, but it is fairly easy to add
    # new ones if needed (e.g: Stark...)

    # %% Functions to calculate broadening HWHM

    def _calc_broadening_HWHM(self):
        """ Calculate broadening HWHM and store in line dataframe (df1)

        Parameters
        ----------

        df: pandas Dataframe
            lines dataframe 

        Returns
        -------

        None 
            Dataframe self.df1 is updated

        Notes
        -----
        
        Called in :py:meth:`radis.lbl.factory.eq_spectrum`, :py:meth:`radis.lbl.factory.non_eq_spectrum`

        Run this method before using `_calc_lineshape`
        """

        # Init variables
        df = self.df1

        if len(df) == 0:
            return  # no lines

        if self.verbose >= 2:
            #            printg('> Calculate broadening HWHM')
            t0 = time()

        if self.input.Tgas is None:
            raise AttributeError(
                "Tgas not defined. Make sure the parent function creates it"
            )
        Tgas = self.input.Tgas
        pressure_mbar = self.input.pressure_mbar
        mole_fraction = self.input.mole_fraction
        # convert from mbar to atm for linebroadening calculation
        pressure_atm = pressure_mbar / 1013.25
        # coefficients tabulation temperature
        Tref = self.input.Tref

        # Get broadenings
        if self._broadening_method == "voigt":
            # Adds hwhm_voigt, hwhm_gauss, hwhm_lorentz:
            self._add_voigt_broadening_HWHM(df, pressure_atm, mole_fraction, Tgas, Tref)
        elif self._broadening_method in ["convolve", "fft"]:
            # Adds hwhm_lorentz:
            self._add_collisional_broadening_HWHM(
                df, pressure_atm, mole_fraction, Tgas, Tref
            )
            # Add hwhm_gauss:
            self._add_doppler_broadening_HWHM(df, Tgas)
        else:
            raise ValueError(
                "Unexpected broadening calculation method: {0}".format(
                    self._broadening_method
                )
            )

        if self.verbose >= 2:
            printg("Calculated broadening HWHM in {0:.2f}s".format(time() - t0))

    def _add_voigt_broadening_HWHM(self, df, pressure_atm, mole_fraction, Tgas, Tref):
        """ Update dataframe with Voigt HWHM

        Returns
        -------

        Note:
            But input pandas Dataframe ``'df'`` is updated with keys: 

            - ``hwhm_voigt``

            - ``hwhm_lorentz``

            - ``hwhm_gauss``

        """

        # Check self broadening is here
        if not "Tdpsel" in list(df.keys()):
            self.warn(
                "Self-broadening temperature coefficient Tdpsel not given in database: used Tdpair instead",
                "MissingSelfBroadeningWarning",
            )
            Tdpsel = None  # if None, voigt_broadening_HWHM uses df.Tdpair
        else:
            Tdpsel = df.Tdpsel

        # Calculate broadening FWHM
        wv, wl, wg = voigt_broadening_HWHM(
            df.airbrd,
            df.selbrd,
            df.Tdpair,
            Tdpsel,
            df.wav,
            df.molar_mass,
            pressure_atm,
            mole_fraction,
            Tgas,
            Tref,
        )

        # Update dataframe
        df["hwhm_voigt"] = wv
        df["hwhm_lorentz"] = wl
        df["hwhm_gauss"] = wg

        return

    def _add_collisional_broadening_HWHM(
        self, df, pressure_atm, mole_fraction, Tgas, Tref
    ):
        """ Update dataframe with collisional HWHM [1]_

        Returns
        -------

        Input pandas Dataframe ``df`` is updated with keys: 

            - hwhm_lorentz

        Notes
        -----

        Temperature and pressure dependent half width

        Hypothesis: we only consider self broadening, and air broadening,
        weighted with mole fractions

        If not ``Tdpsel``, ``Tdpair`` is used

        References
        ----------

        .. [1] `Rothman 1998 (HITRAN 1996) eq (A.12) <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`_

        """

        # Check self broadening is here
        if not "Tdpsel" in list(df.keys()):
            self.warn(
                "Self-broadening temperature coefficient Tdpsel not given in database: used Tdpair instead",
                "MissingSelfBroadeningWarning",
            )
            Tdpsel = None
        else:
            Tdpsel = df.Tdpsel

        # Calculate broadening HWHM
        wl = pressure_broadening_HWHM(
            df.airbrd,
            df.selbrd,
            df.Tdpair,
            Tdpsel,
            pressure_atm,
            mole_fraction,
            Tgas,
            Tref,
        )

        # Update dataframe
        df["hwhm_lorentz"] = wl

        return

    def _add_doppler_broadening_HWHM(self, df, Tgas):
        """ Update dataframe with Gaussian HWHM

        Returns
        -------

        None: input pandas Dataframe 'df' is updated with keys: 

            - ``hwhm_gauss``

        """

        # Calculate broadening HWHM
        wg = doppler_broadening_HWHM(df.wav, df.molar_mass, Tgas)
        # Note @EP: should we use the pressure-shifted wavenumber instead of df.wav?

        # Update dataframe
        df["hwhm_gauss"] = wg

        return

    def _collisional_lineshape(self, dg, wbroad_centered):
        """ Computes collisional broadening over all lines + normalize 
        and raise warnings if an error is detected

        Note that collisional broadening is computed with a different coefficient
        for air broadening and self broadening. In the case of non-air mixtures,
        then the results should be used with caution, or better, the air broadening
        coefficient be replaced with the gas broadening coefficient

        Parameters
        ----------

        dg: pandas Dataframe  (shape N)
            list of lines  (includes `gamma_lb` for broadening calculation)

        wbroad: array       [one per line: shape ``W`` x ``N``]
            wavenumbers (centered on 0)

        Returns
        -------

        pressure_lineshape: pandas Series
            line profile normalized with area = 1

        """

        # Get collisional broadening HWHM
        gamma_lb = dg.hwhm_lorentz

        # Prepare vectorized operations:
        try:  # make it a (1, N) row vector
            gamma_lb = gamma_lb.values.reshape((1, -1))
        except AttributeError:  # probably not a dataframe: assert there is one line only.
            assert type(gamma_lb) is np.float64

        # Calculate broadening for all lines
        # -------
        lineshape = lorentzian_lineshape(wbroad_centered, gamma_lb)

        # Normalize
        # ---------
        # ... 'wbroad_centered' is w_array-(w_shifted_line_center)
        # ... Normalization should not be needed as Lorentzian is normalized already
        # ... but that's only with a good enough wavestep
        # ... Here we compute the integral to check the error that is made
        area = trapz(lineshape.T, x=wbroad_centered.T)
        err = abs((area - 1))
        if self.warnings and self.warnings["CollisionalBroadeningWarning"] != "ignore":
            if (err > self.misc.warning_broadening_threshold).any():
                self.warn(
                    "Large error ({0:.1f}%) ".format(err.max() * 100)
                    + "in pressure broadening. Increase broadening width / reduce wstep. "
                    + "Use .plot_broadening() to visualize each line broadening",
                    "CollisionalBroadeningWarning",
                )
            # Note that although there may be an error here the total broadening
            # is normalized anyway, so the energy is conserved. If we're only
            # looking at a slit-function broadened spectrum (slit>>FWHM) it
            # wont change much. However it does impact multi-slabs calculations
        # ... Renormalize to dampen numeric errors impact
        lineshape /= area

        return lineshape

    def _gaussian_lineshape(self, dg, wbroad_centered):
        """ Computes Doppler (Gaussian) broadening over all lines + normalize 
        and raise warnings if an error is detected

        Parameters
        ----------

        dg: pandas Dataframe
            list of lines   (includes `gamma_db` for broadening calculation)

        wbroad_centered: array       [one per line: shape ``N`` x ``N``]
            wavenumber array (centered on 0, broadening width size)

        Returns
        -------

        gaussian_lineshape: pandas Series
            line profile normalized with area = 1

        """

        # Prepare coefficients, vectorize
        # --------

        # Broadening parameter:
        # ... Doppler broadened half-width:
        # ... Reference: either Wikipedia (in cm-1), or former Sean's code (HITRAN based)
        # ... in CGS, or Laux 2003 (in nm). Both give the same results. Here we
        # ... use the CGS nomenclature (reminder: dg.molar_mass in g/mol)
        gamma_db = dg.hwhm_gauss

        # Prepare vectorized operations:
        try:  # make it a (1,N) row vector
            gamma_db = gamma_db.values.reshape((1, -1))
        except AttributeError:  # probably not a dataframe: assert there is one line only.
            assert type(gamma_db) is np.float64

        # Calculate broadening for all lines
        # -------
        lineshape = gaussian_lineshape(wbroad_centered, gamma_db)

        # Normalize
        # ---------
        # ... normalisation not really needed as the analytical function is normalized
        # ... but that's only with a good enough wavestep
        # ... Here we compute the integral to check the error that is made
        area = trapz(lineshape.T, x=wbroad_centered.T)
        err = abs((area - 1))
        if self.warnings and self.warnings["GaussianBroadeningWarning"] != "ignore":
            # In a "performance" mode (vs "safe" mode), these warnings would be disabled
            if (err > self.misc.warning_broadening_threshold).any():
                self.warn(
                    "Large error ({0:.1f}%) ".format(err.max() * 100)
                    + "in Doppler broadening. Increase broadening width / reduce wstep. "
                    + "Use .plot_broadening() to visualize each line broadening",
                    "GaussianBroadeningWarning",
                )
                # Note that although there may be an error here the total broadening
                # is normalized anyway, so the energy is conserved. If we're only
                # looking at a slit-function broadened spectrum (slit>>FWHM) it
                # wont change much. However it does impact multi-slabs calculations
        lineshape /= area

        return lineshape

    def _voigt_broadening(self, dg, wbroad_centered, jit=True):
        """ Computes voigt broadening over all lines + normalize 

        Uses an approximation of the Voigt profile [1]_, [2]_ that maintains a 
        better accuracy in the far wings.

        Exact for a pure Gaussian and pure Lorentzian

        Parameters
        ----------

        dg: pandas Dataframe    [length ``N``]
            list of lines  (keys includes HWHM half-width half max coefficient
            `gamma_lb` for broadening calculation)

        w_centered: array      [length ``W``]
            wavenumbers (centered on 0)
    
        Other Parameters
        ----------------
        
        jit: boolean
            if ``True``, use just in time compiler. Usually faster when > 10k lines

        Returns
        -------

        lineshape: pandas Series        [shape ``N`` x ``W``]
            line profile  

        References
        --------

        .. [1] `NEQAIR 1996 User Manual, Appendix D <https://ntrs.nasa.gov/search.jsp?R=19970004690>`_

        .. [2] `Whiting 1968 "An empirical approximation to the Voigt profile", JQSRT <https://www.sciencedirect.com/science/article/pii/0022407368900812>`_

        """

        # Prepare coefficients, vectorize
        # --------
        # Get Voigt HWHM
        if not "hwhm_voigt" in list(dg.keys()):
            raise KeyError(
                "hwhm_voigt: Calculate broadening with "
                + "calc_voigt_broadening_HWHM first"
            )
        if not "hwhm_lorentz" in list(dg.keys()):
            raise KeyError(
                "hwhm_lorentz: Calculate broadening with "
                + "calc_voigt_broadening_HWHM first"
            )

        hwhm_lorentz = dg.hwhm_lorentz
        hwhm_voigt = dg.hwhm_voigt

        # Prepare vectorized operations:
        try:  # make it a (1, N) row vector
            hwhm_lorentz = hwhm_lorentz.values.reshape((1, -1))
            #        wg = wg.values.reshape((1,-1))
            hwhm_voigt = hwhm_voigt.values.reshape((1, -1))
        except AttributeError:  # probably not a dataframe: assert there is one line only.
            assert type(hwhm_lorentz) is np.float64
            #        assert type(wg) is np.float64
            assert type(hwhm_voigt) is np.float64

        # Calculate broadening for all lines
        # -------
        lineshape = voigt_lineshape(wbroad_centered, hwhm_lorentz, hwhm_voigt, jit=jit)

        return lineshape

    # %% Function to calculate lineshapes from HWHM

    def _calc_lineshape(self, dg):
        """ Sum over each line (trying to use vectorize operations to be faster)

        Parameters
        ----------

        dg: pandas Dataframe    [length ``N``]
            list of lines  (keys includes HWHM half-width half max coefficient
            `gamma_lb` for broadening calculation)


        Returns
        -------

        line_profile:   (1/cm-1)
            2D array of lineshapes (size B * N, B = lineshape width)

        Notes
        -----

        Implementation:

        - A broadening profile is calculated for each line, using the same
          wavenumber spacing as the calculated data. It is then applied to
          each line.

        - The broadening uses a `nearest` interpolation. It can induce an error
          if the wavenumber spacing is not small enough (I'd recommend ~ 10 wsteps
          per line FWHM). So far there is no automatic check that this criteria
          is applied.  

        Sizes of elements:

        - ``N``: number of lines in database part 'dg'  (only depends on the database)
            (typically: ~ 200 )
        - ``W``: size of wavelength array (the more the more accurate)
            (typically: ~ 10.000)
        - ``B``: restrained window over which we apply the convolution
            (typically: ~ 1000)

        Performance:

        - This function is the bottleneck in the whole Spectrum calculation. 
          Calculation is vectorized for most of the calculations (exp, `outer` ..), 
          but the execution of non vectorized functions, typically `argmin`,  roll` 
          and `convolve`, takes a long time.

        - Profiler test::

            non_eq_emission_spectrum    120s
            | _sum_line_by_line         84s
                | roll                  13s
                | argmin                30s
                | convolve              16s
                | outer                 8s

        See Also
        --------
        
        :py:meth:`~radis.lbl.broadening.BroadenFactory._apply_lineshape`
        
        """
        # TODO automatic wavenumber spacing: ~10 wsteps / FWHM

        if __debug__:
            t0 = time()

        # Init variables
        if self.input.Tgas is None:
            raise AttributeError(
                "Tgas not defined. Make sure the parent function creates it"
            )

        # Generate broadening array (so that it is as large as `broadening_max_width`
        # in cm-1, and keeps the same spacing as the final output wavelength vector)
        wbroad_centered_oneline = self.wbroad_centered  # size (B,)

        shifted_wavenum = dg.shiftwav
        try:  # make it a row vector
            shifted_wavenum = shifted_wavenum.values.reshape((1, -1))
            N = len(dg)
        except AttributeError:  # probably not a dataframe: one line only.
            assert type(shifted_wavenum) is np.float64
            N = 1

        # TODO (performance). Turn into a parabolic mesh? We need high resolution
        # around the peak but low is enough on the wings

        # matrix of absorption (shape W * N)
        wbroad_centered = np.outer(wbroad_centered_oneline, np.ones(N))
        wbroad = wbroad_centered + shifted_wavenum

        if __debug__:
            t1 = time()

        # Calculate lineshape (using precomputed HWHM)
        if self._broadening_method == "voigt":
            jit = True
            line_profile = self._voigt_broadening(dg, wbroad_centered, jit=jit)
        elif self._broadening_method == "convolve":
            # Get pressure and gaussian profiles
            pressure_profile = self._collisional_lineshape(dg, wbroad_centered)
            if __debug__:
                t11 = time()
            gaussian_profile = self._gaussian_lineshape(dg, wbroad_centered)
            if __debug__:
                t12 = time()

            # Convolve and get final line profile:
            line_profile = np.empty_like(pressure_profile)  # size (B, N)
            for i, (x, y) in enumerate(zip(pressure_profile.T, gaussian_profile.T)):
                line_profile[:, i] = np.convolve(x, y, "same")
            line_profile = line_profile / trapz(line_profile.T, x=wbroad.T)  # normalize
            # ... Note that normalization should not be needed as broadening profiles
            # ... are created normalized already. However, we do normalize to reduce
            # ... the impact of any error in line_profiles (due to wstep too big or
            # ... broadening_width too small): at least the energy is conserved, even
            # ... if not perfectly distributed (spectrally). A warning is raised by the
            # ... broadening functions.
        elif self._broadening_method == "fft":
            raise NotImplementedError("FFT")
        else:
            raise ValueError(
                "Unexpected broadening calculation method: {0}".format(
                    self._broadening_method
                )
            )

        if __debug__:
            t2 = time()
            if self.verbose >= 3:
                printg("... Initialized vectors in {0:.1f}s".format(t1 - t0))
                if self._broadening_method == "voigt":
                    printg(
                        "... Calculated Voigt profile (jit={1}) in {0:.1f}s".format(
                            t2 - t1, jit
                        )
                    )
                elif self._broadening_method == "convolve":
                    printg(
                        "... Calculated Lorentzian profile in {0:.1f}s".format(t11 - t1)
                    )
                    printg(
                        "... Calculated Gaussian profile in {0:.1f}s".format(t12 - t11)
                    )
                    printg("... Convolved both profiles in {0:.1f}s".format(t2 - t12))
                elif self._broadening_method == "fft":
                    raise NotImplementedError("FFT")

        return line_profile

    def _calc_lineshape_DLM(self, df):
        """ Generate the lineshape database using the steps defined by the 
        parameters :py:attr:`~radis.lbl.loader.Parameters.dlm_res_L` and 
        :py:attr:`~radis.lbl.loader.Parameters.dlm_res_G`.
        
        Parameters
        ----------
        
        df: pandas DataFrame
            line database
            
        Returns
        -------
        
        line_profile_DLM: dict
            dictionary of Voigt profile template. 
            If ``self._broadening_method == 'fft'``, templates are calculated
            in Fourier space.
    
        wL, wG: array
            Lorentzian and Gaussian FWHM in DLM
            
        wL_dat, wG_dat: array
            Lorentzian and Gaussian FWHM of data lines. 
        
        Reference
        ---------
        
        DLM implemented based on a code snippet from D.v.d.Bekerom.
        See: https://github.com/radis/radis/issues/37

        See Also
        --------
        
        :py:meth:`~radis.lbl.broadening.BroadenFactory._apply_lineshape_DLM`
        
        """

        if __debug__:
            t0 = time()

        # Prepare steps for Lineshape database
        # ------------------------------------

        def lorentzian_step(res_L):
            return (res_L / 0.20) ** 0.5

        def gaussian_step(res_G):
            return (res_G / 0.46) ** 0.5

        def init_w_axis(w_dat, w_step):
            f = 1 + w_step
            N = max(1, int(np.log(w_dat.max() / w_dat.min()) / np.log(f)) + 1) + 1
            return np.min(w_dat) * f ** np.arange(N)

        res_L = self.params.dlm_res_L  # DLM user params
        res_G = self.params.dlm_res_G  # DLM user params

        wL_dat = df.hwhm_lorentz.values * 2  # FWHM
        wG_dat = df.hwhm_gauss.values * 2  # FWHM

        wL_step = lorentzian_step(res_L)
        wG_step = gaussian_step(res_G)

        wL = init_w_axis(wL_dat, wL_step)  # FWHM
        wG = init_w_axis(wG_dat, wG_step)  # FWHM

        # Calculate the Lineshape
        # -----------------------

        line_profile_DLM = {}

        if self._broadening_method == "voigt":
            jit = False  # not enough lines to make the just-in-time FORTRAN compilation useful
            wbroad_centered = self.wbroad_centered

            # Non vectorized loop. Probably slightly slower, but this is not the bottleneck anyway.
            # see commit 6474cb7e on 15/08/2019 for a vectorized version
            for i in range(len(wL)):
                line_profile_DLM[i] = {}
                for j in range(len(wG)):
                    wV_ij = olivero_1977(wG[j], wL[i])  # FWHM
                    lineshape = voigt_lineshape(
                        wbroad_centered, wL[i] / 2, wV_ij / 2, jit=jit
                    )  # FWHM > HWHM
                    line_profile_DLM[i][j] = lineshape

        elif self._broadening_method == "convolve":
            wbroad_centered = self.wbroad_centered

            IL = [
                lorentzian_lineshape(wbroad_centered, wL[i] / 2) for i in range(len(wL))
            ]  # FWHM>HWHM
            IG = [
                gaussian_lineshape(wbroad_centered, wG[i] / 2) for i in range(len(wG))
            ]  # FWHM>HWHM
            # Non vectorized. See Voigt for vectorized.

            # Get all combinations of Voigt lineshapes
            for i in range(len(wL)):
                line_profile_DLM[i] = {}
                for j in range(len(wG)):
                    lineshape = np.convolve(IL[i], IG[j], mode="same")
                    lineshape /= np.trapz(lineshape, x=wbroad_centered)
                    line_profile_DLM[i][j] = lineshape

        elif self._broadening_method == "fft":
            # Unlike real space methods ('convolve', 'voigt'), here we calculate
            # the lineshape on the full spectral range.
            w = self.wavenumber_calc
            wstep = self.params.wstep
            w_lineshape_ft = np.arange(2 * len(w)) / ((2 * len(w)) * wstep)

            IL_FT = [
                lorentzian_FT(w_lineshape_ft, wL[i] / 2) for i in range(len(wL))
            ]  # FWHM>HWHM
            IG_FT = [
                gaussian_FT(w_lineshape_ft, wG[i] / 2) for i in range(len(wG))
            ]  # FWHM>HWHM

            # Get all combinations of Voigt lineshapes (in Fourier space)
            for i in range(len(wL)):
                line_profile_DLM[i] = {}
                for j in range(len(wG)):
                    line_profile_DLM[i][j] = IL_FT[i] * IG_FT[j]

        else:
            raise NotImplementedError(
                "Broadening method with DLM: {0}".format(self._broadening_method)
            )

        if __debug__ and self.verbose >= 3:
            printg(
                "... Precomputed DLM lineshapes ({1}) in {0:.1f}s".format(
                    time() - t0, len(wL) * len(wG)
                )
            )

        return line_profile_DLM, wL, wG, wL_dat, wG_dat

    def plot_broadening(self, i=0, pressure_atm=None, mole_fraction=None, Tgas=None):
        """ just for testing. Recalculate and plot broadening for line of index i

        Parameters
        ----------

        i: int
            line index

        pressure_atm: atm
            if None, defaults to model pressure

        mole_fraction: [0-1]
            if None, defaults to model mole fraction

        Tgas: K
            if None, defaults to model translational temperature

        """
        # TODO #clean: make it a standalone function.

        from publib import set_style, fix_style

        if pressure_atm is None:
            pressure_atm = self.input.pressure_mbar / 1013.25
        if mole_fraction is None:
            mole_fraction = self.input.mole_fraction
        if Tgas is None:
            Tgas = self.input.Tgas

        # Get one line only:
        dg = self.df1.iloc[i]

        wbroad_centered = self.wbroad_centered
        wbroad = wbroad_centered + dg.shiftwav
        #        convolve_profile = self._calc_lineshape(dg)
        # Get Voigt from empirical approximation
        if "hwhm_voigt" in dg:
            voigt_profile = self._voigt_broadening(dg, wbroad_centered, jit=False)
        # Get Voigt from convolution
        pressure_profile = self._collisional_lineshape(dg, wbroad_centered)
        gaussian_profile = self._gaussian_lineshape(dg, wbroad_centered)
        line_profile = np.convolve(pressure_profile, gaussian_profile, "same")
        line_profile /= trapz(line_profile.T, x=wbroad.T)  # normalize

        # Plot!
        set_style("origin")
        plt.figure()
        plt.plot(wbroad, pressure_profile, label="Pressure")
        plt.plot(wbroad, gaussian_profile, label="Doppler")
        plt.plot(wbroad, line_profile, label="Voigt (convolved)", lw=3)
        if "hwhm_voigt" in dg:
            plt.plot(wbroad, voigt_profile, label="Voigt (approximation)")
        plt.xlabel("Wavenumber (cm-1)")
        plt.ylabel("Broadening coefficient")
        plt.title("Line {0}, T={1:.1f}K, P={2:.2f}atm".format(i, Tgas, pressure_atm))
        plt.legend()
        fix_style("origin")

        return

    def _apply_lineshape(self, broadened_param, line_profile, shifted_wavenum):
        """ Multiply `broadened_param` by `line_profile` and project it on the
        correct wavelength given by `shifted_wavenum`

        Parameters
        ----------

        broadened_param: pandas Series (or numpy array)   [size N = number of lines]
            Series to apply lineshape to. Typically linestrength `S` for absorption,
            or `nu * Aul / 4pi * DeltaE` for emission

        line_profile:   (1/cm-1)        2D array of lines_profiles for all lines
                (size B * N, B = width of lineshape)

        shifted_wavenum: (cm-1)     pandas Series (size N = number of lines)
            center wavelength (used to project broaded lineshapes )

        Returns
        -------

        sumoflines: array (size W  = size of output wavenumbers)
            sum of (broadened_param x line_profile)

        Notes
        -----

        Units change during convolution::

            [sumoflines] = [broadened_param] * cm

        See Also
        --------
        
        :py:meth:`~radis.lbl.broadening.BroadenFactory._calc_lineshape`
        
        """

        if __debug__:
            t0 = time()

        #        # Get spectrum range
        wavenumber = self.wavenumber  # final vector of wavenumbers (shape W)
        wavenumber_calc = (
            self.wavenumber_calc
        )  # calculation vector of wavenumbers (shape W + space B on the sides)

        # Vectorize the chunk of lines
        S = broadened_param.reshape((1, -1))
        shifted_wavenum = shifted_wavenum.reshape((1, -1))  # make it a row vector

        # Get broadening array
        wbroad_centered = self.wbroad_centered  # size (B,)
        # index of broadening half width
        iwbroad_half = len(wbroad_centered) // 2

        # Calculate matrix of broadened parameter (for all lines)
        profile_S = line_profile * S

        # ---------------------------
        # Apply line profile

        if __debug__:
            t1 = time()

        # ... First get closest matching line (on the left, and on the right)
        # ... note @dev: wavenumber_calc must be sorted, which it is by construction.
        idcenter_left = (
            np.searchsorted(wavenumber_calc, shifted_wavenum.T, side="left").ravel() - 1
        )
        idcenter_right = np.minimum(idcenter_left + 1, len(wavenumber_calc) - 1)

        # ... Get the fraction of each line distributed to the left and to the right.
        frac_left = (
            shifted_wavenum - wavenumber_calc[idcenter_left]
        ).flatten()  # distance to left
        frac_right = (
            wavenumber_calc[idcenter_right] - shifted_wavenum
        ).flatten()  # distance to right
        dv = frac_left + frac_right
        frac_left, frac_right = frac_right / dv, frac_left / dv

        # ... Initialize array on which to distribute the lineshapes
        sumoflines_calc = zeros_like(wavenumber_calc)

        # Note on performance: it isn't straightforward to vectorize the summation
        # of all lineshapes on the spectral range as some lines may be parly outside
        # the spectral range.
        # to avoid an If / Else condition in the loop, we do a vectorized
        # comparison beforehand and run 3 different loops

        # reminder: wavenumber_calc has size [iwbroad_half+vec_length+iwbroad_half]
        vec_length = len(wavenumber)
        assert len(wavenumber_calc) == vec_length + 2 * iwbroad_half
        boffrangeleft = idcenter_left <= iwbroad_half
        boffrangeright = idcenter_right >= vec_length + iwbroad_half
        binrange = np.ones_like(idcenter_left, dtype=bool) ^ (
            boffrangeleft + boffrangeright
        )

        if __debug__:
            t2 = time()

        #        # Performance for lines below
        #        # ----------
        #        #
        #        # on test case: 6.5k lines x 18.6k grid length
        #        # normal: ~ 36 ms  called 9 times
        #        # with @jit : ~ 200 ms called 9 times (worse!)

        # In range: aggregate both wings
        lines_in = profile_S.T[binrange]
        if len(lines_in) > 0:
            I_low_in_left = idcenter_left[binrange] - iwbroad_half
            I_low_in_right = idcenter_right[binrange] - iwbroad_half
            I_high_in_left = I_low_in_left + 2 * iwbroad_half
            I_high_in_right = I_low_in_right + 2 * iwbroad_half
            for i, (fr_left, fr_right, profS) in enumerate(
                zip(frac_left[binrange], frac_right[binrange], lines_in)
            ):
                sumoflines_calc[I_low_in_left[i] : I_high_in_left[i] + 1] += (
                    fr_left * profS
                )
                sumoflines_calc[I_low_in_right[i] : I_high_in_right[i] + 1] += (
                    fr_right * profS
                )

        # Nomenclature for lines above:
        # - low/high: start/end of a lineshape
        # - left/right: closest spectral grid point on the left/right

        if __debug__:
            t3 = time()

        # Off Range, left : only aggregate the Right wing
        # @dev: the only difference with In range is the extra mask to cut the left wing.
        lines_l = profile_S.T[boffrangeleft]
        if len(lines_l) > 0:
            I_low_l_left = idcenter_left[boffrangeleft] + 1
            I_low_l_right = idcenter_right[boffrangeleft] + 1
            I_high_l_left = I_low_l_left + iwbroad_half - 1
            I_high_l_right = I_low_l_right + iwbroad_half - 1
            for i, (fr_left, fr_right, profS) in enumerate(
                zip(frac_left[boffrangeleft], frac_right[boffrangeleft], lines_l)
            ):
                # cut left wing & peak  with the [iwbroad_half+1:] mask
                sumoflines_calc[I_low_l_left[i] : I_high_l_left[i] + 1] += (
                    fr_left * profS[iwbroad_half + 1 :]
                )
                sumoflines_calc[I_low_l_right[i] : I_high_l_right[i] + 1] += (
                    fr_right * profS[iwbroad_half + 1 :]
                )

        # Off Range, Right : only aggregate the left wing
        lines_r = profile_S.T[boffrangeright]
        if len(lines_r) > 0:
            I_low_r_left = idcenter_left[boffrangeright] - iwbroad_half
            I_low_r_right = idcenter_right[boffrangeright] - iwbroad_half
            I_high_r_left = (
                I_low_r_left + iwbroad_half - 1
            )  # idcenter[boffrangeright]-1
            I_high_r_right = (
                I_low_r_right + iwbroad_half - 1
            )  # idcenter[boffrangeright]-1
            for i, (fr_left, fr_right, profS) in enumerate(
                zip(frac_left[boffrangeright], frac_right[boffrangeright], lines_r)
            ):
                # cut right wing & peak  with the [:iwbroad_half] mask
                sumoflines_calc[I_low_r_left[i] : I_high_r_left[i] + 1] += (
                    fr_left * profS[:iwbroad_half]
                )
                sumoflines_calc[I_low_r_right[i] : I_high_r_right[i] + 1] += (
                    fr_right * profS[:iwbroad_half]
                )

        if __debug__:
            t4 = time()

        # Get valid range (discard wings)
        sumoflines = sumoflines_calc[self.woutrange]

        if __debug__:
            if self.verbose >= 3:
                printg("... Initialized vectors in {0:.1f}s".format(t1 - t0))
                printg(
                    "... Get closest matching line & fraction in {0:.1f}s".format(
                        t2 - t1
                    )
                )
                printg("... Aggregate center lines in {0:.1f}s".format(t3 - t2))
                printg("... Aggregate wing lines in {0:.1f}s".format(t4 - t3))

        return wavenumber, sumoflines

    def _apply_lineshape_DLM(
        self, broadened_param, line_profile_DLM, shifted_wavenum, wL, wG, wL_dat, wG_dat
    ):
        """ Multiply `broadened_param` by `line_profile` and project it on the
        correct wavelength given by `shifted_wavenum`

        Parameters
        ----------

        broadened_param: pandas Series (or numpy array)   [size N = number of lines]
            Series to apply lineshape to. Typically linestrength `S` for absorption,
            or `nu * Aul / 4pi * DeltaE` for emission

        line_profile_DLM:  dict  
            dict of line profiles ::
                
                lineshape = line_profile_DLM[gaussian_index][lorentzian_index]

            If ``self._broadening_method == 'fft'``, templates are given
            in Fourier space.

        shifted_wavenum: (cm-1)     pandas Series (size N = number of lines)
            center wavelength (used to project broaded lineshapes )

        wL: array       (size DL)
            array of all Lorentzian widths in DLM

        wG: array       (size DG)
            array of all Gaussian widths in DLM

        wL_dat: array    (size N)
            FWHM of all lines. Used to lookup the DLM
            
        wG_dat: array    (size N)
            FWHM of all lines. Used to lookup the DLM
            
        Returns
        -------

        sumoflines: array (size W  = size of output wavenumbers)
            sum of (broadened_param x line_profile)

        Notes
        -----

        Units change during convolution::

            [sumoflines] = [broadened_param] * cm
            
        Reference
        ---------
        
        DLM implemented based on a code snippet from D.v.d.Bekerom.
        See: https://github.com/radis/radis/issues/37
            
        See Also
        --------
        
        :py:meth:`~radis.lbl.broadening.BroadenFactory._calc_lineshape_DLM`

        """

        if __debug__:
            t0 = time()

        #        # Get spectrum range
        wavenumber = self.wavenumber  # get vector of wavenumbers (shape W)
        wavenumber_calc = self.wavenumber_calc

        # Vectorize the chunk of lines
        S = broadened_param

        # ---------------------------
        # Apply line profile

        if __debug__:
            t1 = time()

        # ... First get closest matching spectral point  (on the left, and on the right)
        #         ... @dev: np.interp about 30% - 50% faster than np.searchsorted
        iv = np.interp(
            shifted_wavenum, wavenumber_calc, np.arange(len(wavenumber_calc))
        )
        if self._broadening_method == "fft":
            iv += len(wavenumber_calc) // 2  # FFT is done on 2x wavenumber_calc

        iv0 = iv.astype(int)  # size [N]
        iv1 = iv0 + 1

        # DLM: First we calculate the fractional index of the DLM
        #      that corresponds with this line:
        iwL = np.interp(wL_dat, wL, np.arange(len(wL)))
        iwG = np.interp(wG_dat, wG, np.arange(len(wG)))
        iwL0 = iwL.astype(int)  # size [N],   number of values defined by res_L
        iwG0 = iwG.astype(int)  # size [N],   number of values defined by res_G
        iwL1 = iwL0 + 1
        iwG1 = iwG0 + 1

        # DLM : Next calculate how the line is distributed over the 2x2x2 bins we have:
        av = iv - iv0
        awL = (iwL - iwL0) * (wL[iwL1] / wL_dat)
        awG = (iwG - iwG0) * (wG[iwG1] / wG_dat)
        # ... fractions on DLM grid
        awV00 = (1 - awL) * (1 - awG)
        awV10 = awL * (1 - awG)
        awV01 = (1 - awL) * awG
        awV11 = awL * awG

        Iv0 = S * (1 - av)
        Iv1 = S * av

        if __debug__:
            t2 = time()

        # ... Initialize array on which to distribute the lineshapes
        if self._broadening_method in ["voigt", "convolve"]:
            DLM = np.zeros((len(wavenumber_calc), len(wL), len(wG)))
        elif self._broadening_method == "fft":
            DLM = np.zeros((2 * len(wavenumber_calc), len(wL), len(wG)))
        else:
            raise NotImplementedError(self._broadening_method)

        # Distribute all line intensities on the 2x2x2 bins.
        np.add.at(DLM, (iv0, iwL0, iwG0), Iv0 * awV00)
        np.add.at(DLM, (iv1, iwL0, iwG0), Iv1 * awV00)
        np.add.at(DLM, (iv0, iwL1, iwG0), Iv0 * awV10)
        np.add.at(DLM, (iv1, iwL1, iwG0), Iv1 * awV10)
        np.add.at(DLM, (iv0, iwL0, iwG1), Iv0 * awV01)
        np.add.at(DLM, (iv1, iwL0, iwG1), Iv1 * awV01)
        np.add.at(DLM, (iv0, iwL1, iwG1), Iv0 * awV11)
        np.add.at(DLM, (iv1, iwL1, iwG1), Iv1 * awV11)

        # All lines within each bins are convolved with the same lineshape.
        # Let's do it:

        if __debug__:
            t21 = time()

        # For each value from the DLM, retrieve the lineshape and convolve all
        # corresponding lines with it before summing.
        if self._broadening_method in ["voigt", "convolve"]:

            # ... Initialize array on which to distribute the lineshapes
            sumoflines_calc = zeros_like(wavenumber_calc)

            for i in range(len(wL)):
                for j in range(len(wG)):
                    lineshape = line_profile_DLM[i][j]
                    sumoflines_calc += np.convolve(DLM[:, i, j], lineshape, "same")

        elif self._broadening_method == "fft":
            # ... Initialize array in FT space
            Idlm_FT = 1j * np.zeros(len(line_profile_DLM[0][0]))
            for i in range(len(wL)):
                for j in range(len(wG)):
                    lineshape_FT = line_profile_DLM[i][j]
                    Idlm_FT += np.fft.fft(np.fft.fftshift(DLM[:, i, j])) * lineshape_FT
            # Back in real space:
            sumoflines_calc = np.fft.ifftshift(np.fft.ifft(Idlm_FT).real)
            sumoflines_calc = sumoflines_calc[
                len(wavenumber_calc) // 2 : len(wavenumber_calc) // 2
                + len(wavenumber_calc)
            ]
            sumoflines_calc /= self.params.wstep

        else:
            raise NotImplementedError(self._broadening_method)

        if __debug__:
            t3 = time()

        # Get valid range (discard wings)
        sumoflines = sumoflines_calc[self.woutrange]

        if __debug__:
            if self.verbose >= 3:
                printg("... Initialized vectors in {0:.1f}s".format(t1 - t0))
                printg(
                    "... Get closest matching line & fraction in {0:.1f}s".format(
                        t2 - t1
                    )
                )
                printg("... Distribute lines over DLM {0:.1f}s".format(t21 - t2))
                printg(
                    "... Convolve and sum on spectral range {0:.1f}s".format(t3 - t21)
                )

        return wavenumber, sumoflines

    def _broaden_lines(self, df):
        """ Divide over chuncks not to process to many lines in memory at the
        same time (note that this is not where the parallelisation is done: all
        lines are processed on the same core. )
        
        Parameters
        ----------
        
        self: Factory
            contains the ``self.misc.chunksize`` parameter
            
        df: DataFrame
            line dataframe

        See _calc_lineshape for more information
        """
        # --------------------------

        # Reactivate warnings
        reset_warnings(self.warnings)

        # Init arrays
        wavenumber = self.wavenumber
        # Get number of groups for memory splitting
        chunksize = self.misc.chunksize
        # TEMP: @EP use chunksize as the parameter for the broadening optimization strategy:
        # - None
        # - if number: split the number of lines
        # - if 'DLM': use DLM strategy

        try:
            if chunksize is None:

                # Deal with all lines directly (usually faster)
                line_profile = self._calc_lineshape(df)  # usually the bottleneck
                (wavenumber, abscoeff) = self._apply_lineshape(
                    df.S.values, line_profile, df.shiftwav.values
                )

            elif chunksize == "DLM":
                # Use DLM

                line_profile_DLM, wL, wG, wL_dat, wG_dat = self._calc_lineshape_DLM(df)
                (wavenumber, abscoeff) = self._apply_lineshape_DLM(
                    df.S.values,
                    line_profile_DLM,
                    df.shiftwav.values,
                    wL,
                    wG,
                    wL_dat,
                    wG_dat,
                )

            elif is_float(chunksize):
                # Cut lines in smaller bits for better memory handling
                N = int(len(df) * len(wavenumber) / chunksize) + 1
                # Too big may be faster but overload memory.
                # See Performance for more information

                abscoeff = zeros_like(self.wavenumber)

                pb = ProgressBar(N, active=self.verbose)
                for i, (_, dg) in enumerate(df.groupby(arange(len(df)) % N)):
                    line_profile = self._calc_lineshape(dg)
                    (wavenumber, absorption) = self._apply_lineshape(
                        dg.S.values, line_profile, dg.shiftwav.values
                    )
                    abscoeff += absorption
                    pb.update(i)
                pb.done()

            else:
                raise ValueError(
                    "Unexpected value for chunksize: {0}".format(chunksize)
                )

        except MemoryError:
            import traceback

            traceback.print_exc()
            raise MemoryError(
                "Too many lines*wavepoints (see details above). Try to use or reduce the "
                + "chunksize parameter (current={0}{1})".format(
                    chunksize,
                    " so {0:.3e} lines*wavepoints was used".format(
                        len(df) * len(wavenumber)
                    )
                    if chunksize is None
                    else "",
                )
            )

        return wavenumber, abscoeff

    def _broaden_lines_noneq(self, df):
        """ Divide over chuncks not to process to many lines in memory at the
        same time (note that this is not where the parallelisation is done: all
        lines are processed on the same core)

        See _calc_lineshape for more information
        """

        # Reactivate warnings
        reset_warnings(self.warnings)

        # --------------------------

        # Init arrays
        wavenumber = self.wavenumber
        # Get number of groups for memory splitting
        chunksize = self.misc.chunksize

        try:
            if chunksize is None:
                # Deal with all lines directly (usually faster)
                line_profile = self._calc_lineshape(df)  # usually the bottleneck
                (wavenumber, abscoeff) = self._apply_lineshape(
                    df.S.values, line_profile, df.shiftwav.values
                )
                (_, emisscoeff) = self._apply_lineshape(
                    df.Ei.values, line_profile, df.shiftwav.values
                )

            elif chunksize == "DLM":
                # Use DLM

                line_profile_DLM, wL, wG, wL_dat, wG_dat = self._calc_lineshape_DLM(df)
                (wavenumber, abscoeff) = self._apply_lineshape_DLM(
                    df.S.values,
                    line_profile_DLM,
                    df.shiftwav.values,
                    wL,
                    wG,
                    wL_dat,
                    wG_dat,
                )
                (_, emisscoeff) = self._apply_lineshape_DLM(
                    df.Ei.values,
                    line_profile_DLM,
                    df.shiftwav.values,
                    wL,
                    wG,
                    wL_dat,
                    wG_dat,
                )
                # Note @dev: typical results is:
                # >>> abscoeff:
                # ... Precomputed DLM lineshapes in 0.0s
                # ... Initialized vectors in 0.0s
                # ... Get closest matching line & fraction in 0.3s
                # ... Distribute lines over DLM 2.1s
                # ... Convolve and sum on spectral range 0.4s
                # >>> emisscoeff
                # ... Initialized vectors in 0.0s
                # ... Get closest matching line & fraction in 0.2s
                # ... Distribute lines over DLM 2.1s
                # ... Convolve and sum on spectral range 0.3s
                # @EP: #performance.
                # unlike in the non DLM case, the nonequilibruum case here is ~2x
                # the equilibrium case: only the closest matching line is common to the
                # absorption & emission steps. The bottleneck is the distribution
                # of the line over the DLM, which has to be done for both abscoeff & emisscoeff.

            elif is_float(chunksize):
                # Cut lines in smaller bits for better memory handling

                # Get size of numpy array for vectorialization
                N = int(len(df) * len(wavenumber) / chunksize) + 1
                # Too big may be faster but overload memory.
                # See Performance for more information

                abscoeff = zeros_like(self.wavenumber)
                emisscoeff = zeros_like(self.wavenumber)

                pb = ProgressBar(N, active=self.verbose)
                for i, (_, dg) in enumerate(df.groupby(arange(len(df)) % N)):
                    line_profile = self._calc_lineshape(dg)
                    (wavenumber, absorption) = self._apply_lineshape(
                        dg.S.values, line_profile, dg.shiftwav.values
                    )
                    (_, emission) = self._apply_lineshape(
                        dg.Ei.values, line_profile, dg.shiftwav.values
                    )
                    abscoeff += absorption  #
                    emisscoeff += emission
                    pb.update(i)
                pb.done()

            else:
                raise ValueError(
                    "Unexpected value for chunksize: {0}".format(chunksize)
                )

        except MemoryError:
            import traceback

            traceback.print_exc()
            raise MemoryError(
                "See details above. Try to use or reduce the chunksize parameter"
            )

        return wavenumber, abscoeff, emisscoeff

    # %% Generate absorption profile which includes linebroadening factors

    def _calc_broadening(self):
        """
        Loop over all lines, calculate lineshape, and returns the sum of
        absorption coefficient k=S*f over all lines.

        For non-equilibrium, lineshape is calculated once and applied then
        to calculate absorption and emission coefficient.

        Returns
        -------

        abscoeff:  1/(#.cm-2)
            sum of all absorption coefficient k=1/(#.cm-2) for all lines in database
            `df` on the full calculation wavenumber range

        wavenumber: cm-1
            valid calculation wavenumber range

        Notes
        -----
        
        Units:

        - ``abscoeff`` and ``emisscoeff`` still have to be multiplied by the total
          number density (cm-3) to get (cm-1/#) unit.

        Performances:

        This function is the bottleneck of the whole process. However, a fully
        vectorized version is impossible because it takes too much memory. Instead
        I used chunks of variable data length, and ended up parallelizing it.

        - loop version:                 1'53''
        - vectorized, 10 chunks:        1'30''
        - vectorized, 100 chunks:       1'11''
        - vectorized, 1000 chunks:      1'20''
        - vectorized, 10.000 chunks:    2'02''

        - parallel process 8 cpus:      33''

        - fully vectorized w/o roll etc... < 5'

        And then the parallel dispatching overhead becomes comparable with the
        process time, so better not use parallel anymore.


        """

        parallel = self.misc.parallel
        df = self.df1

        if self.verbose >= 2:
            printg(
                "> Calculating line broadening ({0} lines: expect ~ {1:.2f}s on 1 CPU)".format(
                    len(df),
                    self._broadening_time_ruleofthumb
                    * len(df)
                    * len(self.wbroad_centered),
                )
            )
            t0 = time()

        # Just some tests
        try:
            assert len(df.shape) == 2
        except AssertionError:
            warn(
                "Dataframe has only one line. Unexpected behaviour could occur"
                + " because Dataframes will be handled as Series and row/columns"
                + " may be inverted"
            )

        if parallel:
            # Parallel version
            Ngroups = self.misc.Ngroups
            Nprocs = self.misc.Nprocs  # number of processes
            if self.verbose:
                print(
                    "Parallel version on {0} groups, {1} procs".format(Ngroups, Nprocs)
                )
            t1 = time()
            with Pool(Nprocs) as p:  # run on all procs
                # divide all lines in Ngroups
                ret_list = p.map(
                    self._broaden_lines,
                    [group for name, group in df.groupby(arange(len(df)) % Ngroups)],
                )
                # each process returns an array on the whole range, but only
                # including some of the lines
            if self.verbose:
                print("process: {0:.1f}s".format(time() - t1))

            w_arr, a_arr = list(zip(*ret_list))
            wavenumber = w_arr[0]  # they're all the same anyway

            abscoeff = np.array(a_arr).sum(axis=0)  # sum over all lines

        else:
            # Regular version
            (wavenumber, abscoeff) = self._broaden_lines(df)

        if self.verbose >= 2:
            printg("Calculated line broadening in {0:.2f}s".format(time() - t0))

        return wavenumber, abscoeff

    def _calc_broadening_noneq(self):
        """
        Loop over all lines, calculate lineshape, and returns the sum of
        absorption coefficient k=S*f over all lines.

        For non-equilibrium, lineshape is calculated once and applied then
        to calculate absorption and emission coefficient.

        Returns
        -------

        wavenumber: cm-1
            full calculation wavenumber range

        abscoeff:  1/(#.cm-2)
            sum of all absorption coefficient k=1/(#.cm-2) for all lines in database
            `df` on the full calculation wavenumber range

        emisscoeff:  W/sr.cm
            sum of all broadened emission coefficients

        Notes
        -----
        
        Units:

        - Both `abscoeff` and `emisscoeff` still have to be multiplied by the total
          number density (cm-3).

        """

        parallel = self.misc.parallel
        df = self.df1

        if self.verbose >= 2:
            printg(
                "Calculating line broadening ({0:,d} lines: expect ~ {1:.2f}s on 1 CPU)".format(
                    len(df),
                    self._broadening_time_ruleofthumb
                    * len(df)
                    * len(self.wbroad_centered),
                )
            )
            t0 = time()

        # Just some tests
        try:
            assert len(df.shape) == 2
        except AssertionError:
            warn(
                "Dataframe has only one line. Unexpected behaviour could occur"
                + " because Dataframes will be handled as Series and row/columns"
                + " may be inverted"
            )

        if parallel:
            # Parallel version
            Ngroups = cpu_count()
            if self.verbose:
                print(
                    "Parallel version on {0} groups, {1} procs".format(
                        Ngroups, cpu_count()
                    )
                )
            with Pool(cpu_count()) as p:  # run on all procs
                # divide all lines in Ngroups
                ret_list = p.map(
                    self._broaden_lines_noneq,
                    [group for name, group in df.groupby(arange(len(df)) % Ngroups)],
                )
                # each process returns an array on the whole range, but only
                # including some of the lines
            w_arr, a_arr, e_arr = list(zip(*ret_list))
            wavenumber = w_arr[0]
            abscoeff = np.array(a_arr).sum(axis=0)  # sum over all lines
            emisscoeff = np.array(e_arr).sum(axis=0)  # sum over all lines

        else:
            # Regular version
            (wavenumber, abscoeff, emisscoeff) = self._broaden_lines_noneq(df)

        if self.verbose >= 2:
            printg("Calculated line broadening in {0:.2f}s".format(time() - t0))

        return wavenumber, abscoeff, emisscoeff

    # %% Functions to calculate semi-continuum

    def _find_weak_lines(self, weak_rel_intensity_threshold):
        """ Finds weak lines in current line dataframe. 

        Lines are considered weak if recovered by a strong line nearby. These 
        lines are later moved in a semi-continuum, and only strong lines are 
        fully resolved with Voigt broadening (which is costly!)

        To find weak planes, we perform a fast broadening of all lines on a 
        rectangle of same FWHM as the lines. Then, a rough spectrum of linestrengths
        is calculated, and lines linestrength are compared to the rough spectrum 

        This procedure allows to discard many weak lines while preserving the 
        spectrum main features in the less intense parts, what an absolute
        cutoff such as the linestrength 'cutoff' cannot do

        Weak line criteria: "average linestrength S much smaller (parameter alpha) 
        than the approximate spectrum I without its own contribution": 
            
        .. math::
            
            S_{avg} < \\alpha \\times (I - S_{avg})
        

        Returns
        -------

        None
            store weak line status in dataframe as ``self.df1.weak_line``

        """

        # Get inputs
        wavenumber_calc = self.wavenumber_calc  # size W
        wstep = self.params.wstep
        df = self.df1  # lines already scaled with current temperature, size N

        if self.verbose >= 2:
            printg("... classifying lines as weak or strong")
            t0 = time()

        # Get approximate spectral absorption coefficient
        rough_spectrum, S_density_on_grid, line2grid_proj_left = project_lines_on_grid(
            df, wavenumber_calc, wstep
        )

        #     :
        # ~ 1/(#.cm-2)
        # Sizes:
        # - rough_spectrum:     size W
        # - S_density_on_grid:  size N
        # - line2grid_proj:     size N

        # Weak line criteria
        # ... Compare line density (in 1/(#.cm-2) to sum)
        line_is_weak = S_density_on_grid < (
            weak_rel_intensity_threshold
            * (rough_spectrum[line2grid_proj_left] - S_density_on_grid)
        )  # size N

        #        # DEBUG: plot weak lines and strong lines
        #        plt.figure()
        #        plt.plot(wavenumber_calc, rough_spectrum)
        #        for i, w in enumerate(df.shiftwav):
        #            ls = '--' if line_is_weak[i] else '-'
        #            plt.axvline(w, color='k', ls=ls)

        # ... Store weak line label in df
        df["weak_line"] = line_is_weak

        if self.verbose >= 2:
            printg(
                "... {0:,d} lines classified as weak lines ({1:.2f}%) in {2:.1f}s".format(
                    line_is_weak.sum(),
                    line_is_weak.sum() / len(line_is_weak) * 100,
                    time() - t0,
                )
            )

        return

    def _calculate_pseudo_continuum(self, noneq=False):
        """ Find weak lines, add them in pseudo-continuum  (note that pseudo-continuum
        by RADIS definition is actually more of sum of low-resolution lines)
        
        Parameters
        ----------
        
        noneq: bool
            if ``True``, also returns the emisscoeff pseudo continuum (for noneq
            cases). Default ``False``
        

        Returns
        -------

        k_continuum: numpy array    (1/(#.cm-2))
            abscoeff semi-continuum  on wavenumber space

        if noneq:
            
        j_continuum: numpy array  
            emisscoeff semi-continuum  on wavenumber space


        Also:
            self.df1 is updated, lines are removed
            local variables self._Nlines_in_continuum and self._Nlines_calculated
            are created


        Notes
        -----
        
        continuum can be exported in Spectrum is using the ``self.export_continuum``
        boolean. This is not available as a User input but can be edited manually
        in the Factory.
        
        The Weak line characterization [1]_ is only based on abscoeff. For strong 
        nonequilibrium cases there may be lines consired as weak in terms 
        of absorption but not weaks in emission, or the other way around. It should
        be negligible, though, so a dual conditioning (looking at both abscoeff
        and emisscoeff) was not implemented. 
        See :func:`radis.lbl.broadening._find_weak_lines` if you want to change that
        
        Reference
        ---------
        
        .. [1] `RADIS User Guide, RADIS Paper`

        # TODO: export continuum in Spectrum ? (under q['continuum'] ? )

        """

        if self.params.pseudo_continuum_threshold > 0:

            if self.verbose >= 2:
                printg("Calculating pseudo continuum")
            t0 = time()

            # Check inputs
            wavenumber_calc = self.wavenumber_calc
            pseudo_continuum_threshold = self.params.pseudo_continuum_threshold
            wstep = self.params.wstep

            # Calculate rough spectrum, label weak lines
            # ... only guess based on abscoeff. See Notes for noneq case.
            self._find_weak_lines(pseudo_continuum_threshold)

            # Retrieve lines, separate weaks from strong
            df = self.df1
            df_weak_lines = df[df.weak_line]
            df_strong_lines = df[~df.weak_line]

            if not len(df_strong_lines) > 0:
                raise ValueError(
                    "All lines qualified as weak: reduce weak_line_threshold"
                )

            # Calculate continuum
            if noneq:
                k_continuum, j_continuum, _, _, _ = project_lines_on_grid_noneq(
                    df_weak_lines, wavenumber_calc, wstep
                )

                if __debug__:
                    printdbg(
                        "Intensity of k continuum: {0}\n".format(
                            np.trapz(k_continuum, wavenumber_calc)
                        )
                        + "Intensity of lines removed: {0}".format(
                            df_weak_lines.S.sum()
                        )
                    )
                    printdbg(
                        "Intensity of j continuum: {0}\n".format(
                            np.trapz(j_continuum, wavenumber_calc)
                        )
                        + "Intensity of lines removed: {0}".format(
                            df_weak_lines.Ei.sum()
                        )
                    )

            else:
                k_continuum, _, _ = project_lines_on_grid(
                    df_weak_lines, wavenumber_calc, wstep
                )

                if __debug__:
                    printdbg(
                        "Intensity of continuum: {0}\n".format(
                            np.trapz(k_continuum, wavenumber_calc)
                        )
                        + "Intensity of lines removed: {0}".format(
                            df_weak_lines.S.sum()
                        )
                    )

            # Get valid range (discard wings)
            k_continuum = k_continuum[self.woutrange]  # 1/(#.cm-2)
            #            self.continuum_k = k_continuum

            if noneq:
                j_continuum = j_continuum[self.woutrange]  # 1/(#.cm-2)

            # Reduce line dataset to strong lines only
            self.df1 = df_strong_lines

            # Update number of lines

            self._Nlines_in_continuum = len(df_weak_lines)
            self._Nlines_calculated = len(self.df1)

            # Check performances
            time_spent = time() - t0
            # ... Expected broadening time gain (see Rule of Thumb)
            expected_broadening_time_gain = (
                self._broadening_time_ruleofthumb
                * self._Nlines_in_continuum
                * len(self.wbroad_centered)
            )
            if self.verbose >= 2:
                printg(
                    "Calculated pseudo-continuum in {0:.1f}s (expected time saved: {1:.1f}s)".format(
                        time_spent, expected_broadening_time_gain
                    )
                )
                # Add a warning if it looks like it wasnt worth it
            if time_spent > 3 * expected_broadening_time_gain:
                self.warn(
                    "Pseudo-continuum may not be adapted to this kind "
                    + "of spectrum. Time spent on continuum calculation "
                    + "({0:.1f}s) is much longer than expected gain ({1:.1f}s). ".format(
                        time_spent, expected_broadening_time_gain
                    )
                    + "If the calculation takes a lot of time, consider "
                    + "setting pseudo_continuum_threshold=0. If it is fast "
                    + "enough already, you can decrease the linestrength cutoff= "
                    + "parameter to add discarded weak lines to continuum "
                    + "and improve the calculation accuracy.",
                    "PerformanceWarning",
                )

        else:
            k_continuum = None
            self._Nlines_in_continuum = 0
            self._Nlines_calculated = len(self.df1)

            if noneq:
                j_continuum = None

        if noneq:
            return k_continuum, j_continuum
        else:
            return k_continuum

    def _add_pseudo_continuum(self, abscoeff_v, k_continuum):
        """
        Notes
        -----
        
        also used for adding emisscoeff continuum with::
            
            self._add_pseudo_continuum(emisscoeff_v, j_continuum):
        """
        if k_continuum is not None:
            abscoeff_v += k_continuum
        return abscoeff_v


def project_lines_on_grid(df, wavenumber, wstep):
    """ Quickly sums all lines on wavespace grid as rectangles of HWHM corresponding
    to hwhm_voigt and a spectral absorption coefficient value so that linestrength 
    is conserved
    
    i.e. profiles are approximated as a rectangle of width Alpha*FWHM_Voigt, 
    and same linestrength.
    
    Parameters
    ----------
    
    df: pandas Dataframe
        Contains ``shiftwav`` (wavenumbers) and ``S`` (linestrengths) and ``hwhm_voigt`` 
        (Voigt HWHM) size ``N`` (number of lines)
        
    wavenumber: np.array
        spectral grid. Size ``W``. Expected to be regular
        
    wstep: float  (cm-1)
        wavenumber step
        
    Returns
    -------
    
    k_rough_spectrum: np.array
        spectral absorption coefficient for the waverange ``wavenumber``, 
        calculated by assuming a rectangular profile for each line. Size ``W`` 

    S_density_on_grid
        average spectral linestrength intensity of each line (abscoeff ~k), assuming 
        rectangular profile. size ``N`` 

    line2grid_projection
        closest index of the center of each line in ``df`` on the spectral grid 
        ``wavenumber``. Size ``N`` 
        
    """

    shiftwav = df.shiftwav.values  # cm-1  ,   size N (number of lines)
    S = df.S.values  # cm/#  ~   cm-1/(#.cm-2)  ,   size N
    wv = df.hwhm_voigt.values * 2  # HWHM > FWHM

    # ... First get closest matching line (left, and right):
    iwav_on_grid_left = np.searchsorted(wavenumber, shiftwav.T, side="left").ravel() - 1
    iwav_on_grid_right = np.minimum(iwav_on_grid_left + 1, len(wavenumber) - 1)

    # ... Get the fraction of each line distributed to the left and to the right.
    frac_left = (shiftwav - wavenumber[iwav_on_grid_left]).flatten()  # distance to left
    frac_right = (
        wavenumber[iwav_on_grid_right] - shiftwav
    ).flatten()  # distance to right
    dv = frac_left + frac_right
    frac_left, frac_right = frac_right / dv, frac_left / dv

    # ... express FWHM (in nm) in index
    ALPHA = 2  # arbitrary.
    ihwhm_on_grid = np.asarray(ALPHA * wv / 2 // wstep, dtype=np.int64)
    ifwhm_on_grid = ihwhm_on_grid * 2 + 1  # make it odd (conserves center)
    # ... infer min and max index to project lines (sides may be out of range)
    imin_broadened_wav_on_grid_left = iwav_on_grid_left - ihwhm_on_grid
    imax_broadened_wav_on_grid_left = iwav_on_grid_left + ihwhm_on_grid
    imin_broadened_wav_on_grid_right = iwav_on_grid_right - ihwhm_on_grid
    imax_broadened_wav_on_grid_right = iwav_on_grid_right + ihwhm_on_grid

    # Get average intensity, assuming a rectangular profile of width FWHM
    S_density_on_grid = S / (ifwhm_on_grid * wstep)  # ~ 1/(#.cm-2)

    # Cut out of range points
    # ... @dev: you should see imin_broadened_wav_on_grid as a projection array, with
    # ... indexes where Intensity will be summed afterwards.
    # ... here we project the out of range intensities to index -1 and len_grid+1
    # ... (this is faster than )
    len_grid = len(wavenumber)
    imin_broadened_wav_offset_left = imin_broadened_wav_on_grid_left
    imax_broadened_wav_offset_left = imax_broadened_wav_on_grid_left
    imin_broadened_wav_offset_right = imin_broadened_wav_on_grid_right
    imax_broadened_wav_offset_right = imax_broadened_wav_on_grid_right
    imin_broadened_wav_offset_left[imin_broadened_wav_offset_left < 0] = -1
    imin_broadened_wav_offset_right[imin_broadened_wav_offset_right < 0] = -1
    imax_broadened_wav_offset_left[imax_broadened_wav_offset_left > len_grid] = len_grid
    imax_broadened_wav_offset_right[
        imax_broadened_wav_offset_right > len_grid
    ] = len_grid
    imin_broadened_wav_offset_left += 1
    imax_broadened_wav_offset_left += 1
    imin_broadened_wav_offset_right += 1
    imax_broadened_wav_offset_right += 1

    line2grid_projection_left = iwav_on_grid_left  # size N (number of lines)
    line2grid_projection_right = iwav_on_grid_right  # size N (number of lines)
    # ... deal with case where a new point is created (we dont want that)
    # ... just offset it by one unit (we're working with wavenumber_calc anyway,
    # ... these boundary points will be cropped by RADIS at the end of the calculation)
    line2grid_projection_left[line2grid_projection_left == len(wavenumber)] = (
        len(wavenumber) - 1
    )
    line2grid_projection_right[line2grid_projection_right == len(wavenumber)] = (
        len(wavenumber) - 1
    )

    @jit(float64[:](), nopython=True)
    def rough_sum_on_grid():
        """ Sum all lines linestrength density on a spectrum  """
        # Performance
        # ----------
        #
        # Test case:  6.5k lines, 18k grid points
        # - standard: 13.6 ms
        # - with @jit: 0.9 ms
        #
        # Array size
        # ----------
        # rough_spectrum has len_grid+2 because two points are used for out of
        # range intensities. We crop at the end
        rough_spectrum = np.zeros(len_grid + 2)
        for i, (fr_left, fr_right, Iline_density) in enumerate(
            zip(frac_left, frac_right, S_density_on_grid)
        ):
            imin_left = imin_broadened_wav_offset_left[i]
            imax_left = imax_broadened_wav_offset_left[i]
            imin_right = imin_broadened_wav_offset_right[i]
            imax_right = imax_broadened_wav_offset_right[i]
            rough_spectrum[imin_left : imax_left + 1] += fr_left * Iline_density
            rough_spectrum[imin_right : imax_right + 1] += fr_right * Iline_density

        # Nomenclature for lines above:
        # - min/max: start/end of a lineshape
        # - left/right: closest spectral grid point on the left/right

        # crop out of range points
        return rough_spectrum[1:-1]

    # TODO: @dev #performance
    # ... try with k_rough_spectrum.add.at()  ?
    # ... but it probably wont ever be comparable with DLM.

    k_rough_spectrum = rough_sum_on_grid()

    return k_rough_spectrum, S_density_on_grid, line2grid_projection_left


def project_lines_on_grid_noneq(df, wavenumber, wstep):
    """ Quickly sums all lines on wavespace grid as rectangles of HWHM corresponding
    to hwhm_voigt and a spectral absorption coefficient value so that linestrength 
    is conserved
    
    i.e. profiles are approximated as a rectangle of width Alpha*FWHM_Voigt, 
    and same linestrength.
    
    Parameters
    ----------
    
    df: pandas Dataframe
        Contains ``shiftwav`` (wavenumbers) and ``S`` (linestrengths) and ``hwhm_voigt`` 
        (Voigt HWHM) size ``N`` (number of lines)
        
    wavenumber: np.array
        spectral grid. Size ``W``. Expected to be regular
        
    wstep: float  (cm-1)
        wavenumber step
        
    quantity: 'S' or 'Ei'
        use 'S' for Linestrength, 'Ei' for emission integral. Default 'S'
        
    Returns
    -------
    
    k_rough_spectrum: np.array
        spectral absorption coefficient for the waverange ``wavenumber``, 
        calculated by assuming a rectangular profile for each line. Size ``W``
    
    j_rough_spectrum: np.array
        spectral emission coefficient for the waverange ``wavenumber``, 
        calculated by assuming a rectangular profile for each line. Size ``W``
    
    S_density_on_grid
        average spectral linestrength intensity of each line (abscoeff ~k), assuming 
        rectangular profile. size ``N`` 
        
    Ei_density_on_grid
        average spectral emission intensity of each line (emisscoeff ~j), assuming 
        rectangular profile. size ``N`` 
    
    line2grid_projection
        closest index of the center of each line in ``df`` on the spectral grid 
        ``wavenumber``. Size ``N`` 
        
    Notes
    -----
    
    Similar to :py:func:`~radis.lbl.broadening.project_lines_on_grid` except
    that we also calculate the approximate emission intensity (noneq > it cannot 
    be recomputed from the linestrength). 
        
    """

    shiftwav = df.shiftwav.values  # cm-1  ,   size N (number of lines)
    S = df.S.values  # cm/#  ~   cm-1/(#.cm-2)  ,   size N
    Ei = df.Ei.values  # mW/cm3/sr
    wv = df.hwhm_voigt.values * 2  # HWHM > FWHM

    # ... First get closest matching line (left, and right):
    iwav_on_grid_left = np.searchsorted(wavenumber, shiftwav.T, side="left").ravel() - 1
    iwav_on_grid_right = np.minimum(iwav_on_grid_left + 1, len(wavenumber) - 1)

    # ... Get the fraction of each line distributed to the left and to the right.
    frac_left = (shiftwav - wavenumber[iwav_on_grid_left]).flatten()  # distance to left
    frac_right = (
        wavenumber[iwav_on_grid_right] - shiftwav
    ).flatten()  # distance to right
    dv = frac_left + frac_right
    frac_left, frac_right = frac_right / dv, frac_left / dv

    # ... express FWHM (in nm) in index
    ALPHA = 2  # arbitrary.
    ihwhm_on_grid = np.asarray(ALPHA * wv / 2 // wstep, dtype=np.int64)
    ifwhm_on_grid = ihwhm_on_grid * 2 + 1  # make it odd (conserves center)
    # ... infer min and max index to project lines (sides may be out of range)
    imin_broadened_wav_on_grid_left = iwav_on_grid_left - ihwhm_on_grid
    imax_broadened_wav_on_grid_left = iwav_on_grid_left + ihwhm_on_grid
    imin_broadened_wav_on_grid_right = iwav_on_grid_right - ihwhm_on_grid
    imax_broadened_wav_on_grid_right = iwav_on_grid_right + ihwhm_on_grid

    # Get average intensity, assuming a rectangular profile of width FWHM
    S_density_on_grid = S / (ifwhm_on_grid * wstep)  # ~ 1/(#.cm-2)
    Ei_density_on_grid = Ei / (ifwhm_on_grid * wstep)  # mW/cm3/sr/cm-1

    # Cut out of range points
    # ... @dev: you should see imin_broadened_wav_on_grid as a projection array, with
    # ... indexes where Intensity will be summed afterwards.
    # ... here we project the out of range intensities to index -1 and len_grid+1
    # ... (this is faster than )
    len_grid = len(wavenumber)
    imin_broadened_wav_offset_left = imin_broadened_wav_on_grid_left
    imax_broadened_wav_offset_left = imax_broadened_wav_on_grid_left
    imin_broadened_wav_offset_right = imin_broadened_wav_on_grid_right
    imax_broadened_wav_offset_right = imax_broadened_wav_on_grid_right
    imin_broadened_wav_offset_left[imin_broadened_wav_offset_left < 0] = -1
    imin_broadened_wav_offset_right[imin_broadened_wav_offset_right < 0] = -1
    imax_broadened_wav_offset_left[imax_broadened_wav_offset_left > len_grid] = len_grid
    imax_broadened_wav_offset_right[
        imax_broadened_wav_offset_right > len_grid
    ] = len_grid
    imin_broadened_wav_offset_left += 1
    imax_broadened_wav_offset_left += 1
    imin_broadened_wav_offset_right += 1
    imax_broadened_wav_offset_right += 1

    line2grid_projection_left = iwav_on_grid_left  # size N (number of lines)
    line2grid_projection_right = iwav_on_grid_right  # size N (number of lines)
    # ... deal with case where a new point is created (we dont want that)
    # ... just offset it by one unit (we're working with wavenumber_calc anyway,
    # ... these boundary points will be cropped by RADIS at the end of the calculation)
    line2grid_projection_left[line2grid_projection_left == len(wavenumber)] = (
        len(wavenumber) - 1
    )
    line2grid_projection_right[line2grid_projection_right == len(wavenumber)] = (
        len(wavenumber) - 1
    )

    @jit(nopython=True)
    def rough_sum_on_grid():
        """ Sum all lines linestrength density on a spectrum  """
        # Performance
        # ----------
        #
        # Test case:  6.5k lines, 18k grid points
        # - standard: 13.6 ms
        # - with @jit: 0.9 ms
        #
        # Array size
        # ----------
        # rough_spectrum has len_grid+2 because two points are used for out of
        # range intensities. We crop at the end
        k_rough_spectrum = np.zeros(len_grid + 2)
        j_rough_spectrum = np.zeros(len_grid + 2)
        for i in range(len(S_density_on_grid)):
            imin_left = imin_broadened_wav_offset_left[i]
            imax_left = imax_broadened_wav_offset_left[i]
            imin_right = imin_broadened_wav_offset_right[i]
            imax_right = imax_broadened_wav_offset_right[i]
            k_rough_spectrum[imin_left : imax_left + 1] += (
                frac_left[i] * S_density_on_grid[i]
            )
            k_rough_spectrum[imin_right : imax_right + 1] += (
                frac_right[i] * S_density_on_grid[i]
            )
            j_rough_spectrum[imin_left : imax_left + 1] += (
                frac_left[i] * Ei_density_on_grid[i]
            )
            j_rough_spectrum[imin_right : imax_right + 1] += (
                frac_right[i] * Ei_density_on_grid[i]
            )
        # crop out of range points
        return k_rough_spectrum[1:-1], j_rough_spectrum[1:-1]

    k_rough_spectrum, j_rough_spectrum = rough_sum_on_grid()

    return (
        k_rough_spectrum,
        j_rough_spectrum,
        S_density_on_grid,
        Ei_density_on_grid,
        line2grid_projection_left,
    )

    #


if __name__ == "__main__":

    from radis.test.lbl.test_broadening import _run_testcases

    print("Testing broadening.py:", _run_testcases(verbose=True, plot=True))
