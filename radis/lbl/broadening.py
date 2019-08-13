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

PUBLIC FUNCTIONS - BROADENING

- :py:func:`radis.lbl.broadening.gaussian_broadening_FWHM`
- :py:func:`radis.lbl.broadening.gaussian_lineshape`
- :py:func:`radis.lbl.broadening.pressure_broadening_FWHM`
- :py:func:`radis.lbl.broadening.pressure_lineshape`
- :py:func:`radis.lbl.broadening.voigt_broadening_FWHM`
- :py:func:`radis.lbl.broadening.voigt_lineshape`

PRIVATE METHODS - BROADENING
(all computational-heavy functions: calculates all lines broadening,
convolve them, apply them on all calculated range)

- :py:func:`radis.lbl.broadening._whiting`
- :py:func:`radis.lbl.broadening._whiting_jit` : precompiled version
- :py:meth:`radis.lbl.broadening.BroadenFactory._calc_broadening_FWHM`
- :py:meth:`radis.lbl.broadening.BroadenFactory._add_voigt_broadening_FWHM`
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
from numpy import (exp, arange, zeros_like, trapz, pi, sqrt, sin)
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


def gaussian_broadening_FWHM(wav, molar_mass, Tgas):
    ''' Computes Doppler (Gaussian) broadening FWHM over all lines with [1]_, [2]_

    .. math::
        
        2\\frac{w}{c} \\sqrt{\\frac{2N_a k_b T_{gas} \\ln(2)}{M}}
        
    with k and c in CGS

    *generated from the Python formula with* :py:func:`~pytexit.pytexit.py2tex`

    Parameters
    ----------

    wav: array like (nm / cm-1)
        transition waverange  [length N = number of lines]

    molar_mass: array like (g/mol)
        molar mass for isotope of given transition   [length N]

    Tgas: float (K)
        (translational) gas temperature

    Returns
    -------

    fwhm_gauss: numpy array        [shape N]
        gaussian FWHM for all lines

    References
    ----------

    $$ I_{gaussian} = exp(- ln2 (Δw/HWHM)^2 )) $$

    with HWHM full-width half maximum

    .. [1] `Wikipedia <https://en.wikipedia.org/wiki/Doppler_broadening>`_ (in cm-1)

    .. [2] `Laux et al, 2003, "Optical diagnostics of atmospheric pressure air plasmas" <http://iopscience.iop.org/article/10.1088/0963-0252/12/2/301/meta>`_
           (in nm)

    Both give the same results. Here we use the CGS nomenclature 
    (reminder: molar_mass in g/mol)


    '''

    # Broadening parameter:
    # ... Doppler broadened half-width:
    # ... Reference: either Wikipedia (in cm-1), or former Sean's code (HITRAN based)
    # ... in CGS, or Laux 2003 (in nm). Both give the same results. Here we
    # ... use the CGS nomenclature (reminder: dg.molar_mass in g/mol)
    gamma_db = (wav/c_CGS)*sqrt((2*Na*k_b_CGS*Tgas*ln(2)) / molar_mass)  # HWHM (cm-1)

    return 2 * gamma_db

def gaussian_lineshape(dg, w_centered):
    ''' Computes Doppler (Gaussian) lineshape over all lines with [1]_, [2]_

    .. math::
        
        \\frac{1}{\\gamma_{db} \\sqrt{2\\pi}} \\operatorname{exp}\\left(-\\ln(2) 
        \\left(\\frac{w_{centered}}{\\gamma_{db}}\\right)^2\\right)

    *generated from the Python formula with* :py:func:`~pytexit.pytexit.py2tex`

    Parameters
    ----------

    dg: pandas Dataframe    [length N]
        list of lines   (must includes `gamma_db` for broadening calculation)

    w_centered: array       [one per line: shape W x N]
        waverange (nm / cm-1) (centered on 0, broadening width size)

    Tgas: K
        (translational) gas temperature

    Returns
    -------

    lineshape: pandas Series        [shape N x W]
        line profile  

    References
    ----------

    $$ I_{gaussian} = exp(- ln2 (Δw/HWHM)^2 )) $$

    with HWHM full-width half maximum

    .. [1] `Wikipedia <https://en.wikipedia.org/wiki/Doppler_broadening>`_ (in cm-1)

    .. [2] `Laux et al, 2003, "Optical diagnostics of atmospheric pressure air plasmas" <http://iopscience.iop.org/article/10.1088/0963-0252/12/2/301/meta>`_
           (in nm)

    Both give the same results. Here we use the CGS nomenclature 
    (reminder: molar_mass in g/mol)


    '''

    # Prepare coefficients, vectorize
    # --------

    # Broadening parameter:
    # ... Doppler broadened half-width:
    # ... Reference: either Wikipedia (in cm-1), or former Sean's code (HITRAN based)
    # ... in CGS, or Laux 2003 (in nm). Both give the same results. Here we
    # ... use the CGS nomenclature (reminder: dg.molar_mass in g/mol)
    gamma_db = dg.fwhm_gauss / 2     # FWHM > HWHM

    # Prepare vectorized operations:
    try:  # make it a (1,N) row vector
        gamma_db = gamma_db.values.reshape((1, -1))
    except AttributeError:  # probably not a dataframe: assert there is one line only.
        assert type(gamma_db) is np.float64

    # Calculate broadening
    # ------
    lineshape = (1/(gamma_db*sqrt(2*pi)))*exp(-ln(2)*(w_centered/gamma_db)**2)

    return lineshape


def pressure_broadening_FWHM(airbrd, selbrd, Tdpair, Tdpsel, 
                             pressure_atm, mole_fraction, Tgas, Tref):
    ''' Calculates collisional broadening FWHM over all lines by scaling 
    tabulated FWHM for new pressure and mole fractions conditions [1]_

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
        line profile  

    References
    ----------

    .. [1] `Rothman 1998 (HITRAN 1996) eq (A.14) <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`_

    '''

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
        gamma_lb = ((Tref/Tgas)**Tdpair)*((airbrd*pressure_atm*(1-mole_fraction)) +
                                          (selbrd*pressure_atm*mole_fraction))
    else:
        gamma_lb = (((Tref/Tgas)**Tdpair)*(airbrd*pressure_atm*(1-mole_fraction)) +
                    ((Tref/Tgas)**Tdpsel)*(selbrd*pressure_atm*mole_fraction))

    return 2*gamma_lb


def pressure_lineshape(dg, w_centered):
    ''' Computes collisional broadening over all lines [1]_

    .. math:: 
        
        \\frac{1}{\\pi} \\frac{\\gamma_{lb}}{\\gamma_{lb}^2+w_{centered}^2}

    *generated from the Python formula with* :py:func:`~pytexit.pytexit.py2tex`

    Parameters
    ----------

    dg: pandas Dataframe    [length N]
        list of lines. Columns should includes 

        - fwhm_lorentz   (cm-1)
            full-width half max coefficient (FWHM) for pressure broadening 
            calculation)

    w_centered: array       [one per line: shape W x N]
        waverange (nm / cm-1) (centered on 0)

    Returns
    -------

    lineshape: pandas Series        [shape N x W]
        line profile  

    References
    ----------

    .. [1] `Rothman 1998 (HITRAN 1996) eq (A.14) <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`_

    '''

    # Get collisional broadening HWHM
    gamma_lb = dg.fwhm_lorentz / 2        # FWHM > HWHM

    # Prepare vectorized operations:
    try:  # make it a (1, N) row vector
        gamma_lb = gamma_lb.values.reshape((1, -1))
    except AttributeError:  # probably not a dataframe: assert there is one line only.
        assert type(gamma_lb) is np.float64

    # Calculate broadening
    # -------
    lineshape = 1/pi*gamma_lb/((gamma_lb**2)+(w_centered**2))

    return lineshape


def voigt_broadening_FWHM(airbrd, selbrd, Tdpair, Tdpsel, wav, molar_mass, 
                          pressure_atm, mole_fraction, Tgas, Tref):
    ''' Calculate Voigt profile FWHM for lines in df from the Gaussian and 
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

    Environment parameters:

    pressure_atm: float  [atm]

    mole_fraction: float [0-1]

    Tgas: K
        (translational) gas temperature

    Tref: K
        reference temperature at which tabulated HWHM pressure 
        broadening coefficients were tabulated

    Returns
    -------

    wv, wl, wg: numpy array
        Voigt, Lorentz, and Gaussian FWHM    

    References
    ----------

    .. [1] `Olivero 1977 "Empirical fits to the Voigt line width: A brief review" <https://www.sciencedirect.com/science/article/pii/0022407377901613>`_

    .. [2] `Wikipedia <https://en.wikipedia.org/wiki/Doppler_broadening>`_ (in cm-1)

    .. [3] `Rothman 1998 (HITRAN 1996) eq (A.12) <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`_

    '''

    # Collisional broadening HWHM
    # ... Temperature and pressure dependent half width
    # ... Reference: Rothman 1998 (HITRAN 1996) eq (A.12)
    # ... Hypothesis: we only consider self broadening, and air broadening,
    # ... weighted with mole fractions
    # Note: if not Tdpsel in dg Tdpair is used. Lookup parent function
        # | dev note: in that case we simplify the expression by calculation the 
        # | power function once only.
    if Tdpsel is None:
        gamma_lb = ((Tref/Tgas)**Tdpair)*((airbrd*pressure_atm*(1-mole_fraction)) +
                                          (selbrd*pressure_atm*mole_fraction))
    else:
        gamma_lb = (((Tref/Tgas)**Tdpair)*(airbrd*pressure_atm*(1-mole_fraction)) +
                    ((Tref/Tgas)**Tdpsel)*(selbrd*pressure_atm*mole_fraction))

    # Doppler Broadening HWHM:
    # ... Doppler broadened half-width:
    # ... Reference: either HITRAN (in CGS), or Laux 2003 (in nm). Both are equivalent. 
    # ... Here we use the CGS nomenclature (reminder: dg.molar_mass in g/mol)
    gamma_db = (wav/c_CGS)*sqrt((2*Na*k_b_CGS*Tgas*ln(2)) / molar_mass)  # HWHM (cm-1)

    # Calculate broadening
    # -------
    # ... Reference Whiting 1968  Equation (5)
    # ... note that formula is given in wavelength (nm) [doesnt change anything]
    # ... and uses full width half maximum (FWHM)
    wg = 2*gamma_db          # HWHM > FWHM
    wl = 2*gamma_lb          # HWHM > FWHM
#    wv = wl/2 + sqrt((1/4*wl**2+wg**2))
    sd = (wl-wg)/(wl+wg)
    wv = (1 - 0.18121 * (1 - sd**2) - (0.023665 * exp(0.6*sd) + 0.00418 * exp(-1.9*sd))
          * sin(pi * sd))*(wl + wg)

    return wv, wl, wg


def voigt_lineshape(dg, w_centered, jit=True):
    ''' Calculates Voigt lineshape using the approximation of the Voigt profile of
    Whiting [1]_, [2]_ that maintains a good accuracy in the far wings. Exact for a pure 
    Gaussian and pure Lorentzian

    Parameters
    ----------

    dg: pandas Dataframe    [length N]
        list of lines. Columns should include:

        - `fwhm_voigt`: Voigt full-width half max coefficient (HWHM )
            calculated by calc_voigt_broadening_FWHM)

    w_centered: array      [one per line: shape W x N]
        waverange (nm / cm-1) (centered on 0)

    Other Parameters
    ----------------
    
    jit: boolean
        if ``True``, use just in time compiler. Usually faster when > 10k lines

    Returns
    -------

    lineshape: pandas Series        [shape N x W]
        line profile  

    References
    ----------

    .. [1] `NEQAIR 1996 User Manual, Appendix D <https://ntrs.nasa.gov/search.jsp?R=19970004690>`_

    .. [2] `Whiting 1968 "An empirical approximation to the Voigt profile", JQSRT <https://www.sciencedirect.com/science/article/pii/0022407368900812>`_

    '''

    # Prepare coefficients, vectorize
    # --------
    # Get Voigt FWHM
    if not 'fwhm_voigt' in list(dg.keys()):
        raise KeyError('fwhm_voigt: Calculate broadening with ' +
                       'calc_voigt_broadening_FWHM first')
    if not 'fwhm_lorentz' in list(dg.keys()):
        raise KeyError('fwhm_lorentz: Calculate broadening with ' +
                       'calc_voigt_broadening_FWHM first')
    wv = dg.fwhm_voigt
    wl = dg.fwhm_lorentz

    # Prepare vectorized operations:
    try:  # make it a (1, N) row vector
        wl = wl.values.reshape((1, -1))
#        wg = wg.values.reshape((1,-1))
        wv = wv.values.reshape((1, -1))
    except AttributeError:  # probably not a dataframe: assert there is one line only.
        assert type(wl) is np.float64
#        assert type(wg) is np.float64
        assert type(wl) is np.float64

    if jit:
        lineshape = _whiting_jit(wl, wv, w_centered)
    else:
        lineshape = _whiting(wl, wv, w_centered)

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

def _whiting(wl, wv, w_centered):
    ''' 
    Notes
    -----
    
    Performances:
        
    using @jit yield a performance increase from 8.9s down to 5.1s 
    on a 50k lines, 250k wavegrid case (performances.py) 
    '''
    # Calculate some temporary arrays
    # ... fasten up the calculation by 25% (ex: test on 20 cm-1, ~6000 lines:
    # ... 20.5.s > 16.5s) on the total eq_spectrum calculation
    # ... w_wv is typically a (10.001, 1997) array
    w_wv = (w_centered/wv)                  # w_centered can be ~500 Mb
    w_wv_2 = w_wv**2
    wl_wv = wl/wv
    w_wv_225 = np.abs(w_wv)**2.25

    # Calculate!  (>>> this is the performance bottleneck <<< : ~ 2/3 of the time spent
    #              on lineshape equation below + temp array calculation above
    #              In particular exp(...) and ()**2.25 are very expensive <<< )
    # ... Voigt 1st order approximation
    lineshape = ((1-wl_wv)*exp(-2.772*w_wv_2)+wl_wv*1/(1+4*w_wv_2)
                 # ... 2nd order correction
                 + 0.016*(1-wl_wv)*wl_wv*(exp(-0.4*w_wv_225)
                                          - 10/(10+w_wv_225)))
    return lineshape

@jit(float64[:,:](float64[:,:], float64[:,:], float64[:,:]), nopython=True, cache=True) #, parallel=True)
def _whiting_jit(wl, wv, w_centered):
    ''' 
    Notes
    -----
    
    Performances:
        
    using @jit yield a performance increase from 8.9s down to 5.1s 
    on a 50k lines, 250k wavegrid case (performances.py) 
    '''
    # Calculate some temporary arrays
    # ... fasten up the calculation by 25% (ex: test on 20 cm-1, ~6000 lines:
    # ... 20.5.s > 16.5s) on the total eq_spectrum calculation
    # ... w_wv is typically a (10.001, 1997) array
    w_wv = (w_centered/wv)                  # w_centered can be ~500 Mb
    w_wv_2 = w_wv**2
    wl_wv = wl/wv
    w_wv_225 = np.abs(w_wv)**2.25

    # Calculate!  (>>> this is the performance bottleneck <<< : ~ 2/3 of the time spent
    #              on lineshape equation below + temp array calculation above
    #              In particular exp(...) and ()**2.25 are very expensive <<< )
    # ... Voigt 1st order approximation
    lineshape = ((1-wl_wv)*exp(-2.772*w_wv_2)+wl_wv*1/(1+4*w_wv_2)
                 # ... 2nd order correction
                 + 0.016*(1-wl_wv)*wl_wv*(exp(-0.4*w_wv_225)
                                          - 10/(10+w_wv_225)))
    return lineshape

# %% Tools


class BroadenFactory(BaseFactory):
    ''' A class that holds all broadening methods, inherited by 
    :class:`~radis.lbl.factory.SpectrumFactory` eventually 
    
    See Also
    --------

    :class:`~radis.lbl.factory.SpectrumFactory`
    '''

    def __init__(self):

        super(BroadenFactory, self).__init__()

        # Name variables (initialized later in SpectrumFactory)
        self.wbroad_centered = None

        self.wavenumber = None
        self.wavenumber_calc = None
        self.woutrange = None

        self._broadening_method = 'voigt'
        '''str: 'voigt', 'convolve'
        
        Calculates broadening with a direct voigt approximation ('voigt') or
        by convoluting independantly calculated Doppler and collisional
        broadening ('convolve'). First is much faster, 2nd can be used to
        compare results. Not a user available parameter for the moment, but
        you can edit the SpectrumFactory manually::
            
            sf = SpectrumFactory(...)
            sf._broadening_method = 'voigt'  
        '''

        # Predict broadening times (helps trigger warnings for optimization)
        self._broadening_time_ruleofthumb = 1e-7   # s / lines / point
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

    # %% Functions to calculate broadening FWHM

    def _calc_broadening_FWHM(self):
        ''' Calculate broadening FWHM and store in line dataframe (df1)

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
        '''

        # Init variables
        df = self.df1

        if len(df) == 0:
            return  # no lines
        
        if self.verbose >= 2:
            printg('Calculate broadening FWHM')
            t0 = time()

        if self.input.Tgas is None:
            raise AttributeError(
                "Tgas not defined. Make sure the parent function creates it")
        Tgas = self.input.Tgas
        pressure_mbar = self.input.pressure_mbar
        mole_fraction = self.input.mole_fraction
        # convert from mbar to atm for linebroadening calculation
        pressure_atm = pressure_mbar/1013.25
        # coefficients tabulation temperature
        Tref = self.input.Tref

        # Get broadenings
        if self._broadening_method == 'voigt':
            # Adds fwhm_voigt, fwhm_gauss, fwhm_lorentz:
            self._add_voigt_broadening_FWHM(
                df, pressure_atm, mole_fraction, Tgas, Tref)
        elif self._broadening_method == 'convolve':
            # Adds fwhm_lorentz:
            self._add_collisional_broadening_FWHM(
                df, pressure_atm, mole_fraction, Tgas, Tref)
            # Add fwhm_gauss:
            self._add_gaussian_broadening_FWHM(df, Tgas)

        if self.verbose >= 2:
            printg('Calculated broadening FWHM in {0:.1f}s'.format(time()-t0))

    def _add_voigt_broadening_FWHM(self, df, pressure_atm, mole_fraction, Tgas, Tref):
        ''' Update dataframe with Voigt FWHM

        Returns
        -------

        Note:
            But input pandas Dataframe ``'df'`` is updated with keys: 

            - fwhm_voigt

            - fwhm_lorentz

            - fwhm_gauss

        '''

        # Check self broadening is here
        if not 'Tdpsel' in list(df.keys()):
            self.warn('Self-broadening temperature coefficient Tdpsel not given in database: used Tdpair instead',
                      'MissingSelfBroadeningWarning')
            Tdpsel = None     # if None, voigt_broadening_FWHM uses df.Tdpair
        else:
            Tdpsel = df.Tdpsel

        # Calculate broadening FWHM
        wv, wl, wg = voigt_broadening_FWHM(df.airbrd, df.selbrd, df.Tdpair, Tdpsel, 
                                           df.wav, df.molar_mass,
                                           pressure_atm, mole_fraction, Tgas, Tref)
        
        # Update dataframe
        df['fwhm_voigt'] = wv
        df['fwhm_lorentz'] = wl
        df['fwhm_gauss'] = wg

        return

    def _add_collisional_broadening_FWHM(self, df, pressure_atm, mole_fraction, Tgas, Tref):
        ''' Update dataframe with collisional FWHM [1]_

        Returns
        -------

        Input pandas Dataframe 'df' is updated with keys: 

            - fwhm_lorentz

        Notes
        -----

        Temperature and pressure dependent half width

        Hypothesis: we only consider self broadening, and air broadening,
        weighted with mole fractions

        If not Tdpsel in dg use Tdpair

        References
        ----------

        .. [1] `Rothman 1998 (HITRAN 1996) eq (A.12) <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`_

        '''
        
        # Check self broadening is here
        if not 'Tdpsel' in list(df.keys()):
            self.warn('Self-broadening temperature coefficient Tdpsel not given in database: used Tdpair instead',
                      'MissingSelfBroadeningWarning')
            Tdpsel = None
        else:
            Tdpsel = df.Tdpsel

        # Calculate broadening FWHM
        wl = pressure_broadening_FWHM(df.airbrd, df.selbrd, df.Tdpair, Tdpsel,
                                      pressure_atm, mole_fraction, Tgas, Tref)

        # Update dataframe
        df['fwhm_lorentz'] = wl

        return

    def _add_gaussian_broadening_FWHM(self, df, Tgas):
        ''' Update dataframe with Gaussian FWHM

        Returns
        -------

        None: input pandas Dataframe 'df' is updated with keys: 

            - fwhm_gauss

        '''

        # Calculate broadening FWHM
        wg = gaussian_broadening_FWHM(df.wav, df.molar_mass, Tgas)

        # Update dataframe
        df['fwhm_gauss'] = wg

        return

    def _collisional_lineshape(self, dg, wbroad_centered):
        ''' Computes collisional broadening over all lines + normalize 
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

        '''

        # Calculate broadening for all lines
        # -------
        lineshape = pressure_lineshape(dg, wbroad_centered)

        # Normalize
        # ---------
        # ... 'wbroad_centered' is w_array-(w_shifted_line_center)
        # ... Normalization should not be needed as Lorentzian is normalized already
        # ... but that's only with a good enough wavestep
        # ... Here we compute the integral to check the error that is made
        area = trapz(lineshape.T, x=wbroad_centered.T)
        err = abs((area-1))
        if self.warnings and self.warnings['CollisionalBroadeningWarning'] != 'ignore':
            if (err > self.misc.warning_broadening_threshold).any():
                self.warn('Large error ({0:.1f}%) '.format(err.max()*100) +
                          'in pressure broadening. Increase broadening width / reduce wstep. ' +
                          'Use .plot_broadening() to visualize each line broadening',
                          'CollisionalBroadeningWarning')
            # Note that although there may be an error here the total broadening
            # is normalized anyway, so the energy is conserved. If we're only
            # looking at a slit-function broadened spectrum (slit>>FWHM) it
            # wont change much. However it does impact multi-slabs calculations
        # ... Renormalize to dampen numeric errors impact
        lineshape /= area

        return lineshape

    def _gaussian_lineshape(self, dg, wbroad_centered):
        ''' Computes Doppler (Gaussian) broadening over all lines + normalize 
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

        '''

        # Calculate broadening for all lines
        # -------
        lineshape = gaussian_lineshape(dg, wbroad_centered)

        # Normalize
        # ---------
        # ... normalisation not really needed as the analytical function is normalized
        # ... but that's only with a good enough wavestep
        # ... Here we compute the integral to check the error that is made
        area = trapz(lineshape.T, x=wbroad_centered.T)
        err = abs((area-1))
        if self.warnings and self.warnings['GaussianBroadeningWarning'] != 'ignore':
        # In a "performance" mode (vs "safe" mode), these warnings would be disabled
            if (err > self.misc.warning_broadening_threshold).any():
                self.warn('Large error ({0:.1f}%) '.format(err.max()*100) +
                          'in Doppler broadening. Increase broadening width / reduce wstep. ' +
                          'Use .plot_broadening() to visualize each line broadening',
                          'GaussianBroadeningWarning')
                # Note that although there may be an error here the total broadening
                # is normalized anyway, so the energy is conserved. If we're only
                # looking at a slit-function broadened spectrum (slit>>FWHM) it
                # wont change much. However it does impact multi-slabs calculations
        lineshape /= area

        return lineshape

    def _voigt_broadening(self, dg, wbroad_centered, jit=True):
        ''' Computes voigt broadening over all lines + normalize 

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

        '''

        # Calculate broadening for all lines
        # -------
        lineshape = voigt_lineshape(dg, wbroad_centered, jit=jit)

        return lineshape

    # %% Function to calculate lineshapes from FWHM

    def _calc_lineshape(self, dg):
        ''' Sum over each line (trying to use vectorize operations to be faster)

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
          is applied.  # TODO

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

        '''

        # Init variables
        if self.input.Tgas is None:
            raise AttributeError(
                "Tgas not defined. Make sure the parent function creates it")
#        Tgas = self.input.Tgas

#        pressure_mbar = self.input.pressure_mbar
#        mole_fraction = self.input.mole_fraction
#        pressure_atm = pressure_mbar/1013.25     # convert from mbar to atm for linebroadening calculation
#        Tref = self.input.Tref                         # coefficients tabulation temperature

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
        wbroad = wbroad_centered+shifted_wavenum

        # Calculate lineshape (using precomputed FWHM)
        if self._broadening_method == 'voigt':
            line_profile = self._voigt_broadening(dg, wbroad_centered)
        elif self._broadening_method == 'convolve':
            # Get pressure and gaussian profiles
            pressure_profile = self._collisional_lineshape(dg, wbroad_centered)
            gaussian_profile = self._gaussian_lineshape(dg, wbroad_centered)

            # Convolve and get final line profile:
            line_profile = np.empty_like(pressure_profile)     # size (B, N)
            for i, (x, y) in enumerate(zip(pressure_profile.T, gaussian_profile.T)):
                line_profile[:, i] = np.convolve(x, y, 'same')
            line_profile = line_profile / trapz(line_profile.T, x=wbroad.T)  # normalize
            # ... Note that normalization should not be needed as broadening profiles
            # ... are created normalized already. However, we do normalize to reduce
            # ... the impact of any error in line_profiles (due to wstep too big or
            # ... broadening_width too small): at least the energy is conserved, even
            # ... if not perfectly distributed (spectrally). A warning is raised by the
            # ... broadening functions.
        else:
            raise ValueError('Unexpected broadening calculation method: {0}'.format(
                self._broadening_method))

        return line_profile

    def plot_broadening(self, i=0, pressure_atm=None, mole_fraction=None,
                        Tgas=None):
        ''' just for testing. Recalculate and plot broadening for line of index i

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

        '''
        
        from publib import set_style, fix_style

        if pressure_atm is None:
            pressure_atm = self.input.pressure_mbar/1013.25
        if mole_fraction is None:
            mole_fraction = self.input.mole_fraction
        if Tgas is None:
            Tgas = self.input.Tgas

        # Get one line only:
        dg = self.df1.iloc[i]

        wbroad_centered = self.wbroad_centered
        wbroad = wbroad_centered+dg.shiftwav
#        convolve_profile = self._calc_lineshape(dg)
        # Get Voigt from empirical approximation
        voigt_profile = self._voigt_broadening(dg, wbroad_centered, jit=False)
        # Get Voigt from convolution
        pressure_profile = self._collisional_lineshape(dg, wbroad_centered)
        gaussian_profile = self._gaussian_lineshape(dg, wbroad_centered)
        line_profile = np.convolve(pressure_profile, gaussian_profile, 'same')
        line_profile /= trapz(line_profile.T, x=wbroad.T)  # normalize

        # Plot!
        set_style('origin')
        plt.figure()
        plt.plot(wbroad, pressure_profile, label='Pressure')
        plt.plot(wbroad, gaussian_profile, label='Doppler')
        plt.plot(wbroad, line_profile, label='Voigt (convolved)', lw=3)
        plt.plot(wbroad, voigt_profile, label='Voigt (approximation)')
        plt.xlabel('Wavenumber (cm-1)')
        plt.ylabel('Broadening coefficient')
        plt.title('Line {0}, T={1:.1f}K, P={2:.2f}atm'.format(i, Tgas,
                                                              pressure_atm))
        plt.legend()
        fix_style('origin')

        return

    def _apply_lineshape(self, broadened_param, line_profile, shifted_wavenum):
        ''' Multiply `broadened_param` by `line_profile` and project it on the
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

        '''
        
        if __debug__:
            t0 = time()

#        # Get spectrum range
        wavenumber = self.wavenumber  # get vector of wavenumbers (shape W)
        wavenumber_calc = self.wavenumber_calc
        # generate the vector of wavenumbers (shape W + space B on the sides)
        vec_length = len(wavenumber)

        # Vectorize the chunk of lines
        S = broadened_param.values.reshape((1, -1))
        shifted_wavenum = shifted_wavenum.values.reshape((1, -1))  # make it a row vector

        # Get broadening array
        wbroad_centered = self.wbroad_centered  # size (B,)
        # index of broadening half width
        iwbroad_half = len(wbroad_centered)//2

        # Calculate matrix of broadened parameter (for all lines)
        profile_S = line_profile * S

        # ---------------------------
        # Apply line profile

        # ... First get closest matching line:
        idcenter = np.searchsorted(wavenumber_calc, shifted_wavenum.T, side="left").ravel()

        # ... Then distribute the line over the correct location. All of that vectorized
        sumoflines_calc = zeros_like(wavenumber_calc)

        # Note on performance: it isn't easy to vectorize this part as all
        # lines are to be plotted on a different place.

        # to avoid the If / Else condition in the loop, we do a vectorized
        # comparison beforehand and run 3 different loops

        # reminder: wavenumber_calc has size [iwbroad_half+vec_length+iwbroad_half]
        assert len(wavenumber_calc) == vec_length+2*iwbroad_half
        boffrangeleft = (idcenter <= iwbroad_half)
        boffrangeright = (idcenter >= vec_length+iwbroad_half)
        binrange = np.ones_like(idcenter, dtype=bool) ^ (
            boffrangeleft + boffrangeright)

#        # Performance for lines below
#        # ----------
#        #
#        # on test case: 6.5k lines x 18.6k grid length
#        # normal: ~ 36 ms  called 9 times
#        # with @jit : ~ 200 ms called 9 times (worse!)
        
        # Off Range, left : only aggregate the Right wing
        lines_l = profile_S.T[boffrangeleft]
        if len(lines_l) > 0:
            I_low_l = idcenter[boffrangeleft]+1
            # idcenter[boffrangeleft]+iwbroad_half
            I_high_l = I_low_l + iwbroad_half - 1
            for i, profS in enumerate(lines_l):
                # cut left wing + peak
                sumoflines_calc[I_low_l[i]:I_high_l[i] +
                                1] += profS[iwbroad_half+1:]
#                assert(len(absorption_calc[I_low_l[i]:I_high_l[i]+1])==iwbroad_half)
#                assert(len(profS[iwbroad_half+1:])==iwbroad_half)

        # Off Range, Right : only aggregate the left wing
        lines_r = profile_S.T[boffrangeright]
        if len(lines_r) > 0:
            I_low_r = idcenter[boffrangeright]-iwbroad_half
            I_high_r = I_low_r + iwbroad_half - \
                1      # idcenter[boffrangeright]-1
            for i, profS in enumerate(lines_r):
                # cut right wing + peak
                sumoflines_calc[I_low_r[i]:I_high_r[i] +
                                1] += profS[:iwbroad_half]
#                assert(len(profS[:iwbroad_half])==iwbroad_half)

        # In range: aggregate both wings
        lines_i = profile_S.T[binrange]
        if len(lines_i) > 0:
            I_low_i = idcenter[binrange]-iwbroad_half
            # idcenter[binrange]+iwbroad_half
            I_high_i = I_low_i+2*iwbroad_half
            for i, profS in enumerate(lines_i):
                sumoflines_calc[I_low_i[i]:I_high_i[i]+1] += profS
#                assert(len(profS)==2*iwbroad_half+1)

        # Get valid range (discard wings)
        sumoflines = sumoflines_calc[self.woutrange]

        return wavenumber, sumoflines

    def _broaden_lines(self, df):
        ''' Divide over chuncks not to process to many lines in memory at the
        same time (note that this is not where the parallelisation is done: all
        lines are processed on the same core. )
        
        Parameters
        ----------
        
        self: Factory
            contains the ``self.misc.chunksize`` parameter
            
        df: DataFrame
            line dataframe

        See _calc_lineshape for more information
        '''
        # --------------------------

        # Reactivate warnings
        reset_warnings(self.warnings)

        # Init arrays
        wavenumber = self.wavenumber
        # Get number of groups for memory splitting
        chunksize = self.misc.chunksize
        
        
        try:
            if chunksize is None:
                # Deal with all lines directly (usually faster)
                line_profile = self._calc_lineshape(df)
                (wavenumber, abscoeff) = self._apply_lineshape(
                                        df.S, line_profile, df.shiftwav)

            else:
                # Cut lines in smaller bits for better memory handling
                N = int(len(df)*len(wavenumber)/chunksize)+1
                # Too big may be faster but overload memory.
                # See Performance for more information
        
                abscoeff = zeros_like(self.wavenumber)
     
                pb = ProgressBar(N, active=self.verbose)
                for i, (_, dg) in enumerate(df.groupby(arange(len(df)) % N)):
                    line_profile = self._calc_lineshape(dg)
                    (wavenumber, absorption) = self._apply_lineshape(
                        dg.S, line_profile, dg.shiftwav)
                    abscoeff += absorption
                    pb.update(i)
                pb.done()
                
        except MemoryError:
            import traceback
            traceback.print_exc()
            raise MemoryError('See details above. Try to use or reduce the '+\
                              'chunksize parameter (current={0})'.format(chunksize))

        return wavenumber, abscoeff

    def _broaden_lines_noneq(self, df):
        ''' Divide over chuncks not to process to many lines in memory at the
        same time (note that this is not where the parallelisation is done: all
        lines are processed on the same core)

        See _calc_lineshape for more information
        '''

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
                line_profile = self._calc_lineshape(df)
                (wavenumber, abscoeff) = self._apply_lineshape(
                            df.S, line_profile, df.shiftwav)
                (_, emisscoeff) = self._apply_lineshape(
                            df.Ei, line_profile, df.shiftwav)
            
            else:
                # Cut lines in smaller bits for better memory handling
                
                # Get size of numpy array for vectorialization
                N = int(len(df)*len(wavenumber)/chunksize)+1
                # Too big may be faster but overload memory.
                # See Performance for more information
        
                abscoeff = zeros_like(self.wavenumber)
                emisscoeff = zeros_like(self.wavenumber)
        
                pb = ProgressBar(N, active=self.verbose)
                for i, (_, dg) in enumerate(df.groupby(arange(len(df)) % N)):
                    line_profile = self._calc_lineshape(dg)
                    (wavenumber, absorption) = self._apply_lineshape(
                        dg.S, line_profile, dg.shiftwav)
                    (_, emission) = self._apply_lineshape(
                        dg.Ei, line_profile, dg.shiftwav)
                    abscoeff += absorption       #
                    emisscoeff += emission
                    pb.update(i)
                pb.done()

        except MemoryError:
            import traceback
            traceback.print_exc()
            raise MemoryError('See details above. Try to use or reduce the chunksize parameter')

        return wavenumber, abscoeff, emisscoeff

    # %% Generate absorption profile which includes linebroadening factors

    def _calc_broadening(self):
        '''
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


        '''

        parallel = self.misc.parallel
        df = self.df1
        
        if self.verbose>=2:
            printg('Calculating line broadening ({0} lines: expect ~ {1:.1f}s on 1 CPU)'.format(
                    len(df), self._broadening_time_ruleofthumb*len(df)*len(self.wbroad_centered)))
            t0 = time()

        # Just some tests
        try:
            assert len(df.shape) == 2
        except AssertionError:
            warn('Dataframe has only one line. Unexpected behaviour could occur' +
                 ' because Dataframes will be handled as Series and row/columns' +
                 ' may be inverted')

        if parallel:
            # Parallel version
            Ngroups = self.misc.Ngroups
            Nprocs = self.misc.Nprocs  # number of processes
            if self.verbose:
                print('Parallel version on {0} groups, {1} procs'.format(
                    Ngroups, Nprocs))
            t1 = time()
            with Pool(Nprocs) as p:  # run on all procs
                # divide all lines in Ngroups
                ret_list = p.map(self._broaden_lines, [group for name, group in df.groupby(
                    arange(len(df)) % Ngroups)])
                # each process returns an array on the whole range, but only
                # including some of the lines
            if self.verbose:
                print('process: {0:.1f}s'.format(time()-t1))

            w_arr, a_arr = list(zip(*ret_list))
            wavenumber = w_arr[0]  # they're all the same anyway

            abscoeff = np.array(a_arr).sum(axis=0)    # sum over all lines

        else:
            # Regular version
            (wavenumber, abscoeff) = self._broaden_lines(df)

        if self.verbose>=2:
            printg('Calculated line broadening in {0:.1f}s'.format(time()-t0))

        return wavenumber, abscoeff

    def _calc_broadening_noneq(self):
        '''
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

        '''

        parallel = self.misc.parallel
        df = self.df1

        if self.verbose>=2:
            printg('Calculating line broadening ({0:,d} lines: expect ~ {1:.1f}s on 1 CPU)'.format(
                    len(df), self._broadening_time_ruleofthumb*len(df)*len(self.wbroad_centered)))
            t0 = time()

        # Just some tests
        try:
            assert len(df.shape) == 2
        except AssertionError:
            warn('Dataframe has only one line. Unexpected behaviour could occur' +
                 ' because Dataframes will be handled as Series and row/columns' +
                 ' may be inverted')

        if parallel:
            # Parallel version
            Ngroups = cpu_count()
            if self.verbose:
                print('Parallel version on {0} groups, {1} procs'.format(
                    Ngroups, cpu_count()))
            with Pool(cpu_count()) as p:  # run on all procs
                # divide all lines in Ngroups
                ret_list = p.map(self._broaden_lines_noneq, [group for name, group in df.groupby(
                    arange(len(df)) % Ngroups)])
                # each process returns an array on the whole range, but only
                # including some of the lines
            w_arr, a_arr, e_arr = list(zip(*ret_list))
            wavenumber = w_arr[0]
            abscoeff = np.array(a_arr).sum(axis=0)    # sum over all lines
            emisscoeff = np.array(e_arr).sum(axis=0)      # sum over all lines

        else:
            # Regular version
            (wavenumber, abscoeff, emisscoeff) = self._broaden_lines_noneq(df)

        if self.verbose>=2:
            printg('Calculated line broadening in {0:.1f}s'.format(time()-t0))

        return wavenumber, abscoeff, emisscoeff

# %% Functions to calculate semi-continuum

    def _find_weak_lines(self, weak_rel_intensity_threshold):
        ''' Finds weak lines in current line dataframe. 

        Lines are considered weak if recovered by a strong line nearby. These 
        lines are later moved in a semi-continuum, and only strong lines are 
        fully resolved with Voigt broadening (which is costly!)

        To find weak planes, we perform a fast broadening of all lines on a 
        rectangle of same FWHM as the lines. Then, a rough spectrum of linestrengths
        is calculated, and lines linestrength are compared to the rough spectrum 

        This procedure allows to discard many weak lines while preserving the 
        spectrum main features in the less intense parts, what an absolute
        cutoff such as the linestrength 'cutoff' cannot do

        Returns
        -------

        None
            store weak line status in dataframe as ``self.df1.weak_line``

        '''

        # Get inputs
        wavenumber_calc = self.wavenumber_calc     # size W
        wstep = self.params.wstep
        df = self.df1              # lines already scaled with current temperature, size N
        
        if self.verbose>=2:
            printg('... classifying lines as weak or strong')
            t0 = time()

        # Get approximate spectral absorption coefficient
        rough_spectrum, S_density_on_grid, line2grid_proj = project_lines_on_grid(
                                                       df, wavenumber_calc, wstep)
        #     :
        # ~ 1/(#.cm-2)
        # Sizes:
        # - rough_spectrum:     size W
        # - S_density_on_grid:  size N
        # - line2grid_proj:     size N

        # Weak line criteria
        # ... Compare line density (in 1/(#.cm-2) to sum)
        line_is_weak = S_density_on_grid < (
                weak_rel_intensity_threshold * rough_spectrum[line2grid_proj])   # size N

        # ... Store weak line label in df
        df['weak_line'] = line_is_weak

        if self.verbose>=2:
            printg('... {0:,d} lines classified as weak lines ({1:.2f}%) in {2:.1f}s'.format(
                line_is_weak.sum(), line_is_weak.sum()/len(line_is_weak)*100,
                time()-t0))

        return

    def _calculate_pseudo_continuum(self, noneq=False):
        ''' Find weak lines, add them in pseudo-continuum  (note that pseudo-continuum
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

        '''

        if self.params.pseudo_continuum_threshold > 0:
    
            if self.verbose>=2:
                printg('Calculating pseudo continuum')
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
                    'All lines qualified as weak: reduce weak_line_threshold')

            # Calculate continuum
            if noneq:
                k_continuum, j_continuum, _, _, _ = project_lines_on_grid_noneq(
                    df_weak_lines, wavenumber_calc, wstep)

                if __debug__:
                    printdbg('Intensity of k continuum: {0}\n'.format(np.trapz(k_continuum, wavenumber_calc)) +
                             'Intensity of lines removed: {0}'.format(df_weak_lines.S.sum()))
                    printdbg('Intensity of j continuum: {0}\n'.format(np.trapz(j_continuum, wavenumber_calc)) +
                             'Intensity of lines removed: {0}'.format(df_weak_lines.Ei.sum()))

            else:
                k_continuum, _, _ = project_lines_on_grid(
                    df_weak_lines, wavenumber_calc, wstep)
    
                if __debug__:
                    printdbg('Intensity of continuum: {0}\n'.format(np.trapz(k_continuum, wavenumber_calc)) +
                             'Intensity of lines removed: {0}'.format(df_weak_lines.S.sum()))

            # Get valid range (discard wings)
            k_continuum = k_continuum[self.woutrange]   # 1/(#.cm-2)
#            self.continuum_k = k_continuum
            
            if noneq:
                j_continuum = j_continuum[self.woutrange]   # 1/(#.cm-2)

            # Reduce line dataset to strong lines only
            self.df1 = df_strong_lines

        # Update number of lines

            self._Nlines_in_continuum = len(df_weak_lines)
            self._Nlines_calculated = len(self.df1)
    
            # Check performances
            time_spent = time()-t0
            # ... Expected broadening time gain (see Rule of Thumb) 
            expected_broadening_time_gain = (self._broadening_time_ruleofthumb*
                                             self._Nlines_in_continuum*
                                             len(self.wbroad_centered))
            if self.verbose>=2:
                printg('Calculated pseudo-continuum in {0:.1f}s (expected time saved: {1:.1f}s)'.format(
                        time_spent, expected_broadening_time_gain))
                # Add a warning if it looks like it wasnt worth it
            if time_spent > 3*expected_broadening_time_gain:
                self.warn('Pseudo-continuum may not be adapted to this kind '+\
                          'of spectrum. Time spent on continuum calculation '+\
                          '({0:.1f}s) is much longer than expected gain ({1:.1f}s). '.format(
                           time_spent, expected_broadening_time_gain)+\
                          'If the calculation takes a lot of time, consider '+\
                          'setting pseudo_continuum_threshold=0. If it is fast '+\
                          'enough already, you can decrease the linestrength cutoff= '+\
                          'parameter to add discarded weak lines to continuum '+\
                          'and improve the calculation accuracy.',
                          'PerformanceWarning')
                

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
        '''
        Notes
        -----
        
        also used for adding emisscoeff continuum with::
            
            self._add_pseudo_continuum(emisscoeff_v, j_continuum):
        '''
        if k_continuum is not None:
            abscoeff_v += k_continuum
        return abscoeff_v


def project_lines_on_grid(df, wavenumber, wstep):
    ''' Quickly sums all lines on wavespace grid as rectangles of FWHM corresponding
    to fwhm_voigt and a spectral absorption coefficient value so that linestrength 
    is conserved
    
    Parameters
    ----------
    
    df: pandas Dataframe
        Contains ``shiftwav`` (wavenumbers) and ``S`` (linestrengths) and ``wv`` 
        (Voigt FWHM) size ``N`` (number of lines)
        
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
        
    '''

    shiftwav = df.shiftwav.values       # cm-1  ,   size N (number of lines)
    S = df.S.values                     # cm/#  ~   cm-1/(#.cm-2)  ,   size N
    wv = df.fwhm_voigt.values           # FWHM

    # ... First get closest matching line:
    iwav_on_grid = np.searchsorted(wavenumber, shiftwav.T, side="left").ravel()

    # ... express FWHM (in nm) in index 
    ihwhm_on_grid = np.asarray(wv / 2 // wstep, dtype=np.int64)
    ifwhm_on_grid = ihwhm_on_grid * 2 + 1   # make it odd (conserves center)
    # ... infer min and max index to project lines (sides may be out of range)
    imin_broadened_wav_on_grid = iwav_on_grid - ihwhm_on_grid
    imax_broadened_wav_on_grid = iwav_on_grid + ihwhm_on_grid

    # Get average intensity, assuming a rectangular profile
    S_density_on_grid = S/(ifwhm_on_grid*wstep)               # ~ 1/(#.cm-2)

    # Cut out of range points
    # ... you should see imin_broadened_wav_on_grid as a projection array, with
    # ... indexes where Intensity will be summed afterwards. 
    # ... here we project the out of range intensities to index -1 and len_grid+1
    # ... (this is faster than )
    len_grid = len(wavenumber)
    imin_broadened_wav_offset = imin_broadened_wav_on_grid
    imax_broadened_wav_offset = imax_broadened_wav_on_grid
    imin_broadened_wav_offset[imin_broadened_wav_offset < 0] = -1
    imax_broadened_wav_offset[imax_broadened_wav_offset > len_grid] = len_grid
    imin_broadened_wav_offset += 1
    imax_broadened_wav_offset += 1

    line2grid_projection = iwav_on_grid          # size N (number of lines)
    # ... deal with case where a new point is created (we dont want that)
    # ... just offset it by one unit (we're working with wavenumber_calc anyway,
    # ... these boundary points will be cropped by RADIS at the end of the calculation)
    line2grid_projection[line2grid_projection==len(wavenumber)] = len(wavenumber) - 1

    @jit(float64[:](), nopython=True)
    def rough_sum_on_grid():
        ''' Sum all lines linestrength density on a spectrum  '''
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
        rough_spectrum = np.zeros(len_grid+2)
        for i, Iline_density in enumerate(S_density_on_grid):
            imin = imin_broadened_wav_offset[i]
            imax = imax_broadened_wav_offset[i]
            rough_spectrum[imin:imax+1] += Iline_density
        # crop out of range points
        return rough_spectrum[1:-1]

    k_rough_spectrum = rough_sum_on_grid()

    return k_rough_spectrum, S_density_on_grid, line2grid_projection

def project_lines_on_grid_noneq(df, wavenumber, wstep):
    ''' Quickly sums all lines on wavespace grid as rectangles of FWHM corresponding
    to fwhm_voigt and a spectral absorption coefficient value so that linestrength 
    is conserved
    
    Parameters
    ----------
    
    df: pandas Dataframe
        Contains ``shiftwav`` (wavenumbers) and ``S`` (linestrengths) and ``wv`` 
        (Voigt FWHM) size ``N`` (number of lines)
        
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
        
    '''

    shiftwav = df.shiftwav.values       # cm-1  ,   size N (number of lines)
    S = df.S.values                     # cm/#  ~   cm-1/(#.cm-2)  ,   size N
    Ei = df.Ei.values                   # mW/cm3/sr
    wv = df.fwhm_voigt.values           # FWHM

    # ... First get closest matching line:
    iwav_on_grid = np.searchsorted(wavenumber, shiftwav.T, side="left").ravel()

    # ... express FWHM (in nm) in index 
    ihwhm_on_grid = np.asarray(wv / 2 // wstep, dtype=np.int64)
    ifwhm_on_grid = ihwhm_on_grid * 2 + 1   # make it odd (conserves center)
    # ... infer min and max index to project lines (sides may be out of range)
    imin_broadened_wav_on_grid = iwav_on_grid - ihwhm_on_grid
    imax_broadened_wav_on_grid = iwav_on_grid + ihwhm_on_grid

    # Get average intensity, assuming a rectangular profile
    S_density_on_grid = S/(ifwhm_on_grid*wstep)               # ~ 1/(#.cm-2)
    Ei_density_on_grid = Ei/(ifwhm_on_grid*wstep)               # mW/cm3/sr/cm-1

    # Cut out of range points
    # ... you should see imin_broadened_wav_on_grid as a projection array, with
    # ... indexes where Intensity will be summed afterwards. 
    # ... here we project the out of range intensities to index -1 and len_grid+1
    # ... (this is faster than )
    len_grid = len(wavenumber)
    imin_broadened_wav_offset = imin_broadened_wav_on_grid
    imax_broadened_wav_offset = imax_broadened_wav_on_grid
    imin_broadened_wav_offset[imin_broadened_wav_offset < 0] = -1
    imax_broadened_wav_offset[imax_broadened_wav_offset > len_grid] = len_grid
    imin_broadened_wav_offset += 1
    imax_broadened_wav_offset += 1

    line2grid_projection = iwav_on_grid          # size N (number of lines)
    # ... deal with case where a new point is created (we dont want that)
    # ... just offset it by one unit (we're working with wavenumber_calc anyway,
    # ... these boundary points will be cropped by RADIS at the end of the calculation)
    line2grid_projection[line2grid_projection==len(wavenumber)] = len(wavenumber) - 1

    @jit(nopython=True)
    def rough_sum_on_grid():
        ''' Sum all lines linestrength density on a spectrum  '''
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
        k_rough_spectrum = np.zeros(len_grid+2)
        j_rough_spectrum = np.zeros(len_grid+2)
        for i in range(len(S_density_on_grid)):
            imin = imin_broadened_wav_offset[i]
            imax = imax_broadened_wav_offset[i]
            k_rough_spectrum[imin:imax+1] += S_density_on_grid[i]
            j_rough_spectrum[imin:imax+1] += Ei_density_on_grid[i]
        # crop out of range points
        return k_rough_spectrum[1:-1], j_rough_spectrum[1:-1]

    k_rough_spectrum, j_rough_spectrum = rough_sum_on_grid()

    return (k_rough_spectrum, j_rough_spectrum, S_density_on_grid, Ei_density_on_grid, 
            line2grid_projection)

    #
if __name__ == '__main__':

    from radis.test.lbl.test_broadening import _run_testcases
    print('Testing broadening.py:', _run_testcases(verbose=True, plot=True))
