# -*- coding: utf-8 -*-
"""
Summary
-------

Public (front-end) functions to calculate Spectrum with HITRAN / CDSD databanks.
Uses the SpectrumFactory classe from `factory.py`, Spectrum from `spectrum.py`
and line survey from `line_survey.py` 

Routine Listing
---------------

:func:`~radis.lbl.calc.calc_spectrum`

-------------------------------------------------------------------------------

"""

from __future__ import print_function, absolute_import, division, unicode_literals

from radis.lbl.factory import SpectrumFactory
from radis.spectrum import Spectrum, is_spectrum
from radis.phys.blackbody import planck, planck_wn
from radis.phys.convert import nm2cm
from radis.misc.utils import DatabankNotFound
import numpy as np
from numpy import exp, arange, allclose, abs, diff, zeros_like, ones_like
from warnings import warn

# %%


def calc_spectrum(wavenum_min=None,
                  wavenum_max=None,
                  wavelength_min=None,
                  wavelength_max=None,
                  Tgas=None,
                  Tvib=None,
                  Trot=None,
                  pressure=1.01325,
                  overpopulation=None,
                  molecule='',
                  isotope='all',
                  mole_fraction=1,
                  path_length=1,
                  databank='fetch',
                  slit=None,
                  plot=None,
                  medium='air',
                  wstep=0.01,
                  name=None,
                  use_cached=True,
                  **kwargs):
    ''' Multipurpose function to calculate :class:`~radis.spectrum.spectrum.Spectrum`
    under equilibrium, or non-equilibrium, with or without overpopulation. 
    It's a wrapper to :class:`~radis.lbl.factory.SpectrumFactory` class. 
    For advanced used, please refer to the aforementionned class. 

    Parameters
    ----------

    wavenum_min: float (cm-1)
        minimum wavenumber to be processed in cm^-1
    wavenum_max: float (cm-1)
        maximum wavenumber to be processed in cm^-1

    wavelength_min: float (nm)
        minimum wavelength to be processed in nm. Wavelength if in 'air' or 
        'vacuum' depending of the value of the parameter 'medium='

    wavelength_max: float (nm)
        maximum wavelength to be processed in nm. Wavelength if in 'air' or 
        'vacuum' depending of the value of the parameter 'medium='

    Tgas: float (K)
        Gas temperature. If non equilibrium, is used for Ttranslational. 
        Default 300K

    Tvib: float (K)
        Vibrational temperature. If None, equilibrium calculation is run with Tgas

    Trot: float (K)
        Rotational temperature. If None, equilibrium calculation is run with Tgas

    pressure: (bar)
        partial pressure of gas in bar. Default 1.01325 (1 atm)

    overpopulation: dict
        dictionary of overpopulation compared to the given vibrational temperature. 
        Ex with CO2::

        >>> overpopulation = {
        >>>    '(00`0`0)->(00`0`1)': 2.5,
        >>>     '(00`0`1)->(00`0`2)': 1,
        >>>     '(01`1`0)->(01`1`1)': 1,
        >>>     '(01`1`1)->(01`1`2)': 1,
        >>>     }

    molecule: int, str, or ``None``
        molecule id (HITRAN format) or name. If None, the molecule is infered
        from the database files being loaded. It is just good practice to fill
        it. Default ``None``.

    isotope: int, list, str of the form '1,2', or 'all'
        isotope id (sorted by relative density: (eg: 1: CO2-626, 2: CO2-636 for CO2).
        See HITRAN documentation for isotope list for all species. If 'all',
        all isotopes in database are used (this may result in larger computation
        times!). Default 'all'

    mole_fraction: float
        database species mole fraction. Default 1

    path_length: float (cm)
        slab size. Default 1.

    databank: str
        can be either: 

        - the name of a spectral database to use (must be registered 
          in your ``~/.radis``). See :ref:`Configuration file <label_lbl_config_file>`.

        - ``'fetch'``, to fetch automatically from [HITRAN-2016]_ through astroquery. 
        
        .. warning::
            
            [HITRAN-2016]_ is validd for low temperatures (typically < 700 K). For higher
            temperatures you may need [HITEMP-2010]_ 

        Default ``'fetch'``. See :class:`~radis.lbl.loader.DatabankLoader` for more 
        information on line databases, and :data:`~radis.misc.config.DBFORMAT` for 
        your ``~/.radis`` file format 

    slit: float, str, or ``None``
        if float, FWHM of a triangular slit function. If str, path to an 
        experimenta slit function. If None, no slit is applied. Default ``None``.

    plot: str
        any parameter such as 'radiance' (if slit is given), 'radiance_noslit', 
        'absorbance', etc...   Default ``None``

    medium: 'air', vacuum'
        propagating medium: choose whether to return wavelength in air or vacuum.
        Note that the  input wavelength range ("wavenum_min", "wavenum_max") changes. 
        Default 'air'

    Other Parameters
    ----------------

    wstep: cm-1
        Spacing of calculated spectrum. Default 0.01 cm-1

    name: str
        name of the case. If None, a unique ID is generated. Default ``None``

    use_cached: boolean
        use cached files for line database and energy database. Default ``True``

    **kwargs: other inputs forwarded to SpectrumFactory
        See :class:`~radis.lbl.factory.SpectrumFactory` documentation for more details on input

    Example:

    pseudo_continuum_threshold: float
        if not 0, first calculate a rough approximation of the spectrum, then
        moves all lines whose linestrength intensity is less than this threshold
        of the maximum in a semi-continuum. Values above 0.01 can yield significant
        errors, mostly in highly populated areas. 80% of the lines can typically
        be moved in a continuum, resulting in 5 times faster spectra. If 0,
        no semi-continuum is used. Default 0.

    Returns
    -------

    Returns a :class:`~radis.spectrum.spectrum.Spectrum` object::
    
    Use the :meth:`~radis.spectrum.spectrum.Spectrum.get` method to get something 
    among ['radiance', 'radiance_noslit', 'absorbance', etc...]

    Or directly the :meth:`~radis.spectrum.spectrum.Spectrum.plot` method 
    to plot it

    See [1]_ to get an overview of all Spectrum methods

    Notes
    -----

    Calculate line strenghts correcting the CDSD reference one. Then call
    the main routine that sums over all lines


    References
    ----------

    .. [1] RADIS doc: `Spectrum how to? <http://radis.readthedocs.io/en/latest/#the-spectrum-class>`__


    Examples
    --------

    Calculate a CO spectrum from the HITRAN database ::

        s = calc_spectrum(1900, 2300,         # cm-1
                          molecule='CO',
                          isotope='1,2,3',
                          pressure=1.01325,   # bar
                          Tgas=1000, 
                          mole_fraction=0.1, 
                          )
        s.apply_slit(0.5, 'nm')
        s.plot('radiance')
        
    This example uses the :py:meth:`~radis.spectrum.spectrum.Spectrum.apply_slit` 
    and :py:meth:`~radis.spectrum.spectrum.Spectrum.plot` methods. See also
    :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey`:: 
        
        s.line_survey(overlay='radiance')

    Refer to the online :ref:`Examples <label_examples>` for more cases. 

    '''

    # Check inputs

    # ... wavelengths / wavenumbers
    if ((wavelength_min is not None or wavelength_max is not None) and
            (wavenum_min is not None or wavenum_max is not None)):
        raise ValueError("Wavenumber and Wavelength both given... it's time to choose!")

    if (wavenum_min is None and wavenum_max is None):
        assert(wavelength_max is not None)
        assert(wavelength_min is not None)
        wavenum_min = nm2cm(wavelength_max)
        wavenum_max = nm2cm(wavelength_min)
    else:
        assert(wavenum_min is not None)
        assert(wavenum_max is not None)

    # ... temperatures

    if Tgas is None and Trot is None:
        raise ValueError(
            'Choose either Tgas (equilibrium) or Tvib / Trot (non equilibrium)')

    if Tvib is None and Trot is not None or Tvib is not None and Trot is None:
        raise ValueError('Choose both Tvib and Trot')
        
    if databank == 'fetch' and ((Tgas is not None and Tgas > 700) or 
                                (Tvib is not None and Tvib > 700) or 
                                (Trot is not None and Trot > 700)):
        warn("HITRAN is valid for low temperatures (typically < 700 K). "+\
             "For higher temperatures you may need HITEMP or CDSD. See the "+\
             "'databank=' parameter")
        
    # ... others
    if databank is None:
        raise ValueError('Give a databank name')

    if not 'save_memory' in kwargs:
        # no need to save intermediary results as 
         # factory is used once only
        kwargs['save_memory'] = True
       

    def _is_at_equilibrium():
        try:
            assert Tvib is None or Tvib == Tgas
            assert Trot is None or Trot == Tgas
            assert overpopulation is None
            if 'self_absorption' in kwargs:
                assert kwargs['self_absorption']  # == True
            return True
        except AssertionError:
            return False

    _equilibrium = _is_at_equilibrium()

    # Run calculations
    sf = SpectrumFactory(wavenum_min, wavenum_max, medium=medium,
                         molecule=molecule, isotope=isotope,
                         pressure=pressure,
                         wstep=wstep,
                         db_use_cached=use_cached,
                         **kwargs)
    if databank == 'fetch':       # mode to get databank without relying on  Line databases
        if _equilibrium:
            # Get line database from HITRAN
            # and partition functions from HAPI
            sf.fetch_databank(source='astroquery', format='hitran',
                              parfuncfmt='hapi',   # use HAPI partition functions for equilibrium
                              levelsfmt=None,      # no need to load energies
                              )
        else:
            # Also get line database from HITRAN, and calculate partition functions
            # with energy levels from built-in constants (not all molecules
            # are supported!)
            sf.fetch_databank(source='astroquery', format='hitran',
                              parfuncfmt='hapi',   # use HAPI partition functions for equilibrium
                              levelsfmt='radis',     # built-in spectroscopic constants
                              )
    else:   # manual mode: get from files defined in .radis
        sf.load_databank(databank, 
                         load_energies=not _equilibrium   # no need to load/calculate energies at eq.
                         )

    # Use the standard eq_spectrum / non_eq_spectrum functions
    if _equilibrium:
        s = sf.eq_spectrum(Tgas=Tgas, mole_fraction=mole_fraction, path_length=path_length,
                           name=name)
    else:
        s = sf.non_eq_spectrum(Tvib=Tvib, Trot=Trot, Ttrans=Tgas,
                               overpopulation=overpopulation,
                               mole_fraction=mole_fraction,
                               path_length=path_length,
                               name=name)

    if slit is not None:
        s.apply_slit(slit)
    if plot is not None and plot is not False:
        s.plot(plot)
    return s


# --------------------------
if __name__ == '__main__':

    from radis.test.lbl.test_calc import _run_testcases
    print(_run_testcases(verbose=True))
