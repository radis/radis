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
from os.path import exists

# %%


def calc_spectrum(
    wavenum_min=None,
    wavenum_max=None,
    wavelength_min=None,
    wavelength_max=None,
    Tgas=None,
    Tvib=None,
    Trot=None,
    pressure=1.01325,
    overpopulation=None,
    molecule="",
    isotope="all",
    mole_fraction=1,
    path_length=1,
    databank="fetch",
    medium="air",
    wstep=0.01,
    broadening_max_width=10,
    lineshape_optimization="DLM",
    name=None,
    use_cached=True,
    verbose=True,
    **kwargs
):
    """ Multipurpose function to calculate :class:`~radis.spectrum.spectrum.Spectrum`
    under equilibrium, or non-equilibrium, with or without overpopulation. 
    It's a wrapper to :class:`~radis.lbl.factory.SpectrumFactory` class. 
    For advanced used, please refer to the aforementionned class. 

    Parameters
    ----------

    wavenum_min: float(cm^-1) or `~astropy.units.quantity.Quantity`
        minimum wavenumber.
    wavenum_max: float(cm^-1) or `~astropy.units.quantity.Quantity`
        maximum wavenumber.

    wavelength_min: float(nm) or `~astropy.units.quantity.Quantity`
        minimum wavelength. Wavelength in ``'air'`` or 
        ``'vacuum'`` depending of the value of the parameter ``'medium='``
    wavelength_max: float(nm) or `~astropy.units.quantity.Quantity`
        maximum wavelength. Wavelength in ``'air'`` or 
        ``'vacuum'`` depending of the value of the parameter ``'medium='``

    Tgas: float(K) or `~astropy.units.quantity.Quantity`
        Gas temperature. If non equilibrium, is used for Ttranslational. 
        Default ``300`` K

    Tvib: float(K) or `~astropy.units.quantity.Quantity`
        Vibrational temperature. If ``None``, equilibrium calculation is run with Tgas

    Trot: float(K) or `~astropy.units.quantity.Quantity`
        Rotational temperature. If ``None``, equilibrium calculation is run with Tgas

    pressure: float(bar) or `~astropy.units.quantity.Quantity`
        partial pressure of gas. Default ``1.01325`` (1 atm)

    overpopulation: dict
        dictionary of overpopulation compared to the given vibrational temperature. 
        Ex with CO2::

            overpopulation = {
            '(00`0`0)->(00`0`1)': 2.5,
                '(00`0`1)->(00`0`2)': 1,
                '(01`1`0)->(01`1`1)': 1,
                '(01`1`1)->(01`1`2)': 1,
                }

    molecule: int, str, or ``None``
        molecule id (HITRAN format) or name. If ``None``, the molecule can be infered
        from the database files being loaded. See the list of supported molecules 
        in :py:data:`~radis.io.MOLECULES_LIST_EQUILIBRIUM`
        and :py:data:`~radis.io.MOLECULES_LIST_NONEQUILIBRIUM`. 
        Default ``None``. 

    isotope: int, list, str of the form ``'1,2'``, or ``'all'``
        isotope id (sorted by relative density: (eg: 1: CO2-626, 2: CO2-636 for CO2).
        See [HITRAN-2016]_ documentation for isotope list for all species. If ``'all'``,
        all isotopes in database are used (this may result in larger computation
        times!). Default ``'all'``

    mole_fraction: float
        database species mole fraction. Default ``1``

    path_length: float(cm) or `~astropy.units.quantity.Quantity`
        slab size. Default ``1`` cm.

    databank: str
        can be either: 

        - ``'fetch'``, to fetch automatically from [HITRAN-2016]_ through astroquery. 
        
        .. warning::
            
            [HITRAN-2016]_ is valid for low temperatures (typically < 700 K). For higher
            temperatures you may need [HITEMP-2010]_ 

        - the name of a valid database file, in which case the format is inferred. 
          For instance, ``.par`` is recognized as ``hitran/hitemp`` format. 

        - the name of a spectral database registered in your ``~/.radis`` 
          configuration file. This allows to use multiple database files.
          See :ref:`Configuration file <label_lbl_config_file>`.

        Default ``'fetch'``. See :class:`~radis.lbl.loader.DatabankLoader` for more 
        information on line databases, and :data:`~radis.misc.config.DBFORMAT` for 
        your ``~/.radis`` file format 

    medium: ``'air'``, ``'vacuum'``
        propagating medium when giving inputs with ``'wavenum_min'``, ``'wavenum_max'``. 
        Does not change anything when giving inputs in wavenumber. Default ``'air'``

    wstep: float (cm-1)
        Spacing of calculated spectrum. Default ``0.01`` cm-1.

    broadening_max_width: float (cm-1)
        Full width over which to compute the broadening. Large values will create
        a huge performance drop (scales as ~broadening_width^2 without DLM)
        The calculated spectral range is increased (by broadening_max_width/2
        on each side) to take into account overlaps from out-of-range lines.
        Default ``10`` cm-1.

    Other Parameters
    ----------------

    lineshape_optimization: int, ``None``, ``'DLM'``, or ``'auto'``.
        Optimizations for the calculation of the lineshapes:
        
            - If ``None``, all lineshapes are calculated at the same time (can 
              create memory errors). 
            - If ``int``, is given as the ``chunksize`` parameter of 
              :py:class:`~radis.lbl.factory.SpectrumFactory`` to split the line database 
              in several parts so that the number of ``lines * spectral grid points`` is 
              less than ``chunksize`` (reduces memory consumption). Typical values: 
              ``lineshape_optimization=1e6``.
            - If ``'DLM'``, only typical lineshapes are calculated. This can 
              result of speedups of orders of magnitude.  See more about DLM in 
              :ref:`Performance <label_lbl_performance>`. 
              
        Default ``'DLM'``.

    slit: float, str, or ``None``
        if float, FWHM of a triangular slit function. If str, path to an 
        experimental slit function. If None, no slit is applied. Default ``None``.

    plot: str
        any parameter such as 'radiance' (if slit is given), 'radiance_noslit', 
        'absorbance', etc...   Default ``None``

    name: str
        name of the case. If None, a unique ID is generated. Default ``None``

    use_cached: boolean
        use cached files for line database and energy database. Default ``True``

    verbose: boolean, or int
        If ``False``, stays quiet. If ``True``, tells what is going on. 
        If ``>=2``, gives more detailed messages (for instance, details of 
        calculation times). Default ``True``. 

    **kwargs: other inputs forwarded to SpectrumFactory
        For instance: ``warnings``. 
        See :class:`~radis.lbl.factory.SpectrumFactory` documentation for more 
        details on input. 
        For instance:

    pseudo_continuum_threshold: float
        if not 0, first calculate a rough approximation of the spectrum, then
        moves all lines whose linestrength intensity is less than this threshold
        of the maximum in a semi-continuum. Values above 0.01 can yield significant
        errors, mostly in highly populated areas. 80% of the lines can typically
        be moved in a continuum, resulting in 5 times faster spectra. If 0,
        no semi-continuum is used. Default 0.

    Returns
    -------
    
    s: :class:`~radis.spectrum.spectrum.Spectrum`
        Output spectrum.

        Use the :py:meth:`~radis.spectrum.spectrum.Spectrum.get` method to retrieve a 
        spectral quantity (``'radiance'``, ``'radiance_noslit'``, ``'absorbance'``, etc...)

        Or the :py:meth:`~radis.spectrum.spectrum.Spectrum.plot` method to plot it
        directly.

        See [1]_ to get an overview of all Spectrum methods

    References
    ----------

    .. [1] RADIS doc: `Spectrum how to? <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`__


    Examples
    --------

    Calculate a CO spectrum from the HITRAN database::

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

    Refer to the online :ref:`Examples <label_examples>` for more cases, and to 
    the :ref:`Spectrum page <label_spectrum>` for details on post-processing methods. 

    See Also
    --------
    
    :class:`~radis.lbl.factory.SpectrumFactory`, 
    the :ref:`Spectrum page <label_spectrum>`
    """

    # Check inputs

    # ... wavelengths / wavenumbers
    if (wavelength_min is not None or wavelength_max is not None) and (
        wavenum_min is not None or wavenum_max is not None
    ):
        raise ValueError("Wavenumber and Wavelength both given... it's time to choose!")

    if wavenum_min is None and wavenum_max is None:
        assert wavelength_max is not None
        assert wavelength_min is not None
        wavenum_min = nm2cm(wavelength_max)
        wavenum_max = nm2cm(wavelength_min)
    else:
        assert wavenum_min is not None
        assert wavenum_max is not None

    # ... temperatures

    if Tgas is None and Trot is None:
        raise ValueError(
            "Choose either Tgas (equilibrium) or Tvib / Trot (non equilibrium)"
        )

    if Tvib is None and Trot is not None or Tvib is not None and Trot is None:
        raise ValueError("Choose both Tvib and Trot")

    # ... others
    if databank is None:
        raise ValueError("Give a databank name")

    if not "save_memory" in kwargs:
        # no need to save intermediary results as
        # factory is used once only
        kwargs["save_memory"] = True

    if "chunksize" in kwargs:
        raise DeprecationWarning("use lineshape_optimization= instead of chunksize=")

    def _is_at_equilibrium():
        try:
            assert Tvib is None or Tvib == Tgas
            assert Trot is None or Trot == Tgas
            assert overpopulation is None
            if "self_absorption" in kwargs:
                assert kwargs["self_absorption"]  # == True
            return True
        except AssertionError:
            return False

    _equilibrium = _is_at_equilibrium()

    # which columns to keep when loading line database
    if kwargs["save_memory"] >= 2 and _equilibrium:
        drop_columns = "all"
    else:
        drop_columns = "auto"

    # Run calculations
    sf = SpectrumFactory(
        wavenum_min,
        wavenum_max,
        medium=medium,
        molecule=molecule,
        isotope=isotope,
        pressure=pressure,
        wstep=wstep,
        broadening_max_width=broadening_max_width,
        db_use_cached=use_cached,
        verbose=verbose,
        chunksize=lineshape_optimization,  #  if lineshape_optimization != 'auto' else None, #@EP: NotImplemented. DLM use all the time by default
        **kwargs
    )
    if databank == "fetch":  # mode to get databank without relying on  Line databases
        if _equilibrium:
            # Get line database from HITRAN
            # and partition functions from HAPI
            sf.fetch_databank(
                source="astroquery",
                format="hitran",
                parfuncfmt="hapi",  # use HAPI partition functions for equilibrium
                levelsfmt=None,  # no need to load energies
            )
        else:
            # Also get line database from HITRAN, and calculate partition functions
            # with energy levels from built-in constants (not all molecules
            # are supported!)
            sf.fetch_databank(
                source="astroquery",
                format="hitran",
                parfuncfmt="hapi",  # use HAPI partition functions for equilibrium
                levelsfmt="radis",  # built-in spectroscopic constants
            )
    elif exists(databank):
        # Guess format
        if databank.endswith(".par"):
            if verbose:
                print("Infered {0} is a HITRAN file.".format(databank))
            # If non-equilibrium we'll also need to load the energy levels.
            if _equilibrium:
                # Get partition functions from HAPI
                sf.load_databank(
                    path=databank,
                    format="hitran",
                    parfuncfmt="hapi",  # use HAPI partition functions for equilibrium
                    levelsfmt=None,  # no need to load energies
                    drop_columns=drop_columns,
                )
            else:
                # calculate partition functions with energy levels from built-in
                # constants (not all molecules are supported!)
                sf.load_databank(
                    path=databank,
                    format="hitran",
                    parfuncfmt="hapi",  # use HAPI partition functions for equilibrium
                    levelsfmt="radis",  # built-in spectroscopic constants
                    drop_columns=drop_columns,
                )
        else:
            raise ValueError(
                "Couldnt infer the format of the line database: {0}. ".format(databank)
                + "Create a user-defined database in your ~/.radis file "
                + "and define the format there."
            )

    else:  # manual mode: get from user-defined line databases defined in ~/.radis
        sf.load_databank(
            databank,
            load_energies=not _equilibrium,  # no need to load/calculate energies at eq.
            drop_columns=drop_columns,
        )

    #    # Get optimisation strategies
    #    if lineshape_optimization == 'auto':        # NotImplemented: finally we use DLM all the time as default.
    #        if len(sf.df0) > 1e5:
    #            lineshape_optimization = 'DLM'
    #        else:
    #            lineshape_optimization = None
    #        sf.params['chunksize'] = lineshape_optimization

    # Use the standard eq_spectrum / non_eq_spectrum functions
    if _equilibrium:
        s = sf.eq_spectrum(
            Tgas=Tgas, mole_fraction=mole_fraction, path_length=path_length, name=name
        )
    else:
        s = sf.non_eq_spectrum(
            Tvib=Tvib,
            Trot=Trot,
            Ttrans=Tgas,
            overpopulation=overpopulation,
            mole_fraction=mole_fraction,
            path_length=path_length,
            name=name,
        )

    return s


# --------------------------
if __name__ == "__main__":

    from radis.test.lbl.test_calc import _run_testcases

    print(_run_testcases(verbose=True))
