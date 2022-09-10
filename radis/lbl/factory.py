# -*- coding: utf-8 -*-
"""Contains the :py:class:`~radis.lbl.factory.SpectrumFactory` class, which is
the core of the RADIS Line-by-Line module.

Examples
--------

Calculate a CO Spectrum, fetching the lines from HITRAN ::

    # This is how you get a spectrum (see spectrum.py for front-end functions
    # that do just that)
    sf = SpectrumFactory(2125, 2249.9,
                         molecule='CO',
                         isotope=1,
                         cutoff=1e-30,   # for faster calculations. See
                                         # `plot_linestrength_hist` for more details
                         **kwargs)
    sf.fetch_databank()        # autodownload from HITRAN
    s = sf.eq_spectrum(Tgas=300)
    s.plot('abscoeff')          # opacity

    # Here we get some extra informations:
    s.plot('radiance', wunit='nm',
                       Iunit='µW/cm2/sr/nm',  # Iunit is arbitrary. Use whatever makes sense
                       show_points=True)  # show_points to have an idea of the resolution

See :py:mod:`radis.test.lbl.test_factory` for more examples


Routine Listing
---------------
*(the following sections are marked on the slider bar in the Spyder IDE with
the # XXX symbol)*

PUBLIC METHODS

- :meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`             >>> calc equilibrium spectrum
- :meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum`         >>> calc non equilibrium spectrum
- :meth:`~radis.lbl.factory.SpectrumFactory.optically_thin_power`    >>> get total power (equilibrium or non eq)

Most methods are written in inherited class with the following inheritance scheme:

:py:class:`~radis.lbl.loader.DatabankLoader` > :py:class:`~radis.lbl.base.BaseFactory` >
:py:class:`~radis.lbl.broadening.BroadenFactory` > :py:class:`~radis.lbl.bands.BandFactory` >
:py:class:`~radis.lbl.factory.SpectrumFactory`


.. inheritance-diagram:: radis.lbl.factory.SpectrumFactory
   :parts: 1


Notes
-----

Debugging:

Remember one can access all locals variable in a function either by running
%debug after a crash in IPython. Or by adding the line `globals().update(locals())`
at any point


Performance:

Fast version: iterate over chunks of dataframe
Note that we can't use the full dataframe because it's too big and takes too much memory
See :ref:`Performance <label_lbl_performance>` for more details.


for Developers:

- To implement new database formats, see the databases parsers in cdsd.py / hitran.py,
  and the partition function interpolators / calculators methods of SpectrumFactory:
  :py:meth:`~radis.lbl.loader.DatabankLoader._build_partition_function_calculator` and
  :py:meth:`~radis.lbl.loader.DatabankLoader._build_partition_function_interpolator`


----------
"""
from typing import Union
from warnings import warn

import astropy.units as u
import numpy as np
from numpy import arange, exp
from scipy.constants import c
from scipy.optimize import OptimizeResult

from radis import version
from radis.db import MOLECULES_LIST_EQUILIBRIUM, MOLECULES_LIST_NONEQUILIBRIUM
from radis.db.classes import get_molecule, get_molecule_identifier

try:  # Proper import
    from .bands import BandFactory
    from .base import get_wavenumber_range
except ImportError:  # if ran from here
    from radis.lbl.bands import BandFactory
    from radis.lbl.base import get_wavenumber_range

from radis.misc.basics import flatten, is_float, is_range, list_if_float, round_off
from radis.misc.utils import Default
from radis.phys.constants import k_b
from radis.phys.convert import conv2
from radis.phys.units import convert_universal
from radis.phys.units_astropy import convert_and_strip_units
from radis.spectrum.equations import calc_radiance
from radis.spectrum.spectrum import Spectrum

c_cm = c * 100

# %% Main functions


class SpectrumFactory(BandFactory):
    """A class to put together all functions related to loading CDSD / HITRAN
    databases, calculating the broadenings, and summing over all the lines.

    Parameters
    ----------
    wmin, wmax : ``float`` or `~astropy.units.quantity.Quantity`
        a hybrid parameter which can stand for minimum (maximum) wavenumber or minimum
        (maximum) wavelength depending upon the unit accompanying it. If dimensionless,
        ``wunit`` is considered as the accompanying unit.
    wunit: ``'nm'``, ``'cm-1'``
        the unit accompanying wmin and wmax. Can only be passed with wmin
        and wmax. Default is ``"cm-1"``.
    wavenum_min, wavenum_max: ``float(cm^-1)`` or `~astropy.units.quantity.Quantity`
        minimum (maximum) wavenumber to be processed in :math:`cm^{-1}`.
        use astropy.units to specify arbitrary inverse-length units.
    wavelength_min, wavelength_max : ``float(nm)`` or `~astropy.units.quantity.Quantity`
        minimum (maximum) wavelength to be processed in :math:`nm`. This wavelength
        can be in ``'air'`` or ``'vacuum'`` depending on the value of the parameter
        ``medium=``.
        use astropy.units to specify arbitrary length units.
    pressure: ``float(bar)`` or `~astropy.units.quantity.Quantity`
        partial pressure of gas in bar. Default ``1.01325`` (1 atm).
        use astropy.units to specify arbitrary pressure units.
        For example, ``1013.25 * u.mbar``.
    mole_fraction: ``float``  [ 0 - 1]
        species mole fraction. Default ``1``. Note that the rest of the gas
        is considered to be air for collisional broadening.
    path_length: ``float(cm)`` or `~astropy.units.quantity.Quantity`
        path length in cm. Default ``1``.
        use astropy.units to specify arbitrary length units.
    molecule: ``int``, ``str``, or ``None``
        molecule id (HITRAN format) or name. If ``None``, the molecule can be infered
        from the database files being loaded. See the list of supported molecules
        in :py:data:`~radis.db.MOLECULES_LIST_EQUILIBRIUM`
        and :py:data:`~radis.db.MOLECULES_LIST_NONEQUILIBRIUM`.
        Default ``None``.
    isotope: ``int``, ``list``, ``str`` of the form ``'1,2'``, or ``'all'``
        isotope id (sorted by relative density: (eg: 1: CO2-626, 2: CO2-636 for CO2).
        See HITRAN documentation for isotope list for all species. If 'all',
        all isotopes in database are used (this may result in larger computation
        times!). Default ``'all'``
    medium: ``'air'``, ``'vacuum'``
        propagating medium when giving inputs with ``'wavenum_min'``, ``'wavenum_max'``.
        Does not change anything when giving inputs in wavenumber. Default ``'air'``
    diluent: ``str`` or ``dictionary``
            can be a string of a single diluent or a dictionary containing diluent
            name as key and its mole_fraction as value. Default ``air``.

    Other Parameters
    ----------------
    Tref: K
        Reference temperature for calculations (linestrength temperature
        correction). HITRAN database uses 296 Kelvin. Default 296 K
    self_absorption: boolean
        Compute self absorption. If ``False``, spectra are optically thin. Default ``True``.
    truncation: float (:math:`cm^{-1}`)
        Half-width over which to compute the lineshape, i.e. lines are truncated
        on each side after ``truncation`` (:math:`cm^{-1}`) from the line center.
        If ``None``, use no truncation (lineshapes spread on the full spectral range).
        Default is ``300`` :math:`cm^{-1}`

        .. note::
         Large values (> ``50``) can induce a performance drop (computation of lineshape
         typically scale as :math:`~truncation ^2` ). The default ``300`` was
         chosen to maintain a good accuracy, and still exhibit the sub-Lorentzian
         behavior of most lines far (few hundreds :math:`cm^{-1}`) from the line center.
    neighbour_lines: float (:math:`cm^{-1}`)
        The calculated spectral range is increased (by ``neighbour_lines`` cm-1
        on each side) to take into account overlaps from out-of-range lines.
        Default is ``0`` :math:`cm^{-1}`.​
    wstep: float (cm-1) or `'auto'`
        Resolution of wavenumber grid. Default ``0.01`` cm-1.
        If `'auto'`, it is ensured that there
        are slightly more points for each linewidth than the value of ``"GRIDPOINTS_PER_LINEWIDTH_WARN_THRESHOLD"``
        in :py:attr:`radis.config`  (``~/radis.json``)

        .. note::
            wstep = 'auto' is optimized for performances while ensuring accuracy,
            but is still experimental in 0.9.30. Feedback welcome!
    cutoff: float (~ unit of Linestrength: cm-1/(#.cm-2))
        discard linestrengths that are lower that this, to reduce calculation
        times. ``1e-27`` is what is generally used to generate databases such as
        CDSD. If ``0``, no cutoff. Default ``1e-27``.
    parsum_mode: 'full summation', 'tabulation'
        how to compute partition functions, at nonequilibrium or when partition
        function are not already tabulated. ``'full summation'`` : sums over all
        (potentially millions) of rovibrational levels. ``'tabulation'`` :
        builds an on-the-fly tabulation of rovibrational levels (500 - 4000x faster
        and usually accurate within 0.1%). Default ``full summation'``

        .. note::
            parsum_mode= 'tabulation'  is new in 0.9.30, and makes nonequilibrium
            calculations of small spectra extremelly fast. Will become the default
            after 0.9.31.
    pseudo_continuum_threshold: float
        if not ``0``, first calculate a rough approximation of the spectrum, then
        moves all lines whose linestrength intensity is less than this threshold
        of the maximum in a semi-continuum. Values above 0.01 can yield significant
        errors, mostly in highly populated areas. 80% of the lines can typically
        be moved in a continuum, resulting in 5 times faster spectra. If ``0``,
        no semi-continuum is used. Default ``0``.
    save_memory: boolean
        if ``True``, removes databases calculated by intermediate functions (for
        instance, delete the full database once the linestrength cutoff criteria
        was applied). This saves some memory but requires to reload the database
        & recalculate the linestrength for each new parameter. Default ``False``.
    export_populations: ``'vib'``, ``'rovib'``, ``None``
        if not None, store populations in Spectrum. Either store vibrational
        populations ('vib') or rovibrational populations ('rovib'). Default ``None``
    export_lines: boolean
        if ``True``, saves details of all calculated lines in Spectrum. This is
        necessary to later use :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey`,
        but can take some space. Default ``False``.
    chunksize: int, or ``None``
        Splits the lines database in several chuncks during calculation, else
        the multiplication of lines over all spectral range takes too much memory
        and slows the system down. Chunksize let you change the default chunck
        size. If ``None``, all lines are processed directly. Usually faster but
        can create memory problems. Default ``None``
    optimization : ``"simple"``, ``"min-RMS"``, ``None``
        If either ``"simple"`` or ``"min-RMS"`` LDM optimization for lineshape calculation is used:
        - ``"min-RMS"`` : weights optimized by analytical minimization of the RMS-error (See: [Spectral-Synthesis-Algorithm]_)
        - ``"simple"`` : weights equal to their relative position in the grid

        If using the LDM optimization, broadening method is automatically set to ``'fft'``.
        If ``None``, no lineshape interpolation is performed and the lineshape of all lines is calculated.

        Refer to [Spectral-Synthesis-Algorithm]_ for more explanation on the LDM method for lineshape interpolation.

        Default ``"min-RMS"``
    folding_thresh: float
        Folding is a correction procedure thet is applied when the lineshape is calculated with
        the ``fft`` broadening method and the linewidth is comparable to ``wstep``, that prevents
        sinc(v) modulation of the lineshape. Folding continues until the lineshape intensity
        is below ``folding_threshold``. Setting to 1 or higher effectively disables folding correction.

        Range: 0.0 < folding_thresh <= 1.0
        Default: 1e-6
    zero_padding: int
        Zero padding is used in conjunction with the ``fft`` broadening method to prevent circular
        convolution at the cost of performance. When set to -1, padding is set equal to the spectrum length,
        which guarantees a linear convolution.

        Range: 0 <= zero_padding <= len(w), or zero_padding = -1
        Default: -1
    broadening_method: ``"voigt"``, ``"convolve"``, ``"fft"``
        Calculates broadening with a direct voigt approximation ('voigt') or
        by convoluting independantly calculated Doppler and collisional
        broadening ('convolve'). First is much faster, 2nd can be used to
        compare results. This SpectrumFactory parameter can be manually
        adjusted a posteriori with::

            sf = SpectrumFactory(...)
            sf.params.broadening_method = 'voigt'

        Fast fourier transform ``'fft'`` is only available if using the LDM lineshape
        calculation ``optimization``. Because the LDM convolves all lines at the same time,
        and thus operates on large arrays, ``'fft'`` becomes more appropriate than
        convolutions in real space (``'voit'``, ``'convolve'`` )

        By default, use ``"fft"`` for any ``optimization``, and ``"voigt"`` if
        optimization is ``None`` .
    warnings: bool, or one of ``['warn', 'error', 'ignore']``, dict
        If one of ``['warn', 'error', 'ignore']``, set the default behaviour
        for all warnings. Can also be a dictionary to set specific warnings only.
        Example::

            warnings = {'MissingSelfBroadeningWarning':'ignore',
                        'NegativeEnergiesWarning':'ignore',
                        'HighTemperatureWarning':'ignore'}

        See :py:data:`~radis.misc.warning.default_warning_status` for more
        information.
    verbose: boolean, or int
        If ``False``, stays quiet. If ``True``, tells what is going on.
        If ``>=2``, gives more detailed messages (for instance, details of
        calculation times). Default ``True``.

    Examples
    --------

    An example using :class:`~radis.lbl.factory.SpectrumFactory`,
    :meth:`~radis.lbl.loader.DatabankLoader.load_databank`, the
    :class:`~radis.spectrum.spectrum.Spectrum` methods, and
    :py:mod:`~astropy.units` ::

        from radis import SpectrumFactory
        from astropy import units as u
        sf = SpectrumFactory(wavelength_min=4165 * u.nm,
                             wavelength_max=4200 * u.nm,
                             isotope='1,2',
                             truncation=10,  # cm-1
                             optimization=None,
                             medium='vacuum',
                             verbose=1,    # more for more details
                             )
        sf.load_databank('HITRAN-CO2-TEST')        # predefined in ~/radis.json
        s = sf.eq_spectrum(Tgas=300 * u.K, path_length=1 * u.cm)
        s.rescale_path_length(0.01)    # cm
        s.plot('radiance_noslit', Iunit='µW/cm2/sr/nm')

    Refer to the online :ref:`Examples <label_examples>` for more cases.

    .. minigallery:: radis.SpectrumFactory
        :add-heading:

    Notes
    -----

    .. inheritance-diagram:: radis.lbl.factory.SpectrumFactory
       :parts: 1

    High-level wrapper to SpectrumFactory:

    - :func:`~radis.lbl.calc.calc_spectrum`

    Main Methods:

    - :meth:`~radis.lbl.loader.DatabankLoader.load_databank`,
    - :meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`,
    - :meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum`

    For advanced use:

    - :meth:`~radis.lbl.loader.DatabankLoader.fetch_databank`,
    - :meth:`~radis.lbl.loader.DatabankLoader.init_databank`,
    - :meth:`~radis.lbl.loader.DatabankLoader.init_database`,
    - :meth:`~radis.lbl.bands.BandFactory.eq_bands`,
    - :meth:`~radis.lbl.bands.BandFactory.non_eq_bands`

    Inputs and parameters can be accessed a posteriori with :

    - :py:attr:`~radis.lbl.loader.DatabankLoader.input` : physical input
    - :py:attr:`~radis.lbl.loader.DatabankLoader.params` : computational parameters
    - :py:attr:`~radis.lbl.loader.DatabankLoader.misc` : miscallenous parameters (don't change output)

    See Also
    --------
    :func:`~radis.lbl.calc.calc_spectrum`
    """

    # TODO: make it possible to export both 'vib' and 'rovib'

    # TODO
    # -------
    #
    # This code is still designed for one-molecule only: some parameters
    # (e.g.: molar mass used in Doppler broadening) are defined in SpectrumFactory
    # rather than being a column in the dataframe. This should be changed to
    # allow dealing with polymolecules databases

    # TODO
    # store everything in a self.var class instead of self.[] directly

    # TODO : move Tref in load_databank / fetch_databank only

    def __init__(
        self,
        wmin=None,
        wmax=None,
        wunit=Default("cm-1"),
        wavenum_min=None,
        wavenum_max=None,
        wavelength_min=None,
        wavelength_max=None,
        Tref=296,
        pressure=1.01325,
        mole_fraction=1,
        path_length=1,
        wstep=0.01,
        molecule=None,
        isotope="all",
        medium="air",
        truncation=Default(50),
        neighbour_lines=0,
        pseudo_continuum_threshold=0,
        self_absorption=True,
        chunksize=None,
        optimization="simple",
        folding_thresh=1e-6,
        zero_padding=-1,
        broadening_method="voigt",
        cutoff=0,
        parsum_mode="full summation",
        verbose=True,
        warnings=True,
        save_memory=False,
        export_populations=None,
        export_lines=False,
        emulate_gpu=False,
        diluent="air",
        **kwargs,
    ):

        # Initialize BandFactory object
        super(SpectrumFactory, self).__init__()

        # Check inputs (deal with deprecated format)
        if medium not in ["air", "vacuum"]:
            raise ValueError("Wavelength must be one of: 'air', 'vacuum'")
        kwargs0 = kwargs  # kwargs is used to deal with Deprecated names
        if "db_use_cached" in kwargs:
            warn(
                DeprecationWarning(
                    "db_use_cached removed from SpectrumFactory init and moved in load/fetch_databank()"
                )
            )
            kwargs0.pop("db_use_cached")
        if "lvl_use_cached" in kwargs:
            warn(
                DeprecationWarning(
                    "lvl_use_cached removed from SpectrumFactory init and moved in load/fetch_databank()"
                )
            )
            kwargs0.pop("lvl_use_cached")
        if "broadening_max_width" in kwargs:  # changed in 0.9.30
            broadening_max_width = kwargs["broadening_max_width"]
            raise (
                DeprecationWarning(
                    "`broadening_max_width`` (lineshape full-width, also used to compute the effect of neighbour lines) was replaced by `truncation` (lineshape half-width) and `neighbour_lines` (wavenumber range extension on each side). "
                    + f"To keep the current behavior, replace `broadening_max_width={broadening_max_width}` with "
                    + f"`truncation={broadening_max_width/2}, neighbour_lines={broadening_max_width/2}`. "
                    + "We recommended, for most cases: `truncation=300, neighbour_lines=0}`"
                )
            )
        for boolarg in [self_absorption, save_memory, export_lines]:
            if boolarg not in [True, False]:
                raise ValueError(
                    f"Expected boolean parameter. Got `{boolarg.__repr__()}`. Use `True` or `False`"
                )

        if kwargs0 != {}:
            raise TypeError(
                "__init__() got an unexpected keyword argument '{0}'".format(
                    list(kwargs0)[0]
                )
            )

        if not 0 <= pseudo_continuum_threshold < 1:
            raise ValueError("pseudo_continuum_threshold should be in [0-1]")
        if export_populations not in ["vib", "rovib", False, None] and not isinstance(
            export_populations, list
        ):
            raise ValueError(
                "export_populations must be one of 'vib', 'rovib', "
                + "or 'False'. Got '{0}'".format(export_populations)
            )

        # calculate waveranges
        # --------------------

        # Get wavenumber, based on whatever was given as input.
        wavenum_min, wavenum_max, input_wunit = get_wavenumber_range(
            wmin,
            wmax,
            wunit,
            wavenum_min,
            wavenum_max,
            wavelength_min,
            wavelength_max,
            medium,
            return_input_wunit=True,
        )
        # ... Make default Spectrum's output unit consistent with the input waverange
        # ... see https://github.com/radis/radis/issues/456
        # (note: may be overwritten by user after Factory creation)
        self.input_wunit = input_wunit

        # Storing inital value of wstep if wstep != "auto"
        self._wstep = wstep

        # Set default variables from config:
        import radis

        self._sparse_ldm = radis.config[
            "SPARSE_WAVERANGE"
        ]  # target value (can be 'auto'), stored for next calculatinos
        self.params["sparse_ldm"] = radis.config[
            "SPARSE_WAVERANGE"
        ]  # value evaluated at each new spectrum calculation

        # Init variables
        # --------------

        # Get molecule name
        if isinstance(molecule, int):
            molecule == get_molecule(molecule)
        if molecule is not None:
            if (
                molecule
                not in MOLECULES_LIST_EQUILIBRIUM + MOLECULES_LIST_NONEQUILIBRIUM
            ):
                raise ValueError(
                    "Unsupported molecule: {0}.\n".format(molecule)
                    + "Supported molecules are:\n - under equilibrium: {0}".format(
                        MOLECULES_LIST_EQUILIBRIUM
                    )
                    + "\n- under nonequilibrium: {0}".format(
                        MOLECULES_LIST_NONEQUILIBRIUM
                    )
                    + "\n\nNote that RADIS now has ExoMol support, but not all ExoMol molecules are referenced in RADIS. If a molecule is available in ExoMol but does not appear in RADIS yet, please contact the RADIS team or write on https://github.com/radis/radis/issues/319"
                )

        # Store isotope identifier in str format (list wont work in database queries)
        if not isinstance(isotope, str):
            isotope = ",".join([str(k) for k in list_if_float(isotope)])

        # If molecule present in diluent, raise error
        if (isinstance(diluent, str) and diluent == molecule) or (
            isinstance(diluent, dict) and molecule in diluent.keys()
        ):
            raise KeyError(
                "{0} is being called as molecule and diluent, please remove it from diluent.".format(
                    molecule
                )
            )

        # Initialize input conditions
        self.input.wavenum_min = wavenum_min
        self.input.wavenum_max = wavenum_max
        self.input.Tref = convert_and_strip_units(Tref, u.K)
        self.input.pressure_mbar = convert_and_strip_units(pressure, u.bar) * 1e3
        self.input.mole_fraction = mole_fraction

        self.input.path_length = convert_and_strip_units(path_length, u.cm)
        self.input.molecule = (
            molecule  # if None, will be overwritten after reading database
        )
        self.input.state = "X"  # for the moment only ground-state is used
        # (but the code is electronic state aware)
        self.input.isotope = (
            isotope  # if 'all', will be overwritten after reading database
        )
        self.input.self_absorption = self_absorption

        # Initialize computation variables
        self.params.wstep = wstep
        self.params.pseudo_continuum_threshold = pseudo_continuum_threshold
        self.params.diluent = diluent
        self._diluent = None

        if cutoff is None:
            # If None, use no cutoff : https://github.com/radis/radis/pull/259
            cutoff = 0
        self.params.cutoff = cutoff
        self.params.parsum_mode = parsum_mode

        # Time Based variables
        self.verbose = verbose

        if truncation == 0:
            raise ValueError(
                "Lineshape truncation must be > 0. If you want no truncation (compute lineshape on the full spectral range), use `truncation=None`. \nNote (advanced) : no truncation is not physically more accurate. Most molecules exhibit a sub-lorentzian behavior far from the line centers. A truncation at around 40-50 cm-1 is a good choice"
            )
        elif (
            truncation is not None and not isinstance(truncation, Default)
        ) and truncation < 0:
            raise ValueError(
                "Lineshape truncation can't be negative. Truncation must be > 0 or None to compute lineshape on the full spectral range"
            )

        self.misc.export_lines = export_lines
        self.misc.export_populations = export_populations

        if broadening_method == "fft":
            if isinstance(truncation, Default):
                truncation = None
            elif truncation is not None and broadening_method == "fft":
                raise NotImplementedError(
                    "Lines cannot be truncated with `broadening_method='fft'`. Use `broadening_method='voigt'`"
                )
        elif (
            broadening_method == "voigt"
            and truncation is None
            and optimization is not None
        ):
            raise NotImplementedError(
                "Currently `broadening_method='voigt'` doesn't support computation of lineshape on the full spectral range, use `broadening_method='fft'` instead or use a truncation value > 0"
            )

        if isinstance(truncation, Default):
            truncation = truncation.value

        self.params.truncation = self.truncation = truncation  # line truncation
        # self.params.truncation is the input, self.truncation will be the value (different from input if input was None)
        self.params.neighbour_lines = neighbour_lines  # including neighbour lines

        # # reduce neighbour_lines if unnecessary
        # if truncation and truncation < neighbour_lines:
        #     self.warn
        # Define max range on which to load lines (required for fetch_databank / load_databank)
        self.params.wavenum_min_calc = wavenum_min - neighbour_lines
        self.params.wavenum_max_calc = wavenum_max + neighbour_lines
        #  note @dev : if neighbour_lines changes during the SpectrumFactory loop (e.g:  user sets sf.params.neighbour_lines= ;  the database wont be updated.
        # (databank should be reloaded). We store and check later to ensure it is not changed afterwards
        self._neighbour_lines = neighbour_lines

        self.params.broadening_method = broadening_method
        self.params.optimization = optimization
        self.params.folding_thresh = folding_thresh
        self.misc.zero_padding = zero_padding

        # used to split lines into blocks not too big for memory
        self.misc.chunksize = chunksize
        # Other parameters:
        self.save_memory = save_memory
        self.autoupdatedatabase = False  # a boolean to automatically store calculated
        # spectra in a Spectrum database. See init_database
        # for more info
        self.autoretrievedatabase = False  # a boolean to automatically retrieve
        # spectra from database instead of
        # calculating them
        self.SpecDatabase = None  # the database to store spectra. Not to be confused
        # with the databank where lines are stored
        self.database = None  # path to previous database
        # Warnings
        # --------

        # Add cutoff / threshold warning values for the different optimisations
        # (linestrength cutoff, broadening cutoff). These cannot be changed
        # from the Factory input, but can still be modified manually afterwards
        # TODO: replace everything with 'auto' modes.

        # Set default behavior for warnings:
        if isinstance(warnings, dict):
            self.warnings.update(warnings)
        elif warnings in [True, "warn", "warning"]:
            self.warnings["default"] = "warn"
        elif warnings == "error":
            self.warnings["default"] = "error"
        elif warnings in [False, "ignore"]:
            self.warnings["default"] = "ignore"
        else:
            raise ValueError("Unexpected value for warnings: {0}".format(warnings))
        # Set default values for warnings thresholds
        self.misc.warning_linestrength_cutoff = 1e-2
        self.misc.warning_broadening_threshold = 1e-2

        return

    # %% ======================================================================
    # PUBLIC METHODS
    # ------------------------
    # eq_spectrum             >>> calc equilibrium spectrum
    # non_eq_spectrum         >>> calc non equilibrium spectrum
    # eq_bands                >>> returns bands as a list of spectra
    # non_eq_bands            >>> same, with overpopulation factors
    # power                   >>> get total power (equilibrium or non eq)
    # fit_spectrum            >>> fit experimental spectrum
    #
    # XXX =====================================================================

    def eq_spectrum(
        self,
        Tgas,
        mole_fraction=None,
        path_length=None,
        diluent=None,
        pressure=None,
        name=None,
    ) -> Spectrum:
        """Generate a spectrum at equilibrium.

        Parameters
        ----------
        Tgas: float or `~astropy.units.quantity.Quantity`
            Gas temperature (K)
        mole_fraction: float
            database species mole fraction. If None, Factory mole fraction is used.
        diluent: str or dictionary
            can be a string of a single diluent or a dictionary containing diluent
            name as key and its mole_fraction as value
        path_length: float or `~astropy.units.quantity.Quantity`
            slab size (cm). If ``None``, the default Factory
            :py:attr:`~radis.lbl.factory.SpectrumFactor.input.path_length` is used.
        pressure: float or `~astropy.units.quantity.Quantity`
            pressure (bar). If ``None``, the default Factory
            :py:attr:`~radis.lbl.factory.SpectrumFactor.input.pressure` is used.
        name: str
            output Spectrum name (useful in batch)

        Returns
        -------
        s : Spectrum
            Returns a :class:`~radis.spectrum.spectrum.Spectrum` object

                Use the :meth:`~radis.spectrum.spectrum.Spectrum.get` method to get something
                among ``['radiance', 'radiance_noslit', 'absorbance', etc...]``

                Or directly the :meth:`~radis.spectrum.spectrum.Spectrum.plot` method
                to plot it. See [1]_ to get an overview of all Spectrum methods

        Examples
        --------
        ::
            from radis import SpectrumFactory
            sf = SpectrumFactory(
            wavenum_min=2900,
            wavenum_max=3200,
            molecule="OH",
            wstep=0.1,
            )
            sf.fetch_databank("hitemp")

            s1 = sf.eq_spectrum(Tgas=300, path_length=1, pressure=0.1)
            s2 = sf.eq_spectrum(Tgas=500, path_length=1, pressure=0.1)

        .. minigallery:: radis.lbl.SpectrumFactory.eq_spectrum
            :add-heading:

        References
        ----------
        .. [1] RADIS doc: `Spectrum how to? <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`__


        See Also
        --------
        :meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum`
        """
        # %% Preprocessing
        # --------------------------------------------------------------------

        # Check inputs
        if not self.input.self_absorption:
            raise ValueError(
                "Use non_eq_spectrum(Tgas, Tgas) to calculate spectra "
                + "without self_absorption"
            )

        # Convert units
        Tgas = convert_and_strip_units(Tgas, u.K)
        path_length = convert_and_strip_units(path_length, u.cm)
        pressure = convert_and_strip_units(pressure, u.bar)

        # update defaults
        if path_length is not None:
            self.input.path_length = path_length
        if mole_fraction is not None:
            self.input.mole_fraction = mole_fraction
        if pressure is not None:
            self.input.pressure_mbar = pressure * 1e3
        if not is_float(Tgas):
            raise ValueError(
                "Tgas should be float or Astropy unit. Got {0}".format(Tgas)
            )
        self.input.rot_distribution = "boltzmann"  # equilibrium
        self.input.vib_distribution = "boltzmann"  # equilibrium

        # Get temperatures
        self.input.Tgas = Tgas
        self.input.Tvib = Tgas  # just for info
        self.input.Trot = Tgas  # just for info

        # Init variables
        pressure_mbar = self.input.pressure_mbar
        mole_fraction = self.input.mole_fraction
        path_length = self.input.path_length
        verbose = self.verbose

        # New Profiler object
        self._reset_profiler(verbose)

        # Check variables
        self._check_inputs(mole_fraction, max(flatten(Tgas)))

        # Retrieve Spectrum from database if it exists
        if self.autoretrievedatabase:
            s = self._retrieve_from_database()
            if s is not None:
                return s  # exit function

        # %% Start
        # --------------------------------------------------------------------

        self.profiler.start("spectrum_calculation", 1)
        self.profiler.start("spectrum_calc_before_obj", 2)

        if verbose:
            self.print_conditions("Calculating Equilibrium Spectrum")

        # Check database, reset populations, create line dataframe to be scaled
        # --------------------------------------------------------------------
        self._check_line_databank()
        self._reinitialize()  # creates scaled dataframe df1 from df0

        # --------------------------------------------------------------------

        # First calculate the linestrength at given temperature
        self.calc_linestrength_eq(Tgas)  # scales S0 to S (equivalent to S0 in code)
        self._cutoff_linestrength()

        # ----------------------------------------------------------------------

        # Calculate line shift
        self.calc_lineshift()  # scales wav to shiftwav (equivalent to v0)

        # ----------------------------------------------------------------------
        # Line broadening

        # ... generates molefraction for diluents
        self._generate_diluent_molefraction(mole_fraction, diluent)

        # ... calculate broadening  HWHM
        self._calc_broadening_HWHM()

        # ... generates all wstep related entities
        self._generate_wavenumber_arrays()

        # ... find weak lines and calculate semi-continuum (optional)
        I_continuum = self.calculate_pseudo_continuum()
        # ... apply lineshape and get absorption coefficient
        # ... (this is the performance bottleneck)
        wavenumber, abscoeff_v = self._calc_broadening()
        #    :         :
        #   cm-1    1/(#.cm-2)

        # ... add semi-continuum (optional)
        abscoeff_v = self._add_pseudo_continuum(abscoeff_v, I_continuum)
        # Calculate output quantities
        # ----------------------------------------------------------------------

        self.profiler.start("calc_other_spectral_quan", 2)

        # incorporate density of molecules (see equation (A.16) )
        density = mole_fraction * ((pressure_mbar * 100) / (k_b * Tgas)) * 1e-6
        #  :
        # (#/cm3)

        abscoeff = abscoeff_v * density  # cm-1
        # ... # TODO: if the code is extended to multi-species, then density has to be added
        # ... before lineshape broadening (as it would not be constant for all species)

        # get absorbance (technically it's the optical depth `tau`,
        #                absorbance `A` being `A = tau/ln(10)` )
        absorbance = abscoeff * path_length
        # Generate output quantities
        transmittance_noslit = exp(-absorbance)
        emissivity_noslit = 1 - transmittance_noslit
        radiance_noslit = calc_radiance(
            wavenumber, emissivity_noslit, Tgas, unit=self.units["radiance_noslit"]
        )
        assert self.units["abscoeff"] == "cm-1"

        self.profiler.stop(
            "calc_other_spectral_quan", "Calculated other spectral quantities"
        )

        # %% Export
        # --------------------------------------------------------------------

        self.profiler.stop(
            "spectrum_calc_before_obj", "Spectrum calculated (before object generation)"
        )
        self.profiler.start("generate_spectrum_obj", 2)

        # Get conditions
        conditions = self.get_conditions(add_config=True)
        conditions.update(
            {
                "calculation_time": self.profiler.final[list(self.profiler.final)[-1]][
                    "spectrum_calc_before_obj"
                ],
                "lines_calculated": self._Nlines_calculated,
                "lines_cutoff": self._Nlines_cutoff,
                "lines_in_continuum": self._Nlines_in_continuum,
                "thermal_equilibrium": True,
                "diluents": self._diluent,
                "radis_version": version,
                "spectral_points": (
                    self.params.wavenum_max_calc - self.params.wavenum_min_calc
                )
                / self.params.wstep,
                "profiler": dict(self.profiler.final),
            }
        )
        if self.params.optimization != None:
            conditions.update(
                {
                    "NwL": self.NwL,
                    "NwG": self.NwG,
                }
            )
        del self.profiler.final[list(self.profiler.final)[-1]][
            "spectrum_calc_before_obj"
        ]

        # Get populations of levels as calculated in RovibrationalPartitionFunctions
        # ... Populations cannot be calculated at equilibrium (needs energies).
        # ... Use SpectrumFactory.non_eq_spectrum
        populations = None

        # Get lines (intensities + populations)
        lines = self.get_lines()

        # Spectral quantities
        quantities = {
            "wavenumber": wavenumber,
            "abscoeff": abscoeff,
            "absorbance": absorbance,
            "emissivity_noslit": emissivity_noslit,
            "transmittance_noslit": transmittance_noslit,
            "radiance_noslit": radiance_noslit,
        }
        if I_continuum is not None and self._export_continuum:
            quantities.update({"abscoeff_continuum": I_continuum * density})
        conditions["default_output_unit"] = self.input_wunit

        # Store results in Spectrum class
        s = Spectrum(
            quantities=quantities,
            conditions=conditions,
            populations=populations,
            lines=lines,
            units=self.units,
            cond_units=self.cond_units,
            # dont check input (much faster, and Spectrum
            # is freshly baken so probably in a good format
            check_wavespace=False,
            name=name,
            references=dict(self.reftracker),
        )
        # OPTION 2.  Change a posteriori using a Spectrum.method. More universal. Can it be slower?

        # update database if asked so
        if self.autoupdatedatabase:
            self.SpecDatabase.add(s, if_exists_then="increment")
            # Tvib=Trot=Tgas... but this way names in a database
            # generated with eq_spectrum are consistent with names
            # in one generated with non_eq_spectrum

        # Get generation & total calculation time
        self.profiler.stop("generate_spectrum_obj", "Generated Spectrum object")

        #  In the less verbose case, we print the total calculation+generation time:
        self.profiler.stop("spectrum_calculation", "Spectrum calculated")

        return s

    def eq_spectrum_gpu(
        self,
        Tgas,
        mole_fraction=None,
        diluent=None,
        path_length=None,
        pressure=None,
        name=None,
        emulate=False,
    ) -> Spectrum:
        """Generate a spectrum at equilibrium with calculation of lineshapes
        and broadening done on the GPU.

        .. note::
            This method requires CUDA compatible hardware to execute.
            For more information on how to setup your system to run GPU-accelerated methods
            using CUDA and Cython, check :ref:`GPU Spectrum Calculation on RADIS <label_radis_gpu>`

        Parameters
        ----------
        Tgas: float or `~astropy.units.quantity.Quantity`
            Gas temperature (K)
        mole_fraction: float
            database species mole fraction. If None, Factory mole fraction is used.
        path_length: float or `~astropy.units.quantity.Quantity`
            slab size (cm). If ``None``, the default Factory
            :py:attr:`~radis.lbl.factory.SpectrumFactor.input.path_length` is used.
        pressure: float or `~astropy.units.quantity.Quantity`
            pressure (bar). If ``None``, the default Factory
            :py:attr:`~radis.lbl.factory.SpectrumFactor.input.pressure` is used.
        name: str
            output Spectrum name (useful in batch)

        Other Parameters
        ----------------
        emulate: bool
            if ``True``, execute the GPU code on the CPU (useful for development)

        Returns
        -------
        s : Spectrum
            Returns a :class:`~radis.spectrum.spectrum.Spectrum` object

                Use the :meth:`~radis.spectrum.spectrum.Spectrum.get` method to get something
                among ``['radiance', 'radiance_noslit', 'absorbance', etc...]``

                Or directly the :meth:`~radis.spectrum.spectrum.Spectrum.plot` method
                to plot it. See [1]_ to get an overview of all Spectrum methods

        Examples
        --------

        .. minigallery:: radis.lbl.SpectrumFactory.eq_spectrum_gpu
            :add-heading:


        See Also
        --------
        :meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`,
        :meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu_interactive`
        """

        # %% Preprocessing
        # --------------------------------------------------------------------

        # Check inputs
        if not self.input.self_absorption:
            raise ValueError(
                "Use non_eq_spectrum(Tgas, Tgas) to calculate spectra "
                + "without self_absorption"
            )

        # Convert units
        Tgas = convert_and_strip_units(Tgas, u.K)
        path_length = convert_and_strip_units(path_length, u.cm)
        pressure = convert_and_strip_units(pressure, u.bar)

        # update defaults
        if path_length is not None:
            self.input.path_length = path_length
        if mole_fraction is not None:
            self.input.mole_fraction = mole_fraction
        if pressure is not None:
            self.input.pressure_mbar = pressure * 1e3
        if not is_float(Tgas):
            raise ValueError("Tgas should be float.")
        self.input.rot_distribution = "boltzmann"  # equilibrium
        self.input.vib_distribution = "boltzmann"  # equilibrium

        # Get temperatures
        self.input.Tgas = Tgas
        self.input.Tvib = Tgas  # just for info
        self.input.Trot = Tgas  # just for info

        verbose = self.verbose

        # New Profiler object
        self._reset_profiler(verbose)

        # Init variables
        pressure_mbar = self.input.pressure_mbar
        mole_fraction = self.input.mole_fraction
        path_length = self.input.path_length

        # Check variables
        self._check_inputs(mole_fraction, max(flatten(Tgas)))

        # Retrieve Spectrum from database if it exists
        if self.autoretrievedatabase:
            s = self._retrieve_from_database()
            if s is not None:
                return s  # exit function

        ### GET ISOTOPE ABUNDANCE & MOLECULAR MASS ###

        molpar = self.molparam

        try:
            if "id" in self.df0:
                id_set = self.df0[
                    "id"
                ].unique()  # get all the molecules in the dataframe, should ideally be 1 element for GPU
                mol_id = id_set[0]

                assert len(id_set) == 1  # make sure there is only one molecule

            else:
                mol_id = self.df0.attrs["id"]
        except:
            mol_id = get_molecule_identifier(self.input.molecule)

        molecule = get_molecule(mol_id)
        state = self.input.state
        iso_set = self._get_isotope_list(molecule)

        iso_list = list(range(max(iso_set) + 1))  # element 0 is not used

        molarmass_arr = np.empty_like(
            iso_list, dtype=np.float32
        )  # molar mass of each isotope

        Q_interp_list = []
        for iso in iso_list:
            if iso in iso_set:
                params = molpar.df.loc[(mol_id, iso)]
                molarmass_arr[iso] = params.molar_mass
                parsum = self.get_partition_function_interpolator(molecule, iso, state)
                Q_interp_list.append(parsum.at)
            else:
                Q_interp_list.append(lambda T: 1.0)

        molarmass_arr[np.isnan(molarmass_arr)] = 0

        self.profiler.start("spectrum_calculation", 1)
        self.profiler.start("spectrum_calc_before_obj", 2)

        # generate the v_arr
        v_arr = np.arange(
            self.input.wavenum_min,
            self.input.wavenum_max + self.params.wstep,
            self.params.wstep,
        )

        # load the data
        df = self.df0
        v0 = df["wav"].to_numpy(dtype=np.float32)

        if len(iso_set) > 1:
            iso = df["iso"].to_numpy(dtype=np.uint8)
        elif len(iso_set) == 1:
            iso = np.full(len(v0), iso_set[0], dtype=np.uint8)

        da = df["Pshft"].to_numpy(dtype=np.float32)
        El = df["El"].to_numpy(dtype=np.float32)
        na = df["Tdpair"].to_numpy(dtype=np.float32)

        gamma = np.array(
            self._get_lorentzian_broadening(mole_fraction), dtype=np.float32
        )

        self.calc_S0()
        S0 = self.df0["S0"].to_numpy(dtype=np.float32)

        dxG = self.params.dxG
        dxL = self.params.dxL

        _Nlines_calculated = len(v0)

        if verbose >= 2:
            print("Initializing parameters...", end=" ")

        try:
            from radis.lbl.gpu import gpu_init, gpu_iterate
        except (ModuleNotFoundError):
            print("Failed to load GPU module, exiting!")
            exit()

        gpu_init(
            v_arr,
            dxG,
            dxL,
            iso,
            v0,
            da,
            gamma,
            na,
            S0,
            El,
            molarmass_arr,
            Q_interp_list,
            verbose=verbose,
            gpu=(not emulate),
        )

        if verbose >= 2:
            print("Initialization complete!")

        wavenumber = v_arr

        if verbose >= 2:
            print("Calculating spectra...", end=" ")

        abscoeff, transmittance, iter_params = gpu_iterate(
            pressure_mbar * 1e-3,
            Tgas,
            mole_fraction,
            verbose=verbose,
            gpu=(not emulate),
        )
        # Calculate output quantities
        # ----------------------------------------------------------------------

        self.profiler.start("calc_other_spectral_quan", 2)

        # ... # TODO: if the code is extended to multi-species, then density has to be added
        # ... before lineshape broadening (as it would not be constant for all species)

        # get absorbance (technically it's the optical depth `tau`,
        #                absorbance `A` being `A = tau/ln(10)` )
        absorbance = abscoeff * path_length
        # Generate output quantities
        transmittance_noslit = exp(-absorbance)
        emissivity = 1 - transmittance
        emissivity_noslit = 1 - transmittance_noslit
        radiance_noslit = calc_radiance(
            wavenumber, emissivity_noslit, Tgas, unit=self.units["radiance_noslit"]
        )
        assert self.units["abscoeff"] == "cm-1"

        self.profiler.stop(
            "calc_other_spectral_quan", "Calculated other spectral quantities"
        )

        lines = self.get_lines()

        # %% Export
        # --------------------------------------------------------------------

        self.profiler.stop(
            "spectrum_calc_before_obj", "Spectrum calculated (before object generation)"
        )
        self.profiler.start("generate_spectrum_obj", 2)

        # Get lines (intensities + populations)

        conditions = self.get_conditions(add_config=True)
        conditions.update(
            {
                "calculation_time": self.profiler.final[list(self.profiler.final)[-1]][
                    "spectrum_calc_before_obj"
                ],
                "lines_calculated": _Nlines_calculated,
                "thermal_equilibrium": True,
                "diluents": self._diluent,
                "radis_version": version,
                "emulate_gpu": emulate,
                "spectral_points": (
                    self.params.wavenum_max_calc - self.params.wavenum_min_calc
                )
                / self.params.wstep,
                "profiler": dict(self.profiler.final),
            }
        )
        if self.params.optimization != None:
            conditions.update(
                {
                    "NwL": iter_params.N_L,
                    "NwG": iter_params.N_G,
                }
            )
        del self.profiler.final[list(self.profiler.final)[-1]][
            "spectrum_calc_before_obj"
        ]

        # Spectral quantities
        quantities = {
            "wavenumber": wavenumber,
            "abscoeff": abscoeff,
            "absorbance": absorbance,
            "emissivity": emissivity,
            "emissivity_noslit": emissivity_noslit,
            "transmittance_noslit": transmittance_noslit,
            "radiance_noslit": radiance_noslit,
            "transmittance": transmittance,
        }
        conditions["default_output_unit"] = self.input_wunit

        # Store results in Spectrum class
        s = Spectrum(
            quantities=quantities,
            units=self.units,
            conditions=conditions,
            lines=lines,
            cond_units=self.cond_units,
            # dont check input (much faster, and Spectrum
            check_wavespace=False,
            # is freshly baken so probably in a good format
            name=name,
            references=dict(self.reftracker),
        )

        # update database if asked so
        if self.autoupdatedatabase:
            self.SpecDatabase.add(s, if_exists_then="increment")
            # Tvib=Trot=Tgas... but this way names in a database
            # generated with eq_spectrum are consistent with names
            # in one generated with non_eq_spectrum

        # Get generation & total calculation time
        self.profiler.stop("generate_spectrum_obj", "Generated Spectrum object")

        #  In the less verbose case, we print the total calculation+generation time:
        self.profiler.stop("spectrum_calculation", "Spectrum calculated")

        return s

    def eq_spectrum_gpu_interactive(
        self, var="transmittance", slit_FWHM=0.0, plotkwargs={}, *vargs, **kwargs
    ) -> Spectrum:
        """

        Parameters
        ----------
        var : TYPE, optional
            DESCRIPTION. The default is "transmittance".
        slit_FWHM : TYPE, optional
            DESCRIPTION. The default is 0.0.
        *vargs : TYPE
            arguments forwarded to :py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu`
        **kwargs : dict
            arguments forwarded to :py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu`
            In particular, see `emulate=`.
        plotkwargs : dict
            arguments forwarded to :py:meth:`~radis.spectrum.spectrum.Spectrum.plot`

        Returns
        -------
        s : TYPE
            DESCRIPTION.

        Examples
        --------
        ::

            from radis import SpectrumFactory
            from radis.tools.plot_tools import ParamRange


            sf = SpectrumFactory(2200, 2400, # cm-1
                              molecule='CO2',
                              isotope='1,2,3',
                              wstep=0.002,
                              )

            sf.fetch_databank('hitemp')

            s = sf.eq_spectrum_gpu_interactive(Tgas=ParamRange(300.0,2000.0,1200.0), #K
                                           pressure=0.2, #bar
                                           mole_fraction=0.1,
                                           path_length=ParamRange(0,10,2.0), #cm
                                           slit_FWHM=ParamRange(0,1,0), #cm
                                           emulate=False,  # if True, execute the GPU code on the CPU
                                           )

        .. minigallery:: radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu_interactive
            :add-heading:

        """

        import matplotlib.pyplot as plt
        from matplotlib.widgets import Slider

        try:
            from radis.lbl.gpu import gpu_iterate
        except (ModuleNotFoundError):
            print("Failed to load GPU module, exiting!")
            exit()

        self.interactive_params = {}

        for key in kwargs:
            if is_range(kwargs[key]):
                param_range = kwargs[key]
                self.interactive_params[key] = param_range
                param_range.name = key
                kwargs[key] = param_range.valinit

        s = self.eq_spectrum_gpu(*vargs, **kwargs)

        if is_range(slit_FWHM):
            self.interactive_params["slit_FWHM"] = slit_FWHM
            slit_FWHM.name = "slit_FWHM"
            slit_FWHM = slit_FWHM.valinit

        s.apply_slit(slit_FWHM, unit="cm-1")  # to create 'radiance', 'transmittance'
        s.conditions["slit_FWHM"] = slit_FWHM

        # TODO: this should be resolved; there should be only one pressure that is
        # both used as keyword by the user and as input value for spectra.
        s.conditions["pressure"] = s.conditions["pressure_mbar"] * 1e-3

        was_interactive = plt.isinteractive
        plt.ion()

        line = s.plot(var, show=True, **plotkwargs)
        fig = line.figure

        def update_plot(val):
            # TODO : refactor this function and the update() mechanism. Ensure conditions are correct.
            # Update conditions
            # ... at equilibirum, temperatures remain equal :
            s.conditions["Tvib"] = s.conditions["Tgas"]
            s.conditions["Trot"] = s.conditions["Tgas"]
            s.conditions["slit_function"] = s.conditions["slit_FWHM"]

            abscoeff, transmittance, iter_params = gpu_iterate(
                s.conditions["pressure"],
                s.conditions["Tgas"],
                s.conditions["mole_fraction"],
                l=s.conditions["path_length"],
                slit_FWHM=s.conditions["slit_FWHM"],
                verbose=False,
                gpu=(not s.conditions["emulate_gpu"]),
            )

            # This happen inside a Spectrum() method
            for k in list(s._q.keys()):  # reset all quantities
                if k in ["wavespace", "wavelength", "wavenumber"]:
                    pass
                elif k == "abscoeff":
                    s._q["abscoeff"] = abscoeff
                else:
                    del s._q[k]

            _, new_y = s.get(
                var,
                copy=False,  # copy = False saves some time & memory, it's a pointer/reference to the real data, which is fine here as data is just plotted
                wunit=plotkwargs.get("wunit", "default"),
                Iunit=plotkwargs.get("Iunit", "default"),
            )

            line.set_ydata(new_y)
            fig.canvas.draw_idle()

        n_sliders = 0
        for key in self.interactive_params:
            param_range = self.interactive_params[key]
            slider_axis = fig.add_axes([0.25, 0.05 * n_sliders + 0.05, 0.65, 0.03])
            slider = Slider(
                ax=slider_axis,
                label=key,
                valmin=param_range.valmin,
                valmax=param_range.valmax,
                valinit=param_range.valinit,
            )
            self.interactive_params[key].set_widget(slider, s, update_plot)
            n_sliders += 1

        fig.subplots_adjust(bottom=0.05 * n_sliders + 0.15)

        if not was_interactive:
            plt.ioff()

        return s

    def non_eq_spectrum(
        self,
        Tvib,
        Trot,
        Ttrans=None,
        mole_fraction=None,
        path_length=None,
        diluent=None,
        pressure=None,
        vib_distribution="boltzmann",
        rot_distribution="boltzmann",
        overpopulation=None,
        name=None,
    ) -> Spectrum:
        """Calculate emission spectrum in non-equilibrium case. Calculates
        absorption with broadened linestrength and emission with broadened
        Einstein coefficient.

        Parameters
        ----------
        Tvib: float
            vibrational temperature [K]
            can be a tuple of float for the special case of more-than-diatomic
            molecules (e.g: CO2)
        Trot: float
            rotational temperature [K]
        Ttrans: float
            translational temperature [K]. If None, translational temperature is
            taken as rotational temperature (valid at 1 atm for times above ~ 2ns
            which is the RT characteristic time)
        mole_fraction: float
            database species mole fraction. If None, Factory mole fraction is used.
        diluent: str or dictionary
            can be a string of a single diluent or a dictionary containing diluent
            name as key and its mole_fraction as value
        path_length: float or `~astropy.units.quantity.Quantity`
            slab size (cm). If ``None``, the default Factory
            :py:attr:`~radis.lbl.factory.SpectrumFactor.input.path_length` is used.
        pressure: float or `~astropy.units.quantity.Quantity`
            pressure (bar). If ``None``, the default Factory
            :py:attr:`~radis.lbl.factory.SpectrumFactor.input.pressure` is used.

        Other Parameters
        ----------------
        vib_distribution: ``'boltzmann'``, ``'treanor'``
            vibrational distribution
        rot_distribution: ``'boltzmann'``
            rotational distribution
        overpopulation: dict, or ``None``
            add overpopulation factors for given levels::

                {level:overpopulation_factor}
        name: str
            output Spectrum name (useful in batch)

        Returns
        -------
        s : Spectrum
            Returns a :class:`~radis.spectrum.spectrum.Spectrum` object

                Use the :meth:`~radis.spectrum.spectrum.Spectrum.get` method to get something
                among ``['radiance', 'radiance_noslit', 'absorbance', etc...]``

                Or directly the :meth:`~radis.spectrum.spectrum.Spectrum.plot` method
                to plot it. See [1]_ to get an overview of all Spectrum methods

        Examples
        --------
        ::

            from radis import SpectrumFactory
            sf = SpectrumFactory(
            wavenum_min=2000,
            wavenum_max=3000,
            molecule="CO",
            isotope="1,2,3",
            wstep="auto"
            )
            sf.fetch_databank("hitemp", load_columns='noneq')

            s1 = sf.non_eq_spectrum(Tvib=2000, Trot=600, path_length=1, pressure=0.1)
            s2 = sf.non_eq_spectrum(Tvib=2000, Trot=600, path_length=1, pressure=0.1)

        Multi-vibrational temperature. Below we compare non-LTE spectra of CO2 where all
        vibrational temperatures are equal, or where the bending & symmetric modes are in
        equilibrium with rotation ::

            from radis import SpectrumFactory
            sf = SpectrumFactory(
            wavenum_min=2000,
            wavenum_max=3000,
            molecule="CO2",
            isotope="1,2,3",
            )
            sf.fetch_databank("hitemp", load_columns='noneq')
            # nonequilibrium between bending+symmetric and asymmetric modes :
            s1 = sf.non_eq_spectrum(Tvib=(600, 600, 2000), Trot=600, path_length=1, pressure=1)
            # all vibrational temperatures are equal :
            s2 = sf.non_eq_spectrum(Tvib=(2000, 2000, 2000), Trot=600, path_length=1, pressure=1)



        .. minigallery:: radis.lbl.SpectrumFactory.non_eq_spectrum
            :add-heading:

        References
        ----------
        .. [1] RADIS doc: `Spectrum how to? <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`__

        See Also
        --------

        :meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`
        :meth:`~radis.lbl.factory.SpectrumFactory.optically_thin_power`
        """

        # %% Preprocessing
        # --------------------------------------------------------------------
        # Convert units
        Tvib = convert_and_strip_units(Tvib, u.K)
        Trot = convert_and_strip_units(Trot, u.K)
        Ttrans = convert_and_strip_units(Ttrans, u.K)
        path_length = convert_and_strip_units(path_length, u.cm)
        pressure = convert_and_strip_units(pressure, u.bar)

        # check inputs, update defaults
        if path_length is not None:
            self.input.path_length = path_length
        if mole_fraction is not None:
            self.input.mole_fraction = mole_fraction
        if pressure is not None:
            self.input.pressure_mbar = pressure * 1e3
        if isinstance(Tvib, tuple):
            Tvib = tuple([convert_and_strip_units(T, u.K) for T in Tvib])
        elif not is_float(Tvib):
            raise TypeError(
                "Tvib should be float, or tuple (got {0})".format(type(Tvib))
            )
        singleTvibmode = is_float(Tvib)
        if not is_float(Trot):
            raise ValueError("Trot should be float")
        if overpopulation is None:
            overpopulation = {}
        assert vib_distribution in ["boltzmann", "treanor"]
        assert rot_distribution in ["boltzmann"]
        self.input.overpopulation = overpopulation
        self.input.rot_distribution = rot_distribution
        self.input.vib_distribution = vib_distribution
        # Get translational temperature
        Tgas = Ttrans
        if Tgas is None:
            Tgas = Trot  # assuming Ttrans = Trot
        self.input.Tgas = Tgas
        self.input.Tvib = Tvib
        self.input.Trot = Trot

        # Init variables
        path_length = self.input.path_length
        mole_fraction = self.input.mole_fraction
        pressure_mbar = self.input.pressure_mbar
        verbose = self.verbose

        # New Profiler object
        self._reset_profiler(verbose)
        # Check variables
        self._check_inputs(mole_fraction, max(flatten(Tgas, Tvib, Trot)))

        # Retrieve Spectrum from database if it exists
        if self.autoretrievedatabase:
            s = self._retrieve_from_database()
            if s is not None:
                return s  # exit function

        # %% Start
        # --------------------------------------------------------------------

        self.profiler.start("spectrum_calculation", 1)
        self.profiler.start("spectrum_calc_before_obj", 2)
        if verbose:
            self.print_conditions("Calculating Non-Equilibrium Spectrum")

        # Check line database and parameters, reset populations and scaled line dataframe
        # ----------
        self._check_line_databank()

        # add nonequilibrium energies if needed (this may be a bottleneck
        # for a first calculation):
        self._calc_noneq_parameters(vib_distribution, singleTvibmode)
        self._reinitialize()  # creates scaled dataframe df1 from df0

        # ----------------------------------------------------------------------
        # Calculate Populations, Linestrength and Emission Integral
        if singleTvibmode:
            self.calc_populations_noneq(
                Tvib,
                Trot,
                vib_distribution=vib_distribution,
                rot_distribution=rot_distribution,
                overpopulation=overpopulation,
            )
        else:
            self._calc_populations_noneq_multiTvib(
                Tvib,
                Trot,
                vib_distribution=vib_distribution,
                rot_distribution=rot_distribution,
                overpopulation=overpopulation,
            )

        self.calc_linestrength_noneq()
        self.calc_emission_integral()

        # ----------------------------------------------------------------------
        # Cutoff linestrength
        self._cutoff_linestrength()

        # ----------------------------------------------------------------------

        # Calculate lineshift
        self.calc_lineshift()

        # ----------------------------------------------------------------------

        # Line broadening

        # ... generates molefraction for diluents
        self._generate_diluent_molefraction(mole_fraction, diluent)

        # ... calculate broadening  HWHM
        self._calc_broadening_HWHM()

        # ... generates all wstep related entities
        self._generate_wavenumber_arrays()

        # ... find weak lines and calculate semi-continuum (optional)
        k_continuum, j_continuum = self.calculate_pseudo_continuum(noneq=True)

        # ... apply lineshape and get absorption coefficient
        # ... (this is the performance bottleneck)
        wavenumber, abscoeff_v, emisscoeff_v = self._calc_broadening_noneq()
        #    :         :            :
        #   cm-1    1/(#.cm-2)   mW/sr/cm-1

        # ... add semi-continuum (optional)
        abscoeff_v = self._add_pseudo_continuum(abscoeff_v, k_continuum)
        emisscoeff_v = self._add_pseudo_continuum(emisscoeff_v, j_continuum)

        # Calculate output quantities
        # ----------------------------------------------------------------------

        self.profiler.start("calc_other_spectral_quan", 2)

        # incorporate density of molecules (see Rothman 1996 equation (A.16) )
        density = mole_fraction * ((pressure_mbar * 100) / (k_b * Tgas)) * 1e-6
        #  :
        # (#/cm3)

        abscoeff = abscoeff_v * density  # cm-1
        emisscoeff = emisscoeff_v * density  # mW/sr/cm3/cm-1

        # get absorbance (technically it's the optical depth `tau`,
        #                absorbance `A` being `A = tau/ln(10)` )

        # Generate output quantities
        absorbance = abscoeff * path_length  # (adim)
        transmittance_noslit = exp(-absorbance)

        if self.input.self_absorption:
            # Analytical output of computing RTE over a single slab of constant
            # emissivity and absorption coefficient
            b = abscoeff == 0  # optically thin mask
            radiance_noslit = np.zeros_like(emisscoeff)
            radiance_noslit[~b] = (
                emisscoeff[~b] / abscoeff[~b] * (1 - transmittance_noslit[~b])
            )
            radiance_noslit[b] = emisscoeff[b] * path_length
        else:
            # Note that for k -> 0,
            radiance_noslit = emisscoeff * path_length  # (mW/cm2/sr/cm-1)

        # Convert `radiance_noslit` from (mW/sr/cm2/cm-1) to output unit
        radiance_noslit = convert_universal(
            radiance_noslit,
            from_unit="mW/cm2/sr/cm-1",
            to_unit=self.units["radiance_noslit"],
            wavenum=wavenumber,
            per_nm_is_like="mW/cm2/sr/nm",
            per_cm_is_like="mW/cm2/sr/cm-1",
        )
        # Convert 'emisscoeff' from (mW/sr/cm3/cm-1) to output unit
        emisscoeff = convert_universal(
            emisscoeff,
            from_unit="mW/cm3/sr/cm-1",
            to_unit=self.units["emisscoeff"],
            wavenum=wavenumber,
            per_nm_is_like="mW/cm3/sr/nm",
            per_cm_is_like="mW/cm3/sr/cm-1",
        )
        assert self.units["abscoeff"] == "cm-1"

        self.profiler.stop(
            "calc_other_spectral_quan", "Calculated other spectral quantities"
        )

        # Note: emissivity not defined under non equilibrium

        # %% Export
        # ----------------------------------------------------------------------

        self.profiler.stop(
            "spectrum_calc_before_obj", "Spectrum calculated (before object generation)"
        )
        self.profiler.start("generate_spectrum_obj", 2)

        # Get conditions
        conditions = self.get_conditions(add_config=True)
        conditions.update(
            {
                "calculation_time": self.profiler.final[list(self.profiler.final)[-1]][
                    "spectrum_calc_before_obj"
                ],
                "lines_calculated": self._Nlines_calculated,
                "lines_cutoff": self._Nlines_cutoff,
                "lines_in_continuum": self._Nlines_in_continuum,
                "thermal_equilibrium": False,  # dont even try to guess if it's at equilibrium
                "diluents": self._diluent,
                "radis_version": version,
                "spectral_points": (
                    self.params.wavenum_max_calc - self.params.wavenum_min_calc
                )
                / self.params.wstep,
                "profiler": dict(self.profiler.final),
            }
        )
        if self.params.optimization != None:
            conditions.update(
                {
                    "NwL": self.NwL,
                    "NwG": self.NwG,
                }
            )
        del self.profiler.final[list(self.profiler.final)[-1]][
            "spectrum_calc_before_obj"
        ]

        # Get populations of levels as calculated in RovibrationalPartitionFunctions
        populations = self.get_populations(self.misc.export_populations)

        # Get lines (intensities + populations)
        lines = self.get_lines()

        # Spectral quantities
        quantities = {
            "wavenumber": wavenumber,
            "abscoeff": abscoeff,
            "absorbance": absorbance,
            "emisscoeff": emisscoeff,
            "transmittance_noslit": transmittance_noslit,
            "radiance_noslit": radiance_noslit,
        }
        if k_continuum is not None and self._export_continuum:
            quantities.update(
                {
                    "abscoeff_continuum": k_continuum * density,
                    "emisscoeff_continuum": j_continuum * density,
                }
            )
        conditions["default_output_unit"] = self.input_wunit

        # Store results in Spectrum class
        s = Spectrum(
            quantities=quantities,
            conditions=conditions,
            populations=populations,
            lines=lines,
            units=self.units,
            cond_units=self.cond_units,
            # dont check input (much faster, and Spectrum
            check_wavespace=False,
            # is freshly baken so probably in a good format
            name=name,
            references=dict(self.reftracker),
        )

        # update database if asked so
        if self.autoupdatedatabase:
            self.SpecDatabase.add(
                s,
                add_info=["Tvib", "Trot"],
                add_date="%Y%m%d",
                if_exists_then="increment",
            )

        # Get generation & total calculation time
        self.profiler.stop("generate_spectrum_obj", "Generated Spectrum object")

        #  In the less verbose case, we print the total calculation+generation time:
        self.profiler.stop("spectrum_calculation", "Spectrum calculated")

        return s

    def _generate_wavenumber_arrays(self):
        """define wavenumber grid vectors

        `SpectrumFactory.wavenumber` is the output spectral range and
        ``SpectrumFactory.wavenumber_calc`` the spectral range used for calculation, that
        includes neighbour lines within ``neighbour_lines`` distance."""

        import radis

        self.profiler.start("generate_wavenumber_arrays", 2)
        # calculates minimum FWHM of lines
        self._calc_min_width(self.df1)

        # Setting wstep to optimal value and rounding it to a degree 3
        if self._wstep == "auto" or type(self.params.wstep) == list:

            wstep_calc = round_off(
                self.min_width / radis.config["GRIDPOINTS_PER_LINEWIDTH_WARN_THRESHOLD"]
            )

            if type(self.params.wstep) == list:
                self.params.wstep = min(self.params.wstep[1], wstep_calc)
            else:
                self.params.wstep = wstep_calc

            self.warnings["AccuracyWarning"] = "ignore"

        truncation = self.params.truncation
        neighbour_lines = self.params.neighbour_lines

        if truncation and neighbour_lines > truncation:
            self.warn(
                f"Neighbour lines resolved up to {neighbour_lines} cm-1 away from the spectrum. "
                + f"But lines are anyway truncated at {truncation:.2f} cm-1. "
                + f"Choose 'neighbour_lines={truncation:.2f}' to avoid resolving useless lines",
                "PerformanceWarning",
            )

        wavenumber, wavenumber_calc, woutrange = _generate_wavenumber_range(
            self.input.wavenum_min,
            self.input.wavenum_max,
            self.params.wstep,
            neighbour_lines,
        )

        # Generate lineshape array
        if truncation is None:
            # compute lineshape on full range :
            # (note that this means 3x wavenumber_calc will be required when applying lineshapes)
            truncation = wavenumber_calc[-1] - wavenumber_calc[0]

        wbroad_centered = _generate_broadening_range(self.params.wstep, truncation)
        self.truncation = truncation
        # store value for use in lineshape broadening.
        # Note : may be different from self.params.truncation if None was given.

        self.wbroad_centered = wbroad_centered
        self.wavenumber = wavenumber
        self.wavenumber_calc = wavenumber_calc
        self.woutrange = woutrange

        # AccuracyWarning. Check there are enough gridpoints per line.
        self._check_accuracy(self.params.wstep)

        # Set sparse waverange mode

        # Setting wstep to optimal value and rounding it to a degree 3
        if self._sparse_ldm == "auto":
            sparsity = len(wavenumber_calc) / len(self.df1)
            self.params["sparse_ldm"] = (
                sparsity > 1.0
            )  # works ; TODO : set a threshold based on more data
            if self.verbose >= 2:
                print(
                    f"Sparsity (grid points/lines) = {sparsity:.1f}. Set sparse_ldm to {self.params['sparse_ldm']}"
                )

        self.profiler.stop("generate_wavenumber_arrays", "Generated Wavenumber Arrays")

        if radis.config["DEBUG_MODE"]:
            assert (wavenumber_calc[woutrange[0] : woutrange[1]] == wavenumber).all()

        return

    def _generate_diluent_molefraction(self, mole_fraction, diluent):
        from radis.misc.warning import MoleFractionError

        molecule = self.input.molecule
        # Check if diluent is same as molecule or not
        if (isinstance(diluent, str) and diluent == molecule) or (
            isinstance(diluent, dict) and molecule in diluent.keys()
        ):
            raise KeyError(
                "{0} is being called as molecule and diluent, please remove it from diluent.".format(
                    molecule
                )
            )

        # Using Spectrum Factory values of diluents
        if diluent is None:
            # If diluent is air, then add remaining molefraction as 'air'
            if isinstance(self.params.diluent, str):
                diluents = {self.params.diluent: 1 - mole_fraction}
            else:
                diluents = self.params.diluent.copy()
        else:
            if isinstance(diluent, str):
                diluents = {diluent: 1 - mole_fraction}
            else:
                diluents = diluent.copy()

        # Checking mole_fraction of molecule and diluent
        total_mole_fraction = mole_fraction + sum(list(diluents.values()))
        if total_mole_fraction > 1:
            raise MoleFractionError(
                "Total molefraction = {0} of molecule and diluents greater than 1. Please set appropriate molefraction value of molecule and diluents.".format(
                    total_mole_fraction
                )
            )
        elif total_mole_fraction < 1:
            raise MoleFractionError(
                "Total molefraction = {0} of molecule and diluents less than 1. Please set appropriate molefraction value of molecule and diluents".format(
                    total_mole_fraction
                )
            )

        self._diluent = diluents

    def predict_time(self):
        """predict_time(self) uses the input parameters like Spectral Range, Number of lines, wstep,
        truncation to predict the estimated calculation time for the Spectrum
        broadening step(bottleneck step) for the current optimization and broadening_method. The formula
        for predicting time is based on benchmarks performed on various parameters for different optimization,
        broadening_method and deriving its time complexity.

        Benchmarks: https://anandxkumar.github.io/Benchmark_Visualization_GSoC_2021/

        Complexity vs Calculation Time Visualizations
        LBL>Voigt: `LINK <https://public.tableau.com/app/profile/anand.kumar4841/viz/LegacyComplexityvsCalculationTime/Sheet1>`_,
        DIT>Voigt: `LINK <https://public.tableau.com/app/profile/anand.kumar4841/viz/2_096e-07lines_calculated7_185e-091wLwGSpectral_PointslogSpectral_Points/Sheet1>`_,
        DIT>FFT: `LINK <https://public.tableau.com/app/profile/anand.kumar4841/viz/LDMLatestLDMFFTComplexity4_675e-081wLwGSpectralPointslogSpectralPoints/Sheet1>`_

        Returns
        -------
        float: Predicted time in seconds.
        """

        def _is_at_equilibrium():
            try:
                assert self.input.Tvib is None or self.input.Tvib == self.input.Tgas
                assert self.input.Trot is None or self.input.Trot == self.input.Tgas
                assert (
                    self.input.overpopulation is None or self.input.overpopulation == {}
                )
                try:
                    if self.input.self_absorption:
                        assert self.input.self_absorption  # == True
                except KeyError:
                    pass
                return True
            except AssertionError:
                return False

        factor = 1

        if not _is_at_equilibrium():
            factor = 2  #  _apply_broadening_LDM() is called twice

        wstep = self.params.wstep
        n_lines = self.misc.total_lines
        truncation = self.params.truncation
        spectral_points = (
            self.params.wavenum_max_calc - self.params.wavenum_min_calc
        ) / self.params.wstep

        optimization = self.params.optimization
        broadening_method = self.params.broadening_method

        if optimization in ("simple", "min-RMS"):
            NwL = self.NwL
            NwG = self.NwG
            if broadening_method == "voigt":
                estimated_time = (
                    2.096e-07 * n_lines
                    + 7.185e-09
                    * (1 + NwL * NwG)
                    * spectral_points
                    * np.log(spectral_points)
                    * factor
                )
            elif broadening_method == "fft":
                estimated_time = (
                    4.675e-08
                    * (1 + NwL * NwG)
                    * spectral_points
                    * np.log(spectral_points)
                    * factor
                )
            elif broadening_method == "convolve":  # Not benchmarked
                estimated_time = (
                    self._broadening_time_ruleofthumb
                    * len(self.df0)
                    * len(self.wbroad_centered)
                )
            else:
                raise NotImplementedError("broadening_method not implemented")
        elif optimization is None:
            if broadening_method == "voigt":
                estimated_time = 6.6487e-08 * n_lines * truncation / wstep
            elif broadening_method == "convolve":  # Not benchmarked
                estimated_time = (
                    self._broadening_time_ruleofthumb
                    * len(self.df0)
                    * len(self.wbroad_centered)
                )
            else:
                raise NotImplementedError("broadening_method not implemented")
        else:
            raise NotImplementedError("optimization not implemented")

        return estimated_time

    def _get_lorentzian_broadening(self, x):
        """Calculates ratioed self/air broadening based on molefraction"""
        df = self.df0
        # TODO: deal with the case of gamma_self [so we don't forget]
        # TODO (refactor) : move into BaseFactory or BroadenFactory (parent classes)

        gamma_air = df["airbrd"].to_numpy()
        gamma_self = df["selbrd"].to_numpy()
        gamma = x * gamma_self + (1 - x) * gamma_air
        return gamma

    ##    def _get_S0(self, Ia_arr):
    ##        """Returns S0 if it already exists, otherwise computes the value using
    ##        abundance, upper level degeneracy and Einstein's number."""
    ##        df = self.df0
    ##
    ##        # if the column already exists, then return it
    ##        if "S0" in df.columns:
    ##            return df["S0"]
    ##
    ##        ## TO-DO: I don't think 'int' and 'S0' are the same quantity!
    ##        ##elif "int" in df.columns:
    ##        ##    return df["int"]
    ##
    ##        try:
    ##            v0 = df["wav"].to_numpy()
    ##            iso = df["iso"].to_numpy()
    ##            A21 = df["A"].to_numpy()
    ##            Jl = df["jl"].to_numpy()
    ##            DJ = df["branch"].to_numpy()
    ##            Ju = Jl + DJ
    ##            gu = 2 * Ju + 1  # g_up
    ##            S0 = Ia_arr.take(iso) * gu * A21 / (8 * pi * c_cm * v0 ** 2)
    ##            df["S0"] = S0
    ##            return S0
    ##
    ##        except KeyError as err:
    ##            raise KeyError(
    ##                "Could not find wavenumber, Einstein's coefficient, lower state energy or S0 in the dataframe. PLease check the database"
    ##            ) from err

    def optically_thin_power(
        self,
        Tgas=None,
        Tvib=None,
        Trot=None,
        Ttrans=None,
        vib_distribution="boltzmann",
        rot_distribution="boltzmann",
        mole_fraction=None,
        path_length=None,
        unit="mW/cm2/sr",
    ):
        """Calculate total power emitted in equilibrium or non-equilibrium case
        in the optically thin approximation: it sums all emission integral over
        the total spectral range.

        .. warning::

            this is a fast implementation that doesnt take into account
            the contribution of lines outside the given spectral range. It is valid for spectral ranges
            surrounded by no lines, and spectral ranges much broaded than the typical
            line broadening (~ 1-10 cm-1 in the infrared)

        If what you're looking for is an accurate simulation on a narrow spectral range
        you better calculate the spectrum (that does take all of that into account)
        and integrate it with :py:meth:`~radis.spectrum.spectrum.Spectrum.get_power`

        Parameters
        ----------
        Tgas: float
            equilibrium temperature [K]
            If doing a non equilibrium case it should be None. Use Ttrans for
            translational temperature
        Tvib: float
            vibrational temperature [K]
        Trot: float
            rotational temperature [K]
        Ttrans: float
            translational temperature [K]. If None, translational temperature is
            taken as rotational temperature (valid at 1 atm for times above ~ 2ns
            which is the RT characteristic time)
        mole_fraction: float
            database species mole fraction. If None, Factory mole fraction is used.
        path_length: float
            slab size (cm). If None, Factory mole fraction is used.
        unit: str
            output unit. Default ``'mW/cm2/sr'``

        Returns
        -------
        float: Returns total power density in mW/cm2/sr (unless different unit is chosen),
            see ``unit=``.


        See Also
        --------
        :py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`,
        :py:meth:`~radis.spectrum.spectrum.Spectrum.get_power`,
        :py:meth:`~radis.spectrum.spectrum.Spectrum.get_integral`
        """
        # Check inputs

        # ... temperatures

        if Tgas is None and Trot is None:
            raise ValueError(
                "Choose either Tgas (equilibrium) or Tvib / Trot "
                + ". Ttrans (non equilibrium)"
            )

        non_eq_mode = Tgas is None

        if Tvib is None and Trot is not None or Tvib is not None and Trot is None:
            raise ValueError("Choose both Tvib and Trot")

        if non_eq_mode and Ttrans is None:
            Ttrans = Trot

        # update defaults
        if path_length is not None:
            self.input.path_length = path_length
        if mole_fraction is not None:
            self.input.mole_fraction = mole_fraction

        # Get translational temperature
        if Tgas is None:
            Tgas = Ttrans
        self.input.Tgas = Tgas
        self.input.Tvib = Tvib
        self.input.Trot = Trot

        # Init variables
        path_length = self.input.path_length
        mole_fraction = self.input.mole_fraction
        pressure_mbar = self.input.pressure_mbar
        verbose = self.verbose

        # New Profiler object
        self._reset_profiler(verbose)

        self.profiler.start("optically_thin_power_calculation", 1)

        # Make sure database is loaded
        if self.df0 is None:
            if not self.save_memory:
                raise AttributeError("Load databank first (.load_databank())")
            else:
                self._reload_databank()

        if non_eq_mode:
            singleTvibmode = is_float(Tvib)
            # Make sure database has pre-computed non equilibrium quantities
            # (Evib, Erot, etc.)
            self._calc_noneq_parameters(vib_distribution, singleTvibmode)

        # %% Start
        # ----------------------------------------------------------------------

        # Print conditions
        if verbose:
            self.print_conditions(
                "Calculating Radiative Power (optically thin approximation)"
            )

        # Calculate power
        # ---------------------------------------------------

        self._reinitialize()  # creates scaled dataframe df1 from df0

        # ----------------------------------------------------------------------
        # Calculate Populations and Emission Integral
        # (Note: Emission Integral is non canonical quantity, equivalent to
        #  Linestrength for absorption)
        if non_eq_mode:
            self.calc_populations_noneq(
                Tvib,
                Trot,
                vib_distribution=vib_distribution,
                rot_distribution=rot_distribution,
            )
        else:
            self.calc_populations_eq(Tgas)
            self.df1["Aul"] = self.df1.A  # update einstein coefficients
        self.calc_emission_integral()

        #        # ----------------------------------------------------------------------
        #        # Cutoff linestrength  (note that cuting linestrength doesnt make this
        #        # process faster here, but we still give this option to be consistent
        #        # with spectra)
        #        self._cutoff_linestrength()

        # ----------------------------------------------------------------------

        # Sum over all emission integrals (in the valid range)
        b = (self.df1.wav > self.input.wavenum_min) & (
            self.df1.wav < self.input.wavenum_max
        )
        P = self.df1["Ei"][b].sum()  # Ei  >> (mW/sr)

        # incorporate density of molecules (see equation (A.16) )
        density = mole_fraction * ((pressure_mbar * 100) / (k_b * Tgas)) * 1e-6
        Pv = P * density  # (mW/sr/cm3)

        # Optically thin case (no self absorption):
        Ptot = Pv * path_length  # (mW/sr/cm2)

        self.profiler.stop(
            "optically_thin_power_calculation", "Optically thin power calculation"
        )

        return conv2(Ptot, "mW/cm2/sr", unit)

    def fit_spectrum(
        self,
        s_exp,
        model,
        fit_parameters,
        bounds={},
        plot=False,
        solver_options={"maxiter": 300},
        **kwargs,
    ) -> Union[Spectrum, OptimizeResult]:
        """Fit an experimental spectrum with an arbitrary model and an arbitrary
        number of fit parameters.

        Parameters
        ----------
        s_exp : Spectrum
            experimental spectrum. Should have only spectral array only. Use
            :py:meth:`~radis.spectrum.spectrum.Spectrum.take`, e.g::
                sf.fit_spectrum(s_exp.take('transmittance'))
        model : func -> Spectrum
            a line-of-sight model returning a Spectrum. Example :
            :py:func:`~radis.tools.fitting.LTEModel, `:py:func:`~radis.tools.fitting.Tvib12Tvib3Trot_NonLTEModel`
        fit_parameters : dict
            example::

                {fit_parameter:initial_value}
        bounds : dict, optional
            example::

                {fit_parameter:[min, max]}
        fixed_parameters : dict
            fixed parameters given to the model. Example::

                fit_spectrum(fixed_parameters={"vib_distribution":"treanor"})

        Other Parameters
        ----------------
        plot: bool
            if ``True``, plot spectra as they are computed; and plot the convergence of
            the residual. Default ``False``
        solver_options: dict
            parameters forwarded to the solver. More info in `~scipy.optimize.minimize`
            Example::

                {"maxiter": (int)  max number of iteration default ``300``,
                 }
        kwargs: dict
            forwarded to :py:func:`~radis.tools.fitting.fit_spectrum`

        Returns
        -------
        s_best: Spectrum
            best spectrum
        res: OptimizeResults
            output of `~scipy.optimize.minimize`

        Examples
        --------
        See a :ref:`one-temperature fit example <example_one_temperature_fit>`
        and a :ref:`non-LTE fit example <example_multi_temperature_fit>`

        .. minigallery:: radis.lbl.factory.SpectrumFactory.fit_spectrum

        More advanced tools for interactive fitting of multi-dimensional, multi-slabs
        spectra can be found in :py:mod:`fitroom`.

        See Also
        --------
        :py:func:`~radis.tools.fitting.fit_spectrum`,
        :py:func:`~radis.tools.fitting.Tvib12Tvib3Trot_NonLTEModel`,
        :py:mod:`fitroom`

        """
        from radis.tools.fitting import fit_spectrum

        return fit_spectrum(
            self,
            s_exp,
            model,
            fit_parameters,
            bounds=bounds,
            plot=plot,
            solver_options=solver_options,
            **kwargs,
        )

    def print_perf_profile(self, number_format="{:.3f}", precision=16):
        """Prints Profiler output dictionary in a structured manner for
        the last calculated spectrum

        Examples
        --------
        ::

            sf.print_perf_profile()

            # output >>
                spectrum_calculation      0.189s ████████████████
                    check_line_databank              0.000s
                    check_non_eq_param               0.042s ███
                    fetch_energy_5                   0.015s █
                    calc_weight_trans                0.008s
                    reinitialize                     0.002s
                        copy_database                    0.000s
                        memory_usage_warning             0.002s
                        reset_population                 0.000s
                    calc_noneq_population            0.041s ███
                        part_function                    0.035s ██
                        population                       0.006s
                    scaled_non_eq_linestrength       0.005s
                        map_part_func                    0.001s
                        corrected_population_se          0.003s
                    calc_emission_integral           0.006s
                    applied_linestrength_cutoff      0.002s
                    calc_lineshift                   0.001s
                    calc_hwhm                        0.007s
                    generate_wavenumber_arrays       0.001s
                    calc_line_broadening             0.074s ██████
                        precompute_LDM_lineshapes        0.012s
                        LDM_Initialized_vectors          0.000s
                        LDM_closest_matching_line        0.001s
                        LDM_Distribute_lines             0.001s
                        LDM_convolve                     0.060s █████
                        others                           0.001s
                    calc_other_spectral_quan         0.003s
                    generate_spectrum_obj            0.000s
                    others                           -0.016s

        Other Parameters
        ----------------
        precision: int, optional
            total number of blocks. Default 16.

        See Also
        --------

        :py:meth:`~radis.spectrum.spectrum.Spectrum.print_perf_profile`
        """
        from radis.spectrum.utils import print_perf_profile

        profiler = self.profiler.final
        total_time = profiler["spectrum_calculation"]["value"]

        return print_perf_profile(
            profiler, total_time, number_format=number_format, precision=precision
        )

    def generate_perf_profile(self):
        """Generate a visual/interactive performance profile diagram for
        the last calculated spectrum, with ``tuna``.

        Requires a `profiler` key with in Spectrum.conditions, and ``tuna``
        installed.

        .. note::
            :py:meth:`~radis.spectrum.spectrum.Spectrum.print_perf_profile` is an
            ascii-version which does not require ``tuna``.

        Examples
        --------
        ::

            sf = SpectrumFactory(...)
            sf.eq_spectrum(...)
            sf.generate_perf_profile()

        See typical output in https://github.com/radis/radis/pull/325

        .. image:: https://user-images.githubusercontent.com/16088743/128018032-6049be72-1881-46ac-9d7c-1ed89f9c4f42.png
            :alt: https://user-images.githubusercontent.com/16088743/128018032-6049be72-1881-46ac-9d7c-1ed89f9c4f42.png
            :target: https://user-images.githubusercontent.com/16088743/128018032-6049be72-1881-46ac-9d7c-1ed89f9c4f42.png


        .. note::
            You can also profile with `tuna` directly::

                python -m cProfile -o program.prof your_radis_script.py
                tuna your_radis_script.py


        See Also
        --------
        :py:meth:`~radis.spectrum.spectrum.Spectrum.print_perf_profile`
        """
        from radis.spectrum.utils import generate_perf_profile

        profiler = self.profiler.final.copy().copy()
        # Add total calculation time:
        profiler.update({"value": profiler["spectrum_calculation"]["value"]})
        # note: in Spectrum.generate_perf_profile the total time is taken as
        # 'self.conditions["calculation_time"]' which also includes database
        # loading times.

        return generate_perf_profile(profiler)


# %% ======================================================================
# EXTRA FUNCTIONS
# ---------------------
#
# _generate_wavenumber_range
# _generate_broadening_range
#
# XXX =====================================================================


def _generate_wavenumber_range(wavenum_min, wavenum_max, wstep, neighbour_lines):
    """define waverange vectors, with ``wavenumber`` the output spectral range
    and ``wavenumber_calc`` the spectral range used for calculation, that
    includes neighbour lines within ``neighbour_lines`` distance.

    Parameters
    ----------
    wavenum_min, wavenum_max: float
        wavenumber range limits (cm-1)
    wstep: float
        wavenumber step (cm-1)
    neighbour_lines: float
        wavenumber full width of broadening calculation: used to define which
        neighbour lines shall be included in the calculation

    Returns
    -------
    wavenumber: numpy array
        an evenly spaced array between ``wavenum_min`` and ``wavenum_max`` with
        a spacing of ``wstep``
    wavenumber_calc: numpy array
        an evenly spaced array between ``wavenum_min-neighbour_lines`` and
        ``wavenum_max+neighbour_lines`` with a spacing of ``wstep``
    woutrange: (wmin, wmax)
        index to project the full range including neighbour lines `wavenumber_calc`
        on the final range `wavenumber`, i.e. : wavenumber_calc[woutrange[0]:woutrange[1]] = wavenumber
    """
    assert wavenum_min < wavenum_max
    assert wstep > 0

    # Output range
    # generate the final vector of wavenumbers (shape M)
    wavenumber = arange(wavenum_min, wavenum_max + wstep, wstep)

    # generate the calculation vector of wavenumbers (shape M + space on the side)
    # ... Calculation range
    wavenum_min_calc = wavenumber[0] - neighbour_lines  # cm-1
    wavenum_max_calc = wavenumber[-1] + neighbour_lines  # cm-1
    w_out_of_range_left = arange(
        wavenumber[0] - wstep, wavenum_min_calc - wstep, -wstep
    )[::-1]
    w_out_of_range_right = arange(
        wavenumber[-1] + wstep, wavenum_max_calc + wstep, wstep
    )

    # ... deal with rounding errors: 1 side may have 1 more point
    if len(w_out_of_range_left) > len(w_out_of_range_right):
        w_out_of_range_left = w_out_of_range_left[1:]
    elif len(w_out_of_range_left) < len(w_out_of_range_right):
        w_out_of_range_right = w_out_of_range_right[:-1]

    wavenumber_calc = np.hstack((w_out_of_range_left, wavenumber, w_out_of_range_right))
    woutrange = len(w_out_of_range_left), len(w_out_of_range_left) + len(wavenumber)

    assert len(w_out_of_range_left) == len(w_out_of_range_right)
    assert len(wavenumber_calc) == len(wavenumber) + 2 * len(w_out_of_range_left)

    return wavenumber, wavenumber_calc, woutrange


def _generate_broadening_range(wstep, truncation):
    """Generate array on which to compute line broadening.

    Parameters
    ----------
    wstep: float
        wavenumber step (cm-1)
    truncation: float
        wavenumber half-width of broadening calculation: used to define which
        neighbour lines shall be included in the calculation

    Returns
    -------
    wbroad_centered: numpy array
        an evenly spaced array, of odd-parity length, centered on 0, and of width
        ``truncation``
    """

    # create a broadening array, on which lineshape will be calculated.
    # Odd number is important
    wbroad_centered = np.hstack(
        (
            -arange(wstep, truncation + wstep, wstep)[::-1],
            [0],
            arange(wstep, truncation + wstep, wstep),
        )
    )

    assert len(wbroad_centered) % 2 == 1

    return wbroad_centered


# %% Test

# --------------------------
if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/lbl/test_factory.py"]))
