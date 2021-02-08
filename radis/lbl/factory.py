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
import sys
from subprocess import call
from time import time
from warnings import warn

import astropy.units as u
import numpy as np
from numpy import arange, exp
from scipy.constants import N_A, c, k, pi

from radis import get_version
from radis.db import MOLECULES_LIST_EQUILIBRIUM, MOLECULES_LIST_NONEQUILIBRIUM
from radis.db.classes import get_molecule, get_molecule_identifier
from radis.db.molparam import MolParams
from radis.lbl.bands import BandFactory
from radis.lbl.base import get_waverange
from radis.misc import getProjectRoot
from radis.misc.basics import flatten, is_float, list_if_float
from radis.misc.printer import printg
from radis.misc.utils import Default
from radis.phys.constants import k_b
from radis.phys.convert import conv2
from radis.phys.units import convert_emi2nm, convert_rad2nm
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

    wmin, wmax : float or `~astropy.units.quantity.Quantity`
        a hybrid parameter which can stand for minimum (maximum) wavenumber or minimum
        (maximum) wavelength depending upon the unit accompanying it. If dimensionless,
        ``wunit`` is considered as the accompanying unit.
    wunit: string
        the unit accompanying wmin and wmax. Can only be passed with wmin
        and wmax. Default is `cm-1`.
    wavenum_min, wavenum_max: float(cm^-1) or `~astropy.units.quantity.Quantity`
        minimum (maximum) wavenumber to be processed in :math:`cm^{-1}`.
        use astropy.units to specify arbitrary inverse-length units.
    wavelength_min, wavelength_max : float(nm) or `~astropy.units.quantity.Quantity`
        minimum (maximum) wavelength to be processed in :math:`nm`. This wavelength
        can be in ``'air'`` or ``'vacuum'`` depending on the value of the parameter
        ``medium=``.
        use astropy.units to specify arbitrary length units.
    pressure: float(bar) or `~astropy.units.quantity.Quantity`
        partial pressure of gas in bar. Default 1.01325 (1 atm).
        use astropy.units to specify arbitrary pressure units.
        For example, ``1013.25 * u.mbar``.
    mole_fraction: N/D
        species mole fraction. Default 1. Note that the rest of the gas
        is considered to be air for collisional broadening.
    path_length: float(cm) or `~astropy.units.quantity.Quantity`
        path length in cm. Default 1.
        use astropy.units to specify arbitrary length units.
    molecule: int, str, or ``None``
        molecule id (HITRAN format) or name. If ``None``, the molecule can be infered
        from the database files being loaded. See the list of supported molecules
        in :py:data:`~radis.db.MOLECULES_LIST_EQUILIBRIUM`
        and :py:data:`~radis.db.MOLECULES_LIST_NONEQUILIBRIUM`.
        Default ``None``.
    isotope: int, list, str of the form '1,2', or 'all'
        isotope id (sorted by relative density: (eg: 1: CO2-626, 2: CO2-636 for CO2).
        See HITRAN documentation for isotope list for all species. If 'all',
        all isotopes in database are used (this may result in larger computation
        times!). Default 'all'
    medium: ``'air'``, ``'vacuum'``
        propagating medium when giving inputs with ``'wavenum_min'``, ``'wavenum_max'``.
        Does not change anything when giving inputs in wavenumber. Default ``'air'``

    Other Parameters
    ----------------

    *Computation parameters (see :py:attr:`~radis.lbl.loader.DatabankLoader.params`)*

    Tref: K
        Reference temperature for calculations (linestrength temperature
        correction). HITRAN database uses 296 Kelvin. Default 296 K
    self_absorption: boolean
        Compute self absorption. If ``False``, spectra are optically thin. Default ``True``.
    broadening_max_width: float (cm-1)
        Full width over which to compute the broadening. Large values will create
        a huge performance drop (scales as ~broadening_width^2 without DLM)
        The calculated spectral range is increased (by broadening_max_width/2
        on each side) to take into account overlaps from out-of-range lines.
        Default ``10`` cm-1.
    wstep: float (cm-1)
        Spacing of calculated spectrum. Default ``0.01`` cm-1
    cutoff: float (~ unit of Linestrength: cm-1/(#.cm-2))
        discard linestrengths that are lower that this, to reduce calculation
        times. ``1e-27`` is what is generally used to generate databases such as
        CDSD. If ``0``, no cutoff. Default ``1e-27``.
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
        If either ``"simple"`` or ``"min-RMS"`` DLM optimization for lineshape calculation is used:
        - ``"min-RMS"`` : weights optimized by analytical minimization of the RMS-error (See: [Spectral Synthesis Algorithm]_)
        - ``"simple"`` : weights equal to their relative position in the grid

        If using the DLM optimization, broadening method is automatically set to ``'fft'``.
        If ``None``, no lineshape interpolation is performed and the lineshape of all lines is calculated.

        Refer to [Spectral Synthesis Algorithm]_ for more explanation on the DLM method for lineshape interpolation.

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

        Fast fourier transform ``'fft'`` is only available if using the DLM lineshape
        calculation ``optimization``. Because the DLM convolves all lines at the same time,
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
                             broadening_max_width=10,  # cm-1
                             medium='vacuum',
                             verbose=1,    # more for more details
                             )
        sf.load_databank('HITRAN-CO2-TEST')        # predefined in ~/.radis
        s = sf.eq_spectrum(Tgas=300 * u.K, path_length=1 * u.cm)
        s.rescale_path_length(0.01)    # cm
        s.plot('radiance_noslit', Iunit='µW/cm2/sr/nm')

    Refer to the online :ref:`Examples <label_examples>` for more cases.

    .. inheritance-diagram:: radis.lbl.parallel.SpectrumFactory
       :parts: 1

    See Also
    --------

    Alternative:

    :func:`~radis.lbl.calc.calc_spectrum`

    Main Methods:

    :meth:`~radis.lbl.loader.DatabankLoader.load_databank`,
    :meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`,
    :meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum`

    For advanced use:

    :meth:`~radis.lbl.loader.DatabankLoader.fetch_databank`,
    :meth:`~radis.lbl.loader.DatabankLoader.init_databank`,
    :meth:`~radis.lbl.loader.DatabankLoader.init_database`,
    :meth:`~radis.lbl.bands.BandFactory.eq_bands`,
    :meth:`~radis.lbl.bands.BandFactory.non_eq_bands`

    Inputs and parameters can be accessed a posteriori with :

    :py:attr:`~radis.lbl.loader.DatabankLoader.input` : physical input
    :py:attr:`~radis.lbl.loader.DatabankLoader.params` : computational parameters
    :py:attr:`~radis.lbl.loader.DatabankLoader.misc` : miscallenous parameters (don't change output)
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
        broadening_max_width=10,
        pseudo_continuum_threshold=0,
        self_absorption=True,
        chunksize=None,
        optimization="simple",
        folding_thresh=1e-6,
        zero_padding=-1,
        broadening_method=Default("fft"),
        cutoff=1e-27,
        verbose=True,
        warnings=True,
        save_memory=False,
        export_populations=None,
        export_lines=False,
        **kwargs
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
        wavenum_min, wavenum_max = get_waverange(
            wmin,
            wmax,
            wunit,
            wavenum_min,
            wavenum_max,
            wavelength_min,
            wavelength_max,
            medium,
        )

        # calculated range is broader than output waverange to take into account off-range line broadening
        wavenumber, wavenumber_calc = _generate_wavenumber_range(
            wavenum_min, wavenum_max, wstep, broadening_max_width
        )
        wbroad_centered = _generate_broadening_range(wstep, broadening_max_width)
        # Store broadening max width and wstep as hidden variable (to ensure they are not changed afterwards)
        self._wstep = wstep
        self._broadening_max_width = broadening_max_width

        # Get boolean array that extracts the reduced range `wavenumber` from `wavenumber_calc`
        woutrange = np.in1d(wavenumber_calc, wavenumber, assume_unique=True)
        self.wbroad_centered = wbroad_centered
        self.wavenumber = wavenumber
        self.wavenumber_calc = wavenumber_calc
        self.woutrange = woutrange

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
                )

        # Store isotope identifier in str format (list wont work in database queries)
        if not isinstance(isotope, str):
            isotope = ",".join([str(k) for k in list_if_float(isotope)])

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
        self.params.cutoff = cutoff
        self.params.broadening_max_width = broadening_max_width  # line broadening
        self.misc.export_lines = export_lines
        self.misc.export_populations = export_populations
        self.params.wavenum_min_calc = wavenumber_calc[0]
        self.params.wavenum_max_calc = wavenumber_calc[-1]

        # if optimization is ``'simple'`` or ``'min-RMS'``, or None :
        # Adjust default values of broadening method :
        if isinstance(broadening_method, Default):
            if optimization in ("simple", "min-RMS") and broadening_method != "fft":
                if self.verbose >= 3:
                    printg(
                        "LDM algorithm used. Defaulting broadening method from {0} to FFT".format(
                            broadening_method
                        )
                    )
                broadening_method = "fft"
            elif optimization is None and broadening_method != "voigt":
                if self.verbose >= 3:
                    printg(
                        "LDM algorithm not used. Defaulting broadening method from {0} to 'voigt'".format(
                            broadening_method
                        )
                    )
                broadening_method = "voigt"
            else:  # keep default
                broadening_method = broadening_method.value
        self.params.broadening_method = broadening_method
        self.params.optimization = optimization
        self.params.folding_thresh = folding_thresh
        self.params.zero_padding = zero_padding

        # used to split lines into blocks not too big for memory
        self.misc.chunksize = chunksize
        # Other parameters:
        self.verbose = verbose
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
    #
    # XXX =====================================================================

    def eq_spectrum(
        self, Tgas, mole_fraction=None, path_length=None, pressure=None, name=None
    ):
        """Generate a spectrum at equilibrium.

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

        Returns
        -------

        s : Spectrum
            Returns a :class:`~radis.spectrum.spectrum.Spectrum` object

                Use the :meth:`~radis.spectrum.spectrum.Spectrum.get` method to get something
                among ``['radiance', 'radiance_noslit', 'absorbance', etc...]``

                Or directly the :meth:`~radis.spectrum.spectrum.Spectrum.plot` method
                to plot it. See [1]_ to get an overview of all Spectrum methods

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

        # Check variables
        self._check_inputs(mole_fraction, max(flatten(Tgas)))

        # Retrieve Spectrum from database if it exists
        if self.autoretrievedatabase:
            s = self._retrieve_from_database()
            if s is not None:
                return s  # exit function

        # %% Start
        # --------------------------------------------------------------------

        t0 = time()
        if verbose:
            self.print_conditions("Calculating Equilibrium Spectrum")

        # Check database, reset populations, create line dataframe to be scaled
        # --------------------------------------------------------------------
        self._check_line_databank()
        self._reinitialize()  # creates scaled dataframe df1 from df0

        # --------------------------------------------------------------------

        # First calculate the linestrength at given temperature
        self._calc_linestrength_eq(Tgas)  # scales S0 to S (equivalent to S0 in code)
        self._cutoff_linestrength()

        # ----------------------------------------------------------------------

        # Calculate line shift
        self._calc_lineshift()  # scales wav to shiftwav (equivalent to v0)

        # ----------------------------------------------------------------------
        # Line broadening

        # ... calculate broadening  HWHM
        self._calc_broadening_HWHM()

        # ... find weak lines and calculate semi-continuum (optional)
        I_continuum = self._calculate_pseudo_continuum()
        # ... apply lineshape and get absorption coefficient
        # ... (this is the performance bottleneck)
        wavenumber, abscoeff_v = self._calc_broadening()
        #    :         :
        #   cm-1    1/(#.cm-2)

        # ... add semi-continuum (optional)
        abscoeff_v = self._add_pseudo_continuum(abscoeff_v, I_continuum)
        # Calculate output quantities
        # ----------------------------------------------------------------------

        if self.verbose >= 2:
            t1 = time()

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

        if self.verbose >= 2:
            printg(
                "Calculated other spectral quantities in {0:.2f}s".format(time() - t1)
            )

        # %% Export
        # --------------------------------------------------------------------

        t = round(time() - t0, 2)
        if verbose >= 2:
            printg(
                "Spectrum calculated in {0:.2f}s (before object generation)".format(t)
            )
        if self.verbose >= 2:
            t1 = time()

        # Get conditions
        conditions = self.get_conditions()
        conditions.update(
            {
                "calculation_time": t,
                "lines_calculated": self._Nlines_calculated,
                "lines_cutoff": self._Nlines_cutoff,
                "lines_in_continuum": self._Nlines_in_continuum,
                "thermal_equilibrium": True,
                "radis_version": get_version(),
            }
        )

        # Get populations of levels as calculated in RovibrationalPartitionFunctions
        # ... Populations cannot be calculated at equilibrium (needs energies).
        # ... Use SpectrumFactory.non_eq_spectrum
        populations = None

        # Get lines (intensities + populations)
        lines = self.get_lines()

        # Spectral quantities
        quantities = {
            "abscoeff": (wavenumber, abscoeff),
            "absorbance": (wavenumber, absorbance),
            "emissivity_noslit": (wavenumber, emissivity_noslit),
            "transmittance_noslit": (wavenumber, transmittance_noslit),
            # (mW/cm2/sr/nm)
            "radiance_noslit": (wavenumber, radiance_noslit),
        }
        if I_continuum is not None and self._export_continuum:
            quantities.update(
                {"abscoeff_continuum": (wavenumber, I_continuum * density)}
            )

        # Store results in Spectrum class
        s = Spectrum(
            quantities=quantities,
            conditions=conditions,
            populations=populations,
            lines=lines,
            units=self.units,
            cond_units=self.cond_units,
            waveunit=self.params.waveunit,  # cm-1
            # dont check input (much faster, and Spectrum
            warnings=False,
            # is freshly baken so probably in a good format
            name=name,
        )

        # update database if asked so
        if self.autoupdatedatabase:
            self.SpecDatabase.add(s, if_exists_then="increment")
            # Tvib=Trot=Tgas... but this way names in a database
            # generated with eq_spectrum are consistent with names
            # in one generated with non_eq_spectrum

        # Get generation & total calculation time
        if self.verbose >= 2:
            printg("Generated Spectrum object in {0:.2f}s".format(time() - t1))

        #  In the less verbose case, we print the total calculation+generation time:
        t = round(time() - t0, 2)
        if verbose:
            print("Spectrum calculated in {0:.2f}s".format(t))

        return s

    def eq_spectrum_gpu(
        self, Tgas, mole_fraction=None, path_length=None, pressure=None, name=None
    ):

        """Generate a spectrum at equilibrium with calculation of lineshapes
        and broadening done on the GPU.

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

        Returns
        -------
        s : Spectrum
            Returns a :class:`~radis.spectrum.spectrum.Spectrum` object

                Use the :meth:`~radis.spectrum.spectrum.Spectrum.get` method to get something
                among ``['radiance', 'radiance_noslit', 'absorbance', etc...]``

                Or directly the :meth:`~radis.spectrum.spectrum.Spectrum.plot` method
                to plot it. See [1]_ to get an overview of all Spectrum methods

        Notes
        -----

        This method requires CUDA compatible hardware to execute. For more information on how to setup your system to run GPU-accelerated methods using CUDA and Cython, check `GPU Spectrum Calculation on RADIS <https://radis.readthedocs.io/en/latest/lbl/gpu.html>`

        See Also
        --------

        :meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`
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

        molpar = MolParams()

        try:
            id_set = self.df0[
                "id"
            ].unique()  # get all the molecules in the dataframe, should ideally be 1 element for GPU
            mol_id = id_set[0]
            molecule = get_molecule(mol_id)
        except:
            mol_id = get_molecule_identifier(self.input.molecule)
            molecule = get_molecule(mol_id)

        state = self.input.state
        iso_set = self._get_isotope_list(molecule)

        iso_arr = list(range(max(iso_set) + 1))

        Ia_arr = np.empty_like(iso_arr, dtype=np.float32)  # abundance of each isotope
        molarmass_arr = np.empty_like(
            iso_arr, dtype=np.float32
        )  # molar mass of each isotope
        Q_arr = np.empty_like(
            iso_arr, dtype=np.float32
        )  # partitioning function of each isotope
        for iso in iso_arr:
            if iso in iso_set:
                params = molpar.df.loc[(mol_id, iso)]
                Ia_arr[iso] = params.abundance
                molarmass_arr[iso] = params.mol_mass
                Q_arr[iso] = self._get_parsum(molecule, iso, state).at(T=Tgas)

        Ia_arr[np.isnan(Ia_arr)] = 0
        molarmass_arr[np.isnan(molarmass_arr)] = 0

        ### EXPERIMENTAL ###

        project_path = getProjectRoot()
        project_path += "/lbl/py_cuffs/"
        sys.path.insert(1, project_path)

        try:
            import py_cuffs
        except:
            try:

                if verbose >= 2:
                    print("py_cuFFS module not found in directory...")
                    print("Compiling module from source...")

                call(
                    "python setup.py build_ext --inplace",
                    cwd=project_path,
                    shell=True,
                )

                if verbose >= 2:
                    print("Finished compilation...trying to import module again")
                import py_cuffs

                if verbose:
                    print("py_cuFFS imported succesfully!")
            except:
                raise (
                    ModuleNotFoundError(
                        "Failed to load py_cuFFS module, program will exit."
                    )
                )
                exit()

        ### --- ###

        t0 = time()

        # generate the v_arr
        v_arr = np.arange(
            self.input.wavenum_min,
            self.input.wavenum_max + self._wstep,
            self._wstep,
        )

        # load the data
        df = self.df0
        iso = df["iso"].to_numpy(dtype=np.int32)
        v0 = df["wav"].to_numpy(dtype=np.float32)
        da = df["Pshft"].to_numpy(dtype=np.float32)
        El = df["El"].to_numpy(dtype=np.float32)
        na = df["Tdpair"].to_numpy(dtype=np.float32)

        log_2gs = np.array(self._get_log_2gs(), dtype=np.float32)
        log_2vMm = np.array(self._get_log_2vMm(molarmass_arr), dtype=np.float32)
        S0 = np.array(self._get_S0(Ia_arr), dtype=np.float32)

        NwG = 4
        NwL = 8

        _Nlines_calculated = len(v0)

        if verbose >= 2:
            print("Initializing parameters...", end=" ")

        if verbose is False:
            verbose_gpu = 0
        elif verbose is True:
            verbose_gpu = 1
        else:
            verbose_gpu = verbose

        py_cuffs.init(
            v_arr,
            NwG,
            NwL,
            iso,
            v0,
            da,
            log_2gs,
            na,
            log_2vMm,
            S0,
            El,
            Q_arr,
            verbose_gpu,
        )

        if verbose >= 2:
            print("Initialization complete!")

        wavenumber = v_arr

        if verbose >= 2:
            print("Calculating spectra...", end=" ")

        abscoeff = py_cuffs.iterate(
            pressure_mbar * 1e-3,
            Tgas,
            mole_fraction,
            Ia_arr,
            molarmass_arr,
            verbose_gpu,
        )
        # Calculate output quantities
        # ----------------------------------------------------------------------
        if verbose >= 2:
            t1 = time()

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
        if verbose >= 2:
            printg(
                "Calculated other spectral quantities in {0:.2f}s".format(time() - t1)
            )

        lines = self.get_lines()

        # %% Export
        # --------------------------------------------------------------------
        t = round(time() - t0, 2)
        if verbose >= 2:
            t1 = time()
        # Get lines (intensities + populations)

        conditions = self.get_conditions()
        conditions.update(
            {
                "calculation_time": t,
                "lines_calculated": _Nlines_calculated,
                "thermal_equilibrium": True,
                "radis_version": get_version(),
            }
        )

        # Spectral quantities
        quantities = {
            "abscoeff": (wavenumber, abscoeff),
            "absorbance": (wavenumber, absorbance),
            "emissivity_noslit": (wavenumber, emissivity_noslit),
            "transmittance_noslit": (wavenumber, transmittance_noslit),
            # (mW/cm2/sr/nm)
            "radiance_noslit": (wavenumber, radiance_noslit),
        }

        # Store results in Spectrum class
        s = Spectrum(
            quantities=quantities,
            units=self.units,
            conditions=conditions,
            lines=lines,
            cond_units=self.cond_units,
            waveunit=self.params.waveunit,  # cm-1
            # dont check input (much faster, and Spectrum
            warnings=False,
            # is freshly baken so probably in a good format
            name=name,
        )

        # update database if asked so
        if self.autoupdatedatabase:
            self.SpecDatabase.add(s, if_exists_then="increment")
            # Tvib=Trot=Tgas... but this way names in a database
            # generated with eq_spectrum are consistent with names
            # in one generated with non_eq_spectrum

        # Get generation & total calculation time
        if verbose >= 2:
            printg("Generated Spectrum object in {0:.2f}s".format(time() - t1))

        #  In the less verbose case, we print the total calculation+generation time:

        return s

    def non_eq_spectrum(
        self,
        Tvib,
        Trot,
        Ttrans=None,
        mole_fraction=None,
        path_length=None,
        pressure=None,
        vib_distribution="boltzmann",
        rot_distribution="boltzmann",
        overpopulation=None,
        name=None,
    ):
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
            add overpopulation factors for given levels:

            >>> {level:overpopulation_factor}

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

        # Check variables
        self._check_inputs(mole_fraction, max(flatten(Tgas, Tvib, Trot)))

        # Retrieve Spectrum from database if it exists
        if self.autoretrievedatabase:
            s = self._retrieve_from_database()
            if s is not None:
                return s  # exit function

        # %% Start
        # --------------------------------------------------------------------

        t0 = time()
        if verbose:
            self.print_conditions("Calculating Non-Equilibrium Spectrum")

        # Check line database and parameters, reset populations and scaled line dataframe
        # ----------
        self._check_line_databank()
        # add nonequilibrium energies if needed (this may be a bottleneck
        # for a first calculation):
        self._check_noneq_parameters(vib_distribution, singleTvibmode)
        self._reinitialize()  # creates scaled dataframe df1 from df0

        # ----------------------------------------------------------------------
        # Calculate Populations, Linestrength and Emission Integral
        if singleTvibmode:
            self._calc_populations_noneq(
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

        self._calc_linestrength_noneq()
        self._calc_emission_integral()

        # ----------------------------------------------------------------------
        # Cutoff linestrength
        self._cutoff_linestrength()

        # ----------------------------------------------------------------------

        # Calculate lineshift
        self._calc_lineshift()

        # ----------------------------------------------------------------------

        # Line broadening

        # ... calculate broadening  HWHM
        self._calc_broadening_HWHM()

        # ... find weak lines and calculate semi-continuum (optional)
        k_continuum, j_continuum = self._calculate_pseudo_continuum(noneq=True)

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

        if self.verbose >= 2:
            t1 = time()

        # incorporate density of molecules (see Rothman 1996 equation (A.16) )
        density = mole_fraction * ((pressure_mbar * 100) / (k_b * Tgas)) * 1e-6
        #  :
        # (#/cm3)

        abscoeff = abscoeff_v * density  # cm-1
        emisscoeff = emisscoeff_v * density  # mW/sr/cm3/cm-1

        # ... # TODO: if the code is extended to multi-species, then density has to be added
        # ... before lineshape broadening (as it would not be constant for all species)

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
            radiance_noslit = emisscoeff * path_length  # (mW/sr/cm2/cm-1)

        # Convert `radiance_noslit` from (mW/sr/cm2/cm-1) to (mW/sr/cm2/nm)
        radiance_noslit = convert_rad2nm(
            radiance_noslit, wavenumber, "mW/sr/cm2/cm-1", "mW/sr/cm2/nm"
        )
        # Convert 'emisscoeff' from (mW/sr/cm3/cm-1) to (mW/sr/cm3/nm)
        emisscoeff = convert_emi2nm(
            emisscoeff, wavenumber, "mW/sr/cm3/cm-1", "mW/sr/cm3/nm"
        )

        if self.verbose >= 2:
            printg(
                "Calculated other spectral quantities in {0:.2f}s".format(time() - t1)
            )

        # Note: emissivity not defined under non equilibrium

        # %% Export
        # ----------------------------------------------------------------------

        t = round(time() - t0, 2)
        if verbose >= 2:
            printg(
                "Spectrum calculated in {0:.2f}s (before object generation)".format(t)
            )
        if self.verbose >= 2:
            t1 = time()

        # Get conditions
        conditions = self.get_conditions()
        conditions.update(
            {
                "calculation_time": t,
                "lines_calculated": self._Nlines_calculated,
                "lines_cutoff": self._Nlines_cutoff,
                "lines_in_continuum": self._Nlines_in_continuum,
                "thermal_equilibrium": False,  # dont even try to guess if it's at equilibrium
                "radis_version": get_version(),
            }
        )

        # Get populations of levels as calculated in RovibrationalPartitionFunctions
        populations = self.get_populations(self.misc.export_populations)

        # Get lines (intensities + populations)
        lines = self.get_lines()

        # Spectral quantities
        quantities = {
            "abscoeff": (wavenumber, abscoeff),
            "absorbance": (wavenumber, absorbance),
            # (mW/cm3/sr/nm)
            "emisscoeff": (wavenumber, emisscoeff),
            "transmittance_noslit": (wavenumber, transmittance_noslit),
            # (mW/cm2/sr/nm)
            "radiance_noslit": (wavenumber, radiance_noslit),
        }
        if k_continuum is not None and self._export_continuum:
            quantities.update(
                {
                    "abscoeff_continuum": (wavenumber, k_continuum * density),
                    "emisscoeff_continuum": (wavenumber, j_continuum * density),
                }
            )

        # Store results in Spectrum class
        s = Spectrum(
            quantities=quantities,
            conditions=conditions,
            populations=populations,
            lines=lines,
            units=self.units,
            cond_units=self.cond_units,
            waveunit=self.params.waveunit,  # cm-1
            # dont check input (much faster, and Spectrum
            warnings=False,
            # is freshly baken so probably in a good format
            name=name,
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
        if self.verbose >= 2:
            printg("Generated Spectrum object in {0:.2f}s".format(time() - t1))

        #  In the less verbose case, we print the total calculation+generation time:
        t = round(time() - t0, 2)
        if verbose:
            print("Spectrum calculated in {0:.2f}s".format(t))

        return s

    def _get_log_2gs(self):
        """Returns log_2gs if it already exists in the dataframe, otherwise
        computes it using gamma_air."""
        df = self.df0
        # TODO: deal with the case of gamma_self [so we don't forget]

        # if the column already exists, then return
        if "log_2gs" in df.columns:
            return df["log_2gs"]

        try:
            gamma_air = df["airbrd"].to_numpy()
            log_2gs = np.log(2 * gamma_air)
            df["log_2gs"] = log_2gs
            return log_2gs
        except KeyError as err:
            raise KeyError(
                "Cannot find air-broadened half-width or log_2gs in the database... please check the database"
            ) from err

    def _get_log_2vMm(self, molarmass_arr):
        """Returns log_2vMm if it already exists in the dataframe, otherwise
        computes it using the abundance and molar mass for each isotope passed
        in the input."""
        df = self.df0

        # if the column already exists, then return
        if "log_2vMm" in df.columns:
            return df["log_2vMm"]

        try:
            v0 = df["wav"].to_numpy()  # get wavenumber
            iso = df["iso"].to_numpy()  # get isotope
            Mm = molarmass_arr * 1e-3 / N_A
            log_2vMm = np.log(2 * v0) + 0.5 * np.log(
                2 * k * np.log(2) / (c ** 2 * Mm.take(iso))
            )
            df["log_2vMm"] = log_2vMm
            return log_2vMm
        except KeyError as err:
            raise KeyError(
                "Cannot find wavenumber, isotope and/or log_2vMm in the database. Please check the database"
            ) from err

    def _get_S0(self, Ia_arr):
        """Returns S0 if it already exists, otherwise computes the value using
        abundance, gamma_air and Einstein's number."""
        df = self.df0

        # if the column already exists, then return it
        if "S0" in df.columns:
            return df["S0"]

        try:
            v0 = df["wav"].to_numpy()
            iso = df["iso"].to_numpy()
            A21 = df["A"].to_numpy()
            Jl = df["jl"].to_numpy()
            DJ = df["branch"].to_numpy()
            Ju = Jl + DJ
            gu = 2 * Ju + 1  # g_up
            S0 = Ia_arr.take(iso) * gu * A21 / (8 * pi * c_cm * v0 ** 2)
            df["S0"] = S0
            return S0
        except KeyError as err:
            raise KeyError(
                "Could not find wavenumber, Einstein's coefficient, lower state energy or S0 in the dataframe. PLease check the database"
            ) from err

    def optically_thin_power(
        self,
        Tgas=None,
        Tvib=None,
        Trot=None,
        Ttrans=None,
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

        Returns total power density in mW/cm2/sr (unless different unit is chosen),
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

        # Make sure database is loaded
        if self.df0 is None:
            if not self.save_memory:
                raise AttributeError("Load databank first (.load_databank())")
            else:
                self._reload_databank()

        if non_eq_mode:
            # Make sure database has pre-computed non equilibrium quantities
            # (Evib, Erot, etc.)
            try:
                self.df0["Evib"]
            except KeyError:
                self._calc_noneq_parameters()

            try:
                self.df0["Aul"]
            except KeyError:
                self._calc_weighted_trans_moment()
                self._calc_einstein_coefficients()

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
            self._calc_populations_noneq(Tvib, Trot)
        else:
            self._calc_populations_eq(Tgas)
            self.df1["Aul"] = self.df1.A  # update einstein coefficients
        self._calc_emission_integral()

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

        return conv2(Ptot, "mW/cm2/sr", unit)


# %% ======================================================================
# EXTRA FUNCTIONS
# ---------------------
#
# _generate_wavenumber_range
# _generate_broadening_range
#
# XXX =====================================================================


def _generate_wavenumber_range(wavenum_min, wavenum_max, wstep, broadening_max_width):
    """define waverange vectors, with ``wavenumber`` the ouput spectral range
    and ``wavenumber_calc`` the spectral range used for calculation, that
    includes neighbour lines within ``broadening_max_width`` distance.

    Parameters
    ----------

    wavenum_min, wavenum_max: float
        wavenumber range limits (cm-1)

    wstep: float
        wavenumber step (cm-1)

    broadening_max_width: float
        wavenumber full width of broadening calculation: used to define which
        neighbour lines shall be included in the calculation

    Returns
    -------

    wavenumber: numpy array
        an evenly spaced array between ``wavenum_min`` and ``wavenum_max`` with
        a spacing of ``wstep``

    wavenumber_calc: numpy array
        an evenly spaced array between ``wavenum_min-broadening_max_width/2`` and
        ``wavenum_max+broadening_max_width/2`` with a spacing of ``wstep``
    """

    assert wavenum_min < wavenum_max

    # Output range
    # generate the final vector of wavenumbers (shape M)
    wavenumber = arange(wavenum_min, wavenum_max + wstep, wstep)

    # generate the calculation vector of wavenumbers (shape M + space on the side)
    # ... Calculation range
    wavenum_min_calc = wavenumber[0] - broadening_max_width / 2  # cm-1
    wavenum_max_calc = wavenumber[-1] + broadening_max_width / 2  # cm-1
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

    assert len(w_out_of_range_left) == len(w_out_of_range_right)
    assert len(wavenumber_calc) == len(wavenumber) + 2 * len(w_out_of_range_left)

    return wavenumber, wavenumber_calc


def _generate_broadening_range(wstep, broadening_max_width):
    """Generate array on which to compute line broadening.

    Parameters
    ----------

    wstep: float
        wavenumber step (cm-1)

    broadening_max_width: float
        wavenumber full width of broadening calculation: used to define which
        neighbour lines shall be included in the calculation

    Returns
    -------

    wbroad_centered: numpy array
        an evenly spaced array, of odd-parity length, centered on 0, and of width
        ``broadening_max_width``
    """

    # create a broadening array, on which lineshape will be calculated.
    # Odd number is important
    wbroad_centered = np.hstack(
        (
            -arange(wstep, 0.5 * broadening_max_width + wstep, wstep)[::-1],
            [0],
            arange(wstep, 0.5 * broadening_max_width + wstep, wstep),
        )
    )

    assert len(wbroad_centered) % 2 == 1

    return wbroad_centered


# %% Test

# --------------------------
if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/lbl/test_factory.py"]))
