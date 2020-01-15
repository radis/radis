# -*- coding: utf-8 -*-
"""
Contains the :py:class:`~radis.lbl.factory.SpectrumFactory` class, which is 
the core of the RADIS Line-by-Line module. 

Examples
--------

Calculate a CO Spectrum, fetching the lines from HITRAN ::

    # This is how you get a spectrum (see spectrum.py for front-end functions
    # that do just that)
    sf = SpectrumFactory(2125, 2249.9,
                         parallel=False,bplot=False,
                         molecule='CO', 
                         isotope=1,
                         cutoff=1e-30,   # for faster calculations. See
                                         # `plot_linestrength_hist` for more details
                         **kwargs)
    sf.fetch_databank()        # autodownload from HITRAN 
    s = sf.eq_spectrum(Tgas=300)
    s.plot('abscoeff')

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
:py:class:`~radis.lbl.factory.SpectrumFactory` > :py:class:`~radis.lbl.parallel.ParallelFactory`


.. inheritance-diagram:: radis.lbl.parallel.ParallelFactory
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
from __future__ import print_function, absolute_import, division, unicode_literals
from six import string_types
from warnings import warn
from radis.io.hitran import get_molecule
from radis.lbl.bands import BandFactory
from radis.lbl.base import get_waverange
from radis.spectrum.spectrum import Spectrum
from radis.spectrum.equations import calc_radiance
from radis.misc.basics import is_float, list_if_float, flatten
from radis.misc.printer import printg
from radis.phys.convert import conv2
from radis.phys.constants import k_b
from radis.phys.units import convert_rad2nm, convert_emi2nm
from radis.phys.units_astropy import convert_and_strip_units
from radis import get_version
from numpy import exp, arange
from multiprocessing import cpu_count
from time import time
import numpy as np
import astropy.units as u

# %% Main functions


class SpectrumFactory(BandFactory):
    """ A class to put together all functions related to loading CDSD / HITRAN
    databases, calculating the broadenings, and summing over all the lines

    Parameters
    ----------

    wavenum_min: float(cm^-1) or `~astropy.units.quantity.Quantity`
        minimum wavenumber to be processed in cm^-1. 
        use astropy.units to specify arbitrary inverse-length units.

    wavenum_max: float(cm^-1) or `~astropy.units.quantity.Quantity`
        maximum wavenumber to be processed in cm^-1.
        use astropy.units to specify arbitrary inverse-length units.

    wavelength_min: float(nm) or `~astropy.units.quantity.Quantity`
        minimum wavelength to be processed in nm. This wavelength
        can be in ``'air'`` or ``'vacuum'`` depending on the value of the parameter
        ``medium=``.
        use astropy.units to specify arbitrary length units.

    wavelength_max: float(nm) or `~astropy.units.quantity.Quantity`
        maximum wavelength to be processed in nm.
        use astropy.units to specify arbitrary length units.

    Tref: float(K) or `~astropy.units.quantity.Quantity`
        reference temperature for calculations, HITRAN database uses 296 Kelvin
        default: 3400K. 
        use astropy.units to specify arbitrary temperature units.
        For example, ``200 * u.deg_C``.

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
        in :py:data:`~radis.io.MOLECULES_LIST_EQUILIBRIUM`
        and :py:data:`~radis.io.MOLECULES_LIST_NONEQUILIBRIUM`. 
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

    Computation parameters:

    Tref: K
        Reference temperature for calculations (linestrength temperature
        correction). HITRAN database uses 296 Kelvin. Default 296 K

    self_absorption: boolean
        self absorption

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

    bplot: boolean
        plot intermediary results (like slit function generation). Default ``False``.

    save_memory: boolean
        if ``True``, removes databases calculated by intermediate functions (for
        instance, delete the full database once the linestrength cutoff criteria
        was applied). This saves some memory but requires to reload the database
        & recalculate the linestrength for each new parameter. Default ``False``.

    export_populations: ``'vib'``, ``'rovib'``, ``None``
        if not None, store populations in Spectrum. Either store vibrational
        populations ('vib') or rovibrational populations ('rovib'). Default ``None``

    export_lines: boolean
        if ``True``, saves lines in Spectrum. Default ``True``.

    parallel: boolean
        use parallel compution. May only offer a performance improvement if
        one spectrum already takes > 10s to be calculated. Default ``False`` because
        since the latest vectorized version parallelisation is mostly
        counterproductive.

    db_use_cached: boolean, or 'regen'
        use cache for line databases.
        If ``True``, a pandas-readable csv file is generated on first access,
        and later used. This saves on the datatype cast and conversion and
        improves performances a lot. But! ... be sure to delete these files
        to regenerate them if you happen to change the database. If 'regen',
        existing cached files are removed and regenerated. Default ``False``
        From 0.9.16 it is also used to load energy levels from .h5 cache file
        if exist

    lvl_use_cached: boolean, or ``'regen'``, or ``'force'``, or ``None``
        use cache for energy level database
        If None, the same value as ``db_use_cached`` is used.

    Nprocs, Ngroups: int
        parameters used in parallel processing mode. Default ``None``

    chunksize: int, or ``None``
        Splits the lines database in several chuncks during calculation, else
        the multiplication of lines over all spectral range takes too much memory
        and slows the system down. Chunksize let you change the default chunck
        size. If ``None``, all lines are processed directly. Usually faster but
        can create memory problems. Default ``None``
        
        .. note::
            in version 0.9.20 this parameter is temporarily used to accept the ``'DLM'`` 
            argument: in this case, the DLM optimization for lineshape calculation 
            is used. Broadening method is automatically set to ``'fft'``. 
            See :py:attr:`~radis.lbl.broadening.BroadenFactory._broadening_method`.
        
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

    :func:`~radis.lbl.calc.calc_spectrum`, 
    :class:`~radis.lbl.parallel.ParallelFactory`

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
        Nprocs=None,
        Ngroups=None,
        cutoff=1e-27,
        bplot=False,
        parallel=False,
        db_use_cached=True,
        lvl_use_cached=None,
        verbose=True,
        warnings=True,
        save_memory=False,
        export_populations=None,
        export_lines=True,
        **kwargs
    ):

        # Initialize BandFactory object
        super(SpectrumFactory, self).__init__()

        # Check inputs (deal with deprecated format)
        if medium not in ["air", "vacuum"]:
            raise ValueError("Wavelength must be one of: 'air', 'vacuum'")
        kwargs0 = kwargs  # kwargs is used to deal with Deprecated names
        if "use_cached" in kwargs:
            warn(DeprecationWarning("use_cached replaced with db_use_cached"))
            db_use_cached = kwargs0.pop("use_cached")
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
            wavenum_min, wavenum_max, wavelength_min, wavelength_max, medium
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

        # Store isotope identifier in str format (list wont work in database queries)
        if not isinstance(isotope, string_types):
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
        self.params.db_use_cached = db_use_cached
        if lvl_use_cached is None:
            lvl_use_cached = db_use_cached
        self.params.lvl_use_cached = lvl_use_cached
        self.params.pseudo_continuum_threshold = pseudo_continuum_threshold
        self.misc.parallel = parallel
        self.params.cutoff = cutoff
        self.params.broadening_max_width = broadening_max_width  # line broadening
        self.misc.export_lines = export_lines
        self.misc.export_populations = export_populations
        self.params.wavenum_min_calc = wavenumber_calc[0]
        self.params.wavenum_max_calc = wavenumber_calc[-1]

        # in version 0.9.20 the 'chunksize' parameter is temporarily used to accept the ``'DLM'``
        # argument: in this case, the DLM optimization for lineshape calculation
        # is used. Broadening method is automatically set to ``'fft'``.
        # See :py:attr:`~radis.lbl.broadening.BroadenFactory._broadening_method`.
        if chunksize == "DLM":
            self._broadening_method = "fft"
            if self.verbose >= 3:
                print("DLM used. Defaulting broadening method to FFT")
            # TODO: make it a proper parameter in self.misc or self.params

        # used to split lines into blocks not too big for memory
        self.misc.chunksize = chunksize
        if parallel:
            # Set up parameters for //
            if Ngroups is None:
                Ngroups = cpu_count()  # number of processors
            if Nprocs is None:
                Nprocs = cpu_count()
            self.misc.Ngroups = Ngroups
            self.misc.Nprocs = Nprocs
        else:
            if Nprocs is not None:
                print("Choose parallel=True to use Nprocs")
            if Ngroups is not None:
                print("Choose parallel=True to use Ngroups")

        # Other parameters:
        self.bplot = bplot
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
        """ Generate a spectrum at equilibrium

        Parameters
        ----------

        Tgas: float
            Gas temperature (K)

        mole_fraction: float
            database species mole fraction. If None, Factory mole fraction is used.

        path_length: float
            slab size (cm). If None, Factory mole fraction is used.

        pressure: float
            pressure (bar). If None, the default Factory pressure is used.

        name: str
            case name (useful in batch)

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

        Calculate line strenghts correcting the CDSD reference one. Then call
        the main routine that sums over all lines


        References
        ----------

        .. [1] RADIS doc: `Spectrum how to? <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`__


        See Also
        --------

        :meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum`

        """

        try:

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
                    "Tgas should be float. Use ParallelFactory for multiple cases"
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
            self._calc_linestrength_eq(Tgas)
            self._cutoff_linestrength()

            # ----------------------------------------------------------------------

            # Calculate line shift
            self._calc_lineshift()

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
                    "Calculated other spectral quantities in {0:.2f}s".format(
                        time() - t1
                    )
                )

            # %% Export
            # --------------------------------------------------------------------

            t = round(time() - t0, 2)
            if verbose >= 2:
                printg(
                    "Spectrum calculated in {0:.2f}s (before object generation)".format(
                        t
                    )
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
                    "radis_version": get_version(add_git_number=False),
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

        except:
            # An error occured: clean before crashing
            self._clean_temp_file()
            raise

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
        """ Calculate emission spectrum in non-equilibrium case. Calculates
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

        path_length: float
            slab size (cm). If None, Factory mole fraction is used.

        pressure: float
            pressure (bar). If None, the default Factory pressure is used.

        vib_distribution: ``'boltzmann'``, ``'treanor'``
            vibrational distribution

        rot_distribution: ``'boltzmann'``
            rotational distribution

        overpopulation: dict, or ``None``
            add overpopulation factors for given levels:

            >>> {level:overpopulation_factor}

        name: str
            case name (useful in batch)

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

        try:

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
                    + "For parallel processing use ParallelFactory with a "
                    + "list of float or a list of tuple"
                )
            singleTvibmode = is_float(Tvib)
            if not is_float(Trot):
                raise ValueError(
                    "Trot should be float. Use ParallelFactory for multiple cases"
                )
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
            #   cm-1    1/(#.cm-2)   mW/sr/cm_1

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
            emisscoeff = emisscoeff_v * density  # mW/sr/cm3/cm_1

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
                radiance_noslit = emisscoeff * path_length  # (mW/sr/cm2/cm_1)

            # Convert `radiance_noslit` from (mW/sr/cm2/cm_1) to (mW/sr/cm2/nm)
            radiance_noslit = convert_rad2nm(
                radiance_noslit, wavenumber, "mW/sr/cm2/cm_1", "mW/sr/cm2/nm"
            )
            # Convert 'emisscoeff' from (mW/sr/cm3/cm_1) to (mW/sr/cm3/nm)
            emisscoeff = convert_emi2nm(
                emisscoeff, wavenumber, "mW/sr/cm3/cm_1", "mW/sr/cm3/nm"
            )

            if self.verbose >= 2:
                printg(
                    "Calculated other spectral quantities in {0:.2f}s".format(
                        time() - t1
                    )
                )

            # Note: emissivity not defined under non equilibrium

            # %% Export
            # ----------------------------------------------------------------------

            t = round(time() - t0, 2)
            if verbose >= 2:
                printg(
                    "Spectrum calculated in {0:.2f}s (before object generation)".format(
                        t
                    )
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
                    "radis_version": get_version(add_git_number=False),
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

        except:
            # An error occured: clean before crashing
            self._clean_temp_file()
            raise

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
        """ Calculate total power emitted in equilibrium or non-equilibrium case
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

        try:

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

        except:
            # An error occured: clean before crashing
            self._clean_temp_file()
            raise


# %% ======================================================================
# EXTRA FUNCTIONS
# ---------------------
#
# _generate_wavenumber_range
# _generate_broadening_range
#
# XXX =====================================================================


def _generate_wavenumber_range(wavenum_min, wavenum_max, wstep, broadening_max_width):
    """ define waverange vectors, with ``wavenumber`` the ouput spectral range
    and ``wavenumber_calc`` the spectral range used for calculation, that includes
    neighbour lines within ``broadening_max_width`` distance

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
    """ Generate array on which to compute line broadening

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

    print("Testing factory:", pytest.main(["../test/test_factory.py"]))
