.. _label_line_by_line:

=========================
Line-by-line (LBL) module
=========================

This is the core of RADIS: it calculates the spectral densities for a homogeneous
slab of gas, and returns a :py:class:`~radis.spectrum.spectrum.Spectrum` object. 
The detailed RADIS calculation flow chart is given in :ref:`Architecture <label_dev_architecture>`. 

Calculating spectra 
===================

Nonequilibrium Calculations
---------------------------

Nonequilibrium calculations require to know the vibrational and rotational energies of each 
level in order to calculate the nonequilibrium populations. 

You can either let RADIS calculate rovibrational energies
with its built-in :ref:`spectroscopic constants <label_db_spectroscopic_constants>`, 
or supply an energy level database. In the latter case, you need to edit the 
:ref:`Configuration file <label_lbl_config_file>` . 


The Spectrum Factory
--------------------

Most RADIS calculations can be done using the :py:func:`~radis.lbl.calc.calc_spectrum` function. 
Advanced examples require to use the :py:class:`~radis.lbl.factory.SpectrumFactory`
class, which is the core of RADIS line-by-line calculations. 
:py:func:`~radis.lbl.calc.calc_spectrum` is a wrapper to :py:class:`~radis.lbl.factory.SpectrumFactory`
for the simple cases. 

The :py:class:`~radis.lbl.factory.SpectrumFactory` allows you to :

- calculate multiple spectra (batch processing) with a same line database 
- edit the line database manually 
- have access to intermediary calculation variables
- connect to a database of precomputed spectra on your computer

To use the :py:class:`~radis.lbl.factory.SpectrumFactory`, first 
load your own line database with :py:meth:`~radis.lbl.loader.DatabankLoader.load_databank`, 
and then calculate several spectra in batch using 
:py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum` and 
:py:meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum`, 
and :py:mod:`~astropy.units` ::

    import astropy.units as u
    from radis import SpectrumFactory
    sf = SpectrumFactory(wavelength_min=4165 * u.nm, 
                         wavelength_max=4200 * u.nm,
                         path_length=0.1 * u.m,
                         pressure=20 * u.mbar,
                         molecule='CO2',
                         isotope='1,2', 
                         cutoff=1e-25,              # cm/molecule  
                         broadening_max_width=10,   # cm-1
                         )
    sf.load_databank('HITRAN-CO2-TEST')        # this database must be defined in ~/.radis
    s1 = sf.eq_spectrum(Tgas=300 * u.K)
    s2 = sf.eq_spectrum(Tgas=2000 * u.K)
    s3 = sf.non_eq_spectrum(Tvib=2000 * u.K, Trot=300 * u.K)

.. _label_lbl_config_file:

Configuration file
------------------

The ``~/.radis`` configuration file is used to store the list and attributes of the Line databases 
available on your computer. 
Without a configuration file, you are forced to download the corresponding [HITRAN-2016]_ line database 
files everytime, either through the :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank` method 
of :py:class:`~radis.lbl.factory.SpectrumFactory`, or the (default) ``databank='fetch'`` option in 
:py:func:`~radis.lbl.calc.calc_spectrum`

.. note::

    it is also possible to give :py:meth:`~radis.lbl.loader.DatabankLoader.load_databank` the line database path,
    format, and partition function format directly, but this is not recommended and should only be used if for some 
    reason you cannot create a configuration file. 

A ``~/.radis`` is user-dependant, and machine-dependant. It contains a list of database, everyone of which 
is specific to a given molecule. It typically looks like::

    # The CO2 HITEMP 2010 files have been retrieved from ftp://cfa-ftp.harvard.edu/pub/HITEMP-2010/
    # Partition function are that of CDSD-4000, retrived from ftp://ftp.iao.ru/pub/CDSD-4000
    [HITEMP-CO2]
    info = CDSD-HITEMP database, with energy levels calculated from Dunham expansions
    path = 
           PATH_TO\HITEMP-2010\cdsd_hitemp_07
           PATH_TO\HITEMP-2010\cdsd_hitemp_08
           PATH_TO\HITEMP-2010\cdsd_hitemp_09
    format = hitran
    parfunc =  PATH_TO\CDSD-4000\partition_functions.txt
    parfuncfmt = cdsd
    levelsfmt = radis

In the former example, RADIS built-in :ref:`spectroscopic constants <label_db_spectroscopic_constants>` 
are used to calculate the energy levels for CO2. 
It is also possible to use your own Energy level database. For instance::


    # List of databases
    [CDSD-HITEMP-HAMILTONIAN]
    info = CDSD-HITEMP database
    path = 
           PATH_TO\CDSD-HITEMP\cdsd_hitemp_07
           PATH_TO\CDSD-HITEMP\cdsd_hitemp_08
           PATH_TO\CDSD-HITEMP\cdsd_hitemp_09
    format = cdsd
    parfunc = D:\PATH_TO\CDSD-4000\partition_functions.txt
    parfuncfmt = cdsd
    levels_iso1 = D:\PATH_TO\CDSD-4000\626_PJCNn_TvibTrot.levels
    levels_iso2 = D:\PATH_TO\CDSD-4000\636_PJCNn_TvibTrot.levels
    levelsfmt = cdsd
    levelsZPE = 2531.828

The up-to-date format is given in :py:data:`~radis.misc.config.DBFORMAT`:

- ``path`` corresponds to Line databases (here: downloaded from [HITEMP-2010]_) and the ``levels_iso``
  are user generated Energy databases (here: calculated from the [CDSD-4000]_ Hamiltonian on non-distributed code,
  which takes into account non diagonal coupling terms). 

- ``format`` is the databank text file format. It can be one of ``'hitran'`` (for HITRAN / HITEMP 2010), 
  ``'cdsd-hitemp'`` and ``'cdsd-4000'`` for the different CDSD versions (for CO2 only). See full list in 
  :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`. 
  
- ``parfuncfmt``: ``cdsd``, ``hapi`` is the format of the tabulated partition functions used. 
  If ``'hapi'``, then [HAPI]_ is used to retrieve them (valid if your databank is HITRAN data). 
  See full list in :py:data:`~radis.lbl.loader.KNOWN_PARFUNCFORMAT` 
 
- ``parfunc`` is the path to the tabulated partition function to use in in equilibrium calculations 
  (:py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`). If ``parfuncfmt`` is ``'hapi'`` then `parfunc` should be
  the link to the hapi.py file. If not given, then the :py:mod:`~radis.io.hitran.hapi` embedded in RADIS 
  is used (check version)
  
- ``levels_iso#`` are the path to the energy levels to use for each isotope, which are needed for 
  nonequilibrium calculations (:py:meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum`).

- ``levelsfmt`` is the energy levels database format. Typically, ``'radis'``, and various implementation of [CDSD-4000]_ 
  nonequilibrium partitioning of vibrational and rotational energy: ``'cdsd-pc'``, ``'cdsd-pcN'``, ``'cdsd-hamil'``. 
  See full list in :py:data:`~radis.lbl.loader.KNOWN_LVLFORMAT`

  
A default ``~/.radis`` can be generated with :py:func:`~radis.test.utils.setup_test_line_databases`, which 
creates two test databases from fragments of [HITRAN-2016]_ line databases:: 

    from radis.test.utils import setup_test_line_databases
    setup_test_line_databases()
    
which will result in ::


    [HITRAN-CO2-TEST]
    info = HITRAN 2016 database, CO2, 1 main isotope (CO2-626), bandhead: 2380-2398 cm-1 (4165-4200 nm)
    path = PATH_TO\radis\radis\test\files\hitran_co2_626_bandhead_4165_4200nm.par
    format = hitran
    parfuncfmt = hapi
    levelsfmt = radis


    [HITRAN-CO-TEST]
    info = HITRAN 2016 database, CO, 3 main isotopes (CO-26, 36, 28), 2000-2300 cm-1
    path = PATH_TO\radis\radis\test\files\hitran_co_3iso_2000_2300cm.par
    format = hitran
    parfuncfmt = hapi
    levelsfmt = radis


    [HITEMP-CO2-TEST]
    info = HITEMP-2010, CO2, 3 main isotope (CO2-626, 636, 628), 2283.7-2285.1 cm-1
    path = PATH_TO\radis\radis\test\files\cdsd_hitemp_09_fragment.txt
    format = cdsd-hitemp
    parfuncfmt = hapi
    levelsfmt = radis


If you configuration file exists already, the test databases will simply be appended. 
These databases are used in some of the tests cases of RADIS, and the ``~/.radis`` may already contain 
them if you ever started the test suite with::

    cd radis 
    pytest 


Advanced
========

Calculation Flow Chart
----------------------

Refer to :ref:`Architecture <label_dev_architecture>` for an overview of how equilibrium
and nonequilibrium calculations are conducted. 

.. _label_lbl_custom_constants:

Use Custom Spectroscopic constants
----------------------------------

Spectroscopic constants are a property of the RADIS :py:class:`~radis.db.classes.ElectronicState` 
class. All molecules are stored in the :py:class:`~radis.db.molecules.Molecules` dictionary.
You need to update this dictionary before running your calculation in order to use your 
own spectroscopic constants. 

An example of how to use your own spectroscopic constants::

    from radis import calc_spectrum
    from radis.db.molecules import Molecules, ElectronicState

    Molecules['CO2'][1]['X'] = ElectronicState('CO2', isotope=1, state='X', term_symbol='1Σu+',
                                spectroscopic_constants='my_constants.json',  # <<< YOUR FILE HERE 
                                spectroscopic_constants_type='dunham',
                                Ediss=44600,
                                )
    s = calc_spectrum(...)



Vibrational bands
-----------------

To calculate vibrational bands of a given spectrum separately, use the  
:meth:`~radis.lbl.bands.BandFactory.eq_bands` and  :meth:`~radis.lbl.bands.BandFactory.non_eq_bands`
methods. See the :py:func:`~radis.test.lbl.test_bands.test_plot_all_CO2_bandheads` example in 
``radis/test/lbl/test_bands.py`` for more information. 


Connect to a Spectrum Database
------------------------------

In RADIS, the same code can be used to retrieve precomputed spectra if they exist, 
or calculate them and store them if they don't. See :ref:`Precompute Spectra <label_lbl_precompute_spectra>`



.. _label_lbl_performance:

Performance
===========

RADIS is very optimized, making use of C-compiled libraries (NumPy, Numba) for computationally intensive steps, 
and data analysis libraries (Pandas) to handle lines databases efficiently. 
Additionaly, different strategies and parameters are used to improve performances further:

Line Database Reduction Strategies
----------------------------------

By default:

- *linestrength cutoff* : lines with low linestrength are discarded after the new 
  populations are calculated. 
  Parameter: :py:attr:`~radis.lbl.loader.Input.cutoff` 
  (see the default value in the arguments of :py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`)

Additional strategies (deactivated by default):

- *weak lines* (pseudo-continuum): lines which are close to a much stronger line are called weak lines. 
  They are added to a pseudo-continuum and their lineshape is calculated with a simple 
  rectangular approximation.  
  See the default value in the arguments of :py:attr:`~radis.lbl.loader.Parameters.pseudo_continuum_threshold` 
  (see arguments of :py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`)


Lineshape optimizations
-----------------------

Lineshape convolution is usually the performance bottleneck in any line-by-line code. 

Two approaches can be used:

- improve the convolution efficiency. This involves using an efficient convolution algorithm,
  using a reduced convolution kernel, analytical approximations, or multiple spectral grid.
- reduce the number of convolutions (for a given number of lines): this is done using the DLM strategy. 

RADIS implements the two approaches as well as various strategies and parameters 
to calculate the lineshapes efficiently. 

- *broadening width* : lineshapes are calculated on a reduced spectral range. 
  Voigt computation calculation times scale linearly with that parameter. 
  Gaussian x Lorentzian calculation times scale as a square with that parameter. 
  parameters: broadening_max_width

- *Voigt approximation* : Voigt is calculated with an analytical approximation. 
  Parameter : :py:attr:`~radis.lbl.loader.Parameters.broadening_max_width` and 
  default values in the arguments of :py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`. 
  See :py:func:`~radis.lbl.broadening.voigt_lineshape`. 

- *Fortran precompiled* : previous Voigt analytical approximation is 
  precompiled in Fortran to improve performance times. This is always the 
  case and cannot be changed on the user side at the moment. See the source code
  of :py:func:`~radis.lbl.broadening.voigt_lineshape`. 
  
- *Multiple spectral grids* : many LBL codes use different spectral grids to 
  calculate the lineshape wings with a lower resolution. This strategy is not 
  implemented in RADIS. 

- *DLM* :  lines are projected on a Lineshape database to reduce the number of calculated 
  lineshapes from millions to a few dozens.
  With this optimization strategy, the lineshape convolution becomes almost instantaneous 
  and all the other strategies are rendered useless. Projection of all lines on the lineshape 
  database becomes the performance bottleneck.
  parameters: :py:attr:`~radis.lbl.loader.Parameters.dlm_res_L`, 
  :py:attr:`~radis.lbl.loader.Parameters.dlm_res_G`. 
  (this is the default strategy implemented in RADIS)

More details on the parameters below:

Computation parameters
----------------------

If performance is an issue (for instance when calculating polyatomic spectra on large spectral ranges), you 
may want to tweak the computation parameters in :py:func:`~radis.lbl.calc.calc_spectrum` and 
:py:class:`~radis.lbl.factory.SpectrumFactory`. In particular, the parameters that have the highest 
impact on the calculation performances are:

- The ``broadening_max_width``, which defines the spectral range over which the broadening is calculated. 
- The linestrength ``cutoff``, which defines which low intensity lines should be discarded. See 
  :meth:`~radis.lbl.base.BaseFactory.plot_linestrength_hist` to choose a correct cutoff. 
  
Check the [RADIS-2018]_ article for a quantitative assessment of the influence of the different parameters. 

Other strategies are possible, such as calculating the weak lines in a pseudo-continuum. This can 
result in orders of magnitude improvements in computation performances.:

- The ``pseudo_continuum_threshold`` defines which treshold should be used. 

See the :py:func:`~radis.test.lbl.test_broadening.test_abscoeff_continuum` case in ``radis/test/lbl/test_broadening.py`` 
for an example, which can be run with (you will need the CDSD-HITEMP database installed) ::

    pytest radis/test/lbl/test_broadening.py -m "test_abscoeff_continuum"


Database loading
----------------

Line database can be a performance bottleneck, especially for large polyatomic molecules in the [HITEMP-2010]_ 
or [CDSD-4000]_ databases. 
Line database files are automatically cached by RADIS under a ``.h5`` format after they are loaded the first time. 
If you want to deactivate this behaviour, use ``use_cached=False`` in :py:func:`~radis.lbl.calc.calc_spectrum`,
or ``db_use_cached=False, lvl_use_cached=False`` in :py:class:`~radis.lbl.factory.SpectrumFactory`.

If you are downloading the line database from [HITRAN-2016]_ with :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank` 
or the ``databank='fetch'`` option in :py:func:`~radis.lbl.calc.calc_spectrum`, then it is at the moment 
impossible to cache the database. 

You can also use :py:meth:`~radis.lbl.loader.DatabankLoader.init_databank` instead of the default 
:py:meth:`~radis.lbl.loader.DatabankLoader.load_databank`. The former will save the line database parameter,
and only load them if needed. This is useful if used in conjonction with 
:py:meth:`~radis.lbl.loader.DatabankLoader.init_database`, which will retrieve precomputed spectra from 
a database if they exist. 


Manipulate the database
-----------------------

If for any reason, you want to manipulate the line database manually (for instance, keeping only lines emitting 
by a particular level), you need to access the :py:attr:`~radis.lbl.loader.DatabankLoader.df0` attribute of 
:py:class:`~radis.lbl.factory.SpectrumFactory`. 

.. warning::

    never overwrite the ``df0`` attribute, else some metadata may be lost in the process. Only use inplace operations. 
    
For instance::

    sf = SpectrumFactory(
        wavenum_min= 2150.4,
        wavenum_max=2151.4,
        pressure=1,
        isotope=1)
    sf.load_databank('HITRAN-CO-TEST')
    sf.df0.drop(sf.df0[sf.df0.vu!=1].index, inplace=True)   # keep lines emitted by v'=1 only
    sf.eq_spectrum(Tgas=3000, name='vu=1').plot()

:py:attr:`~radis.lbl.loader.DatabankLoader.df0` contains the lines as they are loaded from the database. 
:py:attr:`~radis.lbl.loader.DatabankLoader.df1` is generated during the spectrum calculation, after the 
line database reduction steps, population calculation, and scaling of intensity and broadening parameters 
with the calculated conditions. 

Parallelization
---------------

Two parallelization are built-in RADIS. You can either run several :py:class:`~radis.lbl.factory.SpectrumFactory` 
in parallel. For that, just replace the :py:class:`~radis.lbl.factory.SpectrumFactory` with 
:py:class:`~radis.lbl.parallel.ParallelFactory` in your code, and use lists instead of single values 
for your input parameters. Example::

    from radis import SpectrumFactory
    sf = SpectrumFactory(...)
    sf.init_database(...)              # to store all spectra automatically
    for T in Tlist:
        s = sf.eq_spectrum(T)

Becomes::

    from radis import ParallelFactory
    sf = ParallelFactory(...)
    sf.init_database(...)              # to store all spectra automatically
    sf.eq_spectrum(Tlist)


Another parallelization is possible within one :py:class:`~radis.lbl.factory.SpectrumFactory` instance. 
In that case, the line database is split in different chuncks of lines that are processed independantly. 
See the ``parallel=`` parameter in :py:class:`~radis.lbl.factory.SpectrumFactory`. 

.. warning::
    Because LBL computations are usually more memory-heavy than CPU-heavy, you may not get 
    a lot of improvement by using parallelization. Ensure that your test works. 
    
Parallelized code can be tested against the linear code in `radis/test/lbl/test_parallel.py`, which can be run 
with::

    pytest radis/test/lbl/test_parallel.py 

Profiler
--------

You may want to track where the calculation is taking some time. 
You can set ``verbose=2`` to print the time spent on different operations. Example::

    s = calc_spectrum(1900, 2300,         # cm-1
                      molecule='CO',
                      isotope='1,2,3',
                      pressure=1.01325,   # bar
                      Tvib=1000,          # K
                      Trot=300,           # K
                      mole_fraction=0.1,
                      verbose=2,
                      )

::

    >>> ...
    >>> Fetching vib / rot energies for all 749 transitions
    >>> Fetched energies in 0s
    >>> Calculate weighted transition moment
    >>> Calculated weighted transition moment in 0.0
    >>> Calculating nonequilibrium populations
    >>> sorting lines by vibrational bands
    >>> lines sorted in 0.0s
    >>> Calculated nonequilibrium populations in 0.1s
    >>> scale nonequilibrium linestrength
    >>> scaled nonequilibrium linestrength in 0.0s
    >>> calculated emission integral
    >>> calculated emission integral in 0.0s
    >>> Applying linestrength cutoff
    >>> Applied linestrength cutoff in 0.0s (expected time saved ~ 0.0s)
    >>> Calculating lineshift
    >>> Calculated lineshift in 0.0s
    >>> Calculate broadening FWHM
    >>> Calculated broadening FWHM in 0.0s
    >>> Calculating line broadening (695 lines: expect ~ 0.1s on 1 CPU)
    >>> Calculated line broadening in 0.1s
    >>> process done in 0.4s
    >>> ... 

.. _label_lbl_precompute_spectra:

Precompute Spectra
------------------

See :py:meth:`~radis.lbl.loader.DatabankLoader.init_database`, which is the direct integration 
of :py:class:`~radis.tools.database.SpecDatabase` in a :py:class:`~radis.lbl.factory.SpectrumFactory` 