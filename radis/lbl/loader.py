# -*- coding: utf-8 -*-
""" 

Summary
-------

Module to host the databank loading / database initialisation parts of
SpectrumFactory (and unload the factory.py file). Basically it holds all of the
non-physical machinery, while actual population calculations and line broadening
are still calculated in factory.py

This is done through SpectrumFactory inheritance of the DatabankLoader class
defined here


Routine Listings
----------------

PUBLIC METHODS

- :py:meth:`radis.lbl.loader.DatabankLoader.load_databank`           >>> load line database
- :py:meth:`radis.lbl.loader.DatabankLoader.init_databank`           >>> load loader 
- :py:meth:`radis.lbl.loader.DatabankLoader.fetch_databank`           >>> fetch from HITRAN online
- :py:meth:`radis.lbl.loader.DatabankLoader.init_database`           >>> to interact / generate a SpectrumDatabase
- :py:meth:`radis.lbl.loader.DatabankLoader.get_conditions`
- :py:meth:`radis.lbl.loader.DatabankLoader.get_partition_function_interpolator`
- :py:meth:`radis.lbl.loader.DatabankLoader.get_partition_function_calculator`

PRIVATE METHODS - DATABASE LOADING

- :py:meth:`radis.lbl.loader.DatabankLoader._load_databank`
- :py:meth:`radis.lbl.loader.DatabankLoader._reload_databank`
- :py:meth:`radis.lbl.loader.DatabankLoader._check_line_databank`
- :py:meth:`radis.lbl.loader.DatabankLoader._retrieve_from_database`
- :py:meth:`radis.lbl.loader.DatabankLoader._build_partition_function_interpolator`
- :py:meth:`radis.lbl.loader.DatabankLoader._build_partition_function_calculator`
- :py:meth:`radis.lbl.loader.DatabankLoader._fetch_molecular_parameters`
- :py:meth:`radis.lbl.loader.DatabankLoader._get_temp_file`
- :py:meth:`radis.lbl.loader.DatabankLoader._clean_temp_file`

Most methods are written in inherited class with the following inheritance scheme:
    
:py:class:`~radis.lbl.loader.DatabankLoader` > :py:class:`~radis.lbl.base.BaseFactory` > 
:py:class:`~radis.lbl.broadening.BroadenFactory` > :py:class:`~radis.lbl.bands.BandFactory` > 
:py:class:`~radis.lbl.factory.SpectrumFactory` > :py:class:`~radis.lbl.parallel.ParallelFactory`

.. inheritance-diagram:: radis.lbl.parallel.ParallelFactory
   :parts: 1

Notes
-----

RADIS includes automatic rebuilding of Deprecated cache files + a global variable
to force regenerating them after a given version. See :py:data:`radis.OLDEST_COMPATIBLE_VERSION`

-------------------------------------------------------------------------------

"""
# TODO: on use_cache functions, make a 'clean' / 'reset' option to delete / regenerate
# cache files

# @dev: (on Spyder IDE navigate between sections easily as # XXX makes a reference
# (on the slide bar on the right)

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.db.molparam import MolParams
from radis.io.cdsd import cdsd2df
from radis.io.hitran import (
    hit2df,
    get_molecule,
    parse_global_quanta,
    parse_local_quanta,
)

# from radis.io.hitran import hit2dfTAB
from radis.misc.warning import EmptyDatabaseError
from radis.io.query import fetch_astroquery
from radis.io.tools import drop_object_format_columns, replace_PQR_with_m101
from radis.db.molecules import getMolecule
from radis.levels.partfunc import (
    PartFuncHAPI,
    PartFunc_Dunham,
    RovibParFuncTabulator,
    RovibParFuncCalculator,
)
from radis.levels.partfunc_cdsd import PartFuncCO2_CDSDtab, PartFuncCO2_CDSDcalc
from radis.tools.database import SpecDatabase
from radis.misc.config import getDatabankEntries, printDatabankEntries, getDatabankList
from radis.misc.utils import FileNotFoundError
from radis.misc.basics import compare_dict, compare_lists
from radis.misc.arrays import count_nans
from radis.misc.debug import printdbg
from radis.misc.printer import printg, printr
from radis.misc.log import printwarn
from radis.phys.convert import cm2nm
import os
from os.path import exists, abspath
import pandas as pd
from six import string_types
from radis.misc.warning import warn, default_warning_status
from warnings import catch_warnings, filterwarnings
import numpy as np
from time import time
import gc
from uuid import uuid1
from six.moves import range
import fnmatch
from radis.misc.utils import get_files_from_regex

KNOWN_DBFORMAT = ["hitran", "cdsd-hitemp", "cdsd-4000"]
"""list: Known formats for Line Databases:

- ``'hitran'`` : for HITRAN and HITEMP-2010 
- ``'cdsd-hitemp'`` : CDSD-HITEMP (CO2 only, same lines as HITEMP-2010)
- ``'cdsd-4000'`` : CDSD-4000 (CO2 only)

To install all databases manually see the :ref:`Configuration file <label_lbl_config_file>`
and the :ref:`list of databases <_label_line_databases>` .

See Also
--------

:ref:`Configuration file <label_lbl_config_file>`
"""

KNOWN_LVLFORMAT = ["radis", "cdsd-pc", "cdsd-pcN", "cdsd-hamil", None]
"""list: Known formats for Energy Level Databases (used in non-equilibrium calculations):

- ``'radis'``: energies calculated with Dunham expansions by 
    :class:`~radis.levels.partfunc.PartFunc_Dunham`
- ``'cdsd-pc'``: energies read from precomputed CDSD energies for CO2, with
    ``viblvl=(p,c)`` convention. See :class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDcalc`
- ``'cdsd-pcN'``: energies read from precomputed CDSD energies for CO2, with
    ``viblvl=(p,c,N)`` convention. See :class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDcalc`
- ``'cdsd-hamil'``: energies read from precomputed CDSD energies for CO2, with
    ``viblvl=(p,c,J,N)`` convention, i.e., a each rovibrational level can have a 
    unique vibrational energy (this is needed when taking account Coupling terms)
    See :class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDcalc`
- ``None``: means you can only do Equilibrium calculations.

See Also
--------

:ref:`Configuration file <label_lbl_config_file>`
 """

KNOWN_PARFUNCFORMAT = ["cdsd", "hapi"]
"""list: Known formats for partition function (tabulated files to read), or 'hapi'
to fetch Partition Functions using HITRAN Python interface instead of reading
a tabulated file.

See Also
--------

:ref:`Configuration file <label_lbl_config_file>`
"""

drop_auto_columns_for_dbformat = {
    "hitran": ["ierr", "iref", "lmix", "gp", "gpp"],
    "cdsd-4000": ["wang2"],
    "cdsd-hitemp": ["wang2", "lsrc"],
}
""" dict: drop these columns if using ``drop_columns='auto'`` in load_databank 
Based on the value of ``dbformat=``, some of these columns won't be used.

See Also
--------

- 'hitran': (HITRAN / HITEMP) :data:`~radis.io.hitran.columns_2004`, 
- 'cdsd-hitemp' (CDSD HITEMP): :data:`~radis.io.cdsd.columns_hitemp`, 
- 'cdsd-4000': (CDSD 4000) :data:`~radis.io.cdsd.columns_4000`, 

"""
drop_auto_columns_for_levelsfmt = {
    "radis": [],
    "cdsd-pc": ["v1u", "v2u", "l2u", "v3u", "ru", "v1l", "v2l", "l2l", "v3l", "rl"],
    "cdsd-pcN": ["v1u", "v2u", "l2u", "v3u", "ru", "v1l", "v2l", "l2l", "v3l", "rl"],
    "cdsd-hamil": ["v1u", "v2u", "l2u", "v3u", "ru", "v1l", "v2l", "l2l", "v3l", "rl"],
    None: [],
}
""" dict: drop these columns if using ``drop_columns='auto'`` in load_databank 
Based on the value of ``lvlformat=``, some of these columns won't be used.

See Also
--------

- 'radis': :data:`~radis.io.hitran.columns_2004`, 
- 'cdsd-pc': :data:`~radis.io.hitran.columns_2004`, 
- 'cdsd-pcN' (CDSD-HITEMP): :data:`~radis.io.cdsd.columns_hitemp`, 
- 'cdsd-hamil': :data:`~radis.io.cdsd.columns_4000`, 

"""
drop_all_but_these = [
    "id",
    "iso",
    "wav",
    "int",
    "airbrd",
    "selbrd",
    "Tdpair",
    "Tdpsel",
    "Pshft",
    "El",
]
""" dict: drop all columns but these if using ``drop_columns='all'`` in load_databank 

Note: nonequilibrium calculations wont be possible anymore and it wont be possible
to identify lines with :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey`
    
See Also
--------

- 'hitran': (HITRAN / HITEMP) :data:`~radis.io.hitran.columns_2004`, 
- 'cdsd-hitemp' (CDSD HITEMP): :data:`~radis.io.cdsd.columns_hitemp`, 
- 'cdsd-4000': (CDSD 4000) :data:`~radis.io.cdsd.columns_4000`, 

"""

# Sanity checks
# (make sure all variables are defined everywhere)
assert compare_lists(drop_auto_columns_for_dbformat, KNOWN_DBFORMAT) == 1
assert compare_lists(drop_auto_columns_for_levelsfmt, KNOWN_LVLFORMAT) == 1

# %% Main class

from copy import deepcopy


class ConditionDict(dict):
    """ A class to hold Spectrum calculation input conditions, or computation 
    parameters. Works like a dict except you can also access attribute with::

        v = a.key 

    Also can be copied, deepcopied, and parallelized in multiprocessing
    
    Notes
    -----
    
    for developers: 
        
    Parameters and Input could also have simply derived from the (object) class, 
    but it may have missed some convenients functions implemented for dict. 
    For instance, how to be picked / unpickled. 
    
    See Also
    --------
    
    :class:`~radis.lbl.loader.Input`, 
    :class:`~radis.lbl.loader.Parameter`, 
    """

    def get_params(self):
        """ Returns the variables (and their values) contained in the dictionary, 
        minus some based on their type. Numpy array, dictionaries and pandas DataFrame
        are removed
        
        Tuples are converted to string"""

        # Filter parameters based on type
        def filter_type(params):
            filt_params = {}
            for k, v in params.items():
                # Ignore some
                if type(v) in [np.ndarray, dict, pd.DataFrame] or v is None:
                    continue
                if isinstance(k, string_types):
                    # Also discard all starting with '_'
                    if k.startswith("_"):
                        continue
                # Convert some
                if isinstance(v, tuple):
                    filt_params[k] = ",".join([str(vi) for vi in v])
                # ... or return raw
                else:
                    filt_params[k] = v
            return filt_params

        return filter_type(self)

    # Methods to access keys as attributes
    def __getattr__(self, attr):
        if attr == "__getstate__":
            return dict.__getattr__(attr)
        return self[attr]

    def __setattr__(self, attr, value):
        self[attr] = value

    # Methods needed for Multiprocessing
    def __getstate__(self):
        return dict(self.__dict__)

    def __setstate__(self, state):
        self.__dict__ = state

    def copy(self):
        obj_copy = ConditionDict()
        for k, v in self.items():
            obj_copy[k] = v
        return obj_copy

    def __deepcopy__(self, memo={}):
        obj_copy = ConditionDict()
        for k, v in self.items():
            obj_copy[k] = deepcopy(v)
        return obj_copy


# Below are 3 classes that inherit from ConditionDict:
# - Input: physical input, stored in Spectra
# - Parameters: computational input, stored in Spectra
# - MiscParams: other, purely descriptive parameters, stored in Spectra, but not used
#               when comparing Spectra to precomputed one in SpecDatabases

# class Input(object):
class Input(ConditionDict):
    """ A class to hold Spectrum calculation input conditions. 
    Works like a dict except you can also access attribute with::

        v = a.key 

    Also can be copied, deepcopied, and parallelized in multiprocessing
    """

    #    # hardcode attribute names, to prevent typos and the declaration of unwanted parameters
    #    __slots__ = [
    #         'Tgas', 'Tref', 'Tvib', 'Trot', 'isotope', 'medium', 'mole_fraction',
    #         'molecule', 'overpopulation', 'path_length', 'pressure_mbar', 'rot_distribution',
    #         'self_absorption', 'state', 'vib_distribution', 'wavelength_max',
    #         'wavelength_min', 'wavenum_max', 'wavenum_min']

    def __init__(self):
        super(Input, self).__init__()

        self.Tgas = None  #: float: gas (translational) temperature. Overwritten by SpectrumFactory.eq/noneq_spectrum
        self.Tref = None  #: float: reference temperature for line database.
        self.Tvib = None  #: float: vibrational temperature. Overwritten by SpectrumFactory.eq/noneq_spectrum
        self.Trot = None  #: float: rotational temperature. Overwritten by SpectrumFactory.eq/noneq_spectrum
        self.isotope = None  #: str: isotope list. Can be '1,2,3', etc. or 'all'
        self.mole_fraction = None  #: float: mole fraction
        self.molecule = None  #: str: molecule
        self.overpopulation = None  #: dict: overpopulation
        self.path_length = None  #: float: path length (cm)
        self.pressure_mbar = 1013.25  #: float: pressure (mbar)
        self.rot_distribution = "boltzmann"  #: str: rotational levels distribution
        self.self_absorption = (
            True  #: bool: self absorption (if True, not optically thin)
        )
        self.state = None  #: str: electronic state
        self.vib_distribution = (
            "boltzmann"  #: str: vibrational levels distribution (boltzmann, treanor)
        )
        self.wavenum_max = None  #: str: wavenumber max (cm-1)
        self.wavenum_min = None  #: str: wavenumber min (cm-1)


# class Parameters(object):
class Parameters(ConditionDict):
    """ A class to hold Spectrum calculation computation parameters. Works like 
    a dict except you can also access attribute with::

        v = a.key 

    Also can be copied, deepcopied, and parallelized in multiprocessing
    """

    #    # hardcode attribute names, to prevent typos and the declaration of unwanted parameters
    #    __slots__ = [
    #                 'broadening_max_width', 'chunksize', 'cutoff',
    #                 'db_assumed_sorted', 'db_use_cached', 'dbformat', 'dbpath',
    #                 'export_lines', 'export_populations', 'levelsfmt', 'lvl_use_cached',
    #                 'Ngroups', 'Nprocs', 'parallel', 'parfuncfmt', 'parfuncpath',
    #                 'pseudo_continuum_threshold', 'warning_broadening_threshold',
    #                 'warning_linestrength_cutoff', 'wavenum_max_calc', 'wavenum_min_calc',
    #                 'waveunit', 'wstep']

    def __init__(self):
        super(Parameters, self).__init__()

        # Dev: Init here to be found by autocomplete
        self.broadening_max_width = None  #: float: cutoff for lineshape calculation (cm-1). Overwritten by SpectrumFactory
        self.cutoff = None  #: float: linestrength cutoff (molecule/cm)
        self.db_assumed_sorted = None  #: bool: assume that Line Database is sorted (helps not to parse the whole database)
        self.db_use_cached = (
            None  #: bool: use (and generate) cache files for Line Database
        )
        self.dbformat = None  #: str: format of Line Database. See :data:`~radis.lbl.loader.KNOWN_DBFORMAT`
        self.dbpath = None  #: list: list of filepaths to Line Database
        self.levelsfmt = None  #: str: format of Energy Database. See :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`
        self.lvl_use_cached = (
            None  #: bool: use (and generate) cache files for Energy Database
        )
        self.parfuncfmt = None  #: str: format of tabulated Partition Functions. See #: str: format of Energy Database. See :data:`~radis.lbl.loader.KNOWN_PARFUNCFORMAT`
        self.parfuncpath = None  #: str: filepath to tabulated Partition Functions
        self.pseudo_continuum_threshold = 0  #: float: threshold to assign lines in pseudo continuum. Overwritten in SpectrumFactory
        self.wavenum_max_calc = None  #: float: maximum calculated wavenumber (cm-1) initialized by SpectrumFactory
        self.wavenum_min_calc = None  #: float: minimum calculated wavenumber (cm-1) initialized by SpectrumFactory
        self.waveunit = "cm-1"  #: waverange unit: should be cm-1.
        self.wstep = None  #: float: spectral resolution (cm-1)
        self.dlm_res_L = 0.01  #: float (cm-1): Lorentzian step for DLM lineshape database. Default 0.01 cm-1
        self.dlm_res_G = 0.01  #: float (cm-1): DLM Gaussian step DLM lineshape database. Default 0.01 cm-1
        self.include_neighbouring_lines = True
        """bool: if ``True``, includes the contribution of off-range, neighbouring 
        lines because of lineshape broadening. Default ``True``."""


class MiscParams(ConditionDict):
    """ A class to hold Spectrum calculation descriptive parameters. Unlike 
    :class:`~radis.lbl.loader.Parameters`, these parameters cannot influence the 
    Spectrum output and will not be used when comparing Spectrum with existing, 
    precomputed spectra in :class:`~radis.tools.database.SpecDatabase`
    
    Works like 
    a dict except you can also access attribute with::

        v = a.key 

    Also can be copied, deepcopied, and parallelized in multiprocessing
    """

    def __init__(self):
        super(MiscParams, self).__init__()

        # Dev: Init here to be found by autocomplete
        self.chunksize = None  #: int: divide line database in chunks of lines
        self.export_lines = (
            None  #: bool: export lines in output Spectrum (takes memory!)
        )
        self.export_populations = (
            None  #: bool: export populations in output Spectrum (takes memory!)
        )
        self.Ngroups = None  #: int:
        self.Nprocs = None  #: int:
        self.parallel = None  #: bool: use parallel processing for line broadening step
        self.warning_broadening_threshold = (
            None  #: float: [0-1] raise a warning if the lineshape area is different
        )
        self.warning_linestrength_cutoff = None  #: float [0-1]: raise a warning if the sum of linestrength cut is above that


def format_paths(s):
    """ escape all special characters """
    if s is not None:
        s = s.replace("\\", "/")
    return s


TEMP_FILE_PREFIX = ".radis_"

df_metadata = ["Ia", "molar_mass", "Qref", "Qvib", "Q"]
""" list: metadata of line DataFrames :py:attr:`~radis.lbl.loader.DatabankLoader.df0`, 
:py:attr:`~radis.lbl.loader.DatabankLoader.df1`.
@dev: when having only 1 molecule, 1 isotope, these parameters are 
constant for all rovibrational lines. Thus, it's faster and much more
memory efficient to transport them as attributes of the DataFrame 
rather than columns. The syntax is the same, thus the operations do 
not change, i.e::
    
    k_b / df.molar_mass
    
will work whether molar_mass is a float or a column. 

.. warning:: 
    
    However, in the current Pandas implementation of :py:class:`~pandas.DataFrame`,
    attributes are lost whenever the DataFrame is recreated, transposed,
    pickled.
    
Thus, we use :py:func:`~radis.misc.basics.transfer_metadata` to keep
the attributes after an operation, and :py:func:`~radis.misc.basics.expand_metadata` 
to make them columns before a Serializing operation (ex: multiprocessing)
@dev: all of that is a high-end optimization. Users should not deal 
with internal DataFrames. 

References
----------
https://stackoverflow.com/q/13250499/5622825    
"""


class DatabankLoader(object):
    """

    .. inheritance-diagram:: radis.lbl.factory.SpectrumFactory
       :parts: 1
    
    See Also
    --------

    :class:`~radis.lbl.factory.SpectrumFactory`
    """

    def __init__(self):

        # Name parameters
        # ... defaults values are overwritten by SpectrumFactory input
        # ... values here are just to help autocompletion tools

        self.verbose = True  #: bool, or int: increase verbose level. 0, 1, 2 supported at the moment
        self.save_memory = False  #: bool: if True, tries to save RAM memory (but may take a little for time, saving stuff to files instead of RAM for instance)
        self.parsum_tab = (
            {}
        )  #: dict: store all partition function tabulators, per isotope
        self.parsum_calc = (
            {}
        )  #: dict: store all partition function calculators, per isotope

        # in particular input conditions:
        # ... means all values that can have an impact on the calculated Spectrum
        # ... They will be stored in the Spectrum object
        self.input = Input()
        self.input.isotope = "all"
        self.input.molecule = ""
        self.input.state = ""

        # an computation parameters:
        self.params = Parameters()  # params that can change output (ex: threshold)
        self.misc = MiscParams()  # params that cant (ex: number of CPU, etc.)

        # Setup individual warnings. Value of keys can be:
        # - 'warning' (default: just trigger a warning)
        # - 'error' (raises an error on this warning)
        # - 'ignore'  (do nothing)
        # The key self.warnings['default'] will set the warning behavior for all
        # other warnings
        self.warnings = default_warning_status
        """ dict: Default warnings for SpectrumFactory. See 
        :py:data:`~radis.misc.warnings.default_warning_status`"""

        # Generate unique id for Factory
        self._id = uuid1()

        # Just look up if temp files already exist
        tempfiles = [f for f in os.listdir(".") if f.startswith(TEMP_FILE_PREFIX)]
        if len(tempfiles) > 0:
            self.warn(
                "tempfile already exists: {0} in {1}. ".format(tempfiles, abspath("."))
                + "Consider cleaning if you are not running processes in parallel",
                "default",
            )

        # Init Annotations (Python >3.6) [hints for users]
        try:
            self.load_databank.__annotations__["format"] = KNOWN_DBFORMAT
            self.load_databank.__annotations__["levelsfmt"] = KNOWN_LVLFORMAT
            self.load_databank.__annotations__["parfuncfmt"] = KNOWN_PARFUNCFORMAT
        except:  # old Python version
            pass

        # Variables that will hold the dataframes.
        self.df0 = None
        """pandas DataFrame : initial line database after loading.
                
        If for any reason, you want to manipulate the line database manually (for instance, keeping only lines emitting 
        by a particular level), you need to access the :py:attr:`~radis.lbl.loader.DatabankLoader.df0` attribute of 
        :py:class:`~radis.lbl.factory.SpectrumFactory`. 
        
        .. warning::
        
            never overwrite the ``df0`` attribute, else some metadata may be lost in the process. 
            Only use inplace operations. If reducing the number of lines, add
            a df0.reset_index()
            
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

        See Also
        --------
        
        :py:attr:`~self.radis.lbl.loader.DatabankLoader.df1`
        
        """
        self.df1 = None
        """DataFrame : line database, scaled with populations + linestrength cutoff
        Never edit manually. See all comments about :py:attr:`~self.radis.lbl.loader.DatabankLoader.df0`
        
        See Also
        --------
        
        :py:attr:`~self.radis.lbl.loader.DatabankLoader.df0`
        
        """

        # Temp variable to store databanks information
        self._databank_args = []
        self._databank_kwargs = {}

        self._autoretrieveignoreconditions = []  # HACK. See _retrieve_from_database

    # %% ======================================================================
    # PUBLIC METHODS
    # ------------------------
    # init_databank           >>> init line database (= store reference but load() later)
    # load_databank           >>> load line database
    # init_database           >>> to interact / generate a SpectrumDatabase
    #
    # =========================================================================

    def init_databank(self, *args, **kwargs):
        """ Method to init databank parameters but only load them when needed.
        Databank is reloaded by :meth:`~radis.lbl.loader.DatabankLoader._check_line_databank`

        Same inputs Parameters as :meth:`~radis.lbl.loader.DatabankLoader.load_databank`:


        Parameters
        ----------

        name: a section name specified in your ``~/.radis``
            ``.radis`` has to be created in your HOME (Unix) / User (Windows). If
            not ``None``, all other arguments are discarded.
            Note that all files in database will be loaded and it may takes some
            time. Better limit the database size if you already know what
            range you need. See :ref:`Configuration file <label_lbl_config_file>` and 
            :data:`~radis.misc.config.DBFORMAT` for expected 
            ``~/.radis`` format


        Other Parameters
        ----------------

        path: str, list of str, None
            list of database files, or name of a predefined database in the 
            :ref:`Configuration file <label_lbl_config_file>` (`~/.radis`)
            Accepts wildcards ``*`` to select multiple files
            
        format: ``'hitran'``, ``'cdsd-hitemp'``, ``'cdsd-4000'``, or any of :data:`~radis.lblinit_databank.loader.KNOWN_DBFORMAT`
            database type. ``'hitran'`` for HITRAN/HITEMP, ``'cdsd-hitemp'`` 
            and ``'cdsd-4000'`` for the different CDSD versions. Default ``'hitran'``

        parfuncfmt: ``'hapi'``, ``'cdsd'``, or any of :data:`~radis.lbl.loader.KNOWN_PARFUNCFORMAT`
            format to read tabulated partition function file. If ``hapi``, then
            HAPI (HITRAN Python interface) [1]_ is used to retrieve them (valid if
            your database is HITRAN data). HAPI is embedded into RADIS. Check the
            version.

        parfunc: filename or None
            path to tabulated partition function to use.
            If `parfuncfmt` is `hapi` then `parfunc` should be the link to the
            hapi.py file. If not given, then the hapi.py embedded in RADIS is used (check version)

        levels: dict of str or None
            path to energy levels (needed for non-eq calculations). Format:
            {1:path_to_levels_iso_1, 3:path_to_levels_iso3}. Default ``None``

        levelsfmt: 'cdsd-pc', 'radis' (or any of :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`) or ``None``
            how to read the previous file. Known formats: (see :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`).
            If ``radis``, energies are calculated using the diatomic constants in radis.db database
            if available for given molecule. Look up references there.
            If None, non equilibrium calculations are not possible. Default ``None``.

        db_use_cached: boolean, or ``None``
            if ``True``, a pandas-readable csv file is generated on first access,
            and later used. This saves on the datatype cast and conversion and
            improves performances a lot. But! ... be sure to delete these files
            to regenerate them if you happen to change the database. If ``'regen'``,
            existing cached files are removed and regenerated.
            It is also used to load energy levels from ``.h5`` cache file if exist.
            If ``None``, the value given on Factory creation is used. Default ``None``

        db_assumed_sorted: boolean
            load_databank first reads the first line and check it's relevant.
            This improves database loading times if not all files are required,
            but it assumes database files are sorted in wavenumber!
            Default ``True``
            
        load_energies: boolean
            if ``False``, dont load energy levels. This means that nonequilibrium
            spectra cannot be calculated, but it saves some memory. Default ``True``
    
        include_neighbouring_lines: bool
            ``True``, includes off-range, neighbouring lines that contribute
            because of lineshape broadening. The ``broadening_max_width`` 
            parameter is used to determine the limit. Default ``True``.

        Other arguments are related to how to open the files

        buffer: ``'RAM'``, ``'h5'``, ``'direct'``
            Different modes for loading up database: either directly in 'RAM' mode,
            or in 'h5' mode.

            - 'RAM': is faster but memory hunger
            - 'h5': handles better a bigger database (> 1M lines): slower (up to 3x), but less
              risks of MemoryErrors 
            - 'direct': file is read directly from a single h5 file under key 'df'
              Fastest of all, doesnt check the database validity or format. Use only
              if you have a single, already formatted database file (used by Factory
              when reloading database)
              
            Default ``'RAM'``

        drop_columns: list
            columns names to drop from Line DataFrame after loading the file. 
            Not recommended to use, unless you explicitely want to drop information 
            (for instance if dealing with too large databases). If ``[]``, nothing 
            is dropped. If ``'auto'``, parameters considered unnecessary 
            are dropped. See :data:`~radis.lbl.loader.drop_auto_columns_for_dbformat`
            and :data:`~radis.lbl.loader.drop_auto_columns_for_levelsfmt`. 
            Default ``'auto'``.


        Notes
        -----

        Useful in conjonction with :meth:`~radis.lbl.loader.DatabankLoader.init_database` 
        when dealing with large line databanks when some of the spectra may have 
        been precomputed in a spectrum database (:class:`~radis.tools.database.SpecDatabase`)

        Note that any previously loaded databank is discarded on the method call


        See Also
        --------

        - Download from HITRAN: :meth:`~radis.lbl.loader.DatabankLoader.fetch_databank`
        - Load from local files: :meth:`~radis.lbl.loader.DatabankLoader.load_databank`
        - Reload databank: :meth:`~radis.lbl.loader.DatabankLoader._check_line_databank`

        """

        # Check inputs
        (
            name,
            path,
            dbformat,
            parfunc,
            parfuncfmt,
            levels,
            levelsfmt,
            db_use_cached,
            db_assumed_sorted,
            drop_columns,
            buffer,
            load_energies,
            include_neighbouring_lines,
        ) = self._check_database_params(*args, **kwargs)

        # Store arguments for later. The database will only be loaded if a Spectrum
        # has to be calculated. See Factory.eq_spectrum() and non_eq_spectrum()
        self._databank_args = args
        self._databank_kwargs = kwargs

        # Store the values in factory (for export and comparison with spectrum
        # stored in SpecDatabase)
        self._store_database_params(
            name=name,
            path=path,
            format=dbformat,
            parfunc=parfunc,
            parfuncfmt=parfuncfmt,
            levels=levels,
            levelsfmt=levelsfmt,
            db_use_cached=db_use_cached,
            db_assumed_sorted=db_assumed_sorted,
            load_energies=load_energies,
            include_neighbouring_lines=include_neighbouring_lines,
        )

        # Delete database
        self.df0 = None

    def fetch_databank(
        self,
        source="astroquery",
        format="hitran",
        parfunc=None,
        parfuncfmt="hapi",
        levels=None,
        levelsfmt="radis",
        load_energies=True,
        include_neighbouring_lines=True,
        drop_non_numeric=True,
    ):
        """ Fetch databank with Astroquery [1]_

        Parameters
        ----------

        source: ``'astroquery'``
            where to download database from

        format: ``'hitran'``, ``'cdsd-hitemp'``, ``'cdsd-4000'``, or any of :data:`~radis.lbl.loader.KNOWN_DBFORMAT`
            database type. ``'hitran'`` for HITRAN/HITEMP, ``'cdsd-hitemp'`` 
            and ``'cdsd-4000'`` for the different CDSD versions. Default 'hitran'

        parfuncfmt: ``'cdsd'``, ``'hapi'``, or any of :data:`~radis.lbl.loader.KNOWN_PARFUNCFORMAT`
            format to read tabulated partition function file. If ``hapi``, then
            HAPI (HITRAN Python interface) [2]_ is used to retrieve them (valid if
            your database is HITRAN data). HAPI is embedded into RADIS. Check the
            version.

        parfunc: filename or None
            path to tabulated partition function to use.
            If `parfuncfmt` is `hapi` then `parfunc` should be the link to the
            hapi.py file. If not given, then the hapi.py embedded in RADIS is used (check version)

        levels: dict of str or None
            path to energy levels (needed for non-eq calculations). Format:
            {1:path_to_levels_iso_1, 3:path_to_levels_iso3}. Default ``None``

        levelsfmt: 'cdsd-pc', 'radis' (or any of :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`) or ``None``
            how to read the previous file. Known formats: (see :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`).
            If ``radis``, energies are calculated using the diatomic constants in radis.db database
            if available for given molecule. Look up references there.
            If None, non equilibrium calculations are not possible. Default ``None``.

        load_energies: boolean
            if ``False``, dont load energy levels. This means that nonequilibrium
            spectra cannot be calculated, but it saves some memory. Default ``True``
    
        include_neighbouring_lines: bool
            ``True``, includes off-range, neighbouring lines that contribute
            because of lineshape broadening. The ``broadening_max_width`` 
            parameter is used to determine the limit. Default ``True``.
            
        Other Parameters
        ----------------
        
        drop_non_numeric: boolean
            if ``True``, non numeric columns are dropped. This improves performances, 
            but make sure all the columns you need are converted to numeric formats 
            before hand. Default ``True``. Note that if a cache file is loaded it 
            will be left untouched.

        See Also
        --------

        - Load from local files: :meth:`~radis.lbl.loader.DatabankLoader.load_databank`
        - Load when needed: :meth:`~radis.lbl.loader.DatabankLoader.init_databank`


        References
        ----------

        .. [1] `Astroquery <https://astroquery.readthedocs.io>`_

        .. [2] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_

        """
        # @dev TODO: also add cache file to fetch_databank, similar to load_databank
        # | Should store the waverange, molecule and isotopes in the cache file
        # | metadata to ensures that it is redownloaded if necessary.
        # | see implementation in load_databank.

        # Check inputs
        if source not in ["astroquery"]:
            raise NotImplementedError("source: {0}".format(source))
        if source == "astroquery":
            assert format == "hitran"

        # Get inputs
        dbformat = format

        molecule = self.input.molecule
        if not molecule:
            raise ValueError("Please define `molecule=` so the database can be fetched")

        isotope = self.input.isotope
        if isotope == "all":
            raise ValueError(
                "Please define isotope explicitely (cannot use 'all' with fetch_databank)"
            )
        isotope_list = self._get_isotope_list()

        if include_neighbouring_lines:
            wavenum_min = self.params.wavenum_min_calc
            wavenum_max = self.params.wavenum_max_calc
        else:
            wavenum_min = self.input.wavenum_min
            wavenum_max = self.input.wavenum_max

        # %% Init Line database
        # ---------------------

        frames = []  # lines for all isotopes
        for iso in isotope_list:
            df = fetch_astroquery(
                molecule, iso, wavenum_min, wavenum_max, verbose=self.verbose
            )
            frames.append(df)

        # Merge
        if frames == []:
            raise EmptyDatabaseError("Dataframe is empty")
        else:
            df = pd.concat(frames, ignore_index=True)  # reindex
            if len(df) == 0:
                raise EmptyDatabaseError(
                    "Dataframe is empty on range "
                    + "{0:.2f}-{1:.2f} cm-1".format(wavenum_min, wavenum_max)
                )

        df = parse_local_quanta(df, molecule)
        df = parse_global_quanta(df, molecule)

        # Remove non numerical attributes
        if drop_non_numeric:
            if "branch" in df:
                replace_PQR_with_m101(df)
            df = drop_object_format_columns(df, verbose=self.verbose)

        # Complete database with molecular parameters
        self._fetch_molecular_parameters(df)

        self.df0 = df

        # %% Init Partition functions (with energies)
        # ------------

        self._init_equilibrium_partition_functions(parfunc, parfuncfmt)

        # If energy levels are given, initialize the partition function calculator
        # (necessary for non-equilibrium). If levelsfmt == 'radis' then energies
        # are calculated ab initio from radis internal species database constants
        if load_energies:
            try:
                self._init_rovibrational_energies(levels, levelsfmt)
            except KeyError as err:
                print(err)
                raise KeyError(
                    "Error while fetching rovibrational energies for "
                    + "{0}, iso={1} in RADIS built-in spectroscopic ".format(
                        molecule, isotope
                    )
                    + "constants (see details above). If you only need "
                    + "equilibrium spectra, try using 'load_energies=False' "
                    + "in fetch_databank"
                )

        # %% Store

        # Let's store all params so they can be parsed by "get_conditions()"
        # and saved in output spectra information
        self.params.dbpath = "fetched from " + source
        self.params.dbformat = dbformat
        if levels is not None:
            self.levelspath = ",".join([format_paths(lvl) for lvl in levels.values()])
        else:
            self.levelspath = None
        self.params.levelsfmt = levelsfmt
        self.params.parfuncpath = format_paths(parfunc)
        self.params.parfuncfmt = parfuncfmt

        return

    def load_databank(
        self,
        name=None,
        path=None,
        format=None,
        parfunc=None,
        parfuncfmt=None,
        levels=None,
        levelsfmt=None,
        db_use_cached=None,
        db_assumed_sorted=True,
        load_energies=True,
        include_neighbouring_lines=True,
        drop_columns="auto",
        buffer="RAM",
    ):
        """ Loads databank from shortname in the :ref:`Configuration file <label_lbl_config_file>` 
        (`~/.radis`), or by manually setting all attributes.

        Databank includes:

        - lines
        - partition function & format (tabulated or calculated)
        - (optional) energy levels, format

        It also fetches molecular parameters (molar mass, abundance) for 
        all molecules in database


        Parameters
        ----------

        name: a section name specified in your ``~/.radis``
            ``.radis`` has to be created in your HOME (Unix) / User (Windows). If
            not ``None``, all other arguments are discarded.
            Note that all files in database will be loaded and it may takes some
            time. Better limit the database size if you already know what
            range you need. See :ref:`Configuration file <label_lbl_config_file>` and 
            :data:`~radis.misc.config.DBFORMAT` for expected 
            ``~/.radis`` format


        Other Parameters
        ----------------

        path: str, list of str, None
            list of database files, or name of a predefined database in the 
            :ref:`Configuration file <label_lbl_config_file>` (`~/.radis`)
            Accepts wildcards ``*`` to select multiple files 

        format: ``'hitran'``, ``'cdsd-hitemp'``, ``'cdsd-4000'``, or any of :data:`~radis.lbl.loader.KNOWN_DBFORMAT`
            database type. ``'hitran'`` for HITRAN/HITEMP, ``'cdsd-hitemp'`` 
            and ``'cdsd-4000'`` for the different CDSD versions. Default ``'hitran'``

        parfuncfmt: ``'hapi'``, ``'cdsd'``, or any of :data:`~radis.lbl.loader.KNOWN_PARFUNCFORMAT`
            format to read tabulated partition function file. If ``hapi``, then
            HAPI (HITRAN Python interface) [1]_ is used to retrieve them (valid if
            your database is HITRAN data). HAPI is embedded into RADIS. Check the
            version.

        parfunc: filename or None
            path to tabulated partition function to use.
            If `parfuncfmt` is `hapi` then `parfunc` should be the link to the
            hapi.py file. If not given, then the hapi.py embedded in RADIS is used (check version)

        levels: dict of str or None
            path to energy levels (needed for non-eq calculations). Format:
            {1:path_to_levels_iso_1, 3:path_to_levels_iso3}. Default ``None``

        levelsfmt: 'cdsd-pc', 'radis' (or any of :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`) or ``None``
            how to read the previous file. Known formats: (see :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`).
            If ``radis``, energies are calculated using the diatomic constants in radis.db database
            if available for given molecule. Look up references there.
            If None, non equilibrium calculations are not possible. Default ``None``.

        db_use_cached: boolean, or ``None``
            if ``True``, a pandas-readable csv file is generated on first access,
            and later used. This saves on the datatype cast and conversion and
            improves performances a lot. But! ... be sure to delete these files
            to regenerate them if you happen to change the database. If ``'regen'``,
            existing cached files are removed and regenerated.
            It is also used to load energy levels from ``.h5`` cache file if exist.
            If ``None``, the value given on Factory creation is used. Default ``None``

        db_assumed_sorted: boolean
            load_databank first reads the first line and check it's relevant.
            This improves database loading times if not all files are required,
            but it assumes database files are sorted in wavenumber!
            Default ``True``
            
        load_energies: boolean
            if ``False``, dont load energy levels. This means that nonequilibrium
            spectra cannot be calculated, but it saves some memory. Default ``True``

        include_neighbouring_lines: bool
            ``True``, includes off-range, neighbouring lines that contribute
            because of lineshape broadening. The ``broadening_max_width`` 
            parameter is used to determine the limit. Default ``True``.
            
        Other arguments are related to how to open the files:

        buffer: ``'RAM'``, ``'h5'``, ``'direct'``
            Different modes for loading up database: either directly in 'RAM' mode,
            or in 'h5' mode.

            - 'RAM': is faster but memory hunger
            - 'h5': handles better a bigger database (> 1M lines): slower (up to 3x), but less
              risks of MemoryErrors 
            - 'direct': file is read directly from a single h5 file under key 'df'
              Fastest of all, doesnt check the database validity or format. Use only
              if you have a single, already formatted database file (used by Factory
              when reloading database)
              
            Default ``'RAM'``

        drop_columns: list
            columns names to drop from Line DataFrame after loading the file. 
            Not recommended to use, unless you explicitely want to drop information 
            (for instance if dealing with too large databases). If ``[]``, nothing 
            is dropped. If ``'auto'``, parameters considered useless
            are dropped. See :data:`~radis.lbl.loader.drop_auto_columns_for_dbformat`
            and :data:`~radis.lbl.loader.drop_auto_columns_for_levelsfmt`. 
            If ``'all'``, parameters considered unecessary for equilibrium calculations 
            are dropped, including all information about lines that could be otherwise 
            available in :py:meth:`~radis.spectrum.spectrum.Spectrum` method.
            Warning: nonequilibrium calculations are not possible in this mode. 
            Default ``'auto'``.

        Notes
        -----
        
        Performances of buffer mode:
            
        on the 2Gb CDSD-HITEMP database (1-20), already cached in .h5
        
        - ``'RAM'``: 7.1 s
        - ``'h5'``: 21 s 

        See Also
        --------
        
        - Only load when needed: :meth:`~radis.lbl.loader.DatabankLoader.init_databank`
        - Download from HITRAN: :meth:`~radis.lbl.loader.DatabankLoader.fetch_databank`
        
        :ref:`Configuration file <label_lbl_config_file>` with:
        - all line database formats: :py:data:`~radis.misc.config.DBFORMAT`
        - all energy levels database formats: :py:data:`~radis.misc.config.LVLFORMAT`

        References
        ----------

        .. [1] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_

        """
        # TODO remove tempfile option to clean the code.
        try:

            # %% Check inputs
            # ---------

            (
                name,
                path,
                dbformat,
                parfunc,
                parfuncfmt,
                levels,
                levelsfmt,
                db_use_cached,
                db_assumed_sorted,
                drop_columns,
                buffer,
                load_energies,
                include_neighbouring_lines,
            ) = self._check_database_params(
                name=name,
                path=path,
                format=format,
                parfunc=parfunc,
                parfuncfmt=parfuncfmt,
                levels=levels,
                levelsfmt=levelsfmt,
                db_use_cached=db_use_cached,
                db_assumed_sorted=db_assumed_sorted,
                load_energies=load_energies,
                include_neighbouring_lines=include_neighbouring_lines,
                drop_columns=drop_columns,
                buffer=buffer,
            )
            # Now that we're all set, let's load everything

            # %% Line database
            # ------------
            self.df0 = self._load_databank(
                path,
                dbformat,
                levelsfmt=levelsfmt,
                db_use_cached=db_use_cached,
                db_assumed_sorted=db_assumed_sorted,
                buffer=buffer,
                drop_columns=drop_columns,
                include_neighbouring_lines=include_neighbouring_lines,
            )

            # Check the molecule is what we expected
            if len(set(self.df0.id)) != 1:  # only 1 molecule supported ftm
                raise NotImplementedError(
                    "Only 1 molecule at a time is currently supported "
                    + "in RADIS. Calculate them independently then "
                    + "use MergeSlabs"
                )
            if self.input.molecule not in ["", None]:
                assert self.input.molecule == get_molecule(
                    self.df0.id[0]
                )  # assert molecule is what we expected
            else:
                self.input.molecule = get_molecule(self.df0.id[0])  # get molecule

            # %% Partition functions (with energies)
            # ------------

            self._init_equilibrium_partition_functions(parfunc, parfuncfmt)

            # If energy levels are given, initialize the partition function calculator
            # (necessary for non-equilibrium). If levelsfmt == 'radis' then energies
            # are calculated ab initio from radis internal species database constants
            if load_energies:
                self._init_rovibrational_energies(levels, levelsfmt)

            # %% Store

            # Let's store all params so they can be parsed by "get_conditions()"
            # and saved in output spectra information
            self._store_database_params(
                name=name,
                path=path,
                format=dbformat,
                parfunc=parfunc,
                parfuncfmt=parfuncfmt,
                levels=levels,
                levelsfmt=levelsfmt,
                db_use_cached=db_use_cached,
                db_assumed_sorted=db_assumed_sorted,
                include_neighbouring_lines=include_neighbouring_lines,
            )
            return

        except MemoryError as err:
            # An error occured: clean before crashing
            self._clean_temp_file()
            printr(
                " Error while loading the database. Retry with "
                + "`save_memory=True` option, or `save_memory=2` "
                + "(warning: for equilibrium calculations only)"
            )
            raise err
        except:
            # An error occured: clean before crashing
            self._clean_temp_file()
            raise

    def _check_database_params(
        self,
        name=None,
        path=None,
        format=None,
        parfunc=None,
        parfuncfmt=None,
        levels=None,
        levelsfmt=None,
        db_use_cached=None,
        db_assumed_sorted=True,
        load_energies=True,
        include_neighbouring_lines=True,
        drop_columns="auto",
        buffer="RAM",
    ):
        """ Check that database parameters are valid, in particular that
        paths exist. Loads all parameters if a Database from .radis config file
        was given
        
        Returns
        -------
        
        tuple
            (name, path, dbformat, parfunc, parfuncfmt, levels, levelsfmt,
             db_use_cached, db_assumed_sorted, drop_columns, buffer)
        
        """

        dbformat = format

        # Get database format and path
        # ... either from name (~/.radis config file)
        if name is not None:
            try:
                entries = getDatabankEntries(name)
            except IOError:
                print("There was a problem looking for database name: {0}".format(name))
                raise

            if self.verbose:
                print("Using database: {0}".format(name))
                printDatabankEntries(name)
                print("\n")
            path = entries["path"]
            dbformat = entries["format"]
            if "parfunc" in entries:
                parfunc = entries["parfunc"]
            if "parfuncfmt" in entries:
                parfuncfmt = entries["parfuncfmt"]
            if "levels" in entries:
                levels = entries["levels"]
            if "levelsfmt" in entries:
                levelsfmt = entries["levelsfmt"]
        # ... or directly
        else:
            # make sure at least path and dbformat are given
            if path is None and dbformat is None:
                try:
                    dblist = getDatabankList()
                except IOError:
                    dblist = []
                raise ValueError(
                    "No database name. Please give a path and a dbformat"
                    + ", or use one of the predefined databases in your"
                    + " ~/.radis: {0}".format(",".join(dblist))
                )

        # Check input types are correct
        if isinstance(path, string_types):  # make it a list
            path = get_files_from_regex(path)

        if dbformat not in KNOWN_DBFORMAT:
            # >>>>>>>>>>>
            # Deprecation errors (added in 0.9.21. Remove after 1.0.0)
            if dbformat == "cdsd":
                raise DeprecationWarning(
                    "`cdsd` database format was renamed `cdsd-hitemp` after 0.9.21"
                )
            if dbformat == "cdsd4000":
                raise DeprecationWarning(
                    "`cdsd4000` database format was renamed `cdsd-4000` after 0.9.21"
                )
            # <<<<<<<<<<<
            raise ValueError(
                "Database format ({0}) not in known list: {1}".format(
                    dbformat, KNOWN_DBFORMAT
                )
            )
        if levels is not None and levelsfmt not in KNOWN_LVLFORMAT:
            raise ValueError(
                "Energy level format ({0}) not in known list: {1}".format(
                    levelsfmt, KNOWN_LVLFORMAT
                )
            )
        if parfuncfmt not in [None] + KNOWN_PARFUNCFORMAT:
            raise ValueError(
                "Partition function format ({0}) not in known list: {1}".format(
                    parfuncfmt, KNOWN_PARFUNCFORMAT
                )
            )

        # Check all path exists
        # ... remove empty paths first
        path = [p for p in path if p != ""]
        # ... test paths
        for p in path:
            if not exists(p):
                raise FileNotFoundError("databank lines file: `{0}`".format(p))
        if levels is not None:
            for iso, lvl in levels.items():  # one file per isotope
                if not exists(lvl):
                    raise FileNotFoundError("levels = `{0}`".format(lvl))
        if parfunc is not None:  # all isotopes in same file?
            if not exists(parfunc):
                raise FileNotFoundError("parfunc = `{0}`".format(parfunc))

        # Get cache
        if db_use_cached is None:
            db_use_cached = self.params.db_use_cached

        return (
            name,
            path,
            dbformat,
            parfunc,
            parfuncfmt,
            levels,
            levelsfmt,
            db_use_cached,
            db_assumed_sorted,
            drop_columns,
            buffer,
            load_energies,
            include_neighbouring_lines,
        )

    def _store_database_params(
        self,
        name=None,
        path=None,
        format=None,
        parfunc=None,
        parfuncfmt=None,
        levels=None,
        levelsfmt=None,
        db_use_cached=None,
        db_assumed_sorted=True,
        load_energies=True,
        include_neighbouring_lines=True,
    ):
        """ store all params so they can be parsed by "get_conditions()"
        and saved in output spectra information
        
        Notes
        -----
        
        Only those params stored in self.params will be kept eventually 
        """

        # Store in params if they can change the physical output, else store
        # in misc. See retrieve_from_database() for more information.

        self.params.dbpath = ",".join(
            [format_paths(k) for k in path]
        )  # else it's a nightmare to store
        self.params.dbformat = format
        if levels is not None:
            self.levelspath = ",".join([format_paths(lvl) for lvl in levels.values()])
        else:
            self.levelspath = None
        self.params.levelsfmt = levelsfmt
        self.params.parfuncpath = format_paths(parfunc)
        self.params.parfuncfmt = parfuncfmt
        self.params.db_assumed_sorted = db_assumed_sorted
        self.params.include_neighbouring_lines = include_neighbouring_lines
        self.misc.load_energies = load_energies

    def init_database(
        self,
        path,
        autoretrieve=True,
        autoupdate=True,
        add_info=["Tvib", "Trot"],
        add_date="%Y%m%d",
        compress=True,
    ):
        """ Init a :class:`~radis.tools.database.SpecDatabase` folder in 
        ``path`` to later store our spectra. Spectra can also be automatically 
        retrieved from the database instead of being calculated

        Parameters
        ----------

        path: str
            path to database folder. If it doesnt exist, create it
            Accepts wildcards ``*`` to select multiple files

        autoretrieve: boolean, or ``'force'``
            if ``True``, a database lookup is performed whenever a new spectrum
            is calculated. If the spectrum already exists then it is retrieved
            from the database instead of being calculated. Spectra are considered
            the same if all the stored conditions fit. If set to ``'force'``, an error
            is raised if the spectrum is not found in the database (use it for
            debugging). Default ``True``

        autoupdate: boolean
            if ``True``, all spectra calculated by this Factory are automatically
            exported in database. Default ``True`` (but only if init_database is
            explicitely called by user)

        add_info: list, or ``None``/``False``
            append these parameters and their values if they are in conditions.
            Default ``['Tvib', 'Trot']``

        add_date: str, or ``None``/``False``
            adds date in strftime format to the beginning of the filename.
            Default '%Y%m%d'

        compress: boolean
            if ``True``, Spectrum are read and written in binary format. This is faster,
            and takes less memory space. Default ``True``

        Returns
        -------
        
        db: SpecDatabase
            the database where spectra will be stored or retrieved

        """

        db = SpecDatabase(path, add_info=add_info, add_date=add_date, binary=compress)

        self.SpecDatabase = db
        self.database = format_paths(path)  # just to appear in conditions
        self.autoupdatedatabase = autoupdate
        self.autoretrievedatabase = autoretrieve

        return db

    # %% ======================================================================
    # PRIVATE METHODS - DATABASE LOADING
    # ----------------------------------------------
    # _init_equilibrium_partition_functions
    # _init_rovibrational_energies
    # _load_databank
    # _reload_databank
    # _check_line_databank
    # _retrieve_from_database
    # _get_partition_function_interpolator
    # _get_partition_function_calculator
    # _fetch_molecular_parameters
    # _get_temp_file
    # _clean_temp_file
    #
    # =========================================================================

    def _init_equilibrium_partition_functions(self, parfunc, parfuncfmt):
        """ Initializes equilibrium partition functions in ``self.parsum_tab``

        Parameters
        ----------

        parfuncfmt: 'cdsd', 'hapi' (see :data:`~radis.lbl.loader.KNOWN_PARFUNCFORMAT`)
            format to read tabulated partition function file. If `hapi`, then
            HAPI (HITRAN Python interface) [1]_ is used to retrieve them (valid if
            your database is HITRAN data). HAPI is embedded into RADIS. Check the
            version.

        parfunc: filename or None
            path to tabulated partition function to use.
            If ``parfuncfmt`` is ``hapi`` then ``parfunc`` should be the link to the
            hapi.py file. If not given, then the hapi.py embedded in RADIS is used (check version)

        """

        # Let's get the tabulated partition function (to calculate eq spectra)
        molecule = self.input.molecule
        state = self.input.state
        self.parsum_tab[molecule] = {}
        for iso in self._get_isotope_list():
            self.parsum_tab[molecule][iso] = {}
            ParsumTab = self._build_partition_function_interpolator(
                parfunc, parfuncfmt, self.input.molecule, isotope=iso
            )
            self.parsum_tab[molecule][iso][state] = ParsumTab

    def _init_rovibrational_energies(self, levels, levelsfmt):
        """ Initializes non equilibrium partition (which contain rovibrational
        energies) and store them in ``self.parsum_calc``

        Parameters
        ----------

        levels: dict of str or None
            path to energy levels (needed for non-eq calculations). Format:
            {1:path_to_levels_iso_1, 3:path_to_levels_iso3}. Default ``None``

        levelsfmt: 'cdsd-pc', 'radis' (see :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`) or None
            how to read the previous file. Known formats: :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`
            If ```radis```, energies are calculated using the diatomic constants
            in ``radis.db`` database. If available for given molecule. Look up
            references there. If ``None``, non equilibrium calculations are not
            possible.

        """

        molecule = self.input.molecule
        state = self.input.state
        # only defined for one molecule at the moment (add loops later!)
        self.parsum_calc[molecule] = {}
        # look up if energy levels are defined in an input file:
        if levels is not None:
            for iso, lvl in levels.items():
                self.parsum_calc[molecule][iso] = {}
                ParsumCalc = self._build_partition_function_calculator(
                    lvl, levelsfmt, isotope=iso
                )
                self.parsum_calc[molecule][iso][state] = ParsumCalc
        # energy levels arent specified in a tabulated file, but we can still
        # calculate them directly from Dunham expansions:
        elif levelsfmt == "radis":
            for iso in self._get_isotope_list():
                self.parsum_calc[molecule][iso] = {}
                ParsumCalc = self._build_partition_function_calculator(
                    None, levelsfmt, isotope=iso
                )
                self.parsum_calc[molecule][iso][state] = ParsumCalc

    def _check_line_databank(self):
        """ Make sure database is loaded, loads if it isnt and we have all
        the information needed. Databank has been initialized by 
        :meth:`~radis.lbl.loader.DatabankLoader._init_databank`
        
        """

        # Make sure database is loaded
        if self.df0 is None:
            # Either we're in a save memory mode, i.e, database has been
            # removed from RAM but a HDF5 temp_file was created on disk to
            # be reloaded instantly...
            try:
                if self.save_memory:
                    self._reload_databank()
                    return
            except FileNotFoundError:
                pass  # Deal with it the normal way

            # Now, the databank may have been initialised (references given)
            # but not loaded (because it takes time).
            if self._databank_args != [] or self._databank_kwargs != {}:
                # Load databank
                self.load_databank(*self._databank_args, **self._databank_kwargs)
                # Clean variables
                self._databank_args = []
                self._databank_kwargs = {}
            # Or, just crash
            else:
                raise AttributeError("Load databank first (.load_databank())")

        #        # Reset index
        #        #    (cost ~ 1 ms but is needed if the user manually edited the database
        #        #    in between the load_database() and the calculation command
        self.df0.reset_index(inplace=True, drop=True)  # drop: don't add old index
        # Finally commented: code may crash if users edit the database manually
        # (ex: modify broadening coefficients) and forgot to reset the index,
        # but that's for advanced users anyway. The cost (time+dont know what
        # side effects may occur) is not worth it
        # @EP: reactivated after f014007

        # Check format
        # (it can happen that they changed if database was edited manually. That can break
        # the code during look-up of levels later on.)
        for k in [
            "vu",
            "vl",
            "v1u",
            "v1l",
            "v2u",
            "v2l",
            "l2u",
            "l2l",
            "v3u",
            "v3l",
            "ru",
            "rl",
            "polyu",
            "polyl",
            "wangu",
            "wangl",
            "ranku",
            "rankl",
        ]:
            if (
                k in self.df0.columns
                and self.df0.dtypes[k] != np.int64
                and count_nans(self.df0[k]) == 0
            ):
                self.warn(
                    "Format of column {0} was {1} instead of int. Changed to int".format(
                        k, self.df0.dtypes[k]
                    )
                )
                self.df0[k] = self.df0[k].astype(np.int64)

        # Check metadata ('molar_mass' and 'Ia' are either columns either
        # metadata. They can be lost when transfering object.)
        try:
            self.df0.molar_mass
            self.df0.Ia
        except AttributeError as err:
            raise AttributeError(
                str(err)
                + " : attribute missing in line "
                + "dataframe sf.df0. Make sure you didnt overwrite the line "
                + "dataframe sf.df0 manually. If so, replace any "
                + "`sf.df0=...` line with inplace operations such as "
                + "`sf.df0.drop(..., inplace=True)`. See "
                + "https://stackoverflow.com/q/33103988"
            )

    def _load_databank(
        self,
        database,
        dbformat,
        levelsfmt,
        db_use_cached,
        db_assumed_sorted=False,
        buffer="RAM",
        drop_columns="auto",
        include_neighbouring_lines=True,
    ):
        """ Loads all available database files and keep the relevant one. Returns
        a Pandas dataframe

        Parameters
        ----------

        database: list of str
            list of database files

        db_use_cached: boolean
            if ``True``, a pandas-readable csv file is generated on first access,
            and later used. This saves on the datatype cast and conversion and
            improves performances a lot. But! ... be sure to delete these files
            to regenerate them if you happen to change the database. Default ``False``

        buffer: ``'RAM'``, ``'h5'``, ``'direct'``
            Different modes for loading up a database: either directly in ``'RAM'`` mode,
            or in ``'h5'`` mode.

            - ``'RAM'``: is faster but memory hunger
            - ``'h5'``: handles better a bigger database (> 1M lines): slower (up to 3x), but less
            risks of MemoryErrors 
            - ``'auto'``: choose depending on your number of files in database. If>30 files,
            swith to 'h5'. Former default function, I know switched back to RAM. 
            - ``'direct'``: file is read directly from a single h5 file under key ``'df'``
              Fastest of all, doesnt check the database validity or format. Use only
              if you have a single, already formatted database file (used by Factory
              when reloading database)
              
            Default ``'RAM'``

        drop_columns: list
            columns names to drop from Line DataFrame after loading the file. 
            Not recommended to use, unless you explicitely want to drop information 
            (for instance if dealing with too large databases). If ``[]``, nothing 
            is dropped. If ``'auto'``, parameters considered useless 
            are dropped. See :data:`~radis.lbl.loader.drop_auto_columns_for_dbformat`
            and :data:`~radis.lbl.loader.drop_auto_columns_for_levelsfmt`. 
            If ``'all'``, parameters considered unecessary for equilibrium calculations 
            are dropped, including all information about lines that could be otherwise 
            available in :py:meth:`~radis.spectrum.spectrum.Spectrum` method.
            Warning: nonequilibrium calculations are not possible in this mode. 
            Default ``'auto'``.

        Other Parameters
        ----------------

        include_neighbouring_lines: bool
            ``True``, includes off-range, neighbouring lines that contribute
            because of lineshape broadening. The ``broadening_max_width`` 
            parameter is used to determine the limit. Default ``True``.
            
        Notes
        -----
        
        Performances of buffer mode:
            
        on the 2Gb CDSD-HITEMP database (1-20), already cached in .h5
        
        - ``'RAM'``: 7.1 s
        - ``'h5'``: 21 s 

        """

        # Check inputs
        assert db_use_cached in [True, False, "regen"]
        assert buffer in ["RAM", "h5", "direct"]

        if self.verbose >= 2:
            printg("Loading Line databank")
            t0 = time()

        # Init variables
        verbose = self.verbose
        warnings_default = self.warnings["default"] if self.warnings else False
        if include_neighbouring_lines:
            wavenum_min = self.params.wavenum_min_calc
            wavenum_max = self.params.wavenum_max_calc
        else:
            wavenum_min = self.input.wavenum_min
            wavenum_max = self.input.wavenum_max

        # Check inputs
        if buffer == "direct":
            assert len(database) == 1
            assert database[0].endswith("h5")
        if drop_columns == "auto":
            drop_columns = (
                drop_auto_columns_for_dbformat[dbformat]
                + drop_auto_columns_for_levelsfmt[levelsfmt]
            )

        # subroutine load_and_concat
        # --------------------------------------
        def load_and_concat(files, buffer):
            """ Two modes of storage: either directly in ``'RAM'`` mode, or in ``'h5'``
            mode. ``'RAM'`` is faster but memory hunger, ``'h5'`` handles better
            a bigger database

            Parameters
            ----------

            files: str
                path

            buffer: ``'direct'``, ``'h5'``, ``'RAM'``
                see _load_databank info

            """

            if buffer == "direct":
                df = pd.read_hdf(database[0], key="df")
                df = df.reset_index()
                return df
            elif buffer == "h5":
                tempfile = self._get_temp_file()
            else:
                frames = []

            for i, filename in enumerate(files):
                # first dont load the whole database in unecessary.
                # Check 1st line if we assume it is sorted

                if __debug__:
                    printdbg("Loading {0}/{1}".format(i + 1, len(files)))

                if i % 31 == 30:
                    gc.collect()  # force garbage collection as we may be generating tons of data

                if db_assumed_sorted:
                    # Note on performance: reading the first line of .txt file is still
                    # much faster than reading the whole hdf5 file
                    if dbformat == "cdsd-hitemp":
                        df = cdsd2df(
                            filename,
                            version="hitemp",
                            count=1,
                            cache=False,
                            verbose=False,
                            drop_non_numeric=False,
                        )
                        if df.wav.loc[0] > wavenum_max:
                            if verbose:
                                print(
                                    "Database file {0} > {1:.6f}cm-1: irrelevant and not loaded".format(
                                        filename, wavenum_max
                                    )
                                )
                            continue
                    elif dbformat == "cdsd-4000":
                        df = cdsd2df(
                            filename,
                            version="4000",
                            count=1,
                            cache=False,
                            verbose=False,
                            drop_non_numeric=False,
                        )
                        if df.wav.loc[0] > wavenum_max:
                            if verbose:
                                print(
                                    "Database file {0} > {1:.6f}cm-1: irrelevant and not loaded".format(
                                        filename, wavenum_max
                                    )
                                )
                            continue
                    elif dbformat == "hitran":
                        df = hit2df(
                            filename,
                            count=1,
                            cache=False,
                            verbose=False,
                            drop_non_numeric=False,
                        )
                        if df.wav.loc[0] > wavenum_max:
                            if verbose:
                                print(
                                    "Database file {0} > {1:.6f}cm-1: irrelevant and not loaded".format(
                                        filename, wavenum_max
                                    )
                                )
                            continue
                    else:
                        raise ValueError(
                            "The database format is unknown: {0}".format(dbformat)
                        )

                # Now read all the lines
                # ... this is where the cache files are read/generated.
                if dbformat == "cdsd-hitemp":
                    df = cdsd2df(
                        filename,
                        version="hitemp",
                        cache=db_use_cached,
                        verbose=verbose,
                        drop_non_numeric=True,
                    )
                elif dbformat == "cdsd-4000":
                    df = cdsd2df(
                        filename,
                        version="4000",
                        cache=db_use_cached,
                        verbose=verbose,
                        drop_non_numeric=True,
                    )
                elif dbformat == "hitran":
                    df = hit2df(
                        filename,
                        cache=db_use_cached,
                        verbose=verbose,
                        drop_non_numeric=True,
                    )
                else:
                    raise ValueError("Unknown dbformat: {0}".format(dbformat))

                # Drop the end
                if db_assumed_sorted:
                    if df.wav.iloc[-1] < wavenum_min:
                        if verbose:
                            print(
                                "Database file {0} < {1:.6f}cm-1: irrelevant and rejected".format(
                                    filename, wavenum_min
                                )
                            )
                        continue

                # Drop columns (helps fix some Memory errors)
                dropped = []
                for col in df.columns:
                    if col in drop_columns or (
                        drop_columns == "all" and col not in drop_all_but_these
                    ):
                        del df[col]
                        dropped.append(col)
                if verbose >= 2 and len(dropped) > 0:
                    print("Dropped columns: {0}".format(dropped))

                # Crop to the wavenumber of interest
                df = df[(df.wav >= wavenum_min) & (df.wav <= wavenum_max)]

                if __debug__:
                    if len(df) == 0:
                        printdbg(
                            "File {0} loaded for nothing (out of range)".format(
                                filename
                            )
                        )

                # Select correct isotope(s)
                if self.input.isotope != "all":
                    isotope = [float(k) for k in self.input.isotope.split(",")]
                    df = df[df.iso.isin(isotope)]

                if buffer == "h5":
                    with catch_warnings():  # temporary
                        filterwarnings("ignore", category=DeprecationWarning)
                        # there is a warning "inspect.getargspec() is deprecated,
                        # use inspect.signature() instead" from the pytable library
                        df.to_hdf(tempfile, "df", format="table", append=True)
                        # note that format='table' is slow (5x slower than fixed)
                        # but 'fixed' doesnt allow to append
                else:
                    frames.append(df)

            # Finally: Concatenate all
            if buffer == "h5":
                #                store=pd.HDFStore(tempfile)
                #                df = store.select('df')
                #                store.close()
                try:
                    df = pd.read_hdf(tempfile, "df")
                except KeyError:  # happens if database is empty. A database empty error
                    # will be raised a few lines below
                    df = pd.DataFrame()
                os.remove(tempfile)
                if "Unnamed: 0" in df:
                    # there is no ignore_index option in HDF...
                    # so an index column can be created, so we delete it and regenerate one
                    del df["Unnamed: 0"]
                df = df.reset_index()
            else:
                if frames == []:
                    df = (
                        pd.DataFrame()
                    )  # a database empty error will be raised a few lines below
                else:
                    df = pd.concat(frames, ignore_index=True)  # reindex

            return df

        # end subroutine load_and_concat
        # --------------------------------------

        try:
            df = load_and_concat(database, buffer=buffer)
        except MemoryError:
            # Not enough RAM
            if buffer == "h5":
                # If wasnt 'h5', try again with h5 buffer on disk
                if verbose:
                    print(
                        "Memory error while loading database. Trying through h5 buffer"
                    )
                df = load_and_concat(database, buffer="h5")
            else:
                raise
        finally:  # Clean temp file in any case
            self._clean_temp_file()

        # Final checks

        # ... error in Pandas? Sometimes _metadata is preserved over several runs.
        # ... Clean it here.
        if __debug__ and len(df._metadata) > 0:
            printdbg("df._metadata was not []. Cleaning")
        df._metadata = []

        # ... check database is not empty

        if len(df) == 0:
            msg = (
                "Reference databank "
                + "has 0 lines in range {0:.2f}-{1:.2f}cm-1".format(
                    wavenum_min, wavenum_max
                )
            ) + " ({0:.2f}-{1:.2f}nm) Check your range !".format(
                cm2nm(wavenum_min), cm2nm(wavenum_max)
            )
            raise ValueError(msg)

        maxwavdb = df.wav.max()
        minwavdb = df.wav.min()

        # ... Explicitely write molecule if not given
        if self.input.molecule in [None, ""]:
            id_set = df.id.unique()
            if len(id_set) > 1:
                raise NotImplementedError(
                    "RADIS expects one molecule per run for the "
                    + "moment. Got {0}. Use different runs ".format(id_set)
                    + "and use MergeSlabs(out='transparent' afterwards"
                )
            self.input.molecule = get_molecule(id_set[0])

        # ... explicitely write all isotopes based on isotopes found in the database
        if self.input.isotope == "all":
            self.input.isotope = ",".join(
                [str(k) for k in self._get_isotope_list(df=df)]
            )

        # ... check all requested isotopes are present (only if not using 'all',
        # in which case we assume the user didnt really care about a particular
        # isotope and shouldnt be surprised if it's not there)
        else:
            isotope_list = self._get_isotope_list()

            # check no isotope shows 0 line in this range. Raise an warning if it
            # happens
            for k in isotope_list:
                if not (sum(df.iso == k) > 0):
                    msg = (
                        "Reference databank ({0:.2f}-{1:.2f}cm-1)".format(
                            minwavdb, maxwavdb
                        )
                        + " has 0 lines in range ({0:.2f}-{1:.2f}cm-1)".format(
                            wavenum_min, wavenum_max
                        )
                    ) + " for isotope {0}. Change your range or isotope options".format(
                        k
                    )
                    printwarn(msg, verbose, warnings_default)

        # ... check the database range looks correct
        # ... (i.e, there are some lines on each side of the requested range:
        # ... else, maybe User forgot to add all requested lines in the database '''
        if include_neighbouring_lines:
            broadening = self.params.broadening_max_width
            if minwavdb > wavenum_min + broadening:
                # no lines on left side
                self.warn(
                    "There are no lines in database in range {0:.5f}-{1:.5f}cm-1 ".format(
                        wavenum_min, wavenum_min + broadening
                    )
                    + "to calculate the effect "
                    + "of neighboring lines. Did you add all lines in the database?",
                    "OutOfRangeLinesWarning",
                )
            if maxwavdb < wavenum_max - broadening:
                # no lines on right side
                self.warn(
                    "There are no lines in database in range {0:.5f}-{1:.5f}cm-1 ".format(
                        maxwavdb - broadening, maxwavdb
                    )
                    + "to calculate the effect "
                    + "of neighboring lines. Did you add all lines in the database?",
                    "OutOfRangeLinesWarning",
                )

        # Complete database with molecular parameters
        self._fetch_molecular_parameters(df)

        if self.verbose >= 2:
            printg(
                "Loaded databank in {0:.1f}s ({1:,d} lines)".format(
                    time() - t0, len(df)
                )
            )

        return df

    def _reload_databank(self):
        """ In save_memory mode we're trying to save RAM so reference dataframe (df0)
        will be deleted after scaled database (df1) is created. This makes it
        impossible to calculate another spectrum afterwards, without reloading
        the database: in that case, we have kept the temporary file for some time
        and try to regenerate df0 here """

        path = self._get_temp_file()
        if not exists(path):
            raise FileNotFoundError("temp file not present. Cant reload database")
        dbformat = self.params.dbformat
        db_use_cached = self.params.db_use_cached
        db_assumed_sorted = self.db_assumed_sorted
        levelsfmt = self.params.levelsfmt

        t0 = time()
        self.df0 = self._load_databank(
            [path],
            dbformat,
            levelsfmt,
            db_use_cached=db_use_cached,
            db_assumed_sorted=db_assumed_sorted,
            buffer="direct",
        )

        if __debug__:
            printdbg(
                "Databank reloaded from temporary h5 file in {0:.2f}s".format(
                    time() - t0
                )
            )

    def _get_isotope_list(self, molecule=None, df=None):
        """ Returns list of isotopes for given molecule 
        Parse the Input conditions (fast). If a line database is given, parse the 
        line database instead (slow)
        
        Parameters
        ----------
        
        molecule: str
            molecule 
            
        df: pandas DataFrame, or ``None``
            line database to parse. Default ``None`` 
        """

        if molecule is not None and self.input.molecule != molecule:
            raise ValueError(
                "Expected molecule is {0} according to the inputs, but got {1} ".format(
                    self.input.molecule, molecule
                )
                + "in line database. Check your `molecule=` parameter, or your "
                + "line database."
            )

        if df is None:
            isotope_list = self.input.isotope.split(",")
        else:  # get all isotopes in line database
            isotope_list = set(df.iso)

        return [int(k) for k in isotope_list]

    def _retrieve_from_database(self, ignore_misc=True):
        """ Retrieve a spectrum from a SpecDatabase database, if it matches
        current factory conditions

        Parameters
        ----------

        ignore_misc: boolean
            if ``True``, then all attributes considered as Factory 'descriptive' 
            parameters, as defined in :meth:`~radis.lbl.loader.get_conditions` are ignored when
            comparing the database to current factory conditions. It should
            obviously only be attributes that have no impact on the Spectrum
            produced by the factory. Default ``True`` 

        """

        conditions = {
            k: v
            for (k, v) in self.get_conditions(ignore_misc=ignore_misc).items()
            if k in self.SpecDatabase.conditions()
        }
        conditions = {
            k: v
            for (k, v) in conditions.items()
            if not k in self._autoretrieveignoreconditions
        }

        s = self.SpecDatabase.get(**conditions)
        if len(s) == 1:
            if self.verbose:
                print("Spectrum found in database!")
            return s[0]
        elif len(s) > 1:
            printwarn(
                "Multiple spectra in database match conditions. Picked" + " the first",
                self.verbose,
                self.warnings,
            )
            return s[0]
        else:
            if self.autoretrievedatabase == "force":

                def get_best_match():
                    """ Returns the Spectrum that matches the input conditions
                    better. Used to give a better error message in case no
                    Spectrum was found """
                    best = None
                    score = 0
                    for s in self.SpecDatabase.get():
                        sc = compare_dict(conditions, s.conditions, verbose=False)
                        if sc > score:
                            score = sc
                            best = s
                    return best

                best = get_best_match()
                if best is None:
                    raise ValueError(
                        "No spectrum found in database. Do not use "
                        + "autoretrieve='force'"
                    )
                else:
                    # Print comparison with best
                    print("Differences in best case :")
                    compare_dict(conditions, best.conditions)

                    raise ValueError(
                        "No spectrum found in database that matched "
                        + "given conditions. See best case found above"
                    )
            else:
                # just a print, then calculate
                if self.verbose:
                    print(
                        "No spectrum found in database that "
                        + "matched given conditions."
                    )
                return None

    def _build_partition_function_interpolator(
        self, parfunc, parfuncfmt, molecule, isotope
    ):
        """ Returns an universal partition function object ``parsum`` with the
        following methods defined::

            parsum.at(T)

        Partition functions are interpolated from tabulated values

        """

        if __debug__:
            printdbg(
                "called _build_partition_function_interpolator"
                + "(parfuncfmt={0}, isotope={1})".format(parfuncfmt, isotope)
            )

        isotope = int(isotope)

        # Use HAPI (HITRAN Python interface, integrated in RADIS)
        if parfuncfmt == "hapi":
            parsum = PartFuncHAPI(
                M=molecule, I=isotope, path=parfunc, verbose=self.verbose
            )
        elif parfuncfmt == "cdsd":  # Use tabulated CDSD partition functions
            assert molecule == "CO2"
            parsum = PartFuncCO2_CDSDtab(isotope, parfunc)
        elif parfuncfmt is None:
            # no tabulated partition functions defined. Only non-eq spectra can
            # be calculated if energies are also given
            parsum = None
        else:
            raise ValueError(
                "Unknown format for partition function: {0}".format(parfuncfmt)
            )
            # other formats ?

        return parsum

    def _build_partition_function_calculator(self, levels, levelsfmt, isotope):
        """ Return an universal partition function  object ``parsum`` so that the
        following methods are defined::

            parsum.at(T)
            parsum.at_noneq(Tvib, Trot)

        Partition functions are calculated from energy levels. Populations for
        all levels (independantly of the spectral range) can optionaly be
        calculated with argument ``update_populations=True``  (used to export
        populations of all states in Spectrum object)

        Parameters
        ----------

        levels: str
            energy levels filename

        levelsfmt: str
            energy levels format

        isotope: int
            isotope identifier

        """
        if __debug__:
            printdbg(
                "called _build_partition_function_calculator"
                + "(levelsfmt={0}, isotope={1})".format(levelsfmt, isotope)
            )

        isotope = int(isotope)

        # Forward to the correct partition function calculation engine:

        # ... sum over CDSD levels given by Tashkun / calculated from its Hamiltonian
        if levelsfmt in ["cdsd-pc", "cdsd-pcN", "cdsd-hamil"]:
            parsum = PartFuncCO2_CDSDcalc(
                levels,
                isotope=isotope,
                use_cached=self.params.lvl_use_cached,
                verbose=self.verbose,
                levelsfmt=levelsfmt,
            )

        # calculate energy levels from RADIS Dunham parameters
        elif levelsfmt == "radis":
            state = getMolecule(
                self.input.molecule, isotope, self.input.state, verbose=self.verbose
            )
            parsum = PartFunc_Dunham(
                state, use_cached=self.params.lvl_use_cached, verbose=self.verbose
            )
            # note: use 'levels' (useless here) to specify calculations options
            # for the abinitio calculation ? Like Jmax, etc.

        else:
            raise ValueError("Unknown format for energy levels : {0}".format(levelsfmt))
            # other formats ?

        return parsum

    def _fetch_molecular_parameters(self, df):
        """ Fetch molecular parameters (``molar_mass``, ``abundance``)  from
        Molecular Parameter database

        Parameters
        ----------

        df: pandas Dataframe
            line database with keys ``molar_mass``, ``abundance``

        Returns
        -------
        
        None:
            updates dataframe ``df`` directly

        See Also
        --------

        :class:`~radis.db.molparam.MolParams`
        """

        # order database
        # TODO: replace with attributes of Isotope>ElectronicState objects
        molpar = MolParams()

        if self.verbose >= 2:
            printg("... Fetching molecular parameters for all transitions")
            t0 = time()

        # prefill:

        # Use the fact that isotopes are int, and thus can be considered as
        # index in an array.
        # ... in the following we exploit this to use the np.take function,
        # ... which is amazingly fast
        # ... Read https://stackoverflow.com/a/51388828/5622825 to understand more
        # ... @dev: old versions: see radis <= 0.9.19
        id_set = df.id.unique()

        if len(id_set) == 1:
            id = list(id_set)[0]
            molecule = get_molecule(id)
            iso_set = self._get_isotope_list(molecule)  # df.iso.unique()

            # Shortcut if only 1 molecule & 1 isotope. We attribute molar_mass & abundance
            # as attributes of the line database, instead of columns. Much
            # faster!

            if len(iso_set) == 1:
                params = molpar.df.loc[(id, iso_set[0])]  # fetch all table directly
                df.Ia = params.abundance  # attribute, not column
                df.molar_mass = params.mol_mass  # attribute, not column
            #                # add in metadata so they follow when dataframe is copied/serialized
            #                for k in ['Ia', 'molar_mass']:
            #                    assert k not in df.columns
            #                    if k not in df._metadata:
            #                        df._metadata.append(k)

            # Else, parse for all isotopes. Use np.take that is very fast

            else:

                iso_arr = list(range(max(iso_set) + 1))

                Ia_arr = np.empty_like(iso_arr, dtype=np.float64)
                molarmass_arr = np.empty_like(iso_arr, dtype=np.float64)
                for iso in iso_arr:
                    if iso in iso_set:
                        params = molpar.df.loc[(id, iso)]  # fetch all table directly
                        # ... the trick below is that iso is used as index in the array
                        Ia_arr[iso] = params.abundance
                        molarmass_arr[iso] = params.mol_mass

                df["Ia"] = Ia_arr.take(df.iso)
                df["molar_mass"] = molarmass_arr.take(df.iso)

        elif len(id_set) == 0:

            raise ValueError("No molecule defined in this database.")

        else:
            raise NotImplementedError(
                ">1 molecule. Can use the np.take trick. Need to "
                + "fallback to pandas.map(dict)"
            )
            # TODO: Implement. Read https://stackoverflow.com/a/51388828/5622825 to understand more

        if self.verbose >= 2:
            printg("... Fetched molecular params in {0:.2f}s".format(time() - t0))

        return

    def get_partition_function_interpolator(self, molecule, isotope, elec_state):
        """ Retrieve Partition Function Interpolator

        Parameters
        ----------

        molecule: str

        isotope: int

        elec_state: str
        """

        parsum = self.parsum_tab[molecule][isotope][elec_state]

        # helps IDE find methods
        assert isinstance(parsum, RovibParFuncTabulator)

        return parsum

    def get_partition_function_calculator(self, molecule, isotope, elec_state):
        """ Retrieve Partition Function Calculator

        Parameters
        ----------

        molecule: str

        isotope: int

        elec_state: str
        """

        parsum = self.parsum_calc[molecule][isotope][elec_state]

        # helps IDE find methods
        assert isinstance(parsum, RovibParFuncCalculator)

        return parsum

    def _get_temp_file(self):
        """ Get temp file name (add warnings if temp_file exists already) """

        # Get temp file name
        return TEMP_FILE_PREFIX + str(self._id) + ".h5"

    def _clean_temp_file(self):
        """ Clean the room before leaving  """

        if exists(self._get_temp_file()):
            os.remove(self._get_temp_file())
            if self.verbose:
                print("Cleaned temp file: {0}".format(self._get_temp_file()))

    def get_conditions(self, ignore_misc=False):
        """ Get all parameters defined in the SpectrumFactory

        ignore_misc: boolean
            if ``True``, then all attributes considered as Factory 'descriptive' 
            parameters, as defined in :meth:`~radis.lbl.loader.get_conditions` are ignored when
            comparing the database to current factory conditions. It should
            obviously only be attributes that have no impact on the Spectrum
            produced by the factory. Default ``False`` 

        """

        vardict = self.input.get_params()
        vardict.update(self.params.get_params())
        if not ignore_misc:
            vardict.update(self.misc.get_params())

        return vardict

    def warn(self, message, category="default"):
        """ Trigger a warning, an error or just ignore based on the value defined
        in the :attr:`~radis.lbl.loader.DatabankLoader.warnings` dictionary

        The warnings can thus be deactivated selectively by setting the SpectrumFactory
         :attr:`~radis.lbl.loader.DatabankLoader.warnings` attribute

        Parameters
        ----------

        message: str
            what to print

        category: str
            one of the keys of self.warnings. See :py:attr:`~radis.lbl.loader.DatabankLoader.warnings`

        Notes
        -----

        All warnings in the SpectrumFactory should call to this method rather
        than the default warnings.warn() method, because it allows easier runtime
        modification of how to deal with warnings
        
        See Also
        --------

        :py:attr:`~radis.lbl.loader.DatabankLoader.warnings`

        """

        return warn(message, category=category, status=self.warnings)

    def __del__(self):
        """
        Note: No need to call the parent method. Python does it himself"""

        self._clean_temp_file()


if __name__ == "__main__":

    from radis.test.lbl.test_loader import _run_testcases

    print("Testing loading functions:", _run_testcases())
