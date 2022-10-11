# -*- coding: utf-8 -*-
"""
Summary
-------
Module to host the databank loading / database initialisation parts of
SpectrumFactory. This is done through :py:class:`~radis.lbl.factory.SpectrumFactory`
inheritance of the :py:class:`~radis.lbl.loader.DatabankLoader` class defined here

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

Most methods are written in inherited class with the following inheritance scheme:

:py:class:`~radis.lbl.loader.DatabankLoader` > :py:class:`~radis.lbl.base.BaseFactory` >
:py:class:`~radis.lbl.broadening.BroadenFactory` > :py:class:`~radis.lbl.bands.BandFactory` >
:py:class:`~radis.lbl.factory.SpectrumFactory`
.. inheritance-diagram:: radis.lbl.factory.SpectrumFactory
   :parts: 1

Notes
-----
RADIS includes automatic rebuilding of Deprecated cache files + a global variable
to force regenerating them after a given version. See ``"OLDEST_COMPATIBLE_VERSION"``
key in :py:attr:`radis.config`
-------------------------------------------------------------------------------
"""
# TODO: on use_cache functions, make a 'clean' / 'reset' option to delete / regenerate
# cache files

# @dev: (on Spyder IDE navigate between sections easily as # XXX makes a reference
# (on the slide bar on the right)

import warnings
from copy import deepcopy
from os.path import exists, expanduser, join, splitext
from time import time
from uuid import uuid1

import numpy as np
import pandas as pd

from radis import config
from radis.api.cdsdapi import cdsd2df
from radis.api.hdf5 import hdf2df
from radis.api.hitranapi import hit2df, parse_global_quanta, parse_local_quanta
from radis.api.tools import drop_object_format_columns, replace_PQR_with_m101
from radis.db.classes import get_molecule
from radis.db.molecules import getMolecule
from radis.db.molparam import MOLPARAMS_EXTRA_PATH, MolParams
from radis.db.references import doi
from radis.io.exomol import fetch_exomol
from radis.io.geisa import fetch_geisa
from radis.io.hitemp import fetch_hitemp
from radis.io.hitran import fetch_hitran
from radis.io.query import fetch_astroquery
from radis.levels.partfunc import (
    PartFunc_Dunham,
    PartFuncTIPS,
    RovibParFuncCalculator,
    RovibParFuncTabulator,
)
from radis.levels.partfunc_cdsd import PartFuncCO2_CDSDcalc, PartFuncCO2_CDSDtab
from radis.misc.arrays import count_nans
from radis.misc.basics import compare_dict, compare_lists
from radis.misc.config import getDatabankEntries, getDatabankList, printDatabankEntries
from radis.misc.debug import printdbg
from radis.misc.log import printwarn
from radis.misc.printer import printg
from radis.misc.profiler import Profiler
from radis.misc.utils import get_files_from_regex
from radis.misc.warning import (
    EmptyDatabaseError,
    IrrelevantFileWarning,
    default_warning_status,
    warn,
)
from radis.phys.convert import cm2nm
from radis.tools.database import SpecDatabase
from radis.tools.track_ref import RefTracker

KNOWN_DBFORMAT = [
    "hitran",
    "hitemp",
    "cdsd-hitemp",
    "cdsd-4000",
    "hitemp-radisdb",
    "hdf5-radisdb",
    "geisa",
]
"""list: Known formats for Line Databases:

- ``'hitran'`` : [HITRAN-2020]_ original .par format.
- ``'hitemp'`` : [HITEMP-2010]_ original format (same format as 'hitran').
- ``'cdsd-hitemp'`` : CDSD-HITEMP original format (CO2 only, same lines as HITEMP-2010).
- ``'cdsd-4000'`` : [CDSD-4000]_ original format (CO2 only).
- ``'hitemp-radisdb'`` : HITEMP under RADISDB format (pytables-HDF5 with RADIS column names).
- ``'hdf5-radisdb'`` : arbitrary HDF5 file with RADIS column names.
- ``'geisa'`` : [GEISA-2020]_ original .par format.

To install all databases manually see the :ref:`Configuration file <label_lbl_config_file>`
and the :ref:`list of databases <label_line_databases>` .

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
    "hitran": ["ierr", "iref", "lmix", "gpp"],
    "hitemp": ["ierr", "iref", "lmix", "gpp"],
    "cdsd-4000": ["wang2"],
    "cdsd-hitemp": ["wang2", "lsrc"],
    "hdf5-radisdb": [],
    "hitemp-radisdb": [],
    "geisa": [],
}
""" dict: drop these columns if using ``drop_columns='auto'`` in load_databank
Based on the value of ``dbformat=``, some of these columns won't be used.

See Also
--------
- 'hitran': (HITRAN / HITEMP) :data:`~radis.api.hitranapi.columns_2004`,
- 'cdsd-hitemp' (CDSD HITEMP): :data:`~radis.api.cdsdapi.columns_hitemp`,
- 'cdsd-4000': (CDSD 4000) :data:`~radis.api.cdsdapi.columns_4000`,
- 'geisa': (GEISA 2020) :data:`~radis.api.geisaapi.columns_GEISA`,
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
- 'radis': :data:`~radis.api.hitranapi.columns_2004`,
- 'cdsd-pc': :data:`~radis.api.hitranapi.columns_2004`,
- 'cdsd-pcN' (CDSD-HITEMP): :data:`~radis.api.cdsdapi.columns_hitemp`,
- 'cdsd-hamil': :data:`~radis.api.cdsdapi.columns_4000`,
"""
# TODO @dev : switch from a model where we drop certain useless columns (RADIS==0.9.28)
# to a model where we only-load the required ones initially (if possible with lazy-loading,
# i.e. only load them on demand. See https://github.com/radis/radis/issues/118 )
drop_all_but_these = [
    "id",
    "iso",
    "wav",
    "int",
    "A",
    "airbrd",
    "selbrd",
    "Tdpair",
    "Tdpsel",
    "Pshft",
    "El",
    "gp",
]
""" list: drop all columns but these if using ``drop_columns='all'`` in load_databank
Note: nonequilibrium calculations wont be possible anymore and it wont be possible
to identify lines with :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey`

See Also
--------
- 'hitran': (HITRAN / HITEMP) :data:`~radis.api.hitranapi.columns_2004`,
- 'cdsd-hitemp' (CDSD HITEMP): :data:`~radis.api.cdsdapi.columns_hitemp`,
- 'cdsd-4000': (CDSD 4000) :data:`~radis.api.cdsdapi.columns_4000`,
"""
required_non_eq = [
    "branch",
    "jl",
    "vl",
    "vu",
    "v1l",
    "v2l",
    "l2l",
    "v3l",
    "v1u",
    "v2u",
    "l2u",
    "v3u",
    "polyl",
    "wangl",
    "rankl",
    "polyu",
    "wangu",
    "ranku",
]
"""list: column names required for non-equilibrium calculations.
See load_column= key of fetch_databank() and load_databank() """

broadening_coeff = [
    "gamma_co2",
    "n_co2",
    "gamma_h2o",
    "n_h2o",
    "gamma_he",
    "n_he",
    "gamma_h2",
    "n_h2",
]
"""list: column names required for non-air diluent calculations."""
# TODO refactor : directly go & parse the identifications names (globu, locu, etc.)
# in radis.io.hitran ?
# For the moment we just try to be exhaustive

# @dev: Sanity checks
# (make sure all variables are defined everywhere)
assert compare_lists(drop_auto_columns_for_dbformat, KNOWN_DBFORMAT) == 1
assert compare_lists(drop_auto_columns_for_levelsfmt, KNOWN_LVLFORMAT) == 1

# %% Main class


class ConditionDict(dict):
    """A class to hold Spectrum calculation input conditions
    (:py:class:`~radis.lbl.loader.Input`), computation parameters
    (:py:class:`~radis.lbl.loader.Parameters`), or miscalleneous parameters
    (:py:class:`~radis.lbl.loader.MiscParams`).
    Works like a dict except you can also access attribute with::
        v = a.key   # equivalent to v = a[key]
    Also can be copied, deepcopied, and parallelized in multiprocessing

    Notes
    -----
    for developers:
    Parameters and Input could also have simply derived from the (object) class,
    but it may have missed some convenients functions implemented for dict.
    For instance, how to be picked / unpickled.

    See Also
    --------
    :py:class:`~radis.lbl.loader.Input`,
    :py:class:`~radis.lbl.loader.Parameters`,
    :py:class:`~radis.lbl.loader.MiscParams`
    """

    def get_params(self):
        """Returns the variables (and their values) contained in the
        dictionary, minus some based on their type. Numpy array, dictionaries
        and pandas DataFrame are removed. None is removed in general, except
        for some keys ('cutoff', 'truncation')
        Tuples are converted to string
        """

        # Filter parameters based on type
        def filter_type(params):
            filt_params = {}
            for k, v in params.items():
                # Ignore some
                if type(v) in [np.ndarray, dict, pd.DataFrame]:
                    continue
                if type(v) is None and k not in ["cutoff", "truncation"]:
                    continue
                if isinstance(k, str):
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
        if not attr in self.__slots__:
            raise KeyError(
                f"Undefined attribute `{attr}` for {self.__class__}. Allowed attributes: {self.__slots__}"
            )
        return self[attr]

    def __setattr__(self, attr, value):
        if not attr in self.__slots__:
            raise KeyError(
                f"Undefined attribute `{attr}` for {self.__class__}. Allowed attributes: {self.__slots__}"
            )
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
    """Holds Spectrum calculation input conditions, under the attribute
    :py:attr:`~radis.lbl.loader.DatabankLoader.input` of
    :py:class:`~radis.lbl.factory.SpectrumFactory`.
    Works like a dict except you can also access attribute with::
        v = sf.input.key   # equivalent to v = sf.input[key]

    See Also
    --------
    :py:attr:`~radis.lbl.loader.DatabankLoader.params`,
    :py:attr:`~radis.lbl.loader.DatabankLoader.misc`
    """

    # hardcode attribute names, to prevent typos and the declaration of unwanted parameters
    __slots__ = [
        "Tgas",
        "Tref",
        "Tvib",
        "Trot",
        "isotope",
        "medium",
        "mole_fraction",
        "molecule",
        "overpopulation",
        "path_length",
        "pressure_mbar",
        "rot_distribution",
        "self_absorption",
        "state",
        "vib_distribution",
        "wavelength_max",
        "wavelength_min",
        "wavenum_max",
        "wavenum_min",
    ]

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


# TO-DO: these error estimations are horribly outdated...
def _lorentzian_step(res_L):
    log_pL = np.log((res_L / 0.20) ** 0.5 + 1)
    return log_pL


def _gaussian_step(res_G):
    log_pG = np.log((res_G / 0.46) ** 0.5 + 1)
    return log_pG


# class Parameters(object):
class Parameters(ConditionDict):
    """Holds Spectrum calculation computation parameters, under the attribute
    :py:attr:`~radis.lbl.loader.DatabankLoader.params` of
    :py:class:`~radis.lbl.factory.SpectrumFactory`.
    Works like
    a dict except you can also access attribute with::
        v = sf.params.key    # equivalent to v = sf.params[key]
    Also can be copied, deepcopied, and parallelized in multiprocessing

    See Also
    --------
    :py:attr:`~radis.lbl.loader.DatabankLoader.input`,
    :py:attr:`~radis.lbl.loader.DatabankLoader.misc`
    """

    # hardcode attribute names, to prevent typos and the declaration of unwanted parameters
    __slots__ = [
        # "add_at_used",
        "broadening_method",
        "truncation",
        "neighbour_lines",
        "chunksize",
        "cutoff",
        "db_use_cached",
        "dbformat",
        "dbpath",
        "dxL",
        "dxG",
        "export_lines",
        "export_populations",
        "folding_thresh",
        "include_neighbouring_lines",
        "levelsfmt",
        "lvl_use_cached",
        "optimization",
        "parfuncfmt",
        "parfuncpath",
        "parsum_mode",
        "pseudo_continuum_threshold",
        "sparse_ldm",
        "warning_broadening_threshold",
        "warning_linestrength_cutoff",
        "wavenum_max_calc",
        "wavenum_min_calc",
        "waveunit",
        "wstep",
        "diluent",
    ]

    def __init__(self):
        super(Parameters, self).__init__()

        # Dev: Init here to be found by autocomplete
        self.truncation = None  #: float: cutoff for half-width lineshape calculation (cm-1). Overwritten by SpectrumFactory
        self.neighbour_lines = None  #: float: extra range (cm-1) on each side of the spectrum to account for neighbouring lines. Overwritten by SpectrumFactory
        self.cutoff = None  #: float: linestrength cutoff (molecule/cm)
        self.broadening_method = ""  #: str:``"voigt"``, ``"convolve"``, ``"fft"``
        self.optimization = None  #: str: ``"simple"``, ``"min-RMS"``, ``None``
        self.db_use_cached = (
            None  #: bool: use (and generate) cache files for Line Database
        )
        self.dbformat = None  #: str: format of Line Database. See :data:`~radis.lbl.loader.KNOWN_DBFORMAT`
        self.dbpath = None  #: str: joined list of filepaths to Line Database
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
        self.diluent = {}  # dict: molecule : mole fraction
        self.dxL = _lorentzian_step(
            0.01
        )  #: float : Lorentzian step for LDM lineshape database. Default _lorentzian_step(0.01)
        self.dxG = _gaussian_step(
            0.01
        )  #: float : Gaussian step LDM lineshape database. Default _gaussian_step(0.01)
        # self.add_at_used = None  # use Cython-accelerated code
        self.include_neighbouring_lines = True
        """bool: if ``True``, includes the contribution of off-range, neighbouring
        lines because of lineshape broadening. Default ``True``."""
        self.parsum_mode = "full summation"  #: int : "full summation" or "tabulation"  . calculation mode of parittion function. See :py:class:`~radis.levels.partfunc.RovibParFuncCalculator`
        self.sparse_ldm = "auto"  #: str: "auto", True, False  . Sparse LDM calculation. See :py:meth:`radis.lbl.broadening.BroadenFactory._apply_lineshape_LDM`


class MiscParams(ConditionDict):
    """A class to hold Spectrum calculation descriptive parameters, under the attribute
    :py:attr:`~radis.lbl.loader.DatabankLoader.params` of
    :py:class:`~radis.lbl.factory.SpectrumFactory`.
    Unlike :class:`~radis.lbl.loader.Parameters`, these parameters cannot influence the
    Spectrum output and will not be used when comparing Spectrum with existing,
    precomputed spectra in :class:`~radis.tools.database.SpecDatabase`
    Works like
    a dict except you can also access attribute with::
        v = a.key

    See Also
    --------
    :py:attr:`~radis.lbl.loader.DatabankLoader.input`,
    :py:attr:`~radis.lbl.loader.DatabankLoader.params`,
    """

    # hardcode attribute names, to prevent typos and the declaration of unwanted parameters
    __slots__ = [
        "chunksize",
        "export_lines",
        "export_populations",
        "export_rovib_fraction",
        "load_energies",
        "warning_broadening_threshold",
        "warning_linestrength_cutoff",
        "total_lines",
        "zero_padding",
        "memory_mapping_engine",
        "add_at_used",  # function used in DIT ; a Cython and a pure-Python version exist
    ]

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
        self.export_rovib_fraction = False  #: bool: calculate nu_vib, nu_rot in lines
        self.warning_broadening_threshold = (
            None  #: float: [0-1] raise a warning if the lineshape area is different
        )
        self.warning_linestrength_cutoff = None  #: float [0-1]: raise a warning if the sum of linestrength cut is above that
        self.total_lines = 0  #: int : number of lines in database.
        self.memory_mapping_engine = config[
            "MEMORY_MAPPING_ENGINE"
        ]  # 'pytables', 'vaex', 'feather'

        self.add_at_used = (
            ""  # function used in DIT ; a Cython and a pure-Python version exist
        )


def format_paths(s):
    """escape all special characters."""
    if s is not None:
        s = str(s).replace("\\", "/")
    return s


df_metadata = ["molecule", "iso", "id", "Ia", "molar_mass", "Qref", "Qvib", "Q"]
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
        self.params = Parameters()
        """Computational parameters: :py:class:`~radis.lbl.loader.Parameters`
        they may change the output of calculations (ex: threshold, cutoff, broadening methods, etc.)
        """
        self.misc = MiscParams()
        """Miscelleneous parameters (:py:class:`~radis.lbl.loader.MiscParams`)
        params that cannot change the output of calculations (ex: number of CPU, etc.)
        """
        # Setup individual warnings. Value of keys can be:
        # - 'warning' (default: just trigger a warning)
        # - 'error' (raises an error on this warning)
        # - 'ignore'  (do nothing)
        # The key self.warnings['default'] will set the warning behavior for all
        # other warnings
        self.warnings = default_warning_status.copy()
        """ dict: Default warnings for SpectrumFactory. See
        :py:data:`~radis.misc.warnings.default_warning_status`"""

        # Generate unique id for Factory
        self._id = uuid1()

        # Init Annotations (Python >3.6) [hints for users]
        try:
            self.load_databank.__annotations__["format"] = KNOWN_DBFORMAT
            self.load_databank.__annotations__["levelsfmt"] = KNOWN_LVLFORMAT
            self.load_databank.__annotations__["parfuncfmt"] = KNOWN_PARFUNCFORMAT
        except:  # old Python version
            pass

        # Variables that will hold the dataframes.
        self.df0 = None  # type : pd.DataFrame
        """pd.DataFrame : initial line database after loading.
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
        self.df1 = None  # type : pd.DataFrame
        """pd.DataFrame : line database, scaled with populations + linestrength cutoff
        Never edit manually. See all comments about :py:attr:`~self.radis.lbl.loader.DatabankLoader.df0`

        See Also
        --------
        :py:attr:`~self.radis.lbl.loader.DatabankLoader.df0`
        """

        # Temp variable to store databanks information
        self._databank_args = []
        self._databank_kwargs = {}

        self._autoretrieveignoreconditions = []  # HACK. See _retrieve_from_database

        # Molecular parameters
        self.molparam = MolParams(extra_file_json=MOLPARAMS_EXTRA_PATH)
        """MolParam: contains information about molar mass; isotopic abundance.

        See :py:class:`~radis.db.molparam.MolParams`"""
        # TODO @dev : Refactor : turn it into a Dictinoary? (easier to store as JSON Etc.)

        # Profiler
        self.profiler = None

    def _reset_profiler(self, verbose):
        """Reset :py:class:`~radis.misc.profiler.Profiler`

        See Also
        --------
        :py:func:`radis.lbl.factory.SpectrumFactory.print_perf_profile"""

        self.profiler = Profiler(verbose)

    def _reset_references(self):
        """Reset :py:class:`~radis.tools.track_refs.RefTracker`"""

        # Track bibliography references
        self.reftracker = RefTracker()
        # ... init with RADIS itself:
        self.reftracker.add(doi["RADIS-2018"], "calculation")

    # %% ======================================================================
    # PUBLIC METHODS
    # ------------------------
    # init_databank           >>> init line database (= store reference but load() later)
    # load_databank           >>> load line database
    # init_database           >>> to interact / generate a SpectrumDatabase
    #
    # =========================================================================

    def init_databank(self, *args, **kwargs):
        """Method to init databank parameters but only load them when needed.

        Databank is reloaded by :py:meth:`~radis.lbl.loader.DatabankLoader._check_line_databank`
        Same inputs Parameters as :meth:`~radis.lbl.loader.DatabankLoader.load_databank`:

        Parameters
        ----------
        name: a section name specified in your ``~/radis.json``
            ``.radis`` has to be created in your HOME (Unix) / User (Windows). If
            not ``None``, all other arguments are discarded.
            Note that all files in database will be loaded and it may takes some
            time. Better limit the database size if you already know what
            range you need. See :ref:`Configuration file <label_lbl_config_file>` and
            :data:`~radis.misc.config.DBFORMAT` for expected
            ``~/radis.json`` format

        Other Parameters
        ----------------
        path: str, list of str, None
            list of database files, or name of a predefined database in the
            :ref:`Configuration file <label_lbl_config_file>` (`~/radis.json`)
            Accepts wildcards ``*`` to select multiple files
        format: ``'hitran'``, ``'cdsd-hitemp'``, ``'cdsd-4000'``, or any of :data:`~radis.lblinit_databank.loader.KNOWN_DBFORMAT`
            database type. ``'hitran'`` for HITRAN/HITEMP, ``'cdsd-hitemp'``
            and ``'cdsd-4000'`` for the different CDSD versions. Default ``'hitran'``
        parfuncfmt: ``'hapi'``, ``'cdsd'``, or any of :data:`~radis.lbl.loader.KNOWN_PARFUNCFORMAT`
            format to read tabulated partition function file. If ``hapi``, then
            HAPI (HITRAN Python interface) [1]_ is used to retrieve them (valid if
            your database is HITRAN data). HAPI is embedded into RADIS. Check the
            version. If partfuncfmt is None then ``hapi`` is used. Default ``hapi``.
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
            If ``None``, non equilibrium calculations are not possible. Default ``'radis'``.
        db_use_cached: boolean, or ``None``
            if ``True``, a pandas-readable csv file is generated on first access,
            and later used. This saves on the datatype cast and conversion and
            improves performances a lot. But! ... be sure to delete these files
            to regenerate them if you happen to change the database. If ``'regen'``,
            existing cached files are removed and regenerated.
            It is also used to load energy levels from ``.h5`` cache file if exist.
            If ``None``, the value given on Factory creation is used. Default ``None``
        load_energies: boolean
            if ``False``, dont load energy levels. This means that nonequilibrium
            spectra cannot be calculated, but it saves some memory. Default ``True``
        include_neighbouring_lines: bool
            ``True``, includes off-range, neighbouring lines that contribute
            because of lineshape broadening. The ``neighbour_lines``
            parameter is used to determine the limit. Default ``True``.
        drop_columns: list
            columns names to drop from Line DataFrame after loading the file.
            Not recommended to use, unless you explicitely want to drop information
            (for instance if dealing with too large databases). If ``[]``, nothing
            is dropped. If ``'auto'``, parameters considered unnecessary
            are dropped. See :data:`~radis.lbl.loader.drop_auto_columns_for_dbformat`
            and :data:`~radis.lbl.loader.drop_auto_columns_for_levelsfmt`.
            Default ``'auto'``.
        load_columns: list, ``'all'``, ``'equilibrium'``, ``'noneq'``
            columns names to load.
            If ``'equilibrium'``, only load the columns required for equilibrium
            calculations. If ``'noneq'``, also load the columns required for
            non-LTE calculations. See :data:`~radis.lbl.loader.drop_all_but_these`.
            If ``'all'``, load everything. Note that for performances, it is
            better to load only certain columsn rather than loading them all
            and dropping them with ``drop_columns``.
            Default ``'equilibrium'``.

            .. warning::
                if using ``'equilibrium'``, not all parameters will be available
                for a Spectrum :py:func:`~radis.spectrum.spectrum.Spectrum.line_survey`.

        *Other arguments are related to how to open the files*

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
        # TODO : refactor drop_columns/load_columns

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
            lvl_use_cached,
            drop_columns,
            load_columns,
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
            lvl_use_cached=lvl_use_cached,
            load_energies=load_energies,
            include_neighbouring_lines=include_neighbouring_lines,
        )

        # Delete database
        self.df0 = None  # type : pd.DataFrame
        self._reset_references()  # bibliographic references

    def columns_list_to_load(self, load_columns_type):
        # Which columns to load
        if load_columns_type == "equilibrium":
            columns = list(drop_all_but_these)
        elif load_columns_type == "noneq":
            columns = list(set(drop_all_but_these) | set(required_non_eq))
        elif load_columns_type == "diluent":
            columns = list(broadening_coeff)
        else:
            raise ValueError(
                f"Expected a list or 'all' for `load_columns`, got `load_columns={load_columns_type}"
            )
        return columns

    def fetch_databank(
        self,
        source="hitran",
        database="default",
        parfunc=None,
        parfuncfmt="hapi",
        levels=None,
        levelsfmt="radis",
        load_energies=False,
        include_neighbouring_lines=True,
        parse_local_global_quanta=True,
        drop_non_numeric=True,
        db_use_cached=True,
        lvl_use_cached=True,
        memory_mapping_engine="default",
        load_columns="equilibrium",
        parallel=True,
        extra_params=None,
    ):
        """Fetch the latest files from [HITRAN-2020]_, [HITEMP-2010]_ (or newer),
        [ExoMol-2020]_  or [GEISA-2020] , and store them locally in memory-mapping
        formats for extremelly fast access.

        Parameters
        ----------
        source: ``'hitran'``, ``'hitemp'``, ``'exomol'``, ``'geisa'``
            which database to use.
        database: ``'full'``, ``'range'``, name of an ExoMol database, or ``'default'``
            if fetching from HITRAN, ``'full'`` download the full database and register
            it, ``'range'`` download only the lines in the range of the molecule.

            .. note::
                ``'range'`` will be faster, but will require a new download each time
                you'll change the range. ``'full'`` is slower and takes more memory, but
                will be downloaded only once.

            Default is ``'full'``.

            If fetching from HITEMP, only ``'full'`` is available.

            if fetching from ''`exomol`'', use this parameter to choose which database
            to use. Keep ``'default'`` to use the recommended one. See all available databases
            with :py:func:`radis.io.exomol.get_exomol_database_list`

            By default, databases are download in ``~/.radisdb``.
            Can be changed in ``radis.config["DEFAULT_DOWNLOAD_PATH"]`` or in
            ``~/radis.json`` config file


        Other Parameters
        ----------------
        parfuncfmt: ``'cdsd'``, ``'hapi'``, ``'exomol'``, or any of :data:`~radis.lbl.loader.KNOWN_PARFUNCFORMAT`
            format to read tabulated partition function file. If ``'hapi'``, then
            [HAPI]_ (HITRAN Python interface) is used to retrieve [TIPS-2020]_
            tabulated partition functions.
            If ``'exomol'`` then partition functions are downloaded from ExoMol.
            Default ``'hapi'``.
        parfunc: filename or None
            path to a tabulated partition function file to use.
        levels: dict of str or None
            path to energy levels (needed for non-eq calculations). Format::

                {1:path_to_levels_iso_1, 3:path_to_levels_iso3}.

            Default ``None``
        levelsfmt: ``'cdsd-pc'``, ``'radis'`` (or any of :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`) or ``None``
            how to read the previous file. Known formats: (see :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`).
            If ``radis``, energies are calculated using the diatomic constants in radis.db database
            if available for given molecule. Look up references there.
            If ``None``, non equilibrium calculations are not possible. Default ``'radis'``.
        load_energies: boolean
            if ``False``, dont load energy levels. This means that nonequilibrium
            spectra cannot be calculated, but it saves some memory. Default ``False``
        include_neighbouring_lines: bool
            if ``True``, includes off-range, neighbouring lines that contribute
            because of lineshape broadening. The ``neighbour_lines``
            parameter is used to determine the limit. Default ``True``.
        parse_local_global_quanta: bool, or ``'auto'``
            if ``True``, parses the HITRAN/HITEMP 'glob' and 'loc' columns to extract
            quanta identifying the lines. Required for nonequilibrium calculations,
            or to use :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey`,
            but takes up more space.
        drop_non_numeric: boolean
            if ``True``, non numeric columns are dropped. This improves performances,
            but make sure all the columns you need are converted to numeric formats
            before hand. Default ``True``. Note that if a cache file is loaded it
            will be left untouched.
        db_use_cached: bool, or ``'regen'``
            use cached
        memory_mapping_engine: ``'pytables'``, ``'vaex'``, ``'feather'``
            which library to use to read HDF5 files (they are incompatible: ``'pytables'`` is
            row-major while ``'vaex'`` is column-major) or other memory-mapping formats
            If ``'default'``, use the value from ~/radis.json `["MEMORY_MAPPING_ENGINE"]`
        parallel: bool
            if ``True``, uses joblib.parallel to load database with multiple processes
            (works only for HITEMP files)
        load_columns: list, ``'all'``, ``'equilibrium'``, ``'noneq'``, ``diluent``,
            columns names to load.
            If ``'equilibrium'``, only load the columns required for equilibrium
            calculations. If ``'noneq'``, also load the columns required for
            non-LTE calculations. See :data:`~radis.lbl.loader.drop_all_but_these`.
            If ``'all'``, load everything. Note that for performances, it is
            better to load only certain columsn rather than loading them all
            and dropping them with ``drop_columns``.
            If ``diluent`` then all additional columns required for calculating spectrum
            in that diluent is loaded.
            Default ``'equilibrium'``.

            .. warning::
                if using ``'equilibrium'``, not all parameters will be available
                for a Spectrum :py:func:`~radis.spectrum.spectrum.Spectrum.line_survey`.
                If you are calculating equilibrium (LTE) spectra, it is recommended to
                use ``'equilibrium'``. If you are calculating non-LTE spectra, it is
                recommended to use ``'noneq'``.

        Notes
        -----
        HITRAN is fetched with Astroquery [1]_ or [HAPI]_,  and HITEMP with
        :py:func:`~radis.io.hitemp.fetch_hitemp`
        HITEMP files are generated in a ~/.radisdb database.

        See Also
        --------
        :meth:`~radis.lbl.loader.DatabankLoader.load_databank`,
        :meth:`~radis.lbl.loader.DatabankLoader.init_databank`

        References
        ----------
        .. [1] `Astroquery <https://astroquery.readthedocs.io>`_
        """
        # @dev TODO: also add cache file to fetch_databank, similar to load_databank
        # | Should store the waverange, molecule and isotopes in the cache file
        # | metadata to ensures that it is redownloaded if necessary.
        # | see implementation in load_databank.

        # Check inputs
        if source == "astroquery":
            warnings.warn(
                DeprecationWarning(
                    "source='astroquery' replaced with source='hitran' in 0.9.28"
                )
            )
            source = "hitran"
        if source not in ["hitran", "hitemp", "exomol", "geisa"]:
            raise NotImplementedError("source: {0}".format(source))
        if source == "hitran":
            dbformat = "hitran"
            if database == "default":
                database = "full"
        elif source == "hitemp":
            dbformat = (
                "hitemp-radisdb"  # downloaded in RADIS local databases ~/.radisdb
            )
            if database == "default":
                database = "full"
        elif source == "exomol":
            dbformat = (
                "exomol-radisdb"  # downloaded in RADIS local databases ~/.radisdb
            )
        elif source == "geisa":
            dbformat = "geisa"
            if database == "default":
                database = "full"

        local_databases = config["DEFAULT_DOWNLOAD_PATH"]

        if [parfuncfmt, source].count("exomol") == 1:
            self.warn(
                f"Using lines from {source} but partition functions from {parfuncfmt}"
                + "for consistency we recommend using lines and partition functions from the same database",
                "AccuracyWarning",
            )
        if memory_mapping_engine == "default":
            memory_mapping_engine = self.misc.memory_mapping_engine

        # Get inputs
        molecule = self.input.molecule
        isotope = self.input.isotope
        if not molecule:
            raise ValueError("Please define `molecule=` so the database can be fetched")

        if include_neighbouring_lines:
            wavenum_min = self.params.wavenum_min_calc
            wavenum_max = self.params.wavenum_max_calc
        else:
            wavenum_min = self.input.wavenum_min
            wavenum_max = self.input.wavenum_max

        # Let's store all params so they can be parsed by "get_conditions()"
        # and saved in output spectra information
        self.params.dbformat = dbformat
        self.misc.load_energies = load_energies
        self.levels = levels
        if levels is not None:
            self.levelspath = ",".join([format_paths(lvl) for lvl in levels.values()])
        else:
            self.levelspath = None
        self.params.levelsfmt = levelsfmt
        self.params.parfuncpath = format_paths(parfunc)
        self.params.parfuncfmt = parfuncfmt
        self.params.db_use_cached = db_use_cached
        self.params.lvl_use_cached = lvl_use_cached

        # Which columns to load
        columns = []
        if "all" in load_columns:
            columns = None  # see fetch_hitemp, fetch_hitran, etc.
        elif isinstance(load_columns, str) and load_columns in ["equilibrium", "noneq"]:
            columns = self.columns_list_to_load(load_columns)
        elif load_columns == "diluent":
            raise ValueError(
                "Please use diluent along with 'equilibrium' or 'noneq' in a list like ['diluent','noneq']"
            )

        elif isinstance(load_columns, list) and "all" not in load_columns:
            for load_columns_type in load_columns:
                if load_columns_type in ["equilibrium", "noneq", "diluent"]:
                    for col in self.columns_list_to_load(load_columns_type):
                        columns.append(col)
                elif load_columns_type in list(
                    set(drop_all_but_these)
                    | set(required_non_eq)
                    | set(broadening_coeff)
                ):
                    columns.append(load_columns_type)
                else:
                    raise ValueError("invalid column name provided")
            columns = list(set(columns))

        # %% Init Line database
        # ---------------------
        self._reset_references()  # bibliographic references

        if source == "hitran":
            self.reftracker.add(doi["HITRAN-2020"], "line database")  # [HITRAN-2020]_

            if database == "full":
                self.reftracker.add(doi["HAPI"], "data retrieval")  # [HAPI]_

                if memory_mapping_engine == "auto":
                    memory_mapping_engine = "vaex"

                if isotope == "all":
                    isotope_list = None
                else:
                    isotope_list = ",".join([str(k) for k in self._get_isotope_list()])

                df, local_paths = fetch_hitran(
                    molecule,
                    isotope=isotope_list,
                    local_databases=join(local_databases, "hitran"),
                    load_wavenum_min=wavenum_min,
                    load_wavenum_max=wavenum_max,
                    columns=columns,
                    cache=db_use_cached,
                    verbose=self.verbose,
                    return_local_path=True,
                    engine=memory_mapping_engine,
                    parallel=parallel,
                    extra_params=extra_params,
                )
                self.params.dbpath = ",".join(local_paths)

                # ... explicitely write all isotopes based on isotopes found in the database
                if isotope == "all":
                    self.input.isotope = ",".join(
                        [str(k) for k in self._get_isotope_list(df=df)]
                    )

            elif database == "range":
                self.reftracker.add(
                    doi["Astroquery"], "data retrieval"
                )  # [Astroquery]_

                # Query one isotope at a time
                if isotope == "all":
                    raise ValueError(
                        "Please define isotope explicitely (cannot use 'all' with fetch_databank('hitran'))"
                    )
                isotope_list = self._get_isotope_list()

                frames = []  # lines for all isotopes
                for iso in isotope_list:
                    df = fetch_astroquery(
                        molecule, iso, wavenum_min, wavenum_max, verbose=self.verbose
                    )
                    if len(df) > 0:
                        frames.append(df)
                    else:
                        self.warn(
                            "No line for isotope n°{}".format(iso),
                            "EmptyDatabaseWarning",
                            level=2,
                        )

                # Merge
                if frames == []:
                    raise EmptyDatabaseError(
                        f"{molecule} has no lines on range "
                        + "{0:.2f}-{1:.2f} cm-1".format(wavenum_min, wavenum_max)
                    )
                if len(frames) > 1:
                    # Note @dev : may be faster/less memory hungry to keep lines separated for each isotope. TODO : test both versions
                    for df in frames:
                        assert "iso" in df.columns
                    df = pd.concat(frames, ignore_index=True)  # reindex
                else:
                    df = frames[0]
                self.params.dbpath = "fetched from hitran"
            else:
                raise ValueError(
                    f"Got `database={database}`. When fetching HITRAN, choose `database='full'` to download all database (once for all) or `database='range'` to download only the lines in the current range."
                )

        elif source == "hitemp":
            self.reftracker.add(doi["HITEMP-2010"], "line database")  # [HITEMP-2010]_

            if memory_mapping_engine == "auto":
                memory_mapping_engine = "vaex"

            if database != "full":
                raise ValueError(
                    f"Got `database={database}`. When fetching HITEMP, only the `database='full'` option is available."
                )

            # Download, setup local databases, and fetch (use existing if possible)

            if isotope == "all":
                isotope_list = None
            else:
                isotope_list = ",".join([str(k) for k in self._get_isotope_list()])

            df, local_paths = fetch_hitemp(
                molecule,
                isotope=isotope_list,
                local_databases=join(local_databases, "hitemp"),
                load_wavenum_min=wavenum_min,
                load_wavenum_max=wavenum_max,
                columns=columns,
                cache=db_use_cached,
                verbose=self.verbose,
                return_local_path=True,
                engine=memory_mapping_engine,
                parallel=parallel,
            )
            self.params.dbpath = ",".join(local_paths)

            # ... explicitely write all isotopes based on isotopes found in the database
            if isotope == "all":
                self.input.isotope = ",".join(
                    [str(k) for k in self._get_isotope_list(df=df)]
                )

        elif source == "exomol":
            self.reftracker.add(doi["ExoMol-2020"], "line database")  # [ExoMol-2020]

            if memory_mapping_engine == "auto":
                memory_mapping_engine = "vaex"

            if database in ["full", "range"]:
                raise ValueError(
                    f"Got `database={database}`. When fetching ExoMol, use the `database=` key to retrieve a specific database. Use `database='default'` to get the recommended database. See more informatino in radis.io.fetch_exomol()"
                )

            # Download, setup local databases, and fetch (use existing if possible)
            if memory_mapping_engine not in ["vaex", "feather", "pytables"]:
                raise NotImplementedError(
                    f"{memory_mapping_engine} with ExoMol files. Define radis.config['MEMORY_MAPPING_ENGINE'] = 'vaex' or 'feather'"
                )

            if isotope == "all":
                raise ValueError(
                    "Please define isotope explicitely (cannot use 'all' with fetch_databank('exomol'))"
                )
            isotope_list = self._get_isotope_list()

            local_paths = []
            frames = []  # lines for all isotopes
            partition_function_exomol = {
                molecule: {}
            }  # partition function tabulators for all isotpes
            for iso in isotope_list:
                df, local_path, Z_exomol = fetch_exomol(
                    molecule,
                    database=database,
                    isotope=iso,
                    local_databases=join(local_databases, "exomol"),
                    load_wavenum_min=wavenum_min,
                    load_wavenum_max=wavenum_max,
                    columns=columns,
                    cache=db_use_cached,
                    verbose=self.verbose,
                    return_local_path=True,
                    return_partition_function=True,
                    engine=memory_mapping_engine,
                )
                # @dev refactor : have a DatabaseClass from which we load lines and partition functions
                if len(df) > 0:
                    frames.append(df)
                local_paths.append(local_path)
                partition_function_exomol[molecule][iso] = Z_exomol

            # Merge
            if frames == []:
                raise EmptyDatabaseError(
                    f"{molecule} has no lines on range "
                    + "{0:.2f}-{1:.2f} cm-1".format(wavenum_min, wavenum_max)
                )
            if len(frames) > 1:
                # Note @dev : may be faster/less memory hungry to keep lines separated for each isotope. TODO : test both versions
                for df in frames:
                    if "iso" not in df.columns:
                        assert "iso" in df.attrs
                        df["iso"] = df.attrs["iso"]
                # Keep attributes:
                from radis.misc.basics import intersect

                attrs = frames[0].attrs
                for df in frames[1:]:
                    attrs = intersect(attrs, df.attrs)
                del attrs["iso"]  # added as a column (different for each line)
                # Merge:
                df = pd.concat(frames, ignore_index=True)  # reindex
                df.attrs = attrs
                self.params.dbpath = ",".join(local_paths)
            else:
                df = frames[0]
                self.params.dbpath = local_paths[0]

        elif source == "geisa":

            self.reftracker.add(doi["GEISA-2020"], "line database")

            if memory_mapping_engine == "auto":
                memory_mapping_engine = "vaex"

            if database != "full":
                raise ValueError(
                    f"Got `database={database}`. When fetching GEISA, only the `database='full'` option is available."
                )

            # Download, setup local databases, and fetch (use existing if possible)

            if isotope == "all":
                isotope_list = None
            else:
                isotope_list = ",".join([str(k) for k in self._get_isotope_list()])

            df, local_paths = fetch_geisa(
                molecule,
                isotope=isotope_list,
                local_databases=join(local_databases, "geisa"),
                load_wavenum_min=wavenum_min,
                load_wavenum_max=wavenum_max,
                columns=columns,
                cache=db_use_cached,
                verbose=self.verbose,
                return_local_path=True,
                engine=memory_mapping_engine,
                parallel=parallel,
            )
            self.params.dbpath = ",".join(local_paths)

            # ... explicitely write all isotopes based on isotopes found in the database
            if isotope == "all":
                self.input.isotope = ",".join(
                    [str(k) for k in self._get_isotope_list(df=df)]
                )

        else:
            raise NotImplementedError("source: {0}".format(source))

        if len(df) == 0:
            raise EmptyDatabaseError(
                f"{molecule} has no lines on range "
                + "{0:.2f}-{1:.2f} cm-1".format(wavenum_min, wavenum_max)
            )

        # Always sort line database by wavenumber (required to SPARSE_WAVERANGE mode)
        df.sort_values("wav", ignore_index=True, inplace=True)

        # Post-processing of the line database
        # (note : this is now done in 'fetch_hitemp' before saving to the disk)
        # spectroscopic quantum numbers will be needed for nonequilibrium calculations, and line survey.
        if parse_local_global_quanta and "locu" in df and source != "geisa":
            df = parse_local_quanta(df, molecule, verbose=self.verbose)
        if (
            parse_local_global_quanta and "globu" in df and source != "geisa"
        ):  # spectroscopic quantum numbers will be needed for nonequilibrium calculations :
            df = parse_global_quanta(df, molecule, verbose=self.verbose)

        # Remove non numerical attributes
        if drop_non_numeric:
            if "branch" in df:
                replace_PQR_with_m101(df)
            df = drop_object_format_columns(df, verbose=self.verbose)

        self.df0 = df  # type : pd.DataFrame
        self.misc.total_lines = len(df)  # will be stored in Spectrum metadata

        # %% Init Partition functions (with energies)
        # ------------

        if parfuncfmt == "exomol":
            self._init_equilibrium_partition_functions(
                parfunc,
                parfuncfmt,
                predefined_partition_functions=partition_function_exomol,
            )
        else:
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

        self._remove_unecessary_columns(df)

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
        db_use_cached=True,
        lvl_use_cached=True,
        load_energies=False,
        include_neighbouring_lines=True,
        drop_columns="auto",
        load_columns="equilibrium",
    ):
        """Loads databank from shortname in the :ref:`Configuration file.
        <label_lbl_config_file>` (`~/radis.json`), or by manually setting all
        attributes.

        Databank includes:
        - lines
        - partition function & format (tabulated or calculated)
        - (optional) energy levels, format

        Parameters
        ----------
        name: a section name specified in your ``~/radis.json``
            ``.radis`` has to be created in your HOME (Unix) / User (Windows). If
            not ``None``, all other arguments are discarded.
            Note that all files in database will be loaded and it may takes some
            time. Better limit the database size if you already know what
            range you need. See :ref:`Configuration file <label_lbl_config_file>` and
            :data:`~radis.misc.config.DBFORMAT` for expected
            ``~/radis.json`` format

        Other Parameters
        ----------------
        path: str, list of str, None
            list of database files, or name of a predefined database in the
            :ref:`Configuration file <label_lbl_config_file>` (`~/radis.json`)
            Accepts wildcards ``*`` to select multiple files
        format: ``'hitran'``, ``'cdsd-hitemp'``, ``'cdsd-4000'``, or any of :data:`~radis.lbl.loader.KNOWN_DBFORMAT`
            database type. ``'hitran'`` for HITRAN/HITEMP, ``'cdsd-hitemp'``
            and ``'cdsd-4000'`` for the different CDSD versions. Default ``'hitran'``
        parfuncfmt: ``'hapi'``, ``'cdsd'``, or any of :data:`~radis.lbl.loader.KNOWN_PARFUNCFORMAT`
            format to read tabulated partition function file. If ``hapi``, then
            HAPI (HITRAN Python interface) [1]_ is used to retrieve them (valid if
            your database is HITRAN data). HAPI is embedded into RADIS. Check the
            version. If partfuncfmt is None then ``hapi`` is used. Default ``hapi``.
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
            If ``None``, non equilibrium calculations are not possible. Default ``'radis'``.
        db_use_cached: boolean, or ``None``
            if ``True``, a pandas-readable csv file is generated on first access,
            and later used. This saves on the datatype cast and conversion and
            improves performances a lot. But! ... be sure to delete these files
            to regenerate them if you happen to change the database. If ``'regen'``,
            existing cached files are removed and regenerated.
            It is also used to load energy levels from ``.h5`` cache file if exist.
            If ``None``, the value given on Factory creation is used. Default ``True``
        load_energies: boolean
            if ``False``, dont load energy levels. This means that nonequilibrium
            spectra cannot be calculated, but it saves some memory. Default ``True``
        include_neighbouring_lines: bool
            ``True``, includes off-range, neighbouring lines that contribute
            because of lineshape broadening. The ``neighbour_lines``
            parameter is used to determine the limit. Default ``True``.
        *Other arguments are related to how to open the files:*
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
        load_columns: list, ``'all'``, ``'equilibrium'``, ``'noneq'``
            columns names to load.
            If ``'equilibrium'``, only load the columns required for equilibrium
            calculations. If ``'noneq'``, also load the columns required for
            non-LTE calculations. See :data:`~radis.lbl.loader.drop_all_but_these`.
            If ``'all'``, load everything. Note that for performances, it is
            better to load only certain columsn rather than loading them all
            and dropping them with ``drop_columns``.
            Default ``'equilibrium'``.

            .. warning::
                if using ``'equilibrium'``, not all parameters will be available
                for a Spectrum :py:func:`~radis.spectrum.spectrum.Spectrum.line_survey`.


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
        # %% Check inputs
        # ---------

        # use radis default for calculations non equilibrium calculations
        # if the levelsfmt is not specified in the databank
        if levelsfmt is None:
            levelsfmt = "radis"

        (
            name,
            path,
            dbformat,
            parfunc,
            parfuncfmt,
            levels,
            levelsfmt,
            db_use_cached,
            lvl_use_cached,
            drop_columns,
            load_columns,
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
            lvl_use_cached=lvl_use_cached,
            load_energies=load_energies,
            drop_columns=drop_columns,
            load_columns=load_columns,
            include_neighbouring_lines=include_neighbouring_lines,
        )
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
            load_energies=load_energies,
            lvl_use_cached=lvl_use_cached,
            include_neighbouring_lines=include_neighbouring_lines,
        )
        # Now that we're all set, let's load everything

        # %% Line database
        # ------------
        self._reset_references()  # bibliographic references

        self.df0 = self._load_databank(
            path,
            dbformat,
            levelsfmt=levelsfmt,
            db_use_cached=db_use_cached,
            drop_columns=drop_columns,
            load_columns=load_columns,
            include_neighbouring_lines=include_neighbouring_lines,
        )
        self.misc.total_lines = len(self.df0)  # will be stored in Spectrum metadata

        # Check the molecule is what we expected
        if self.input.molecule not in ["", None]:
            assert self.input.molecule == get_molecule(
                self.df0.attrs["id"]
            )  # assert molecule is what we expected
        else:
            self.input.molecule = get_molecule(self.df0.attrs["id"])  # get molecule

        # %% Partition functions (with energies)
        # ------------

        self._init_equilibrium_partition_functions(parfunc, parfuncfmt)

        # If energy levels are given, initialize the partition function calculator
        # (necessary for non-equilibrium). If levelsfmt == 'radis' then energies
        # are calculated ab initio from radis internal species database constants
        if load_energies:
            self._init_rovibrational_energies(levels, levelsfmt)

        return

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
        lvl_use_cached=None,
        load_energies=False,
        drop_columns="auto",
        load_columns="required",
        include_neighbouring_lines=True,
    ):
        """Check that database parameters are valid, in particular that paths
        exist. Loads all parameters if a Database from radis.json config file was
        given.

        Returns
        -------
        tuple
            (name, path, dbformat, parfunc, parfuncfmt, levels, levelsfmt,
             db_use_cached, load_energies, include_neighbouring_lines, drop_columns)
        """

        dbformat = format

        # Get database format and path
        # ... either from name (~/radis.json config file)
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
                    + " ~/radis.json: {0}".format(",".join(dblist))
                )

        # Check database format
        if dbformat not in KNOWN_DBFORMAT:
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

        # Line database path
        if isinstance(path, str):  # make it a list
            path = [path]

        # ... Parse all paths and read wildcards
        path_list = [expanduser(p) for p in path]
        new_paths = []
        for pathrg in path_list:
            path = get_files_from_regex(pathrg)

            # Ensure that `path` does not contain the cached dataset files in
            # case a wildcard input is given by the user. For instance, if the
            # given input is "cdsd_hitemp_09_frag*", path should not contain both
            # "cdsd_hitemp_09_fragment.txt" and "cdsd_hitemp_09_fragment.h5".

            # Reference: https://github.com/radis/radis/issues/121

            filtered_path = [fname for fname in path]
            for fname in path:
                for likely_fname_cache in [
                    splitext(fname)[0] + ".h5",
                    fname + ".h5",
                    splitext(fname)[0] + ".hdf5",
                    fname + ".hdf5",
                ]:
                    if likely_fname_cache in path and likely_fname_cache != fname:
                        filtered_path.remove(likely_fname_cache)
            new_paths += filtered_path

            # Raise errors if no file / print which files were selected.
            if len(filtered_path) == 0:
                self.warn(f"Path `{pathrg}` match no file", "DatabaseNotFoundError")
            if self.verbose >= 3 and pathrg != filtered_path:
                printg(
                    f"Regex `{pathrg}` match {len(filtered_path)} files: {filtered_path}"
                )

        path = new_paths

        # ... Check all path exists
        # ... remove empty paths first
        path = [p for p in path if p != ""]
        # ... test paths
        for p in path:
            if not exists(p):
                raise FileNotFoundError("databank lines file: `{0}`".format(p))

        # Energy levels and partition functions
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
            lvl_use_cached,
            drop_columns,
            load_columns,
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
        lvl_use_cached=None,
        load_energies=False,
        include_neighbouring_lines=True,
    ):
        """store all params so they can be parsed by "get_conditions()" and
        saved in output spectra information.

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
        self.levels = levels
        if levels is not None:
            self.levelspath = ",".join([format_paths(lvl) for lvl in levels.values()])
        else:
            self.levelspath = None
        self.params.levelsfmt = levelsfmt
        self.params.parfuncpath = format_paths(parfunc)
        self.params.parfuncfmt = parfuncfmt
        self.params.include_neighbouring_lines = include_neighbouring_lines
        self.params.db_use_cached = db_use_cached
        self.params.lvl_use_cached = lvl_use_cached
        self.misc.load_energies = load_energies

    def init_database(
        self,
        path,
        autoretrieve=True,
        autoupdate=True,
        add_info=["Tvib", "Trot"],
        add_date="%Y%m%d",
        compress=True,
        **kwargs,
    ):
        """Init a :class:`~radis.tools.database.SpecDatabase` folder in
        ``path`` to later store our spectra. Spectra can also be automatically
        retrieved from the database instead of being calculated.

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
        compress: boolean, or 2
            if ``True``, Spectrum are read and written in binary format. This is faster,
            and takes less memory space. Default ``True``.
            If ``2``, additionaly remove all redundant quantities.

        Other Parameters
        ----------------
        **kwargs: **dict
            arguments sent to :py:class:`~radis.tools.database.SpecDatabase` initialization.

        Returns
        -------
        db: SpecDatabase
            the database where spectra will be stored or retrieved


        .. minigallery:: radis.lbl.loader.DatabankLoader.init_database
        """

        db = SpecDatabase(
            path, add_info=add_info, add_date=add_date, binary=compress, **kwargs
        )

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
    # get_partition_function_interpolator
    # get_partition_function_calculator
    #
    # =========================================================================

    def _init_equilibrium_partition_functions(
        self, parfunc, parfuncfmt, predefined_partition_functions={}
    ):
        """Initializes equilibrium partition functions in ``self.parsum_tab``

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

        Other Parameters
        ----------------
        predefined_partition_functions: dict
            ::
                {molecule: {isotope: PartitionFunctionTabulator object}}
        """

        # Let's get the tabulated partition function (to calculate eq spectra)
        molecule = self.input.molecule
        state = self.input.state
        self.parsum_tab[molecule] = {}
        for iso in self._get_isotope_list():
            self.parsum_tab[molecule][iso] = {}
            ParsumTab = self._build_partition_function_interpolator(
                parfunc,
                parfuncfmt,
                self.input.molecule,
                isotope=iso,
                predefined_partition_functions=predefined_partition_functions,
            )
            self.parsum_tab[molecule][iso][state] = ParsumTab

    def _init_rovibrational_energies(self, levels, levelsfmt):
        """Initializes non equilibrium partition (which contain rovibrational
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
                    lvl,
                    levelsfmt,
                    isotope=iso,
                    parsum_mode=self.params.parsum_mode,
                )
                self.parsum_calc[molecule][iso][state] = ParsumCalc
        # energy levels arent specified in a tabulated file, but we can still
        # calculate them directly from Dunham expansions:
        elif levelsfmt == "radis":
            for iso in self._get_isotope_list():
                self.parsum_calc[molecule][iso] = {}
                ParsumCalc = self._build_partition_function_calculator(
                    None,
                    levelsfmt,
                    isotope=iso,
                    parsum_mode=self.params.parsum_mode,
                )
                self.parsum_calc[molecule][iso][state] = ParsumCalc

    def _check_line_databank(self):
        """Make sure database is loaded, loads if it isnt and we have all the
        information needed.
        Databank has been initialized by
        :meth:`~radis.lbl.loader.DatabankLoader.init_databank`
        """
        self.profiler.start("check_line_databank", 2)
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

        self.profiler.stop("check_line_databank", "Check line databank")

    def _load_databank(
        self,
        database,
        dbformat,
        levelsfmt,
        db_use_cached,
        drop_columns,
        load_columns,
        include_neighbouring_lines=True,
    ) -> pd.DataFrame:

        """Loads all available database files and keep the relevant one.
        Returns a Pandas dataframe.

        Parameters
        ----------
        database: list of str
            list of database files
        db_use_cached: boolean, or ``'regen'``
            if ``True``, a pandas-readable csv file is generated on first access,
            and later used. This saves on the datatype cast and conversion and
            improves performances a lot. But! ... be sure to delete these files
            to regenerate them if you happen to change the database. Default ``False``
            If ``'regen'`` regenerate existing cache files.
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
        load_columns: list, ``'all'``, ``'equilibrium'``, ``'noneq'``
            columns names to load.
            If ``'equilibrium'``, only load the columns required for equilibrium
            calculations. If ``'noneq'``, also load the columns required for
            non-LTE calculations. See :data:`~radis.lbl.loader.drop_all_but_these`.
            If ``'all'``, load everything. Note that for performances, it is
            better to load only certain columsn rather than loading them all
            and dropping them with ``drop_columns``.
            Default ``'equilibrium'``.

            .. warning::
                if using ``'equilibrium'``, not all parameters will be available
                for a Spectrum :py:func:`~radis.spectrum.spectrum.Spectrum.line_survey`.


        Other Parameters
        ----------------
        include_neighbouring_lines: bool
            ``True``, includes off-range, neighbouring lines that contribute
            because of lineshape broadening. The ``neighbour_lines``
            parameter is used to determine the limit. Default ``True``.
        """
        # Check inputs
        assert db_use_cached in [True, False, "regen", "force"]

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

        if drop_columns == "auto":
            drop_columns = (
                drop_auto_columns_for_dbformat[dbformat]
                + drop_auto_columns_for_levelsfmt[levelsfmt]
            )

        # which columns to load
        columns = []
        if "all" in load_columns:
            columns = None  # see fetch_hitemp, fetch_hitran, etc.
        elif isinstance(load_columns, str) and load_columns in ["equilibrium", "noneq"]:
            columns = self.columns_list_to_load(load_columns)
        elif load_columns == "diluent":
            raise ValueError(
                "Please use diluent along with 'equilibrium' or 'noneq' in a list like ['diluent','noneq']"
            )

        elif isinstance(load_columns, list) and "all" not in load_columns:
            for load_columns_type in load_columns:
                if load_columns_type in ["equilibrium", "noneq", "diluent"]:
                    for col in self.columns_list_to_load(load_columns_type):
                        columns.append(col)
                elif load_columns_type in list(
                    set(drop_all_but_these)
                    | set(required_non_eq)
                    | set(broadening_coeff)
                ):
                    columns.append(load_columns_type)
                else:
                    raise ValueError("invalid column name provided")
            columns = list(set(columns))

        # subroutine load_and_concat
        # --------------------------------------
        def load_and_concat(files):
            """Contatenate many files in RAM
            Parameters
            ----------
            files: list of str
                elist of path to database files ::
                    [PATH/TO/01_1000-1150_HITEMP2010.par,
                     PATH/TO/01_1150-1300_HITEMP2010.par,
                     PATH/TO/01_1300-1500_HITEMP2010.par]
            """

            frames = []

            for i, filename in enumerate(files):

                if __debug__:
                    printdbg("Loading {0}/{1}".format(i + 1, len(files)))

                # Read all the lines
                # ... this is where the cache files are read/generated.
                try:
                    if dbformat in ["cdsd-hitemp", "cdsd-4000"]:
                        if dbformat == "cdsd-4000":
                            self.reftracker.add(
                                doi["CDSD-4000"], "line database"
                            )  # [CDSD-4000]_
                        if dbformat == "cdsd-hitemp":
                            self.warn(
                                "Missing doi for CDSD-HITEMP. Use HITEMP-2010?",
                                "MissingReferenceWarning",
                            )
                        df = cdsd2df(
                            filename,
                            version="hitemp" if dbformat == "cdsd-hitemp" else "4000",
                            cache=db_use_cached,
                            # load_columns=columns,  # not possible with "pytables-fixed"
                            verbose=verbose,
                            drop_non_numeric=True,
                            load_wavenum_min=wavenum_min,
                            load_wavenum_max=wavenum_max,
                            engine="pytables",
                        )
                        # TODO: implement load_columns
                    elif dbformat in ["hitran", "hitemp"]:
                        if dbformat == "hitran":
                            self.reftracker.add(
                                doi["HITRAN-2020"], "line database"
                            )  # [HITRAN-2020]_
                        if dbformat == "hitemp":
                            self.reftracker.add(
                                doi["HITEMP-2010"], "line database"
                            )  # [HITEMP-2010]_
                        df = hit2df(
                            filename,
                            cache=db_use_cached,
                            # load_columns=columns,  # not possible with "pytables-fixed"
                            verbose=verbose,
                            drop_non_numeric=True,
                            load_wavenum_min=wavenum_min,
                            load_wavenum_max=wavenum_max,
                            engine="pytables",
                        )
                    elif dbformat in ["hdf5-radisdb", "hitemp-radisdb"]:
                        if dbformat == "hitemp-radisdb":
                            self.reftracker.add(
                                doi["HITEMP-2010"], "line database"
                            )  # [HITEMP-2010]_
                        if dbformat == "hdf5-radisdb":
                            self.warn(
                                f"Missing doi reference for database used {filename}",
                                "MissingReferenceWarning",
                            )
                        df = hdf2df(
                            filename,
                            columns=columns,
                            # cache=db_use_cached,
                            verbose=verbose,
                            # drop_non_numeric=True,
                            isotope=self.input.isotope
                            if self.input.isotope != "all"
                            else None,
                            load_wavenum_min=wavenum_min,
                            load_wavenum_max=wavenum_max,
                            engine="guess",
                        )
                    elif dbformat in ["exomol"]:
                        # self.reftracker.add("10.1016/j.jqsrt.2020.107228", "line database")  # [ExoMol-2020]
                        raise NotImplementedError("use fetch_databank('exomol')")

                    else:
                        raise ValueError("Unknown dbformat: {0}".format(dbformat))
                except IrrelevantFileWarning as err:
                    if db_use_cached == "force":
                        raise
                    else:
                        # Irrelevant file, just print and continue.
                        if verbose >= 2:
                            printg(str(err))
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
                # TODO : is it still needed since we use load_only_wavenum_above ?
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

                frames.append(df)

            # Finally: Concatenate all
            if frames == []:
                df = (
                    pd.DataFrame()
                )  # a database empty error will be raised a few lines below
            else:
                df = pd.concat(frames, ignore_index=True)  # reindex

            return df

        # end subroutine load_and_concat
        # --------------------------------------

        df = load_and_concat(database)

        # Final checks

        # ... error in Pandas? Sometimes _metadata is preserved over several runs.
        # ... Clean it here.
        if __debug__ and len(df.attrs) > 0:
            printdbg("df.attrs was not []. Cleaning")
        df.attrs = []

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
            raise EmptyDatabaseError(msg)

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
            neighbour_lines = self.params.neighbour_lines
            if neighbour_lines > 0 and minwavdb > wavenum_min + neighbour_lines:
                # no lines on left side
                self.warn(
                    "There are no lines in database in range {0:.5f}-{1:.5f}cm-1 ".format(
                        wavenum_min, wavenum_min + neighbour_lines
                    )
                    + "to calculate the effect "
                    + "of neighboring lines. Did you add all lines in the database?",
                    "OutOfRangeLinesWarning",
                )
            if neighbour_lines > 0 and maxwavdb < wavenum_max - neighbour_lines:
                # no lines on right side
                self.warn(
                    "There are no lines in database in range {0:.5f}-{1:.5f}cm-1 ".format(
                        maxwavdb - neighbour_lines, maxwavdb
                    )
                    + "to calculate the effect "
                    + "of neighboring lines. Did you add all lines in the database?",
                    "OutOfRangeLinesWarning",
                )

        if self.verbose >= 2:
            printg(
                "Loaded databank in {0:.1f}s ({1:,d} lines)".format(
                    time() - t0, len(df)
                )
            )

        self._remove_unecessary_columns(df)

        return df

    def _remove_unecessary_columns(self, df):
        """Remove unecessary columns and add values as attributes

        Returns
        -------
        None: DataFrame updated inplace
        """

        # Discard molecule column if unique
        if "id" in df.columns:
            id_set = df.id.unique()
            if len(id_set) != 1:  # only 1 molecule supported ftm
                raise NotImplementedError(
                    "Only 1 molecule at a time is currently supported "
                    + "in SpectrumFactory. Use radis.calc_spectrum, which "
                    + "calculates them independently then use MergeSlabs"
                )

            df.drop("id", axis=1, inplace=True)
            df_metadata.append("id")
            df.attrs["id"] = id_set[0]
        else:
            assert "id" in df.attrs or "molecule" in df.attrs

        if "iso" in df.columns:
            isotope_set = df.iso.unique()

            if len(isotope_set) == 1:
                df.drop("iso", axis=1, inplace=True)
                df_metadata.append("iso")
                df.attrs["iso"] = isotope_set[0]
        else:
            assert "iso" in df.attrs

    def _get_isotope_list(self, molecule=None, df=None):
        """Returns list of isotopes for given molecule Parse the Input
        conditions (fast). If a line database is given, parse the line database
        instead (slow)

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
        """Retrieve a spectrum from a SpecDatabase database, if it matches
        current factory conditions.

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
            for (k, v) in self.get_conditions(
                ignore_misc=ignore_misc, add_config=True
            ).items()
            if k in self.SpecDatabase.conditions()
        }
        conditions = {
            k: v
            for (k, v) in conditions.items()
            if k not in self._autoretrieveignoreconditions
        }

        if "wstep" in conditions and conditions["wstep"] == "auto":
            if "GRIDPOINTS_PER_LINEWIDTH_WARN_THRESHOLD" in conditions:
                # ignore wstep='auto' but still make sure that 'GRIDPOINTS_PER_LINEWIDTH_WARN_THRESHOLD' is the same
                del conditions["wstep"]
                # TODO Refactor. Include this in a large Parameter metaclass where
                # 'GRIDPOINTS_PER_LINEWIDTH_WARN_THRESHOLD' is a parameter from which
                # depends 'wstep' evaluation

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
                    """Returns the Spectrum that matches the input conditions
                    better.
                    Used to give a better error message in case no
                    Spectrum was found
                    """
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
                    print("Differences between us (left) and best case (right):")
                    compare_dict(conditions, best.conditions)

                    raise ValueError(
                        "No spectrum found in database that matched "
                        + "given conditions. See best case found above. You can also add parameters in SpectrumFactory._ignoreautoretrieveconditions to ignore them"
                    )
            else:
                # just a print, then calculate
                if self.verbose:
                    print(
                        "No spectrum found in database that "
                        + "matched given conditions. You can also add parameters in SpectrumFactory._ignoreautoretrieveconditions to ignore them"
                    )
                return None

    def _build_partition_function_interpolator(
        self, parfunc, parfuncfmt, molecule, isotope, predefined_partition_functions={}
    ):
        """Returns an universal partition function object ``parsum`` with the
        following methods defined::
            parsum.at(T)

        Partition functions are interpolated from tabulated values

        Other Parameters
        ----------------
        predefined_partition_functions: dict
            ::
                {molecule: {isotope: PartitionFunctionTabulator object}}
        """

        if __debug__:
            printdbg(
                "called _build_partition_function_interpolator"
                + "(parfuncfmt={0}, isotope={1})".format(parfuncfmt, isotope)
            )

        isotope = int(isotope)

        if parfuncfmt in ["hapi", "tips"] or parfuncfmt is None:
            assert len(predefined_partition_functions) == 0
            self.reftracker.add(doi["TIPS-2020"], "partition function")
            self.reftracker.add(doi["HAPI"], "partition function")
            # Use TIPS-2017 through HAPI (HITRAN Python interface, integrated in RADIS)
            # no tabulated partition functions defined. Only non-eq spectra can
            # be calculated if energies are also given
            parsum = PartFuncTIPS(
                M=molecule, I=isotope, path=parfunc, verbose=self.verbose
            )
        elif parfuncfmt == "cdsd":  # Use tabulated CDSD partition functions
            self.reftracker.add(doi["CDSD-4000"], "partition function")
            assert len(predefined_partition_functions) == 0
            assert molecule == "CO2"
            parsum = PartFuncCO2_CDSDtab(isotope, parfunc)
        elif parfuncfmt == "exomol":
            self.reftracker.add(doi["ExoMol-2020"], "partition function")
            # Just read dictionary of predefined partition function
            assert len(predefined_partition_functions) > 0
            parsum = predefined_partition_functions[molecule][isotope]
        else:
            raise ValueError(
                "Unknown format for partition function: {0}".format(parfuncfmt)
            )
            # other formats ?

        return parsum

    def _build_partition_function_calculator(
        self, levels, levelsfmt, isotope, parsum_mode="full summation"
    ):
        """Return an universal partition function  object ``parsum`` so that
        the following methods are defined::
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

        Other Parameters
        ----------------
        parsum_mode: 'full summation', 'tabulation'
            calculation mode. ``'tabulation'`` is much faster but not all possible
            distributions are implemented. See ``mode`` in
            :py:class:`~radis.levels.partfunc.RovibParFuncCalculator`
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
            self.reftracker.add(doi["CDSD-4000"], "rovibrational energies")
            parsum = PartFuncCO2_CDSDcalc(
                levels,
                isotope=isotope,
                use_cached=self.params.lvl_use_cached,
                verbose=self.verbose,
                levelsfmt=levelsfmt,
                mode=parsum_mode,
            )

        # calculate energy levels from RADIS Dunham parameters
        elif levelsfmt == "radis":
            self.reftracker.add(doi["RADIS-2018"], "rovibrational energies")
            state = getMolecule(
                self.input.molecule, isotope, self.input.state, verbose=self.verbose
            )
            if state.doi is not None:
                self.reftracker.add(state.doi, "spectroscopic constants")
            parsum = PartFunc_Dunham(
                state,
                use_cached=self.params.lvl_use_cached,
                verbose=self.verbose,
                mode=parsum_mode,
            )
            # note: use 'levels' (useless here) to specify calculations options
            # for the abinitio calculation ? Like Jmax, etc.

        else:
            raise ValueError("Unknown format for energy levels : {0}".format(levelsfmt))
            # other formats ?

        return parsum

    def get_abundance(self, molecule, isotope):
        """Get isotopic abundance

        Parameters
        ----------
        molecule: str
        isotope: int, or list
            isotope number, sorted in terrestrial abundance

        Examples
        --------
        Use it from SpectrumFactory::

            sf.get_abundance("H2O", 1)
            sf.get_abundance("CH4", [1,2,3])

        .. minigallery:: radis.lbl.loader.DatabankLoader.get_abundance
            :add-heading:

        See Also
        --------
        :py:meth:`~radis.lbl.loader.DatabankLoader.set_abundance`
        """

        if isinstance(molecule, str):
            from radis.db.classes import get_molecule_identifier

            molecule = get_molecule_identifier(molecule)

        if isinstance(isotope, int):
            return self.molparam.df.loc[(molecule, isotope)].abundance
        elif isinstance(isotope, list):
            return np.array(
                [self.molparam.df.loc[(molecule, iso)].abundance for iso in isotope]
            )
        else:
            raise ValueError(isotope)

    def set_abundance(self, molecule, isotope, abundance):
        """Set isotopic abundance

        Parameters
        ----------
        molecule: str
        isotope: int, or list
            isotope number, sorted in terrestrial abundance
        abundance: float, or list

        Examples
        --------
        ::

            from radis import SpectrumFactory

            sf = SpectrumFactory(
                2284.2,
                2284.6,
                wstep=0.001,  # cm-1
                pressure=20 * 1e-3,  # bar
                mole_fraction=400e-6,
                molecule="CO2",
                isotope="1,2",
                verbose=False
            )
            sf.load_databank("HITEMP-CO2-TEST")
            print("Abundance of CO2[1,2]", sf.get_abundance("CO2", [1, 2]))
            sf.eq_spectrum(2000).plot("abscoeff")

            #%% Set the abundance of CO2(626) to 0.8; and the abundance of CO2(636) to 0.2 (arbitrary):
            sf.set_abundance("CO2", [1, 2], [0.8, 0.2])
            print("New abundance of CO2[1,2]", sf.get_abundance("CO2", [1, 2]))

            sf.eq_spectrum(2000).plot("abscoeff", nfig="same")

        .. minigallery:: radis.lbl.loader.DatabankLoader.set_abundance
            :add-heading

        See Also
        --------
        :py:meth:`~radis.lbl.loader.DatabankLoader.get_abundance`

        """

        if isinstance(molecule, str):
            from radis.db.classes import get_molecule_identifier

            molecule = get_molecule_identifier(molecule)

        self.molparam.terrestrial_abundances = False

        if isinstance(isotope, int):
            self.molparam.df.loc[(molecule, isotope), "abundance"] = abundance
        elif isinstance(isotope, list):
            assert len(isotope) == len(abundance)
            self.molparam.df.loc[(molecule, isotope), "abundance"] = abundance
        else:
            raise ValueError(isotope)

    def get_partition_function_interpolator(self, molecule, isotope, elec_state):
        """Retrieve Partition Function Interpolator.

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
        """Retrieve Partition Function Calculator.

        Parameters
        ----------
        molecule: str
        isotope: int
        elec_state: str
        """

        parsum = self.get_partition_function_molecule(molecule)[isotope][elec_state]

        # helps IDE find methods
        assert isinstance(parsum, RovibParFuncCalculator)

        # Update partition sum calculation mode (if has been reset by user)
        if parsum.mode != self.params.parsum_mode:
            parsum.mode = self.params.parsum_mode

        return parsum

    def get_partition_function_molecule(self, molecule):
        """Retrieve Partition Function for Molecule.

        Parameters
        ----------
        molecule: str
        """
        try:
            parsum = self.parsum_calc[molecule]
        except KeyError as err:
            raise KeyError(
                "Error while Retrieving Partition Function of Molecule!"
                + " Load the energies levels with SpectrumFactory.load_databank"
                + "('path', load_energies=True). If using SpectrumFactory.fetch_databank()"
                + " consider adding arguement load_energies=True"
            ) from err

        return parsum

    def get_conditions(self, ignore_misc=False, add_config=False):
        """Get all parameters defined in the SpectrumFactory.

        Other Parameters
        ----------------
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

        if add_config:
            import radis
            from radis.spectrum.utils import CONFIG_PARAMS

            vardict.update(
                {k: v for k, v in radis.config.items() if k in CONFIG_PARAMS}
            )

        return vardict

    def warn(self, message, category="default", level=0):
        """Trigger a warning, an error or just ignore based on the value
        defined in the :attr:`~radis.lbl.loader.DatabankLoader.warnings`
        dictionary.
        The warnings can thus be deactivated selectively by setting the SpectrumFactory
         :attr:`~radis.lbl.loader.DatabankLoader.warnings` attribute

        Parameters
        ----------
        message: str
            what to print
        category: str
            one of the keys of self.warnings. See :py:attr:`~radis.lbl.loader.DatabankLoader.warnings`
        level: int
            warning level. Only print warnings when verbose level is higher
            than the warning levels. i.e., warnings of level 1 appear only
            if ``verbose==True``, warnings of level 2 appear only
            for ``verbose>=2``, etc..  Warnings of level 0 appear only the time.
            Default ``0``

        Examples
        --------
        ::
            if not ((df.Erotu > tol).all() and (df.Erotl > tol).all()):
                self.warn(
                    "There are negative rotational energies in the database",
                    "NegativeEnergiesWarning",
                )

        Notes
        -----
        All warnings in the SpectrumFactory should call to this method rather
        than the default warnings.warn() method, because it allows easier runtime
        modification of how to deal with warnings

        See Also
        --------
        :py:attr:`~radis.lbl.loader.DatabankLoader.warnings`
        """
        if level > self.verbose:
            return
        else:
            return warn(message, category=category, status=self.warnings)


if __name__ == "__main__":

    from radis.test.lbl.test_loader import _run_testcases

    print("Testing loading functions:", _run_testcases())
