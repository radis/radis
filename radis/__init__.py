# -*- coding: utf-8 -*-
"""

::

                *(((((((
                 ((((((((((((              ,(((((
                 ((((((((((((((((/   *((((((((((*
                  ((((((((((((((((( ((((((((((((
                      (((((((( (((((((((((((
                         *
                       @@  *@@       ..  /@@
                  @@&  @@  *@@       @@  /@@  @@%
              @@  @@&  @@  *@@  @@&  @@  /@@  @@%  @@
              @@  @@&  @@  *@@  @@&  @@  /@@  @@%  @@
              @@  @@&  @@  *@@  @@&  @@  /@@  @@%  @@  (@
         ,@   @@  @@&  @@  *@@  @@&  @@  /@@  @@%  @@
         @@   @@  @@&  @@  ,.
                                    ,%&&&&&&&&&&&&&&&&&&&
          &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
           &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
             &&&&&&&&&&&&&&&&@@@@@@&@@@&&&@@@&&&&&&&&
               &&&&&&&&&&&&&&&@@@@@@&&&&&&&&&&&&&&&
                 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                   &&&&&&&&&&&&&&&&&&&&&&&&&&&.
                       &&&&&&&&&&&&&&&&&&&
                               .**.
                                &&&,
                                 &&

"""

import importlib
import os

from .misc.utils import Chdir as _chdir
from .misc.utils import getProjectRoot

# %% Config files

# @dev: refactor in progress.
# So far there are config files in ~/radis.json (for databanks), global variables
# here, and a radis/config.json file.
# Everything should be merged in a user JSON file ~/radis.json (json) overriding
# the default one.
# # molecular parameters for non-HITRAN species    should also be removed #TODO


class _LazyConfig(dict):
    """Lazy-loading wrapper for RADIS global config."""

    def __init__(self):
        super().__init__()
        self._is_loaded = False

    def _ensure_loaded(self):
        if self._is_loaded:
            return
        from .misc.config import get_config

        super().update(get_config())
        self._is_loaded = True

    def __getitem__(self, key):
        self._ensure_loaded()
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        self._ensure_loaded()
        return super().__setitem__(key, value)

    def __contains__(self, key):
        self._ensure_loaded()
        return super().__contains__(key)

    def __iter__(self):
        self._ensure_loaded()
        return super().__iter__()

    def __len__(self):
        self._ensure_loaded()
        return super().__len__()

    def get(self, key, default=None):
        self._ensure_loaded()
        return super().get(key, default)

    def items(self):
        self._ensure_loaded()
        return super().items()

    def keys(self):
        self._ensure_loaded()
        return super().keys()

    def values(self):
        self._ensure_loaded()
        return super().values()

    def update(self, *args, **kwargs):
        self._ensure_loaded()
        return super().update(*args, **kwargs)

    def copy(self):
        self._ensure_loaded()
        return super().copy()

    def __repr__(self):
        self._ensure_loaded()
        return super().__repr__()


config = _LazyConfig()
"""dict: RADIS configuration parameters

Parameters
----------
"DEBUG_MODE":False

    bool: change this at runtime with::

        import radis
        radis.config["DEBUG_MODE"] = True

    Use the :py:func:`~radis.misc.debug.printdbg` function in ``radis.misc``, typically with::

        if __debug__: printdbg(...)

    so that printdbg are removed by the Python preprocessor when running in
    optimize mode::

        python -O *.py

"ALLOW_OVERWRITE": False

    bool: Whether to allow RADIS to overwrite values of a databank entry that is already registered.
    If True, allow RADIS to download files and make changes to the local database folder even when the databank is already registered, and then re-register it with the values for the new files, overwriting the previously registered values. This could occur when e.g. running `fetch_databank` with the `cache` argument set to "regen" or the registered "paths" not including the path to the databank itself or some of the registered files missing from the disk, resulting in files being re-downloaded and the "download_date" (and possibly "paths") value being updated.
    If False, no changes are made to the local database folder when a databank is already registered, and an exception is raised in cases where re-registering would have been required, so e.g. in the example above, this setting would prevent the "download_date" from being updated and hence an exception would be raised instead.
    False by default, so RADIS only ever inserts new databank entries and doesn't modify existing ones, hence a databank entry is only ever modified by the user and responsibility is on them to ensure it is adequate for the code to run as RADIS wouldn't automatically take steps to rectify inadequate entries.

"AUTO_UPDATE_SPEC": False

    bool: experimental feature
    used to autoupdate .spec files to the latest format, by simply saving
    them again once they're loaded and fixed.
    Warning! Better have a copy of your files before that, or a way to regenerate
    them.

    Example : Add to the top of your script (once is enough!)::

        import radis
        radis.config["AUTO_UPDATE_SPEC"] = True

    You can also see: :py:func:`~radis.tools.database._update_to_latest_format`


"AUTO_UPDATE_DATABASE": False
	bool: experimental feature
	used to autoupdate .h5 files to the latest format, by simply saving
	them again once they're loaded and fixed.
	Warning! Better have a copy of your files before that, or a way to regenerate
	them.

	Example. Add to the top of your script (once is enough!)::

	    import radis
	    radis.AUTO_UPDATE_DATABASE = True

	See Also
	--------
	:py:func:`~radis.api.hdf5.hdf2df`

"MISSING_BROAD_COEF" : False

    str: Accepted values: False (default) and "air". If "air", missing boradening coefficients are replaced by those of air.

    See Also
    --------
    :py:func:`~radis.lbl.broadening._calc_broadening_HWHM`

"OLDEST_COMPATIBLE_VERSION": "0.9.1"
    str: forces to regenerate cache files that were created in a previous version

    See Also
    --------
    :py:func:`~radis.api.cache_files.load_h5_cache_file`

"GRIDPOINTS_PER_LINEWIDTH_WARN_THRESHOLD": 3
    float: to determine the optimal value
    of wstep using minimum FWHM value of spectrum.
    Makes sure there are enough gridpoints per line.

    See more in :  :py:meth:`radis.lbl.broadening.BroadenFactory._check_accuracy`


"GRIDPOINTS_PER_LINEWIDTH_ERROR_THRESHOLD": 1
    float: to determine the minimum feasible value
    of wstep using minimum FWHM value of spectrum.
    Makes sure there are enough gridpoints per line.

    See more in :     :py:meth:`radis.lbl.broadening.BroadenFactory._check_accuracy`

"SPARSE_WAVERANGE": True
    bool: if True, allow special optimizations to improve performances for
    large spectra that are sparse, i.e there are no lines within ``truncation``
    distance. In that case, spectral calculations can be an order of magnitude
    faster, and result in less memory use.
    This optimization may result in a small but unnecessary overhead
    when there are lines everywhere in the spectral range considered. In such
    cases, you may deactivate it by setting ``radis.config['SPARSE_WAVERANGE'] = False``
    Default ``True``

    See more in  :py:meth:`radis.lbl.broadening.BroadenFactory._apply_lineshape_LDM`

"RESAMPLING_TOLERANCE_THRESHOLD" 5e-3
    an error if raises if areas do not match by a value above this threshold,
    during resampling. See :py:meth:`~radis.spectrum.spectrum.Spectrum.resample`

    Default ``5e-3``


Notes
-----

Default values are read from the `radis/default_radis.json <https://github.com/radis/radis/blob/develop/radis/default_radis.json>`__ file.

All values are overridden at runtime by the keys in the user JSON file ``~/radis.json (json)``
(in particular, the list of databases)

See more in the :ref:`Configuration file <label_lbl_config_file>` documentation.
"""
# TODO : Refactor in progress.


# %% Version


def get_version(verbose=False, add_git_number=True):
    """Reads `__version.txt__
    <https://github.com/radis/radis/blob/master/radis/__version__.txt>`__ and
    retrieve version number. If ``add_git_number``, also appends Git commit
    number if we're on a gitted session.

    Examples
    --------

    ::

        import radis
        print(radis.get_version())
        >>> '0.9.17'
    """

    # First get version
    with open(os.path.join(getProjectRoot(), "__version__.txt")) as version_file:
        version = version_file.read().strip()

    # Now get git info
    if add_git_number:
        import subprocess
        import sys

        cd = _chdir(os.path.dirname(__file__))
        try:
            label = subprocess.check_output("git describe", stderr=subprocess.DEVNULL)
        except:
            if verbose:
                print(f"couldnt get git version: {sys.exc_info()[1]}")
            # probably not a git session. drop
        else:
            commit = label.decode().strip().split("-")[-1]
            version = version + "-" + commit
        finally:
            cd.__del__()

    return version


__version__ = get_version(add_git_number=False)
version = get_version(add_git_number=False)


# %% Global namespace

__all__ = [
    "config",
    "version",
    "__version__",
    "cdsd2df",
    "hit2df",
    "MOLECULES_LIST_EQUILIBRIUM",
    "MOLECULES_LIST_NONEQUILIBRIUM",
    "Molecules",
    "getMolecule",
    "fetch_hitemp",
    "fetch_hitran",
    "fetch_exomol",
    "fetch_astroquery",
    "fetch_geisa",
    "LevelsList",
    "SpectrumFactory",
    "calc_spectrum",
    "spectrum_test",
    "MergeSlabs",
    "SerialSlabs",
    "sPlanck",
    "planck",
    "planck_wn",
    "get_diff",
    "get_distance",
    "get_ratio",
    "get_residual",
    "get_residual_integral",
    "plot_diff",
    "calculated_spectrum",
    "experimental_spectrum",
    "transmittance_spectrum",
    "PerfectAbsorber",
    "Radiance",
    "Radiance_noslit",
    "Transmittance",
    "Transmittance_noslit",
    "Spectrum",
    "SpecDatabase",
    "load_spec",
    "plot_spec",
    "save",
    "get_eq_mole_fraction",
    "plot_slit",
    "get_effective_FWHM",
    "get_FWHM",
]

_EXPORT_TO_MODULE = {
    "cdsd2df": "api",
    "hit2df": "api",
    "MOLECULES_LIST_EQUILIBRIUM": "db",
    "MOLECULES_LIST_NONEQUILIBRIUM": "db",
    "Molecules": "db",
    "getMolecule": "db",
    "fetch_hitemp": "io",
    "fetch_hitran": "io",
    "fetch_exomol": "io",
    "fetch_astroquery": "io",
    "fetch_geisa": "io",
    "LevelsList": "lbl",
    "SpectrumFactory": "lbl",
    "calc_spectrum": "lbl",
    "spectrum_test": "lbl",
    "MergeSlabs": "los",
    "SerialSlabs": "los",
    "sPlanck": "phys",
    "planck": "phys",
    "planck_wn": "phys",
    "get_diff": "spectrum",
    "get_distance": "spectrum",
    "get_ratio": "spectrum",
    "get_residual": "spectrum",
    "get_residual_integral": "spectrum",
    "plot_diff": "spectrum",
    "calculated_spectrum": "spectrum",
    "experimental_spectrum": "spectrum",
    "transmittance_spectrum": "spectrum",
    "PerfectAbsorber": "spectrum",
    "Radiance": "spectrum",
    "Radiance_noslit": "spectrum",
    "Transmittance": "spectrum",
    "Transmittance_noslit": "spectrum",
    "Spectrum": "spectrum",
    "SpecDatabase": "tools",
    "load_spec": "tools",
    "plot_spec": "tools",
    "save": "tools",
    "get_eq_mole_fraction": "tools",
    "plot_slit": "tools",
    "get_effective_FWHM": "tools",
    "get_FWHM": "tools",
}

_SUBMODULES = {
    "api",
    "db",
    "gpu",
    "io",
    "lbl",
    "levels",
    "los",
    "misc",
    "phys",
    "spectrum",
    "tools",
}


def _load_submodule(name):
    module = importlib.import_module(f".{name}", __name__)
    globals()[name] = module
    return module


def __getattr__(name):
    if name in _SUBMODULES:
        return _load_submodule(name)

    module_name = _EXPORT_TO_MODULE.get(name)
    if module_name is not None:
        module = _load_submodule(module_name)
        value = getattr(module, name)
        globals()[name] = value
        return value

    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__():
    return sorted(set(globals()) | set(__all__) | _SUBMODULES)
