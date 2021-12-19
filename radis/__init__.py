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

import os

from .misc.config import get_config
from .misc.utils import Chdir as _chdir
from .misc.utils import getProjectRoot

# %% Config files

# @dev: refactor in progress.
# So far there are config files in ~/radis.json (for databanks), global variables
# here, and a radis/config.json file.
# Everything should be merged in a user JSON file ~/radis.json (json) overriding
# the default one.

config = get_config()
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

"AUTO_UPDATE_SPEC": False

    bool: experimental feature
    used to autoupdate .spec files to the latest format, by simply saving
    them again once they're loaded and fixed.
    Warning! Better have a copy of your files before that, or a way to regenerate
    them.

    Example : Add to the top of your script (once is enough!)::

        import radis
        radis.config["AUTO_UPDATE_SPEC"] = True

    See Also
    --------
    :py:func:`~radis.tools.database._update_to_latest_format`


"AUTO_UPDATE_DATABASE": False
	bool: experimental feature
	used to autoupdate .h5 files to the latest format, by simply saving
	them again once they're loaded and fixed.
	Warning! Better have a copy of your files before that, or a way to regenerate
	them.

	Examples
	--------

	Add to the top of your script (once is enough!)::

	    import radis
	    radis.AUTO_UPDATE_DATABASE = True

	See Also
	--------
	:py:func:`~radis.io.hdf5.hdf2df`


"OLDEST_COMPATIBLE_VERSION": "0.9.1"
    str: forces to regenerate cache files that were created in a previous version

    See Also
    --------
    :py:func:`~radis.io.cache_files.load_h5_cache_file`


"USE_CYTHON": True
    bool: try to use Cython functions when possible

    See Also
    --------
    :py:func:`~radis.misc.arrays.add_at`


"GRIDPOINTS_PER_LINEWIDTH_WARN_THRESHOLD": 3
    float: to determine the optimal value
    of wstep using minimum FWHM value of spectrum.
    Makes sure there are enough gridpoints per line.

    See Also
    --------
    :py:meth:`~radis.lbl.broadening.BroadenFactory._check_accuracy`


"GRIDPOINTS_PER_LINEWIDTH_ERROR_THRESHOLD": 1
    float: to determine the minimum feasible value
    of wstep using minimum FWHM value of spectrum.
    Makes sure there are enough gridpoints per line.

    See Also
    --------
    :py:meth:`~radis.lbl.broadening.BroadenFactory._check_accuracy`

"SPARSE_WAVERANGE": True
    bool: if True, allow special optimizations to improve performances for
    large spectra that are sparse, i.e there are no lines within ``truncation``
    distance. In that case, spectral calculations can be an order of magnitude
    faster, and result in less memory use.
    This optimization may result in a small but unnecessary overhead
    when there are lines everywhere in the spectral range considered. In such
    cases, you may deactivate it by setting ``radis.config['SPARSE_WAVERANGE'] = False``
    Default ``True``

    See Also
    --------
    :py:meth:`~radis.lbl.broadening.BroadenFactory._apply_lineshape_LDM`



Notes
-----

Default values are read from the ``radis/config.json`` file.

All values are overriden at runtime by the keys in the user JSON file ``~/radis.json (json)``
(in particular, the list of databases)
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
                print("couldnt get git version: {0}".format(sys.exc_info()[1]))
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
]

# prevent cyclic importants:
from . import db, io, lbl, los, misc, phys, spectrum, tools
from .db import *  # database of molecules
from .io import *  # input / output
from .lbl import *  # line-by-line module
from .levels import *  # rovibrational energies and partition functions
from .los import *  # line-of-sight module
from .phys import *  # conversion functions, blackbody objects
from .spectrum import *  # Spectrum object
from .test import *  # test
from .tools import *  # slit, database, line survey, etc.

__all__.extend(db.__all__)
__all__.extend(io.__all__)
__all__.extend(lbl.__all__)
__all__.extend(los.__all__)
__all__.extend(phys.__all__)
__all__.extend(spectrum.__all__)
__all__.extend(tools.__all__)
