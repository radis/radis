# -*- coding: utf-8 -*-
"""
Default values for global Radis parameters

"""
# to be added in the new radis.json configuration file


DEBUG_MODE = False
"""bool: change this at runtime with::

    import radis
    radis.DEBUG_MODE = True

Use the :py:func:`~radis.misc.debug.printdbg` function in ``radis.misc``, typically with::

    if __debug__: printdbg(...)

so that printdbg are removed by the Python preprocessor when running in
optimize mode::

    python -O *.py
"""

AUTO_UPDATE_SPEC = False
"""bool: experimental feature
used to autoupdate .spec files to the latest format, by simply saving
them again once they're loaded and fixed.
Warning! Better have a copy of your files before that, or a way to regenerate
them.

Examples
--------

Add to the top of your script (once is enough!)::

    import radis
    radis.AUTO_UPDATE_SPEC = True

See Also
--------
:py:func:`~radis.tools.database._update_to_latest_format`
"""

OLDEST_COMPATIBLE_VERSION = "0.9.1"
"""str: forces to regenerate cache files that were created in a previous version

See Also
--------
:py:func:`~radis.io.cache_files.load_h5_cache_file`
"""

USE_CYTHON = True
"""bool: try to use Cython functions when possible

See Also
--------
:py:func:`~radis.misc.arrays.add_at`
"""

WARN_THRESHOLD = 3
"""float: to determine the optimal value
of wstep using minimum FWHM value of spectrum.
Makes sure there are enough gridpoints per line.

See Also
--------
:py:func:`~radis.lbl.broadening._check_accuarcy`
"""
