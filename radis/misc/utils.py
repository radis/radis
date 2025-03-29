# -*- coding: utf-8 -*-
"""

@author: Erwan

-------------------------------------------------------------------------------


"""


import importlib
import inspect
import os
import sys
from fnmatch import translate
from os.path import basename, dirname, join
from re import IGNORECASE, compile


def getProjectRoot():
    """Return the full path of the project root."""

    return dirname(dirname(__file__))


def import_from_module(module, name):
    """Import object 'name' from module 'module' raises AttributeError if name
    doesnt exist.

    Parameters
    ----------

    module, name: str
        module path, object name
    """
    impmodule = importlib.import_module(module)
    return getattr(impmodule, name)


class Chdir:
    """because we need to change directory to get into the RADIS folder to find
    the Git files and version, when imported from another program. This class
    then ensures we get back in the correct directory.

    Examples
    --------

    Do::

        cd = Chdir(os.path.dirname(__file__))
        try:
            (...)
        except:
            (...)
        finally:
            (...)
            cd.__del__()
    """

    def __init__(self, newPath):
        self.savedPath = os.getcwd()
        os.chdir(newPath)

    def __del__(self):
        os.chdir(self.savedPath)


# %%
# ==============================================================================
# Function arguments
# ==============================================================================


class Default:
    """Contains a value. Used to know whether a function argument equal to its
    default value was explicitly given by the user or not. This allows to
    prevent user errors.

    Examples
    --------

    Check if a value is Default::

        from radis.misc.utils import Default
        a = Default("42")
        isinstance(a, Default)
        >>> True

        a.value
        >>> 42
    """

    def __repr__(self):
        return "{}".format(self.value)

    def __init__(self, value):
        self.value = value


def getarglist(function):
    """Get list of arguments in a function.

    See https://stackoverflow.com/a/41188411/5622825
    """

    if sys.version_info[0] == 2:
        from inspect import getargspec

        return getargspec(function).args

    else:
        from inspect import signature

        return list(signature(function).parameters)

def calculate_auto_chunksize(nlines, ngridpoints, available_mem=None):
    """Auto-determine chunksize to balance speed/memory.
    Defaults to using 70% of available RAM.
    """
    if available_mem is None:
        import psutil
        available_mem = psutil.virtual_memory().available / 1e9  # GB

    # Empirical memory usage (1e6 lines × 1e5 gridpoints ≈ 1GB)
    MEM_PER_PAIR = 1e-11  # GB per (line×gridpoint)
    
    safe_mem = available_mem * 0.7
    chunksize = int(safe_mem / (ngridpoints * MEM_PER_PAIR))
    
    # Constrain between 1k lines and full dataset
    return max(1000, min(chunksize, nlines))

def get_default_arg(func, arg):
    """Get default value of argument ``arg`` in function ``func``

    Adapted from https://stackoverflow.com/questions/12627118/get-a-function-arguments-default-value
    """
    signature = inspect.signature(func)
    items = dict(signature.parameters.items())
    if not arg in items:
        raise ValueError("Function {0} has no argument `{1}`".format(func, arg))
    elif items[arg].default is inspect.Parameter.empty:
        raise ValueError("No default value for argument `{0}` in {1}".format(arg, func))
    else:
        return items[arg].default


# %% Other stuff


class DatabankNotFound(FileNotFoundError):
    """Used when a line database is not found in radis.rc."""

    pass


# %%
# ==============================================================================
# Optional packages
# ==============================================================================


class NotInstalled(object):
    """A class to deal with optional packages. Will raise an error only if the
    package is used (but not if imported only)

    Examples
    --------

    ::

        some function = NotInstalled("you should install C drivers")
        a = some_function   # no error
        a()             # will raise the NotInstalled error and display the message

    """

    def __init__(self, name, info=""):
        self.__name = name
        self.__info = info

    def __getattr__(self, item):
        raise ImportError(
            "The {0} package is required to use this "
            "feature. {1}".format(self.__name, self.__info)
        )

    def __call__(self, *args, **kwargs):
        raise ImportError(
            "The {0} package is required to use this "
            "feature. {1}".format(self.__name, self.__info)
        )


not_installed_vaex_args = (
    "vaex",
    "You must install Vaex to use these features. Vaex is a fast, "
    + "memory-mapped DataFrame library. However is not available yet on latest Python versions. "
    + "Use Pytables (slower) as an alternative in your Radis.json config file. To use Pytables, set "
    + '"MEMORY_MAPPING_ENGINE": "pytables" and "DATAFRAME_ENGINE": "pandas"',
)
not_installed_nvidia_args = (
    "nvidia-cufft",
    "Nvidia was not installed on your computer. `nvidia-cufft` is"
    + "not available on Mac OS. ALternative solutions are being developped "
    + "with Vulkan, see https://github.com/radis/radis/pull/624. (latest update Oct. 2024). ",
)


def get_files_from_regex(path):
    """Returns a list of absolute paths of all the files whose names match the
    input regular expression."""
    directory_name = dirname(path)
    if directory_name == "":
        directory_name = "."
    regex = basename(path)

    file_names = []

    pattern = compile(translate(regex), IGNORECASE)

    for file in os.listdir(directory_name):
        if pattern.fullmatch(file):
            file_names.append(join(directory_name, file))

    return file_names


# %% Test


def _test(*args, **kwargs):

    print("Project root:", getProjectRoot())

    return True


if __name__ == "__main__":
    _test()
