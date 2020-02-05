# -*- coding: utf-8 -*-
"""
Tools to deal with HDF5 cache files
HDF5 cache files are used to cache Energy Database files, and Line Database
files, and yield a much faster access time. 

Routine Listing
---------------


- :func:`~radis.misc.cache_file.check_cache_file`
- :func:`~radis.misc.cache_file.check_not_deprecated`
- :func:`~radis.misc.cache_file.save_to_hdf`

See Also
--------

:py:func:`~radis.io.hitran.hit2df`, 
:py:func:`~radis.io.cdsd.cdsd2df`

-------------------------------------------------------------------------------


"""
# TODO: start using .feather format. Faster than HDF5 for pandas Dataframe,
# and can be multithreaded.
# see:
# https://gist.github.com/gansanay/4514ec731da1a40d8811a2b3c313f836
# and pd.read_feather(file, nthreads=3)

from __future__ import absolute_import, print_function, division

# Note: don't import unicode_literals because it breaks the df.to_hdf of
# save_to_hdf because of a stupid unicode/str error in Python 2.7
import os
import h5py
import radis
from warnings import warn
from os.path import exists
from radis import OLDEST_COMPATIBLE_VERSION
from radis.misc.basics import compare_dict, is_float
from radis.misc.printer import printr, printm
import pandas as pd
from packaging.version import parse


class DeprecatedFileError(DeprecationWarning):
    pass


"""str: forces to regenerate cache files that were created in a previous version"""

# Just make sure LAST_BACKWARD_COMPATIBLE_VERSION is valid
assert radis.__version__ >= OLDEST_COMPATIBLE_VERSION

# Utils


def load_h5_cache_file(
    cachefile,
    use_cached,
    metadata,
    current_version,
    last_compatible_version=OLDEST_COMPATIBLE_VERSION,
    verbose=True,
):
    """ Function to load a h5 cache file
    
    Parameters
    ----------
    
    cachefile: str
        cache file path
    
    use_cached: str
        use cache file if value is not ``False``:
        
        - if ``True``, use (and generate if doesnt exist) cache file.
        - if ``'regen'``, delete cache file (if exists) so it is regenerated
        - if ``'force'``, use cache file and raises an error if it doesnt exist
        
        if using the cache file, check if the file is deprecated. If it 
        is deprecated, regenerate the file unless ``'force'`` was used
        (in that case, raise an error)
        
    metadata: dict
        values are compared to cache file attributes. If they dont match,
        the file is considered deprecated. See ``use_cached`` to know
        how to handle deprecated files
        
    current_version: str
        version is compared to cache file version (part of attributes). 
        If current version is superior, a simple warning is triggered.
        
    last_compatible_version: str
        if file version is inferior to this, file is considered deprecated. 
        See ``use_cached`` to know how to handle deprecated files. 
        Default :data:`~radis.OLDEST_COMPATIBLE_VERSION`. 
        
    Returns
    -------
    
    df: pandas DataFrame, or None
        None if no cache file was found, or if it was deleted
        
        
    See Also
    --------
    
    :data:`~radis.OLDEST_COMPATIBLE_VERSION`
        
    """

    # 1. know if we have to load the file
    if not use_cached:
        return None
    elif use_cached == "regen" and exists(cachefile):
        os.remove(cachefile)
        if verbose:
            printm("Deleted h5 cache file : {0}".format(cachefile))
        return None

    # 2. check the file is here
    if not exists(cachefile):
        if use_cached == "force":
            raise ValueError("Cache file {0} doesnt exist".format(cachefile))
        else:
            return None  # File doesn't exist. It's okay.

    # 3. read file attributes to know if it's deprecated
    try:
        check_not_deprecated(
            cachefile,
            metadata,
            current_version=current_version,
            last_compatible_version=last_compatible_version,
        )
    # ... if deprecated, raise an error only if 'force'
    except DeprecatedFileError as err:
        if use_cached == "force":
            raise err
        else:
            if verbose:
                printr(
                    "File {0} deprecated:\n{1}\nDeleting it!".format(
                        cachefile, str(err)
                    )
                )
            os.remove(cachefile)
            return None
    # 4. File is not not deprecated: read the content.
    else:
        df = None
        if verbose >= 2:
            printm("Reading cache file ({0})".format(cachefile))
        try:
            df = pd.read_hdf(cachefile, "df")
        except KeyError as err:  # An error happened during file reading.
            # Fail safe by deleting cache file (unless we explicitely wanted it
            # with 'force')
            if use_cached == "force":
                raise
            else:
                if verbose:
                    printr(
                        "An error happened during cache file reading "
                        + "{0}:\n{1}\n".format(cachefile, str(err))
                        + "Deleting cache file to regenerate it"
                    )
                os.remove(cachefile)
                df = None

    return df


def get_cache_file(fcache, verbose=True):
    """ Load HDF5 cache file 
    
    Parameters
    ----------
    
    fcache: str
        file name
    
    Other Parameters
    ----------------
    
    verbose: bool
        If >=2, also warns if non numeric values are present (it would make 
        calculations slower)

    Notes
    -----
    
    we could start using FEATHER format instead. See notes in cache_files.py 
    """

    # Load file

    df = pd.read_hdf(fcache, "df")

    # Check file

    # ... 'object' columns slow everything down (not fixed format strings!)
    if verbose >= 2:
        _warn_if_object_columns(df, fcache)

    return df


def check_cache_file(fcache, use_cached=True, metadata={}, verbose=True):
    """ Quick function that check status of cache file generated by RADIS:

    Parameters
    ----------
    
    fcache: str
        cache file name
        
    use_cached: ``True``, ``False``, ``'force'``, ``'regen'``
        see notes below. Default ``True``.
        
    metadata: dict
        attributes to check
        
    verbose: boolean
        print stuff
    
    Notes
    -----
    
    The function first checks the existence of ``fcache``. What is does depends
    on the value of ``use_cached``:

    - if ``True``, check it exists. 
    - if ``'regen'``, delete cache file to regenerate it later. 
    - if ``'force'``, raise an error if file doesnt exist. 

    Then look if it is deprecated (we just look at the attributes, the file 
    is never fully read). Deprecation is done by :py:func:`~radis.misc.cache_files.check_not_deprecated`
    comparing the ``metadata=`` content.

    - if deprecated, deletes it to regenerate later unless 'force' was used
    
    See Also
    --------
    
    :py:func:`~radis.misc.cache_files.check_not_deprecated`
    """

    # Test existence of file:
    if not use_cached:
        return  # we dont want a cache file, no need to test it
    elif use_cached == "regen":
        if exists(fcache):
            os.remove(fcache)
            if verbose:
                print(("Deleted h5 cache file : {0}".format(fcache)))
    elif use_cached == "force":
        if not exists(fcache):
            raise ValueError("Cache file {0} doesnt exist".format(fcache))
    else:  # use_cached == True
        pass  # just use the file as is

    # If file is still here, test if it is deprecated:
    # (we just read the attributes, the file is never fully read)
    if exists(fcache):
        if verbose:
            print(("Using cache file: {0}".format(fcache)))
        try:
            check_not_deprecated(
                fcache,
                metadata=metadata,
                current_version=radis.__version__,
                last_compatible_version=OLDEST_COMPATIBLE_VERSION,
            )
        except DeprecatedFileError as err:
            if use_cached == "force":
                raise
            else:  # delete file to regenerate it in the end of the script
                if verbose:
                    printr(
                        "File {0} deprecated:\n{1}\nDeleting it!".format(
                            fcache, str(err)
                        )
                    )
                os.remove(fcache)

    return


def check_not_deprecated(
    file,
    metadata,
    current_version=None,
    last_compatible_version=OLDEST_COMPATIBLE_VERSION,
):
    """ Make sure cache file is not deprecated: checks that ``metadata`` is the same,
    and that the version under which the file was generated is valid.

    Parameters
    ----------

    file: str
        a `` .h5`` cache file for Energy Levels 

    metadata: dict
        list of variables used to create the file. If the values dont match,
        an error is raised. 

    current_version: str, or ``None``
        current version number. If the file was generated in a previous version
        a warning is raised. If ``None``, current version is read from 
        :data:`radis.__version__`.

    last_backward_compatible_version: str
        If the file was generated in a non-compatible version, an error is raised.
        Default :py:data:`~radis.OLDEST_COMPATIBLE_VERSION`

    """

    # Get attributes (metadata+version)
    hf = h5py.File(file, "r")
    try:
        attrs = dict(hf.attrs)
    except OSError:
        attrs = {}
    finally:
        hf.close()

    # Raise an error if version is not found
    try:
        file_version = attrs.pop("version")
    except KeyError:
        raise DeprecatedFileError(
            "File {0} has been generated in a deprecated ".format(file)
            + "version. Delete it to regenerate it on next run"
        )

    # Get current version
    if current_version is None:
        current_version = radis.__version__

    # If file version is anterior to a major change
    # ... Update here versions afterwhich Deprecated cache file is not safe
    # ... (example: a key name was changed)
    if parse(file_version) < parse(last_compatible_version):
        raise DeprecatedFileError(
            "File {0} has been generated in a deprecated ".format(file)
            + "version ({0}). Oldest compatible version is {1}. ".format(
                file_version, last_compatible_version
            )
            + "Delete the file to regenerate it on next run"
        )

    # If file version is outdated: Warning, but no error
    if parse(current_version) > parse(file_version):
        warn(
            DeprecationWarning(
                "File {0} has been generated in ".format(file)
                + "a deprecated version ({0}) compared to current ({1})".format(
                    file_version, current_version
                )
                + ". Delete it to regenerate it on next run"
            )
        )
        out = False
    elif parse(current_version) == parse(file_version):
        out = True
    else:
        raise ValueError(
            "Cache file ({0}) generated with a future version ({1} > {2})? ".format(
                file, file_version, current_version
            )
            + "Do you own a DeLorean? Delete the file manually if you understand what happened"
        )

    # Compare metadata
    metadata = _h5_compatible(metadata)
    out, compare_string = compare_dict(
        metadata, attrs, verbose=False, return_string=True
    )
    if out != 1:
        raise DeprecatedFileError(
            "Metadata in file {0} dont match ".format(file)
            + "expected values. See comparison below:"
            + "\n\tExpected\tFile\n{0}".format(compare_string)
        )

    return out


def _h5_compatible(a_dict):
    """ Make dictionary ``a_dict`` h5 compatible """
    out = {}
    for k, v in a_dict.items():
        if v is None:
            continue  # dont store None
        elif is_float(v):
            out[k] = v
        else:
            out[k] = str(v)  # convert to str
    return out


def _warn_if_object_columns(df, fname):
    """  'object' columns slow everything down (not fixed format strings!) """
    objects = [k for k, v in df.dtypes.items() if v == object]
    if len(objects) > 0:
        warn(
            "Dataframe in {0} contains `object` format columns: {1}. ".format(
                fname, objects
            )
            + "Operations will be slower. Try to convert them to numeric."
        )


def save_to_hdf(
    df, fname, metadata, version=None, key="df", overwrite=True, verbose=True
):
    """ Save energy levels to HDF5 file. Add metadata and version 

    Parameters
    ----------

    df: a pandas DataFrame
        data will be stored in the key ``'df'`` 

    fname: str
        ``.h5`` file where to store. 

    metadata: dict
         dictionary of values that were used to generate the DataFrame. Metadata
         will be asked again on file load to ensure it hasnt changed. ``None`` 
         values will not be stored.

    version: str, or ``None``
        file version. If ``None``, the current :data:`radis.__version__` is used. 
        On file loading, a warning will be raised if the current version is 
        posterior, or an error if the file version is set to be uncompatible.
        See :py:data:`~radis.OLDEST_COMPATIBLE_VERSION`

    key: str
        dataset name. Default ``'df'`` 

    overwrite: boolean
        if ``True``, overwrites file. Else, raise an error if it exists.

    Other Parameters
    ----------------
    
    verbose: bool
        If >=2, also warns if non numeric values are present (it would make 
        calculations slower)

    Notes
    -----

    ``None`` values are not stored

    """

    # Check file
    assert fname.endswith(".h5")
    assert "version" not in metadata
    # ... 'object' columns slow everything down (not fixed format strings!)
    if verbose >= 2:
        _warn_if_object_columns(df, fname)

    # Update metadata format
    metadata = _h5_compatible(metadata)

    # Overwrite file
    if exists(fname) and not overwrite:
        raise ValueError("File exist: {0}".format(fname))
    hf = h5py.File(fname, "w")

    try:

        # Start by adding version
        if version is None:
            version = radis.__version__
        hf.attrs["version"] = version

        # Add metadata
        for k, v in metadata.items():
            hf.attrs[k] = v

    except:
        raise

    finally:
        hf.close()

    # now export dataframe
    df.to_hdf(fname, key, format="fixed", mode="a", complevel=1, complib="blosc")


#    df.to_hdf(fname, 'df', format='fixed', mode='a')  # would be 10-20% faster, but take 2x more space


def filter_metadata(arguments, discard_variables=["self", "verbose"]):
    """ Filter arguments (created with  ``locals()`` at the beginning
    of the script) to extract metadata. Metadata is stored as attributes
    in the cached file:

    - remove variables in ``discard_variables``
    - remove variables that start with ``'_'``
    - remove varibles whose value is ``None``

    Parameters
    ----------

    arguments: dict
        list of local variables. For instance::

            arguments = locals()

    discard_variables: list of str
        variable names to discard

    Returns
    -------

    metadata: dict
        a (new) dictionary built from arguments by removing ``discard_variables``
        and variables starting with ``'_'``

    Examples
    --------

    How to get only function argument::

        def some_function(*args):
            metadata = locals()     # stores only function arguments because it's the first line

            ...

            metadata = filter_metadata(metadata)
            save_to_hdf(df, fname, metadata=metadata)

            ...

    """

    metadata = {k: v for (k, v) in arguments.items() if not k.startswith("_")}
    metadata = {k: v for (k, v) in metadata.items() if k not in discard_variables}
    metadata = {k: v for (k, v) in metadata.items() if v is not None}

    return metadata
