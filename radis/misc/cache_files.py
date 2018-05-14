# -*- coding: utf-8 -*-
"""
Tools to deal with cache files
"""

import h5py
import radis
from warnings import warn
from os.path import exists
from radis.misc.basics import compare_dict, is_float


class DeprecatedFileError(DeprecationWarning):
    pass


LAST_COMPATIBLE_VERSION = '0.1.18'
'''str: forces to regenerate cache files that were created in a previous version'''

# Just make sure LAST_BACKWARD_COMPATIBLE_VERSION is valid
assert radis.__version__ >= LAST_COMPATIBLE_VERSION

# Utils


def check_not_deprecated(file, metadata, current_version=None, last_compatible_version=LAST_COMPATIBLE_VERSION):
    ''' Make sure cache file is not deprecated: checks that ``metadata`` is the same,
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
        Default :data:`~radis.misc.cache_files.LAST_COMPATIBLE_VERSION`

    '''

    # Get attributes (metadata+version)
    hf = h5py.File(file, 'r')
    try:
        attrs = dict(hf.attrs)
    except OSError:
        attrs = {}
    finally:
        hf.close()

    # Raise an error if version is not found
    try:
        file_version = attrs.pop('version')
    except KeyError:
        raise DeprecatedFileError('File {0} has been generated in a deprecated '.format(file) +
                                  'version. Delete it to regenerate it on next run')

    # Get current version
    if current_version is None:
        current_version = radis.__version__

    # If file version is anterior to a major change
    # ... Update here versions afterwhich Deprecated cache file is not safe
    # ... (example: a key name was changed)
    if file_version < last_compatible_version:
        raise DeprecatedFileError('File {0} has been generated in a deprecated '.format(file) +
                                  'version ({0}). Last compatible version is {1}. '.format(
            file_version, last_compatible_version) +
            'Delete the file to regenerate it on next run')

    # If file version is outdated: Warning, but no error
    if current_version > file_version:
        warn(DeprecationWarning('File {0} has been generated in '.format(file) +
                                'a deprecated version ({0}) compared to current ({1})'.format(
            file_version, current_version) +
            '. Delete it to regenerate it on next run'))
        out = False
    elif current_version == file_version:
        out = True
    else:
        raise ValueError('File generated with a future version?')

    # Compare metadata
    metadata = _h5_compatible(metadata)
    out, compare_string = compare_dict(
        metadata, attrs, verbose=False, return_string=True)
    if out != 1:
        raise DeprecatedFileError('Metadata in file {0} dont match '.format(file) +
                                  'expected values. See comparison below:' +
                                  '\n\tExpected\tFile\n{0}'.format(compare_string))

    return out


def _h5_compatible(a_dict):
    ''' Make dictionary ``a_dict`` h5 compatible '''
    out = {}
    for k, v in a_dict.items():
        if v is None:
            continue   # dont store None
        elif is_float(v):
            out[k] = v
        else:
            out[k] = str(v)   # convert to str
    return out


def save_to_hdf(df, fname, metadata, version=None, key='df', overwrite=True):
    ''' Save energy levels to HDF5 file. Add metadata and version 

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
        See :data:`~radis.misc.cache_files.LAST_COMPATIBLE_VERSION`

    key: str
        dataset name. Default ``'df'`` 

    overwrite: boolean
        if ``True``, overwrites file. Else, raise an error if it exists.

    Notes
    -----

    ``None`` values are not stored

    '''

    assert fname.endswith('.h5')
    assert 'version' not in metadata

    # Update metadata format
    metadata = _h5_compatible(metadata)

    # Overwrite file
    if exists(fname) and not overwrite:
        raise ValueError('File exist: {0}'.format(fname))
    hf = h5py.File(fname, 'w')

    try:

        # Start by adding version
        if version is None:
            version = radis.__version__
        hf.attrs['version'] = version

        # Add metadata
        for k, v in metadata.items():
            hf.attrs[k] = v

    except:
        raise

    finally:
        hf.close()

    # now export dataframe
    df.to_hdf(fname, key, format='fixed', mode='a',
              complevel=1, complib='blosc')
#    df.to_hdf(fname, 'df', format='fixed', mode='a')  # would be 10-20% faster, but take 2x more space


def filter_metadata(arguments, discard_variables=['self', 'verbose']):
    ''' Filter arguments (created with  ``locals()`` at the beginning
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

    '''

    metadata = {k: v for (k, v) in arguments.items() if not k.startswith('_')}
    metadata = {k: v for (k, v) in metadata.items()
                if k not in discard_variables}
    metadata = {k: v for (k, v) in metadata.items() if v is not None}

    return metadata
