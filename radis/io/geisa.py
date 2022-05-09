# -*- coding: utf-8 -*-
"""
Summary
-----------------------

HITRAN database parser

-----------------------


"""


# import os
# import sys
import time
from collections import OrderedDict
from os.path import exists, getmtime

import radis
from radis.io.cache_files import cache_file_name, load_h5_cache_file, save_to_hdf
from radis.io.tools import drop_object_format_columns, parse_hitran_file

# %% Parsing functions

# General case : GEISA 2020
# Provided by Thibault Delahaye
# fmt: off
columns_GEISA = OrderedDict(
    [
        # name  # format  # type  # description  # unit
        ("wav", ("a12", float, "vacuum wavenumber", "cm-1")),
        ("int", ("a11", float, "intensity at 296K", "cm-1/(molecule/cm-2)")),
        ("airbrd", ("a6", float, "air-broadened half-width at 296K", "cm-1.atm-1")),
        ("El", ("a10", float, "lower-state energy", "cm-1")),
        ("globu", ("a25", str, "electronic and vibrational global upper quanta", "")),
        ("globl", ("a25", str, "electronic and vibrational global lower quanta", "")),
        ("locu", ("a15", str, "electronic and vibrational local upper quanta", "")),
        ("locl", ("a15", str, "electronic and vibrational local lower quanta", "")),
        ("Tdpgair", ("a4", float, "temperature-dependance exponent for Gamma air", "")),
        ("iso", ("a3", int, "GEISA isotope number", "")),
        ("mol", ("a3", int, "GEISA molecular number", "")),
        ("id", ("a3", str, "Internal GEISA code for the data identification", "")),
        ("idH", ("a2", int, "Hitran molecular number", "")),
        ("isoH", ("a1", int, "Hitran isotope number", "")),
        ("A", ("a10", float, "Einstein A coefficient", "s-1")),
        ("selbrd", ("a7", float, "self-broadened half-width at 296K", "cm-1.atm-1")),
        ("Pshft", ("a9", float, "air pressure-induced line shift at 296K", "cm-1.atm-1")),
        ("Tdppair", ("a6", float, "temperature-dependance exponent for air pressure-induced line shift", "")),
        ("ierrA", ("a10", float, "estimated accuracy on the line position", "cm-1")),
        ("ierrB", ("a11", float, "estimated accuracy on the intensity of the line", "cm-1/(molecule/cm-2)")),
        ("ierrC", ("a6", float, "estimated accuracy on the air collision halfwidth", "cm-1.atm-1")),
        ("ierrF", ("a4", float, "estimated accuracy on the temperature dependence coefficient of the air-broadening halfwidth", "")),
        ("ierrO", ("a9", float, "estimated accuracy on the air pressure shift of the line transition at 296K", "cm-1.atm-1")),
        ("ierrR", ("a6", float, "estimated accuracy on the temperature dependence coefficient of the air pressure shift", "")),
        ("ierrN", ("a7", float, "estimated accuracy on the self-broadened at 296K", "cm-1.atm-1")),
        ("Tdpgself", ("a4", float, "temperature-dependance exponent for self-broadening halfwidth", "")),
        ("ierrS", ("a4", float, "estimated accuracy on the temperature dependence coefficient of the self-broadening halfwidth", "")),
        ("Pshfts", ("a8", float, "self pressure-induced line shift at 296K", "cm-1.atm-1")),
        ("ierrT", ("a8", float, "estimated accuracy on the self-pressure shift of the line transition at 296K", "cm-1.atm-1")),
        ("Tdppself", ("a4", float, "temperature-dependance exponent for self pressure-induced line shift", "")),
        ("ierrU", ("a4", float, "estimated accuracy on the temperature dependence coefficient of the self pressure shift", "")),
    ]
)
""" OrderedDict: parsing order of GEISA2020 format """

# The columns required for equilibrium calculations only
equilibrium_columns = OrderedDict(
    [
        # name  # format  # type  # description  # unit
        ("wav", ("a12", float, "vacuum wavenumber", "cm-1")),
        ("int", ("a11", float, "intensity at 296K", "cm-1/(molecule/cm-2)")),
        ("airbrd", ("a6", float, "air-broadened half-width at 296K", "cm-1.atm-1")),
        ("iso", ("a3", int, "GEISA isotope number", "")),
        ("mol", ("a3", int, "GEISA molecular number", "")),
        ("id", ("a3", str, "Internal GEISA code for the data identification", "")),
        ("selbrd", ("a7", float, "self-broadened half-width at 296K", "cm-1.atm-1")),
        ("Pshft", ("a9", float, "air pressure-induced line shift at 296K", "cm-1.atm-1")),
        ("Pshfts", ("a8", float, "self pressure-induced line shift at 296K", "cm-1.atm-1")),
    ]
)
# fmt: on


def gei2df(
    fname,
    cache=True,
    load_columns=None,
    verbose=True,
    drop_non_numeric=True,
    load_wavenum_min=None,
    load_wavenum_max=None,
    engine="pytables",
):
    """Convert a GEISA [1]_ file to a Pandas dataframe.
    Parameters
    ----------
    fname: str
        GEISA file name.
    cache: boolean, or 'regen'
        if ``True``, a pandas-readable HDF5 file is generated on first access,
        and later used. This saves on the datatype cast and conversion and
        improves performances a lot (but changes in the database are not
        taken into account). If ``False``, no database is used. If 'regen', temp
        file are reconstructed. Default ``True``.
    load_columns: list
        columns to load. If ``None``, loads everything.
        .. note::
            this is only relevant when loading from a cache file. To generate
            the cache file, all columns are loaded anyway.
    Other Parameters
    ----------------
    drop_non_numeric: boolean
        if ``True``, non numeric columns are dropped. This improves performances,
        but make sure all the columns you need are converted to numeric formats
        before hand. Default ``True``. Note that if a cache file is loaded it
        will be left untouched.
    load_wavenum_min, load_wavenum_max: float
        if not ``'None'``, only load the cached file if it contains data for
        wavenumbers above/below the specified value. See :py:func`~radis.io.cache_files.load_h5_cache_file`.
        Default ``'None'``.
    engine: 'pytables', 'vaex'
        format for Hdf5 cache file, `pytables` by default.
    Returns
    -------
    df: pandas Dataframe
        dataframe containing all lines and parameters.
    Notes
    -----
    GEISA Database 2020 release can be downloaded from [2]_
    References
    ----------
    .. [1] `The 2020 edition of the GEISA spectroscopic database, Thibault Delahaye et al., 2021 <https://www.sciencedirect.com/science/article/abs/pii/S0022285221000928>`_
    .. [2] `GEISA Database 2020 release <https://geisa.aeris-data.fr/interactive-access/?db=2020&info=ftp>`_
    See Also
    --------
    :func:`~radis.io.hitran.hit2df`
    :func:`~radis.io.cdsd.cdsd2df`
    """

    # Try to access last modification time of original file
    metadata = {}
    metadata["last_modification"] = time.ctime(getmtime(fname))
    # Make sure the wavenums are legit
    if load_wavenum_min and load_wavenum_max:
        assert load_wavenum_min < load_wavenum_max
    # Verbose
    if verbose >= 2:
        print("Opening file {0}, cache={1})".format(fname, cache))
        print("Last Modification time: {0}".format(metadata["last_modification"]))

    parse_columns = columns_GEISA

    # Attempt to use cache file
    fcache = cache_file_name(fname, engine=engine)
    if cache and exists(fcache):
        relevant_if_metadata_above = (
            {"wavenum_max": load_wavenum_min} if load_wavenum_min else {}
        )  # not relevant if wavenum_max of file is < wavenum min required
        relevant_if_metadata_below = (
            {"wavenum_min": load_wavenum_max} if load_wavenum_max else {}
        )  # not relevant if wavenum_min of file is > wavenum max required
        df = load_h5_cache_file(
            fcache,
            cache,
            columns=load_columns,
            valid_if_metadata_is=metadata,
            relevant_if_metadata_above=relevant_if_metadata_above,
            relevant_if_metadata_below=relevant_if_metadata_below,
            current_version=radis.__version__,
            last_compatible_version=radis.config["OLDEST_COMPATIBLE_VERSION"],
            verbose=verbose,
            engine=engine,
        )
        if df is not None:
            return df

    # If cache files are not found, commence reading of full file
    df = parse_hitran_file(fname, parse_columns, is_geisa=True)

    # Remove non numerical attributes
    if drop_non_numeric:
        # The function replace_PQR_with_m101 below will be developed
        # in later updates for GEISA's non-equilibrium calculations
        # replace_PQR_with_m101(df)

        # Basically ripping out 4 quanta columns
        df = drop_object_format_columns(df, verbose=verbose)

    # Generate cache file for later use
    if cache:
        new_metadata = {
            # Last modification time of the original file :
            "last_modification": time.ctime(getmtime(fname)),
            "wavenum_min": df.wav.min(),
            "wavenum_max": df.wav.max(),
        }
        if verbose:
            print(
                "Generating cache file {0} with metadata :\n{1}".format(
                    fcache, new_metadata
                )
            )
        try:
            save_to_hdf(
                df,
                fcache,
                metadata=new_metadata,
                version=radis.__version__,
                key="df",
                overwrite=True,
                verbose=verbose,
                engine=engine,
            )
        except PermissionError:
            if verbose:
                print(
                    "An error occured in cache file generation. Lookup access rights."
                )
            pass

    return df


# Test cases for
if __name__ == "__main__":
    from radis.test.test_io import test_geisa

    print("Testing GEISA: ")
    test_geisa()
