# -*- coding: utf-8 -*-
"""
Summary
-------

HITRAN database parser

The quantum formats follow the supplementary material of Gordon et al.,
`The HITRAN2020 molecular spectroscopic database (2022)
<https://doi.org/10.1016/j.jqsrt.2021.107949>`_, with the newer HITRAN2024
formats where they differ (Gordon et al., 2026,
`HITRAN2024 <https://doi.org/10.1016/j.jqsrt.2026.109807>`_).


Routine Listing
---------------

- :func:`~radis.api.hitranapi.hit2df`
- :func:`~radis.api.hitranapi.parse_local_quanta`
- :func:`~radis.api.hitranapi.parse_global_quanta`
- :func:`~radis.api.hitranapi.hitranxsc`




"""

import os
import sys

# from radis.test.utils import getTestFile
import time
import warnings
from collections import OrderedDict
from os.path import abspath, exists, expanduser, getmtime, join, split

import pandas as pd
from numpy import int64

from radis.db.classes import (  # get_molecule_identifier,
    HITRAN_CLASS1,
    HITRAN_CLASS2,
    HITRAN_CLASS3,
    HITRAN_CLASS4,
    HITRAN_CLASS5,
    HITRAN_CLASS6,
    HITRAN_CLASS7,
    HITRAN_CLASS8,
    HITRAN_CLASS9,
    HITRAN_GROUP1,
    HITRAN_GROUP2,
    HITRAN_GROUP3,
    HITRAN_GROUP4,
    HITRAN_GROUP5,
    HITRAN_GROUP6,
    HITRAN_GROUP7,
    get_molecule,
)
from radis.misc.config import getDatabankEntries
from radis.misc.warning import DatabaseAlreadyExists

try:
    from .cache_files import load_h5_cache_file, save_to_hdf
    from .dbmanager import DatabaseManager
    from .hdf5 import DataFileManager
    from .tools import (
        drop_object_format_columns,
        parse_hitran_file,
        replace_PQR_with_m101,
    )
except ImportError:
    if __name__ == "__main__":  # running from this file, as a script
        from radis.api.cache_files import load_h5_cache_file, save_to_hdf
        from radis.api.dbmanager import DatabaseManager
        from radis.api.hdf5 import DataFileManager
        from radis.api.tools import (
            drop_object_format_columns,
            parse_hitran_file,
            replace_PQR_with_m101,
        )
    else:
        raise

# %% Parsing functions

# General case : HITRAN 2004
# fmt: off
columns_2004 = OrderedDict(
    [
        # name    # format # type  # description                                 # unit
        ("id", ("a2", int, "Molecular number", "")),
        ("iso", ("a1", int, "isotope number", "")),
        ("wav", ("a12", float, "vacuum wavenumber", "cm-1")),
        ("int", ("a10", float, "intensity at 296K", "cm-1/(molecule/cm-2)")),
        ("A", ("a10", float, "Einstein A coefficient", "s-1")),
        ("airbrd", ("a5", float, "air-broadened half-width at 296K", "cm-1.atm-1")),
        ("selbrd", ("a5", float, "self-broadened half-width at 296K", "cm-1.atm-1")),
        ("El", ("a10", float, "lower-state energy", "cm-1")),
        ("Tdpair", ("a4", float, "temperature-dependance exponent for Gamma air", "")),
        ("Pshft", ("a8", float, "air pressure-induced line shift at 296K", "cm-1.atm-1")),
        ("globu", ("a15", str, "electronic and vibrational global upper quanta", "")),
        ("globl", ("a15", str, "electronic and vibrational global lower quanta", "")),
        ("locu", ("a15", str, "electronic and vibrational local upper quanta", "")),
        ("locl", ("a15", str, "electronic and vibrational local lower quanta", "")),
        ("ierr", ("a6", str, "ordered list of indices corresponding to uncertainty estimates of transition parameters", "")),
        ("iref", ("a12", str, "ordered list of reference identifiers for transition parameters", "")),
        ("lmix", ("a1", str, "flag indicating the presence of additional data and code relating to line-mixing", "")),
        ("gp", ("a7", float, "upper state degeneracy", "")),
        ("gpp", ("a7", float, "lower state degeneracy", "")),
    ]
)
""" OrderedDict: parsing order of HITRAN 2004 format """
# fmt: on


def get_hitran_isotopologues(molecule, max_isotopologue=9):
    """Isotopologues of ``molecule`` available in HITRAN, read from HAPI's
    isotopologue table.

    Isotopologues above ``max_isotopologue`` (ex: 10-12 for CO2) are left out:
    they have negligible abundance and were never downloaded.

    Examples
    --------
    ::

        get_hitran_isotopologues("CO")
        >>> [1, 2, 3, 4, 5, 6]
    """
    from hapi import ISO

    from radis.db.classes import get_molecule_identifier

    mol_id = get_molecule_identifier(molecule)
    return [iso for (m, iso) in sorted(ISO) if m == mol_id and iso <= max_isotopologue]


PARAMETER_GROUPS_HITRAN = {
    "par_line": "PARLIST_DOTPAR",
    "id": "PARLIST_ID",
    "standard": "PARLIST_STANDARD",
    "labels": "PARLIST_LABELS",
    "voigt": "PARLIST_VOIGT_ALL",
    "ht": "PARLIST_HT_ALL",
}


def cast_to_int64_with_missing_values(dg, keys, dataframe_type="pandas"):
    """replace missing values of int64 columns with -1"""
    for c in keys:
        if dataframe_type == "pandas":
            dtype = dg.dtypes[c]
        elif dataframe_type == "vaex":
            dtype = dg.data_type(c)

        if dtype != int64:
            # replace empty strings by -1, e.g. HCN
            # Warning: -1 may be a valid non-equilibirum quantum number for some
            # molecules, e.g. H2O, see https://github.com/radis/radis/issues/280#issuecomment-896120510
            if dataframe_type == "pandas":
                dg[c] = dg[c].replace(r"^\s+$", -1, regex=True)
                dg[c] = dg[c].fillna(-1).astype(int64)  # replace nans with -1
            elif dataframe_type == "vaex":
                dg[c] = dg[c].str.replace(
                    r"^\s+$", "-1", regex=True
                )  # replace empty strings by -1
            else:
                raise NotImplementedError(dataframe_type)


def cast_all_to_int64(df, keys):
    """Convert integer strings to int64. Edited with AI assistance."""
    for key in keys:
        df[key] = df[key].apply(
            lambda x: int(x) if x.strip().lstrip("-").isdigit() else x
        )


def extract_columns(df, extracted_values, columns):
    """extracts column from extracted_values"""
    for column in columns:
        df[column] = extracted_values.apply(lambda x: x.get(column))

    cast_all_to_int64(df, columns)

    for column in columns:
        df[column] = df.evaluate(df[column])


def hit2df(
    fname,
    cache=True,
    verbose=True,
    drop_non_numeric=True,
    load_wavenum_min=None,
    load_wavenum_max=None,
    engine="pytables",
    output="pandas",
    parse_quanta=True,
    cache_directory_path=None,
    fast_parsing=True,
):
    """Convert a HITRAN/HITEMP [1]_ file to a Pandas dataframe

    Parameters
    ----------
    fname: str
        HITRAN-HITEMP file name
    cache: boolean, or ``'regen'`` or ``'force'``
        if ``True``, a pandas-readable HDF5 file is generated on first access,
        and later used. This saves on the datatype cast and conversion and
        improves performances a lot (but changes in the database are not
        taken into account). If False, no database is used. If ``'regen'``, temp
        file are reconstructed. Default ``True``.

    Other Parameters
    ----------------
    drop_non_numeric: boolean
        if ``True``, non numeric columns are dropped. This improves performances,
        but make sure all the columns you need are converted to numeric formats
        before hand. Default ``True``. Note that if a cache file is loaded it
        will be left untouched.
    load_wavenum_min, load_wavenum_max: float
        if not ``'None'``, only load the cached file if it contains data for
        wavenumbers above/below the specified value. See :py:func`~radis.api.cache_files.load_h5_cache_file`.
        Default ``'None'``.
    engine: 'pytables', 'vaex'
        format for Hdf5 cache file. Default `pytables`
    parse_quanta: bool
        if ``True``, parse local & global quanta (required to identify lines
        for non-LTE calculations ; but sometimes lines are not labelled.)
    fast_parsing: bool
        Retained for compatibility. Fast quantum parsing is temporarily
        disabled because it mishandles unassigned quantum numbers.
    output : str
        output format of data as pandas Dataformat or vaex Dataformat
    cache_directory_path : str or None, optional
        Directory to store/read cache files. If None, use the directory of `fname`.

    Returns
    -------
    df: pandas Dataframe or Vaex Dataframe
        dataframe containing all lines and parameters


    References
    ----------

    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__



    Notes
    -----

    Performances: see CDSD-HITEMP parser
    Quantum formats are described in :mod:`radis.api.hitranapi`.
    Edited with AI assistance.


    See Also
    --------

    :func:`~radis.api.cdsdapi.cdsd2df`
    """
    if engine == "pytables" and output == "vaex":
        raise ValueError(
            "Incompatible options: engine='pytables' cannot be used with output='vaex'. "
            "Please choose either:\n"
            "- engine='pytables' with output='pandas', or\n"
            "- engine='vaex' with output='pandas', or\n"
            "- engine='vaex' with output='vaex' (slower)."
        )

    metadata = {}
    # Last modification time of the original file :
    metadata["last_modification"] = time.ctime(getmtime(fname))
    if verbose >= 2:
        print(f"Opening file {fname} (cache={cache})")
        print(f"Last modification time: {metadata['last_modification']}")
    if load_wavenum_min and load_wavenum_max:
        assert load_wavenum_min < load_wavenum_max

    columns = columns_2004

    # Use cache file if possible
    if cache_directory_path:
        base_filename = os.path.basename(fname)
        possible_cache_file = os.path.join(cache_directory_path, base_filename)
        fcache = DataFileManager(engine).cache_file(possible_cache_file)
    else:
        fcache = DataFileManager(engine).cache_file(
            fname
        )  # Use default cache directory

    if cache and exists(fcache):
        relevant_if_metadata_above = (
            {"wavenum_max": load_wavenum_min} if load_wavenum_min else {}
        )  # not relevant if wavenum_max of file is < wavenum min required
        relevant_if_metadata_below = (
            {"wavenum_min": load_wavenum_max} if load_wavenum_max else {}
        )  # not relevant if wavenum_min of file is > wavenum max required
        from radis import __version__, config

        df = load_h5_cache_file(
            fcache,
            cache,
            valid_if_metadata_is=metadata,
            relevant_if_metadata_above=relevant_if_metadata_above,
            relevant_if_metadata_below=relevant_if_metadata_below,
            current_version=__version__,
            last_compatible_version=config["OLDEST_COMPATIBLE_VERSION"],
            verbose=verbose,
            engine=engine,
        )
        if df is not None:
            return df

    #  %% Start reading the full file

    # Detect the molecule by reading the start of the file
    try:
        with open(fname) as f:
            mol = get_molecule(int(f.read(2)))
    except UnicodeDecodeError as err:
        raise ValueError(
            f"You're trying to read a binary file {fname} instead of an HITRAN file"
        ) from err
    df = parse_hitran_file(fname, columns, output=output, molecule=mol)
    df = post_process_hitran_data(
        df,
        molecule=mol,
        dataframe_type=output,
        parse_quanta=parse_quanta,
        fast_parsing=fast_parsing,
    )
    # cached file mode but cached file doesn't exist yet (else we had returned)
    if cache:
        new_metadata = {
            # Last modification time of the original file :
            "last_modification": time.ctime(getmtime(fname)),
            "wavenum_min": df.wav.min(),
            "wavenum_max": df.wav.max(),
        }
        if verbose:
            print(f"Generating cache file {fcache} with metadata :\n{new_metadata}")
        from radis import __version__

        try:
            save_to_hdf(
                df,
                fcache,
                metadata=new_metadata,
                version=__version__,
                overwrite=True,
                verbose=verbose,
                engine=engine,
            )
        except PermissionError:
            if verbose:
                print(sys.exc_info())
                print(
                    "An error occurred in cache file generation. Lookup access rights"
                )
            pass

    # TODO : get only wavenum above/below 'load_wavenum_min', 'load_wavenum_max'
    # by parsing df.wav.   Completely irrelevant files are discarded in 'load_h5_cache_file'
    # but files that have partly relevant lines are fully loaded.
    # Note : cache file is generated with the full line list.

    return df


def post_process_hitran_data(
    df,
    molecule,
    verbose=True,
    drop_non_numeric=True,
    parse_quanta=True,
    add_HITRAN_uncertainty_code=False,
    dataframe_type="pandas",
    fast_parsing=True,
):
    """Parsing non-equilibrium parameters in HITRAN/HITEMP [1]_ file to and return final Pandas Dataframe

    Parameters
    ----------
    df: pandas Dataframe
      dataframe containing generic parameters

    molecule: str
       molecule name

    Other Parameters
    ----------------
    drop_non_numeric: boolean
        if ``True``, non numeric columns are dropped. This improves performances,
        but make sure all the columns you need are converted to numeric formats
        before hand. Default ``True``. Note that if a cache file is loaded it
        will be left untouched.
    parse_quanta: bool
        if ``True``, parse local & global quanta (required to identify lines
        for non-LTE calculations ; but sometimes lines are not labelled.)
    add_HITRAN_uncertainty_code: bool
        if ``True``, a column which contains HITRAN uncertainty code is converted to integer and not dropped.
    dataframe_type: str
        pandas or vaex
    fast_parsing: bool
        Retained for compatibility. Fast quantum parsing is temporarily
        disabled because it mishandles unassigned quantum numbers.

    Returns
    -------
    df: pandas Dataframe or vaex Dataframe
        dataframe containing all lines and parameters


    References
    ----------

    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__



    Notes
    -----

    Performances: see CDSD-HITEMP parser
    Quantum formats are described in :mod:`radis.api.hitranapi`.
    Edited with AI assistance.


    See Also
    --------

    :func:`~radis.io.cdsd.cdsd2df`
    """
    nmol = len(df["id"].unique())
    if nmol == 0:
        raise ValueError("Databank looks empty")
    elif nmol != 1:
        # Crash, give explicitly error messages
        try:
            if dataframe_type == "pandas":
                secondline = df.iloc[1]
            elif dataframe_type == "vaex":
                secondline = df.values[1]
            else:
                raise NotImplementedError(dataframe_type)
        except IndexError:
            secondline = ""
        raise ValueError(
            f"Multiple molecules in database ({nmol} : {[get_molecule(idi) for idi in df['id'].unique()]}). Current "
            + "spectral code only computes 1 species at the time. Use MergeSlabs. "
            + "Verify the parsing was correct by looking at the first row below: "
            + f"\n{df.iloc[0]}"
            + "\n----------------\nand the second row "
            + f"below: \n{secondline}"
        )

    if parse_quanta:
        # Add local quanta attributes, based on the HITRAN group
        try:
            df = parse_local_quanta(
                df,
                molecule,
                verbose=verbose,
                dataframe_type=dataframe_type,
                fast_parsing=fast_parsing,
            )
        except ValueError as err:
            # Empty strings (unlabelled lines) have been reported for HITEMP2010-H2O.
            # In this case, do not parse (makes non-equilibrium calculations impossible).
            # see https://github.com/radis/radis/issues/211
            if verbose:
                print(str(err))
                print("-" * 10)
                print(
                    f"Impossible to parse local quanta in {molecule}, probably an unlabelled line. Ignoring, but nonequilibrium calculations will not be possible. See details above."
                )

        # Add global quanta attributes, based on the HITRAN class
        try:
            df = parse_global_quanta(
                df,
                molecule,
                verbose=verbose,
                dataframe_type=dataframe_type,
                fast_parsing=fast_parsing,
            )
        except ValueError as err:
            # Empty strings (unlabelled lines) have been reported for HITEMP2010-H2O.
            # In this case, do not parse (makes non-equilibrium calculations impossible).
            # see https://github.com/radis/radis/issues/211
            if verbose:
                print(str(err))
                print("-" * 10)
                print(
                    f"Impossible to parse global quanta in {molecule}, probably an unlabelled line. Ignoring, but nonequilibrium calculations will not be possible. See details above."
                )

    # Remove non numerical attributes
    if drop_non_numeric:
        if "branch" in df:
            replace_PQR_with_m101(df)
        ierr = None
        if ("ierr" in df) and add_HITRAN_uncertainty_code:
            df["ierr"] = df["ierr"].astype(int64)
            ierr = df["ierr"]
        if dataframe_type != "vaex":
            df = drop_object_format_columns(df, verbose=verbose)
            if add_HITRAN_uncertainty_code and ierr is not None:
                df["ierr"] = ierr

    return df


# %% Hitran global quanta classes
def _parse_HITRAN_class1_fast_parsing(df, verbose=True, dataframe_type="pandas"):
    r"""Diatomic molecules: CO, HF, HCl, HBr, HI, N2, NO+


    Parameters
    ----------
    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe

    Notes
    -----
    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    Disabled in dispatch because blank assignments are not preserved.

    HITRAN syntax [1]_ :

    >>>       v
    >>>  13x I2

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """
    # Define slice positions for upper and lower vib states
    _GLOBU_SLICES = {
        "vu": (13, 15),
    }
    _GLOBL_SLICES = {
        "vl": (13, 15),
    }

    if dataframe_type == "vaex":
        # Ensure string type
        df["globu"] = df["globu"].astype(str)
        df["globl"] = df["globl"].astype(str)

        # Slice and convert in-place
        for name, (i0, i1) in _GLOBU_SLICES.items():
            df[name] = df["globu"].str.slice(i0, i1).str.strip().astype("int64")
        for name, (i0, i1) in _GLOBL_SLICES.items():
            df[name] = df["globl"].str.slice(i0, i1).str.strip().astype("int64")

        # Drop originals
        df.drop(["globu", "globl"], inplace=True)
        return df

    else:
        # pandas path: cast to str first
        for name, (i0, i1) in _GLOBU_SLICES.items():
            series = (
                df["globu"].astype(str).str.slice(i0, i1).str.strip().replace("", "0")
            )
            df[name] = series.astype("int64")
        for name, (i0, i1) in _GLOBL_SLICES.items():
            series = (
                df["globl"].astype(str).str.slice(i0, i1).str.strip().replace("", "0")
            )
            df[name] = series.astype("int64")

        # Drop originals
        df.drop(columns=["globu", "globl"], inplace=True)
        return df


def _parse_HITRAN_class1(df, verbose=True, dataframe_type="pandas"):
    r"""Diatomic molecules: CO, HF, HCl, HBr, HI, N2, NO+


    Parameters
    ----------
    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe

    Notes
    -----
    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    HITRAN syntax [1]_ :

    >>>       v
    >>>  13x I2

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """

    if dataframe_type == "vaex":
        df = _parse_HITRAN_global_fields(df, [("v", 13, 15, True)], dataframe_type)
        df.drop("globu", inplace=True)
        df.drop("globl", inplace=True)

        return df
    elif dataframe_type == "pandas":
        # 1. Parse
        dgu = (
            df["globu"].astype(str).str.extract(r"[ ]{13}(?P<vu>[\d ]{2})", expand=True)
        )
        dgl = (
            df["globl"].astype(str).str.extract(r"[ ]{13}(?P<vl>[\d ]{2})", expand=True)
        )

        # 2. Convert to numeric
        cast_to_int64_with_missing_values(dgu, ["vu"], dataframe_type=dataframe_type)
        cast_to_int64_with_missing_values(dgl, ["vl"], dataframe_type=dataframe_type)

        # 3. Clean
        del df["globu"]
        del df["globl"]

        return pd.concat([df, dgu, dgl], axis=1)
    else:
        raise NotImplementedError(dataframe_type)


def _parse_HITRAN_class2b(df, verbose=True, dataframe_type="pandas"):
    r"""Parse linear triatomic class in HITRAN [1]_: N2O, OCS, HCN

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe

    Notes
    -----
    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    HITRAN syntax:

    >>>     v1 v2 l2 v3
    >>>  7x I2 I2 I2 I2

    Note: I2 in regexp: [\d ]{2}

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    if dataframe_type == "pandas":
        # 1. Parse
        dgu = (
            df["globu"]
            .astype(str)
            .str.extract(
                r"[ ]{7}(?P<v1u>[-\d ]{2})(?P<v2u>[-\d ]{2})(?P<l2u>[-\d ]{2})(?P<v3u>[-\d ]{2})",
                expand=True,
            )
        )
        dgl = (
            df["globl"]
            .astype(str)
            .str.extract(
                r"[ ]{7}(?P<v1l>[-\d ]{2})(?P<v2l>[-\d ]{2})(?P<l2l>[-\d ]{2})(?P<v3l>[-\d ]{2})",
                expand=True,
            )
        )

        # 2. Convert to numeric
        cast_to_int64_with_missing_values(dgu, ["v1u", "v2u", "l2u", "v3u"])
        cast_to_int64_with_missing_values(dgl, ["v1l", "v2l", "l2l", "v3l"])

        # 3. Clean
        del df["globu"]
        del df["globl"]

        return pd.concat([df, dgu, dgl], axis=1)
    elif dataframe_type == "vaex":
        df = _parse_HITRAN_global_fields(
            df,
            [
                ("v1", 7, 9, True),
                ("v2", 9, 11, True),
                ("l2", 11, 13, True),
                ("v3", 13, 15, True),
            ],
            dataframe_type,
        )
        df.drop("globu", inplace=True)
        df.drop("globl", inplace=True)

        return df
    else:
        raise NotImplementedError(dataframe_type)


def _parse_HITRAN_class2a(df, verbose=True, dataframe_type="pandas"):
    r"""Parse linear triatomic with large Fermi resonance in HITRAN [1]_: CO2

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe

    Notes
    -----
    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    HITRAN syntax:

    >>>     v1 v2 l2 v3 r
    >>>  6x I2 I2 I2 I2 I1

    Note: I2 in regexp: [\d ]{2}

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """

    if dataframe_type == "pandas":
        # 1. Parse
        dgu = (
            df["globu"]
            .astype(str)
            .str.extract(
                r"[ ]{6}(?P<v1u>[-\d ]{2})(?P<v2u>[-\d ]{2})(?P<l2u>[-\d ]{2})(?P<v3u>[-\d ]{2})(?P<ru>[\d ])",
                expand=True,
            )
        )
        dgl = (
            df["globl"]
            .astype(str)
            .str.extract(
                r"[ ]{6}(?P<v1l>[-\d ]{2})(?P<v2l>[-\d ]{2})(?P<l2l>[-\d ]{2})(?P<v3l>[-\d ]{2})(?P<rl>[\d ])",
                expand=True,
            )
        )

        # 2. Convert to numeric
        cast_to_int64_with_missing_values(dgu, ["v1u", "v2u", "l2u", "v3u", "ru"])
        cast_to_int64_with_missing_values(dgl, ["v1l", "v2l", "l2l", "v3l", "rl"])

        # 3. Clean
        del df["globu"]
        del df["globl"]

        return pd.concat([df, dgu, dgl], axis=1)
    elif dataframe_type == "vaex":
        # Define slices for extracting values from globu and globl
        _GLOBU_SLICES = {
            "v1u": (6, 8),
            "v2u": (8, 10),
            "l2u": (10, 12),
            "v3u": (12, 14),
            "ru": (14, 15),
        }
        _GLOBL_SLICES = {
            "v1l": (6, 8),
            "v2l": (8, 10),
            "l2l": (10, 12),
            "v3l": (12, 14),
            "rl": (14, 15),
        }

        for name, (i0, i1) in _GLOBU_SLICES.items():
            df[name] = df["globu"].str.slice(i0, i1).str.strip()

        for name, (i0, i1) in _GLOBL_SLICES.items():
            df[name] = df["globl"].str.slice(i0, i1).str.strip()

        numeric_cols = list(_GLOBU_SLICES.keys()) + list(_GLOBL_SLICES.keys())
        for col in numeric_cols:
            df[col] = df[col].str.replace(r"^\s*$", "-1", regex=True).astype("int64")

        # 3. Clean
        df.drop("globu", inplace=True)
        df.drop("globl", inplace=True)

        return df


def _parse_HITRAN_class3_fast_parsing(df, verbose=True, dataframe_type="pandas"):
    r"""Parse non-linear triatomic in HITRAN [1]_: H2O, O3, SO2, NO2, HOCl, H2S, HO2, HOBr

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe
    Notes
    -----

    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    Disabled in dispatch because blank assignments are not preserved.

    HITRAN syntax:

    >>>     v1 v2 v3
    >>>  9x I2 I2 I2

    Note: I2 in regexp: [\d ]{2}

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    _GLOBU_SLICES = {
        "v1u": (9, 11),
        "v2u": (11, 13),
        "v3u": (13, 15),
    }

    _GLOBL_SLICES = {
        "v1l": (9, 11),
        "v2l": (11, 13),
        "v3l": (13, 15),
    }
    if dataframe_type == "vaex":
        # Vaex string slicing and assignment
        for name, (i0, i1) in _GLOBU_SLICES.items():
            df[name] = df["globu"].str.slice(i0, i1).str.strip().astype("int64")
        for name, (i0, i1) in _GLOBL_SLICES.items():
            df[name] = df["globl"].str.slice(i0, i1).str.strip().astype("int64")
        df.drop("globu", inplace=True)
        df.drop("globl", inplace=True)
        return df
    elif dataframe_type == "pandas":
        for name, (i0, i1) in _GLOBU_SLICES.items():
            series = df["globu"].str.slice(i0, i1).str.strip().replace("", "0")
            df[name] = series.astype("int64")
        for name, (i0, i1) in _GLOBL_SLICES.items():
            series = df["globl"].str.slice(i0, i1).str.strip().replace("", "0")
            df[name] = series.astype("int64")
        df.drop(columns=["globu", "globl"], inplace=True)
        return df
    else:
        raise NotImplementedError(dataframe_type)


def _parse_HITRAN_class3(df, verbose=True, dataframe_type="pandas"):
    r"""Parse non-linear triatomic in HITRAN [1]_: H2O, O3, SO2, NO2, HOCl, H2S, HO2, HOBr

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe
    Notes
    -----
    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    HITRAN syntax:

    >>>     v1 v2 v3
    >>>  9x I2 I2 I2

    Note: I2 in regexp: [\d ]{2}

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    if dataframe_type == "pandas":
        # 1. Parse
        dgu = (
            df["globu"]
            .astype(str)
            .str.extract(
                #        '[ ]{9}(?P<v1u>[\d ]{2})(?P<v2u>[\d ]{2})(?P<v3u>[\d ]{2})',
                r"[ ]{9}(?P<v1u>[\-\d ]{2})(?P<v2u>[\-\d ]{2})(?P<v3u>[\-\d ]{2})",
                expand=True,
            )
        )
        dgl = (
            df["globl"]
            .astype(str)
            .str.extract(
                #        '[ ]{9}(?P<v1l>[\d ]{2})(?P<v2l>[\d ]{2})(?P<v3l>[\d ]{2})',
                r"[ ]{9}(?P<v1l>[\-\d ]{2})(?P<v2l>[\-\d ]{2})(?P<v3l>[\-\d ]{2})",
                expand=True,
            )
        )
        # ... note @EP: in HITRAN H2O files, for iso=2, vibrational levels are
        # ... somehow negative. The regex above is adapted to catch negation signs with \-

        # 2. Convert to numeric
        cast_to_int64_with_missing_values(
            dgu,
            [
                "v1u",
                "v2u",
                "v3u",
            ],
        )
        cast_to_int64_with_missing_values(dgl, ["v1l", "v2l", "v3l"])

        # 3. Clean
        del df["globu"]
        del df["globl"]

        return pd.concat([df, dgu, dgl], axis=1)
    elif dataframe_type == "vaex":
        df = _parse_HITRAN_global_fields(
            df,
            [("v1", 9, 11, True), ("v2", 11, 13, True), ("v3", 13, 15, True)],
            dataframe_type,
        )
        df.drop("globu", inplace=True)
        df.drop("globl", inplace=True)

        return df
    else:
        raise NotImplementedError(dataframe_type)


def _parse_HITRAN_global_fields(df, fields, dataframe_type):
    """Parse fixed-width global fields, retaining the original state labels.

    Each field is (name, start, end, numeric). Ammonia's single-character
    assignments are strings because HITRAN also uses special characters for
    partially assigned states. Edited with AI assistance.
    """
    if dataframe_type not in ("pandas", "vaex"):
        raise NotImplementedError(dataframe_type)
    parsed = {}
    for source, suffix in (("globu", "u"), ("globl", "l")):
        for name, start, end, numeric in fields:
            values = df[source].str.slice(start, end).str.strip()
            if numeric:
                if dataframe_type == "pandas":
                    values = values.replace("", "-1").fillna("-1").astype("int64")
                else:
                    values = values.str.replace(r"^\s*$", "-1", regex=True).astype(
                        "int64"
                    )
            parsed[name + suffix] = values
    for name, values in parsed.items():
        df[name] = values
    return df


def _parse_HITRAN_class1b(df, dataframe_type="pandas"):
    """Parse HITRAN2024 class 1b, updating HITRAN2020. Edited with AI assistance."""
    return _parse_HITRAN_global_fields(
        df,
        [
            ("electronic_state", 6, 8, False),
            ("omega", 8, 11, False),
            ("v", 13, 15, True),
        ],
        dataframe_type,
    )


def _parse_HITRAN_class4(df, mol, dataframe_type="pandas"):
    """Parse HITRAN2024 class 4, updating HITRAN2020. Edited with AI assistance."""
    if mol == "NH3":
        fields = [
            ("v1", 1, 2, False),
            ("v2", 2, 3, False),
            ("v3", 3, 4, False),
            ("v4", 4, 5, False),
            ("l3", 6, 7, False),
            ("l4", 7, 8, False),
            ("l", 9, 10, False),
            ("vib_symmetry", 11, 15, False),
        ]
    else:
        fields = [
            ("v1", 5, 7, True),
            ("v2", 7, 9, True),
            ("v3", 9, 11, True),
            ("v4", 11, 13, True),
            ("vib_symmetry", 13, 15, False),
        ]
    return _parse_HITRAN_global_fields(df, fields, dataframe_type)


def _parse_HITRAN_class5(df, mol, dataframe_type="pandas"):
    """Parse HITRAN2024 class 5, updating HITRAN2020. Edited with AI assistance."""
    if mol == "C2H2":
        fields = [
            ("v1", 1, 2, True),
            ("v2", 2, 3, True),
            ("v3", 3, 4, True),
            ("v4", 4, 6, True),
            ("v5", 6, 8, True),
            ("l4", 8, 10, True),
            ("l5", 10, 12, True),
            ("vib_parity", 12, 13, False),
            ("vib_symmetry", 14, 15, False),
        ]
    elif mol == "C4H2":
        fields = [(f"v{i}", i, i + 1, True) for i in range(1, 10)]
        fields += [("vib_parity", 11, 12, False), ("vib_symmetry", 13, 15, False)]
    elif mol == "HC3N":
        fields = [(f"v{i}", i + 1, i + 2, True) for i in range(1, 8)]
        fields += [("l5", 9, 11, True), ("l6", 11, 13, True), ("l7", 13, 15, True)]
    elif mol == "C2N2":
        fields = [(f"v{i}", 2 * (i - 1), 2 * i, True) for i in range(1, 6)]
        fields += [
            ("l", 10, 12, True),
            ("vib_parity", 12, 13, False),
            ("r", 13, 14, True),
            ("vib_symmetry", 14, 15, False),
        ]
    else:
        raise ValueError(f"Unknown HITRAN class 5 molecule {mol}")
    return _parse_HITRAN_global_fields(df, fields, dataframe_type)


def _parse_HITRAN_class6(df, mol, dataframe_type="pandas"):
    """Parse HITRAN2024 class 6, updating HITRAN2020. Edited with AI assistance."""
    fields = [("v1", 3, 5, True), ("v2", 5, 7, True), ("v3", 7, 9, True)]
    if mol == "H2O2":
        fields += [("n", 9, 10, True), ("tau", 10, 11, True)]
    else:
        fields += [("v4", 9, 11, True)]
    fields += [("v5", 11, 13, True), ("v6", 13, 15, True)]
    return _parse_HITRAN_global_fields(df, fields, dataframe_type)


def _parse_HITRAN_class7(df, mol, dataframe_type="pandas"):
    """Parse HITRAN2024 class 7, updating HITRAN2020. Edited with AI assistance."""
    if mol == "H3+":
        fields = [
            ("v1", 5, 7, True),
            ("v2", 7, 9, True),
            ("l2", 9, 11, True),
            ("total_symmetry", 12, 15, False),
        ]
    else:
        fields = [
            ("v1", 0, 2, True),
            ("v2", 2, 4, True),
            ("v3", 4, 6, True),
            ("l3", 6, 8, True),
            ("v4", 8, 10, True),
            ("l4", 10, 12, True),
            ("vib_symmetry", 12, 15, False),
        ]
    return _parse_HITRAN_global_fields(df, fields, dataframe_type)


def _parse_HITRAN_class8(df, dataframe_type="pandas"):
    """Parse HITRAN2024 class 8, updating HITRAN2020. Edited with AI assistance."""
    return _parse_HITRAN_global_fields(
        df,
        [
            ("v1", 3, 5, True),
            ("v2", 5, 7, True),
            ("v3", 7, 9, True),
            ("v4", 9, 11, True),
            ("multiplicity", 11, 13, False),
            ("vib_symmetry", 13, 15, False),
        ],
        dataframe_type,
    )


def _parse_HITRAN_class9(df, dataframe_type="pandas"):
    """Parse HITRAN2024 class 9, updating HITRAN2020. Edited with AI assistance."""
    return _parse_HITRAN_global_fields(df, [("vib_band", 0, 15, False)], dataframe_type)


# %% HITRAN Local quanta


def _parse_HITRAN_group1_fast_parsing(df, verbose=True, dataframe_type="pandas"):
    r"""Parse asymmetric rotors (:py:attr:`~radis.db.classes.HITRAN_GROUP1` ):
    H2O, O3, SO2, NO2, HNO3, H2CO, HOCl, H2O2, COF2, H2S, HO2, HCOOH, ClONO2, HOBr, C2H4

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe


    Notes
    -----

    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    Disabled in dispatch because blank assignments are not preserved.

    HITRAN syntax: [1]_

    References
    ----------

    .. [1] `Table 4 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """
    _LOCU_SLICES = {
        "ju": (0, 3),
        "Kau": (3, 6),
        "Kcu": (6, 9),
        "Fu": (9, 14),
        "symu": (14, 15),
    }
    _LOCL_SLICES = {
        "jl": (0, 3),
        "Kal": (3, 6),
        "Kcl": (6, 9),
        "Fl": (9, 14),
        "syml": (14, 15),
    }
    if dataframe_type == "vaex":
        # Use vaex string slicing and assignment
        for name, (i0, i1) in _LOCU_SLICES.items():

            if name in ("ju", "Kau", "Kcu"):
                df[name] = df["locu"].str.slice(i0, i1).str.strip().astype("int64")
            else:
                df[name] = df["locu"].str.slice(i0, i1).str.strip()
        for name, (i0, i1) in _LOCL_SLICES.items():
            # series = df["locl"].str.slice(i0, i1).str.strip().str.replace("", "0")
            if name in ("jl", "Kal", "Kcl"):
                df[name] = df["locl"].str.slice(i0, i1).str.strip().astype("int64")
            else:
                df[name] = df["locl"].str.slice(i0, i1).str.strip()
        df.drop("locu", inplace=True)
        df.drop("locl", inplace=True)
        return df
    elif dataframe_type == "pandas":
        # str.slice + astype in one go
        for name, (i0, i1) in _LOCU_SLICES.items():
            series = df["locu"].str.slice(i0, i1).str.strip().replace("", "0")
            if name in ("ju", "Kau", "Kcu"):
                df[name] = series.astype("int64")
            else:
                df[name] = series

        for name, (i0, i1) in _LOCL_SLICES.items():
            series = df["locl"].str.slice(i0, i1).str.strip().replace("", "0")
            if name in ("jl", "Kal", "Kcl"):
                df[name] = series.astype("int64")
            else:
                df[name] = series

        df.drop(columns=["locu", "locl"], inplace=True)
        return df
    else:
        raise NotImplementedError(dataframe_type)


def _parse_HITRAN_group1(df, verbose=True, dataframe_type="pandas"):
    r"""Parse asymmetric rotors (:py:attr:`~radis.db.classes.HITRAN_GROUP1` ):
    H2O, O3, SO2, NO2, HNO3, H2CO, HOCl, H2O2, COF2, H2S, HO2, HCOOH, ClONO2, HOBr, C2H4

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe


    Notes
    -----
    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    HITRAN syntax: [1]_


    References
    ----------

    .. [1] `Table 4 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """

    if dataframe_type == "vaex":
        df = _parse_HITRAN_local_states(
            df,
            [
                ("j", 0, 3, "int64"),
                ("Ka", 3, 6, "int64"),
                ("Kc", 6, 9, "int64"),
                ("F", 9, 14, "str"),
                ("sym", 14, 15, "str"),
            ],
            dataframe_type,
        )
        df.drop("locu", inplace=True)
        df.drop("locl", inplace=True)

        return df
    elif dataframe_type == "pandas":
        # 1. Parse

        # Ref [1] : locu
        # --------------
        # J'  | Ka' | Kc' | F'  | Sym'
        # I3  | I3  | I3  | A5  | A1

        dgu = (
            df["locu"]
            .astype(str)
            .str.extract(
                r"(?P<ju>[\d ]{3})(?P<Kau>[\-\d ]{3})(?P<Kcu>[\-\d ]{3})(?P<Fu>.{5})(?P<symu>.)",
                expand=True,
            )
        )
        # Ref [1] : locl
        # --------------
        # J'' | Ka''| Kc''| F'' | Sym''
        # I3  | I3  | I3  | A5  | A1
        dgl = (
            df["locl"]
            .astype(str)
            .str.extract(
                r"(?P<jl>[\d ]{3})(?P<Kal>[\-\d ]{3})(?P<Kcl>[\-\d ]{3})(?P<Fl>.{5})(?P<syml>.)",
                expand=True,
            )
        )
        # ... note @EP: in HITRAN H2O files, for iso=2, the Kau, Kcu can somehow
        # ... be negative. The regex above is adapted to catch negation signs with \-
        # 2. Convert to numeric
        cast_to_int64_with_missing_values(dgu, ["ju", "Kau", "Kcu"])
        cast_to_int64_with_missing_values(dgl, ["jl", "Kal", "Kcl"])

        # 3. Clean
        del df["locu"]
        del df["locl"]

        return pd.concat([df, dgu, dgl], axis=1)
    else:
        raise NotImplementedError(dataframe_type)


def _parse_HITRAN_group2_fast_parsing(df, verbose=True, dataframe_type="pandas"):
    r"""Parse diatomic and linear molecules (:py:attr:`~radis.db.classes.HITRAN_GROUP2` ):
    CO2, N2O, CO, HF, HCl, HBr, HI, OCS, N2, HCN, C2H2, NO+

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe


    Notes
    -----

    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    Disabled in dispatch because blank assignments are not preserved.

    HITRAN syntax: [1]


    References
    ----------

    .. [1] `Table 4 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    # Define slicing positions
    _LOCU_SLICES = {
        "magnetic": (0, 1),
        "Fu": (10, 15),
    }
    _LOCL_SLICES = {
        "branch": (5, 6),
        "jl": (6, 9),
        "syml": (9, 10),
        "Fl": (10, 15),
    }

    if dataframe_type == "vaex":
        # Use vaex string slicing and assignment
        for name, (i0, i1) in _LOCU_SLICES.items():
            df[name] = df["locu"].str.slice(i0, i1).str.strip()
        for name, (i0, i1) in _LOCL_SLICES.items():
            series = df["locl"].str.slice(i0, i1).str.strip()
            if name == "jl":
                df[name] = series.astype("int64")
            else:
                df[name] = series
        df.drop("locu", inplace=True)
        df.drop("locl", inplace=True)
        return df
    elif dataframe_type == "pandas":
        # pandas: str.slice + replace + astype
        for name, (i0, i1) in _LOCU_SLICES.items():
            series = df["locu"].str.slice(i0, i1).str.strip()
            if name != "magnetic":
                series = series.replace("", "0")
            df[name] = series

        for name, (i0, i1) in _LOCL_SLICES.items():
            series = df["locl"].str.slice(i0, i1).str.strip().replace("", "0")
            if name == "jl":
                df[name] = series.astype("int64")
            else:
                df[name] = series

        df.drop(columns=["locu", "locl"], inplace=True)
        return df
    else:
        raise NotImplementedError(dataframe_type)


def _parse_HITRAN_group2(df, verbose=True, dataframe_type="pandas"):
    r"""Parse diatomic and linear molecules (:py:attr:`~radis.db.classes.HITRAN_GROUP2` ):
    CO2, N2O, CO, HF, HCl, HBr, HI, OCS, N2, HCN, C2H2, NO+

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database
    dataframe_type : str
        pandas or vaex

    Returns
    -------
        pandas Dataframe or Vaex Dataframe


    Notes
    -----
    HITRAN2020 format (Gordon et al., 2022), with HITRAN2024 updates; see
    :mod:`radis.api.hitranapi`. Edited with AI assistance.

    HITRAN syntax: [1]


    References
    ----------

    .. [1] `Table 4 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """

    # 1. Parse

    # Ref [1] : locu
    # --------------
    #     | F'  |
    # 10X | A5  |
    if dataframe_type == "pandas":
        dgu = (
            df["locu"]
            .astype(str)
            .str.extract(r"(?P<magnetic>.)[ ]{9}(?P<Fu>.{5})", expand=True)
        )

        # Ref [1] : locl
        # --------------
        #     | Br  | J'' | Sym''| F'' |
        # 5X  | A1  | I3  | A1   | A5  |
        dgl = (
            df["locl"]
            .astype(str)
            .str.extract(
                r"[ ]{5}(?P<branch>[\S ])(?P<jl>[\d ]{3})(?P<syml>.)(?P<Fl>.{5})",
                expand=True,
            )
        )

        # 2. Convert to numeric

        # dgl['jl'] = dgl.jl.apply(pd.to_numeric)
        cast_to_int64_with_missing_values(dgl, ["jl"], dataframe_type=dataframe_type)
        dgl["branch"] = dgl["branch"].str.strip()

        # 3. Clean
        del df["locu"]
        del df["locl"]

        return pd.concat([df, dgu, dgl], axis=1)
    elif dataframe_type == "vaex":
        # 1.Parse
        _LOCU_SLICES = {
            "magnetic": (0, 1),
            "Fu": (10, 15),
        }

        _LOCL_SLICES = {
            "branch": (5, 6),
            "jl": (6, 9),
            "syml": (9, 10),
            "Fl": (10, 15),
        }

        for name, (i0, i1) in _LOCU_SLICES.items():
            df[name] = df["locu"].str.slice(i0, i1).str.strip()

        for name, (i0, i1) in _LOCL_SLICES.items():
            df[name] = df["locl"].str.slice(i0, i1).str.strip()

        df["jl"] = df["jl"].str.replace(r"^\s*$", "-1", regex=True).astype("int64")

        # 3. Clean
        df.drop("locu", inplace=True)
        df.drop("locl", inplace=True)

        return df
    else:
        raise NotImplementedError(dataframe_type)


#                # 10X included in Fu
#               'Fu',     ('a15',  str,  'upper state total angular momentum including nuclear spin'  ,''         )),(
#               #'locl',   ('a15',  str,   'electronic and vibrational local lower quanta'  ,''                      )),(
#               'branch', ('a6',   str,     'O, P, Q, R, S branch symbol'                  ,''                      )),(
#               'jl',     ('a3',   int,    'lower state rotational quantum number'         ,''                      )),(
#               'sym',    ('a1',   str,     'symmetry'                                     ,''                      )),(
#               'Fl',     ('a5',   str,   'lower state total angular momentum including nuclear spin', ''         )),(


def _parse_HITRAN_local_fields(df, fields, dataframe_type):
    """Parse fixed-width local fields. Edited with AI assistance."""
    if dataframe_type not in ("pandas", "vaex"):
        raise NotImplementedError(dataframe_type)
    parsed = {}
    for source, name, start, end, dtype in fields:
        values = df[source].str.slice(start, end).str.strip()
        if dtype != "str":
            if dataframe_type == "pandas":
                values = values.replace("", "-1").fillna("-1").astype(dtype)
            else:
                values = values.str.replace(r"^\s*$", "-1", regex=True).astype(dtype)
        parsed[name] = values
    for name, values in parsed.items():
        df[name] = values
    return df


def _parse_HITRAN_local_states(df, fields, dataframe_type):
    """Parse both local states. Edited with AI assistance."""
    return _parse_HITRAN_local_fields(
        df,
        [
            (source, name + suffix, start, end, dtype)
            for source, suffix in (("locu", "u"), ("locl", "l"))
            for name, start, end, dtype in fields
        ],
        dataframe_type,
    )


def _parse_HITRAN_group2b(df, dataframe_type="pandas"):
    """Parse HITRAN2024 group 2b, updating HITRAN2020. Edited with AI assistance."""
    fields = [(f"l{i}", 2 * (i - 6), 2 * (i - 5), "str") for i in range(6, 10)]
    df = _parse_HITRAN_local_states(df, fields, dataframe_type)
    return _parse_HITRAN_local_fields(
        df,
        [
            ("locl", "branch", 9, 10, "str"),
            ("locl", "jl", 10, 13, "int64"),
            ("locl", "syml", 13, 14, "str"),
        ],
        dataframe_type,
    )


def _parse_HITRAN_group3(df, dataframe_type="pandas"):
    """Parse HITRAN2024 group 3, updating HITRAN2020. Edited with AI assistance."""
    return _parse_HITRAN_local_states(
        df,
        [
            ("j", 2, 5, "int64"),
            ("rot_symmetry", 5, 7, "str"),
            ("alpha", 7, 10, "int64"),
            ("F", 10, 15, "str"),
        ],
        dataframe_type,
    )


def _parse_HITRAN_group4(df, mol, dataframe_type="pandas"):
    """Parse HITRAN2024 group 4, updating HITRAN2020. Edited with AI assistance."""
    if mol == "NH3":
        fields = [
            ("j", 0, 2, "int64"),
            ("K", 2, 5, "int64"),
            ("rot_l", 5, 7, "int64"),
            ("rot_symmetry", 8, 11, "str"),
            ("total_symmetry", 11, 14, "str"),
        ]
    elif mol == "C2H6":
        fields = [
            ("j", 0, 3, "int64"),
            ("K", 3, 6, "int64"),
            ("rot_l", 6, 8, "int64"),
            ("sym", 8, 11, "str"),
            ("F", 11, 15, "str"),
        ]
    else:
        fields = [
            ("j", 0, 3, "int64"),
            ("K", 3, 6, "int64"),
            ("rot_l", 6, 8, "int64"),
            ("rot_symmetry", 8, 10, "str"),
            ("sym", 10, 11, "str"),
            ("F", 11, 15, "str"),
        ]
    return _parse_HITRAN_local_states(df, fields, dataframe_type)


def _parse_HITRAN_group5(df, mol, dataframe_type="pandas"):
    """Parse HITRAN2024 group 5, updating HITRAN2020. Edited with AI assistance."""
    if mol == "H3+":
        fields = [
            ("j", 0, 3, "int64"),
            ("sym", 3, 5, "str"),
            ("G", 5, 8, "int64"),
            ("U", 8, 10, "str"),
            ("K", 10, 13, "str"),
        ]
    else:
        fields = [
            ("j", 0, 3, "int64"),
            ("K", 3, 6, "int64"),
            ("rot_symmetry", 8, 11, "str"),
            ("total_symmetry", 11, 14, "str"),
        ]
    return _parse_HITRAN_local_states(df, fields, dataframe_type)


def _parse_HITRAN_group6(df, dataframe_type="pandas"):
    """Parse HITRAN2024 group 6, updating HITRAN2020. Edited with AI assistance."""
    return _parse_HITRAN_local_fields(
        df,
        [
            ("locu", "Fu", 10, 15, "str"),
            ("locl", "branch_N", 1, 2, "str"),
            ("locl", "Nl", 2, 5, "int64"),
            ("locl", "branch_J", 5, 6, "str"),
            ("locl", "jl", 6, 9, "int64"),
            ("locl", "Fl", 9, 14, "str"),
            ("locl", "transition_type", 14, 15, "str"),
        ],
        dataframe_type,
    )


def _parse_HITRAN_group7(df, mol, dataframe_type="pandas"):
    """Parse HITRAN2024 group 7, updating HITRAN2020. Edited with AI assistance."""
    if mol == "OH":
        fields = [
            ("locl", "branch_NJ", 1, 3, "str"),
            ("locl", "jl", 3, 8, "float64"),
            ("locl", "symmetry_pair", 8, 10, "str"),
        ]
    else:
        fields = [
            ("locu", "magnetic", 0, 1, "str"),
            ("locl", "branch_NJ", 2, 4, "str"),
            ("locl", "jl", 4, 9, "float64"),
            ("locl", "syml", 9, 10, "str"),
        ]
    fields += [("locu", "Fu", 10, 15, "str"), ("locl", "Fl", 10, 15, "str")]
    return _parse_HITRAN_local_fields(df, fields, dataframe_type)


# %% Reading function


def parse_local_quanta(
    df, mol, verbose=True, dataframe_type="pandas", fast_parsing=True
):
    r"""
    Parameters
    ----------

    df: pandas Dataframe

    mol: str
        molecule name

    Notes
    -----
    Uses the HITRAN2024 Table S2 local-state formats. Newly supported groups
    retain ``locu`` and ``locl``; blank numeric assignments are -1. Local
    angular momenta are named ``rot_l`` to distinguish them from global
    vibrational angular momenta. Open-shell branches are kept as separate
    strings rather than interpreted as closed-shell P/Q/R branches. String
    labels are removed by pandas post-processing when ``drop_non_numeric`` is
    true; parsing these fields does not itself enable non-LTE calculations.
    These formats update HITRAN2020 (Gordon et al., 2022); see
    :mod:`radis.api.hitranapi`. Edited with AI assistance. ``fast_parsing`` is
    retained for compatibility but currently has no effect.
    """

    # Fast paths are disabled until they preserve unassigned quantum numbers.
    if mol in HITRAN_GROUP1:
        df = _parse_HITRAN_group1(df, verbose=verbose, dataframe_type=dataframe_type)

    elif mol == "C4H2":
        df = _parse_HITRAN_group2b(df, dataframe_type=dataframe_type)
    elif mol in HITRAN_GROUP2:
        df = _parse_HITRAN_group2(df, verbose=verbose, dataframe_type=dataframe_type)
    elif mol in HITRAN_GROUP3:
        df = _parse_HITRAN_group3(df, dataframe_type=dataframe_type)
    elif mol in HITRAN_GROUP4:
        df = _parse_HITRAN_group4(df, mol, dataframe_type=dataframe_type)
    elif mol in HITRAN_GROUP5:
        df = _parse_HITRAN_group5(df, mol, dataframe_type=dataframe_type)
    elif mol in HITRAN_GROUP6:
        df = _parse_HITRAN_group6(df, dataframe_type=dataframe_type)
    elif mol in HITRAN_GROUP7:
        df = _parse_HITRAN_group7(df, mol, dataframe_type=dataframe_type)
    else:
        raise ValueError(f"Unknown group for molecule {mol}. Cant parse local quanta")

    return df


def parse_global_quanta(
    df, mol, verbose=True, dataframe_type="pandas", fast_parsing=True
):
    r"""

    Parameters
    ----------

    df: pandas Dataframe

    mol: str
        molecule name

    Notes
    -----
    Uses the nine-class global-state formats in Table S1 of the HITRAN2024
    supplementary material. Previously unsupported formats retain ``globu``
    and ``globl`` alongside the parsed columns. In those formats, blank
    numeric assignments are represented by -1; NH3 assignments are strings
    to preserve partial labels. String-valued fields are dropped when
    ``post_process_hitran_data`` is called with ``drop_non_numeric=True``.
    These formats update HITRAN2020 (Gordon et al., 2022); see
    :mod:`radis.api.hitranapi`. Edited with AI assistance. ``fast_parsing`` is
    retained for compatibility but currently has no effect.
    """
    if mol in HITRAN_CLASS1:
        if mol in ("O2", "NO", "OH", "ClO", "SO", "S2"):
            df = _parse_HITRAN_class1b(df, dataframe_type=dataframe_type)
        else:
            df = _parse_HITRAN_class1(
                df, verbose=verbose, dataframe_type=dataframe_type
            )
    elif mol in HITRAN_CLASS2:
        if mol == "CO2":
            df = _parse_HITRAN_class2a(
                df, verbose=verbose, dataframe_type=dataframe_type
            )
        else:
            df = _parse_HITRAN_class2b(
                df, verbose=verbose, dataframe_type=dataframe_type
            )
    elif mol in HITRAN_CLASS3:
        df = _parse_HITRAN_class3(df, verbose=verbose, dataframe_type=dataframe_type)
    elif mol in HITRAN_CLASS4:
        df = _parse_HITRAN_class4(df, mol, dataframe_type=dataframe_type)
    elif mol in HITRAN_CLASS5:
        df = _parse_HITRAN_class5(df, mol, dataframe_type=dataframe_type)
    elif mol in HITRAN_CLASS6:
        df = _parse_HITRAN_class6(df, mol, dataframe_type=dataframe_type)
    elif mol in HITRAN_CLASS7:
        df = _parse_HITRAN_class7(df, mol, dataframe_type=dataframe_type)
    elif mol in HITRAN_CLASS8:
        df = _parse_HITRAN_class8(df, dataframe_type=dataframe_type)
    elif mol in HITRAN_CLASS9:
        df = _parse_HITRAN_class9(df, dataframe_type=dataframe_type)
    else:
        raise ValueError(f"Unknown class for molecule {mol}. Cant parse global quanta")

    return df


# %%


class HITRANDatabaseManager(DatabaseManager):
    def __init__(
        self,
        name,
        molecule,
        local_databases,
        engine="default",
        extra_params=None,
        verbose=True,
        parallel=True,
    ):
        """
        See Also
        --------
        HITRANDatabaseManager is compatible with Exojax :py:class:`exojax.spec.api.Mdbhitran`
        """
        super().__init__(
            name,
            molecule,
            local_databases,
            engine=engine,
            extra_params=extra_params,
            verbose=verbose,
            parallel=parallel,
        )
        self.downloadable = True
        self.base_url = None
        self.Nlines = None
        self.wmin = None
        self.wmax = None

        self.actual_file = None

    def get_filenames(self):
        if self.engine == "vaex":
            return [join(self.local_databases, f"{self.molecule}.hdf5")]
        elif self.engine == "pytables":
            return [join(self.local_databases, f"{self.molecule}.h5")]
        else:
            raise NotImplementedError()

    def download_and_parse(
        self,
        local_file,
        cache=True,
        parse_quanta=True,
        add_HITRAN_uncertainty_code=False,
    ):
        """Download from HITRAN and parse into ``local_file``.
        Also add metadata

        Overwrites :py:meth:`radis.api.dbmanager.DatabaseManager.download_and_parse`
        which downloads from a list of URL, because here we use [HAPI]_ to
        download the files.

        Parameters
        ----------
        opener: an opener with an .open() command
        gfile : file handler. Filename: for info"""

        from hapi import LOCAL_TABLE_CACHE, db_begin, fetch

        from radis.db.classes import get_molecule_identifier

        if isinstance(local_file, list):
            assert (
                len(local_file) == 1
            )  # fetch_hitran stores all lines of a given molecule in one file
            local_file = local_file[0]

        wmin = 0
        wmax = 1e10

        def download_all_hitran_isotopes(molecule, directory, extra_params):
            """Download all isotopologues of the given molecule, one request each.

            .. warning::
                this won't download higher isotopologues (ex : isotopologue
                10-11-12 for CO2). Neglected for the moment, they're irrelevant
                for most calculations anyway
            """
            directory = abspath(expanduser(directory))
            max_fetch_retries = 3
            # create temp folder :
            from radis.misc.basics import make_folders

            make_folders(*split(directory))

            # Use tqdm for isotope progress
            import io
            import sys

            from tqdm import tqdm

            # Suppress HAPI's verbose output (version banner, directory messages, etc.)
            old_stdout = sys.stdout
            sys.stdout = io.StringIO()
            try:
                db_begin(directory)
            finally:
                sys.stdout = old_stdout

            isotope_list = []
            data_file_list = []
            header_file_list = []
            uncertainty_params = ["ierr", "iref"] if add_HITRAN_uncertainty_code else []
            if extra_params == "all":
                fetch_kwargs = dict(
                    ParameterGroups=[*PARAMETER_GROUPS_HITRAN],
                    Parameters=uncertainty_params,
                )
            elif extra_params is None:
                fetch_kwargs = dict(Parameters=uncertainty_params)
            else:
                raise ValueError("extra_params can only be 'all' or None ")

            def is_complete(data_file):
                """A complete HITRAN response ends on a record boundary, in both
                the fixed-width ``.par`` and the comma-separated formats."""
                with open(data_file, "rb") as f:
                    f.seek(-1, os.SEEK_END)
                    return f.read(1) == b"\n"

            for iso in tqdm(
                get_hitran_isotopologues(molecule),
                desc="Downloading isotopologues",
                disable=not self.verbose,
            ):
                file = f"{molecule}_{iso}"
                data_file = join(directory, file + ".data")
                if exists(data_file):
                    if cache != "regen":  # else remove without printing message
                        from radis.misc.printer import printr

                        printr(f"File already exist: {data_file}. Deleting it.`")
                    os.remove(data_file)

                retry_delay = 1
                truncated = False
                for attempt in range(max_fetch_retries):
                    # Suppress HAPI's verbose output (65536 bytes written...)
                    old_stdout = sys.stdout
                    sys.stdout = io.StringIO()
                    try:
                        fetch(
                            file,
                            get_molecule_identifier(molecule),
                            iso,
                            wmin,
                            wmax,
                            **fetch_kwargs,
                        )
                    finally:
                        sys.stdout = old_stdout

                    ### We test if the download went well ###
                    df = pd.DataFrame(LOCAL_TABLE_CACHE[file]["data"])
                    truncated = False
                    if len(df["molec_id"]) == 0:
                        reason = "The database looks empty."
                    elif not is_complete(data_file):
                        # HAPI doesn't check that the transfer completed: a
                        # truncated file still parses, silently giving an
                        # incomplete line list
                        truncated = True
                        reason = "The response was cut short."
                    else:
                        break  # everything worked well
                    warning_msg = (
                        f"HITRAN fetch failed. {reason} "
                        f"Waiting {retry_delay} seconds and trying again (attempt {attempt + 1}/{max_fetch_retries})."
                    )
                    if attempt == max_fetch_retries - 1:
                        warning_msg += "\nLAST ATTEMPT"
                    warnings.warn(warning_msg, UserWarning, stacklevel=2)
                    time.sleep(retry_delay)
                    retry_delay *= 2  # exponential

                if truncated:  # never return a partial line list silently
                    raise ConnectionError(
                        f"Download of {file} from HITRAN was truncated every time "
                        f"({max_fetch_retries} attempts). The line list would be incomplete."
                    )

                isotope_list.append(iso)
                data_file_list.append(file + ".data")
                header_file_list.append(file + ".header")
            return isotope_list, data_file_list, header_file_list

        molecule = self.molecule
        wmin_final = 100000
        wmax_final = -1

        # create database in a subfolder to isolate molecules from one-another
        # (HAPI doesn't check and may mix molecules --> see failure at https://app.travis-ci.com/github/radis/radis/jobs/548126303#L2676)
        tempdir = join(self.tempdir, molecule)
        extra_params = self.extra_params

        # Use HAPI only to download the files, then we'll parse them with RADIS's
        # parsers, and convert to RADIS's fast HDF5 file formats.

        max_retries = 2
        retry_count = 0

        while retry_count < max_retries:
            try:
                (
                    isotope_list,
                    data_file_list,
                    header_file_list,
                ) = download_all_hitran_isotopes(molecule, tempdir, extra_params)
                break  # Success, exit retry loop
            except ValueError as e:
                if "invalid literal for int() with base 10: '\\x00\\x00'" in str(e):
                    retry_count += 1
                    if retry_count < max_retries:
                        warnings.warn(
                            f"ValueError encountered during download. Retrying (attempt {retry_count}/{max_retries})...",
                            UserWarning,
                            stacklevel=2,
                        )
                    else:
                        raise  # Re-raise on final attempt
                else:
                    raise  # Re-raise if it's a different ValueError

        writer = self.get_datafile_manager()

        # Create HDF5 cache file for all isotopes
        Nlines = 0
        for iso, data_file in zip(isotope_list, data_file_list):
            df = pd.DataFrame(LOCAL_TABLE_CACHE[data_file.split(".")[0]]["data"])
            if add_HITRAN_uncertainty_code and "ierr" not in df and "par_line" in df:
                df["ierr"] = df["par_line"].str.slice(127, 133).str.strip()
                df["iref"] = df["par_line"].str.slice(133, 145).str.strip()
            df.rename(
                columns={
                    "molec_id": "id",
                    "local_iso_id": "iso",
                    "nu": "wav",
                    "sw": "int",
                    "a": "A",
                    "gamma_air": "airbrd",
                    "gamma_self": "selbrd",
                    "elower": "El",
                    "n_air": "Tdpair",
                    "delta_air": "Pshft",
                    "global_upper_quanta": "globu",
                    "global_lower_quanta": "globl",
                    "local_upper_quanta": "locu",
                    "local_lower_quanta": "locl",
                    "gp": "gp",
                    "gpp": "gpp",
                },
                inplace=True,
            )
            df = post_process_hitran_data(
                df,
                molecule=molecule,
                parse_quanta=parse_quanta,
                add_HITRAN_uncertainty_code=add_HITRAN_uncertainty_code,
            )
            if add_HITRAN_uncertainty_code and "ierr" not in df.columns:
                raise KeyError(
                    "HITRAN uncertainty column 'ierr' was requested but not loaded"
                )

            wmin_final = min(wmin_final, df.wav.min())
            wmax_final = max(wmax_final, df.wav.max())
            Nlines += len(df)

            writer.write(
                local_file, df, append=True
            )  # create temporary files if required

        # Open all HDF5 cache files and export in a single file with Vaex
        writer.combine_temp_batch_files(
            local_file, sort_values="wav"
        )  # used for vaex mode only
        # Note: by construction, in Pytables mode the database is not sorted
        # by 'wav' but by isotope

        self.wmin = wmin_final
        self.wmax = wmax_final

        # Add metadata
        from radis import __version__

        writer.add_metadata(
            local_file,
            {
                "wavenumber_min": self.wmin,
                "wavenumber_max": self.wmax,
                "download_date": self.get_today(),
                "download_url": "downloaded by HAPI, parsed & store with RADIS",
                "total_lines": Nlines,
                "version": __version__,
            },
        )

        # # clean downloaded files  TODO
        # for file in data_file_list + header_file_list:
        #     os.remove(join(self.local_databases, "downloads", file))

    def register(self, download_files):
        """register in ~/radis.json"""

        if self.is_registered():
            dict_entries = getDatabankEntries(
                self.name
            )  # just update previous register details
        else:

            dict_entries = {}

        if download_files:
            files = [self.actual_file]
            if self.wmin is None or self.wmax is None:
                print(
                    "Somehow wmin and wmax was not given for this database. Reading from the files"
                )
                ##  fix:
                # (can happen if database was downloaded & parsed, but registration failed a first time)
                df_full = self.load(
                    local_files=files,
                    columns=["wav"],
                    within=[],
                    lower_bound=[],
                    upper_bound=[],
                )
                self.wmin = df_full.wav.min()
                self.wmax = df_full.wav.max()
                print(
                    f"Somehow wmin and wmax was not given for this database. Read {self.wmin}, {self.wmax} directly from the files"
                )

        info = f"HITRAN {self.molecule} lines ({self.wmin:.1f}-{self.wmax:.1f} cm-1) with TIPS-2021 (through HAPI) for partition functions"
        dict_entries.update(
            {
                "path": files,
                "format": "hdf5-radisdb",
                "info": info,
                "wavenumber_min": self.wmin,
                "wavenumber_max": self.wmax,
                "download_date": self.get_today(),
            }
        )

        # Add energy level calculation
        from radis.db import MOLECULES_LIST_NONEQUILIBRIUM

        if self.molecule in MOLECULES_LIST_NONEQUILIBRIUM:
            dict_entries[
                "info"
            ] += " and RADIS spectroscopic constants for rovibrational energies (nonequilibrium)"
            dict_entries["levelsfmt"] = "radis"

        try:
            super().register(dict_entries)
        except DatabaseAlreadyExists as e:
            raise Exception(
                'If you want RADIS to overwrite the existing entry for a registered databank, set the config option "ALLOW_OVERWRITE" to True.'
            ) from e


# %%
def hitranxsc(hitranXSC):
    """Parse Hitran Cross-section files manually downloaded from https://hitran.org/xsc/
    Returns a dictionary

    Example
    -------
    Read and plot a local acetone cross-section file:

        datafile = 'CH3COCH3_233.4_375.2_700.0-1780.0_13.xsc'

        data = hitran_crosssection(datafile)

        # %% Plot cross section:
        import matplotlib.pyplot as plt
        plt.plot(data['wavenumber'], data['spectrum'])
        plt.xlabel('Wavenumber (cm-1)')
        plt.ylabel('Cross section (cm2/molecule)')
        # add title with molecule name, pressure and temperature:
        plt.title(data['name'] + ' at ' + str(data['P']) + ' Torr and ' + str(data['T']) + ' K')
        plt.show()


    Notes
    -----

    Code adapted from https://fr.mathworks.com/matlabcentral/fileexchange/74716-load-hitran-absorption-cross-section
    by Mohammadamir Ghaderi

    Under licence conditions ::

        Copyright (c) 2020, Mohammadamir Ghaderi
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

        * Redistributions of source code must retain the above copyright notice, this
        list of conditions and the following disclaimer.

        * Redistributions in binary form must reproduce the above copyright notice,
        this list of conditions and the following disclaimer in the documentation
        and/or other materials provided with the distribution
        * Neither the name of  nor the names of its
        contributors may be used to endorse or promote products derived from this
        software without specific prior written permission.
        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
        AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
        IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
        FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
        DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
        SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
        CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
        OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
        OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



    """
    import numpy as np

    with open(hitranXSC, "r") as file:
        line = file.readline()
        info = line.split()
        data = {}
        data["molec"] = info[0]  # molecule
        min_nu = float(info[1])  # minimum wavenumber (cm-1)
        max_nu = float(info[2])  # maximum wavenumber (cm-1)
        N = int(info[3])  # number of datapoints
        data["T"] = float(info[4])  # Temperature (K)
        data["P"] = float(info[5])  # Pressure (Torr)
        data["InstrumentResolution"] = float(info[7])
        data["name"] = info[8]  # Familiar name of molecule
        if len(info) > 10:
            data["Broadener"] = info[10]  # Broadening species
        data["wavenumber"] = np.linspace(min_nu, max_nu, N)  # Wavenumber range (cm-1)

    # Read all lines of file (except from header) and flatten into an array; knowing that all rows aren't the same length
    with open(hitranXSC, "r") as file:
        lines = file.readlines()[1:]
        lines = [line.split() for line in lines]
        lines = [item for sublist in lines for item in sublist]
        data["spectrum"] = np.array(lines, dtype=np.float64)
        # remove zeros from end of file
        while data["spectrum"][-1] == 0:
            data["spectrum"] = data["spectrum"][:-1]
    return data


# ======================================================
# %% Test


if __name__ == "__main__":

    from radis.test.io.test_hitran_cdsd import _run_testcases

    print("Testing HITRAN parsing: ", _run_testcases())
    from radis.test.io.test_query import _run_testcases

    print("Testing HITRAN fetch: ", _run_testcases())
