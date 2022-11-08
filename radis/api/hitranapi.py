# -*- coding: utf-8 -*-
"""
Summary
-------

HITRAN database parser


Routine Listing
---------------

- :func:`~radis.api.hitranapi.hit2df`
- :func:`~radis.api.hitranapi.parse_local_quanta`
- :func:`~radis.api.hitranapi.parse_global_quanta`


-------------------------------------------------------------------------------


"""


import os
import sys

# from radis.test.utils import getTestFile
import time
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
    HITRAN_CLASS10,
    HITRAN_GROUP1,
    HITRAN_GROUP2,
    HITRAN_GROUP3,
    HITRAN_GROUP4,
    HITRAN_GROUP5,
    HITRAN_GROUP6,
    get_molecule,
)

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


PARAMETER_GROUPS_HITRAN = {
    "par_line": "PARLIST_DOTPAR",
    "id": "PARLIST_ID",
    "standard": "PARLIST_STANDARD",
    "labels": "PARLIST_LABELS",
    "voigt": "PARLIST_VOIGT_ALL",
    "ht": "PARLIST_HT_ALL",
}


def cast_to_int64_with_missing_values(dg, keys):
    """replace missing values of int64 columns with -1"""
    for c in keys:
        if dg.dtypes[c] != int64:
            dg[c].replace(
                r"^\s+$", -1, regex=True, inplace=True
            )  # replace empty strings by -1, e.g. HCN
            # Warning: -1 may be a valid non-equilibirum quantum number for some
            # molecules, e.g. H2O, see https://github.com/radis/radis/issues/280#issuecomment-896120510
            dg[c] = dg[c].fillna(-1).astype(int64)  # replace nans with -1


def hit2df(
    fname,
    cache=True,
    verbose=True,
    drop_non_numeric=True,
    load_wavenum_min=None,
    load_wavenum_max=None,
    engine="pytables",
    parse_quanta=True,
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

    Returns
    -------
    df: pandas Dataframe
        dataframe containing all lines and parameters


    References
    ----------

    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__



    Notes
    -----

    Performances: see CDSD-HITEMP parser


    See Also
    --------

    :func:`~radis.api.cdsdapi.cdsd2df`
    """
    metadata = {}
    # Last modification time of the original file :
    metadata["last_modification"] = time.ctime(getmtime(fname))
    if verbose >= 2:
        print("Opening file {0} (cache={1})".format(fname, cache))
        print("Last modification time: {0}".format(metadata["last_modification"]))
    if load_wavenum_min and load_wavenum_max:
        assert load_wavenum_min < load_wavenum_max

    columns = columns_2004

    # Use cache file if possible
    fcache = DataFileManager(engine).cache_file(fname)
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
            "You're trying to read a binary file {0} ".format(fname)
            + "instead of an HITRAN file"
        ) from err

    df = parse_hitran_file(fname, columns)

    df = post_process_hitran_data(
        df,
        molecule=mol,
        parse_quanta=parse_quanta,
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
            print(
                "Generating cache file {0} with metadata :\n{1}".format(
                    fcache, new_metadata
                )
            )
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
                print("An error occured in cache file generation. Lookup access rights")
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
):
    """Parsing non-equilibrum parameters in HITRAN/HITEMP [1]_ file to and return final Pandas Dataframe

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

    Returns
    -------
    df: pandas Dataframe
        dataframe containing all lines and parameters


    References
    ----------

    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__



    Notes
    -----

    Performances: see CDSD-HITEMP parser


    See Also
    --------

    :func:`~radis.io.cdsd.cdsd2df`
    """

    # %% Post processing

    # assert one molecule per database only. Else the groupbase data reading
    # above doesnt make sense
    nmol = len(df["id"].unique())
    if nmol == 0:
        raise ValueError("Databank looks empty")
    elif nmol != 1:
        # Crash, give explicity error messages
        try:
            secondline = df.iloc[1]
        except IndexError:
            secondline = ""
        raise ValueError(
            "Multiple molecules in database ({0} : {1}). Current ".format(
                nmol, [get_molecule(idi) for idi in df["id"].unique()]
            )
            + "spectral code only computes 1 species at the time. Use MergeSlabs. "
            + "Verify the parsing was correct by looking at the first row below: "
            + "\n{0}".format(df.iloc[0])
            + "\n----------------\nand the second row "
            + "below: \n{0}".format(secondline)
        )

    if parse_quanta:
        # Add local quanta attributes, based on the HITRAN group
        try:
            df = parse_local_quanta(df, molecule, verbose=verbose)
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
            df = parse_global_quanta(df, molecule, verbose=verbose)
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
        df = drop_object_format_columns(df, verbose=verbose)

    return df


# %% Hitran global quanta classes


def _parse_HITRAN_class1(df, verbose=True):
    r"""Diatomic molecules: CO, HF, HCl, HBr, HI, N2, NO+


    Parameters
    ----------
    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----
    HITRAN syntax [1]_ :

    >>>       v
    >>>  13x I2

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """

    # 1. Parse
    dgu = df["globu"].astype(str).str.extract(r"[ ]{13}(?P<vu>[\d ]{2})", expand=True)
    dgl = df["globl"].astype(str).str.extract(r"[ ]{13}(?P<vl>[\d ]{2})", expand=True)

    # 2. Convert to numeric
    cast_to_int64_with_missing_values(dgu, ["vu"])
    cast_to_int64_with_missing_values(dgl, ["vl"])

    # 3. Clean
    del df["globu"]
    del df["globl"]

    return pd.concat([df, dgu, dgl], axis=1)


def _parse_HITRAN_class2(df, verbose=True):
    r"""Diatomic molecules with different electronic levels: O2


    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

    HITRAN syntax [1]_ :

    >>>      X   v
    >>>  12x A1 I2

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    if verbose > 2:
        print(
            f"parse_global_quanta not implemented for molecules of HITRAN class 2 ({HITRAN_CLASS2}). Non-LTE calculations will not be possible."
        )
    return df


def _parse_HITRAN_class3(df, verbose=True):
    r"""Diatomic molecules with doublet-Pi electronic state: NO, OH, ClO


    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -------

    HITRAN syntax [1]_:

    >>>      X i     v1
    >>>  7x A1 A3 2x I2

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    if verbose > 2:
        print(
            "parse_global_quanta not implemented for molecules of HITRAN class 3 ({HITRAN_CLASS3}). Non-LTE calculations will not be possible."
        )
    return df


def _parse_HITRAN_class4(df, verbose=True):
    r"""Parse linear triatomic class in HITRAN [1]_: N2O, OCS, HCN

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database

    Notes
    -----

    HITRAN syntax:

    >>>     v1 v2 l2 v3
    >>>  7x I2 I2 I2 I2

    Note: I2 in regexp: [\d ]{2}

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """

    # 1. Parse
    dgu = (
        df["globu"]
        .astype(str)
        .str.extract(
            r"[ ]{7}(?P<v1u>[\d ]{2})(?P<v2u>[\d ]{2})(?P<l2u>[\d ]{2})(?P<v3u>[\d ]{2})",
            expand=True,
        )
    )
    dgl = (
        df["globl"]
        .astype(str)
        .str.extract(
            r"[ ]{7}(?P<v1l>[\d ]{2})(?P<v2l>[\d ]{2})(?P<l2l>[\d ]{2})(?P<v3l>[\d ]{2})",
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


def _parse_HITRAN_class5(df, verbose=True):
    r"""Parse linear triatomic with large Fermi resonance in HITRAN [1]_: CO2

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database

    Notes
    -----

    HITRAN syntax:

    >>>     v1 v2 l2 v3 r
    >>>  6x I2 I2 I2 I2 I1

    Note: I2 in regexp: [\d ]{2}

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """

    # 1. Parse
    dgu = (
        df["globu"]
        .astype(str)
        .str.extract(
            r"[ ]{6}(?P<v1u>[\d ]{2})(?P<v2u>[\d ]{2})(?P<l2u>[\d ]{2})(?P<v3u>[\d ]{2})(?P<ru>\d)",
            expand=True,
        )
    )
    dgl = (
        df["globl"]
        .astype(str)
        .str.extract(
            r"[ ]{6}(?P<v1l>[\d ]{2})(?P<v2l>[\d ]{2})(?P<l2l>[\d ]{2})(?P<v3l>[\d ]{2})(?P<rl>\d)",
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


def _parse_HITRAN_class6(df, verbose=True):
    r"""Parse non-linear triatomic in HITRAN [1]_: H2O, O3, SO2, NO2, HOCl, H2S, HO2, HOBr

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database

    Notes
    -----

    HITRAN syntax:

    >>>     v1 v2 v3
    >>>  9x I2 I2 I2

    Note: I2 in regexp: [\d ]{2}

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """

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


def _parse_HITRAN_class7(df, verbose=True):
    r"""Parse linear tetratomic in HITRAN [1]_: C2H2

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database

    Notes
    -----

    HITRAN syntax:

    >>>

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    if verbose > 2:
        print(
            "parse_global_quanta not implemented for molecules of HITRAN class 7 ({HITRAN_CLASS7}). Non-LTE calculations will not be possible."
        )
    return df


def _parse_HITRAN_class8(df, verbose=True):
    r"""Pyramidal tetratomic in HITRAN [1]_: NH3, PH3


    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

    HITRAN syntax:

    >>>

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    if verbose > 2:
        print(
            "parse_global_quanta not implemented for molecules of HITRAN class 8 ({HITRAN_CLASS8}). Non-LTE calculations will not be possible."
        )
    return df


def _parse_HITRAN_class9(df, verbose=True):
    r"""Non-linear tetratomic in HITRAN [1]_: H2CO, H2O2, COF2


    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

    HITRAN syntax:

    >>>

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    if verbose > 2:
        print(
            "parse_global_quanta not implemented for molecules of HITRAN class 9 ({HITRAN_CLASS9}). Non-LTE calculations will not be possible."
        )
    return df


def _parse_HITRAN_class10(df, verbose=True):
    r"""Pentatomic or greater polyatomic in HITRAN [1]_


    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

    HITRAN syntax:

    >>>

    References
    ----------

    .. [1] `Table 3 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__

    """
    if verbose > 2:
        print(
            "parse_global_quanta not implemented for molecules of HITRAN class 10. Non-LTE calculations will not be possible."
        )
    return df


# %% HITRAN Local quanta


def _parse_HITRAN_group1(df, verbose=True):
    r"""Parse asymmetric rotors (:py:attr:`~radis.db.classes.HITRAN_GROUP1` ):
    H2O, O3, SO2, NO2, HNO3, H2CO, HOCl, H2O2, COF2, H2S, HO2, HCOOH, ClONO2, HOBr, C2H4

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

    HITRAN syntax: [1]_


    References
    ----------

    .. [1] `Table 4 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """

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


def _parse_HITRAN_group2(df, verbose=True):
    r"""Parse diatomic and linear molecules (:py:attr:`~radis.db.classes.HITRAN_GROUP2` ):
    CO2, N2O, CO, HF, HCl, HBr, HI, OCS, N2, HCN, C2H2, NO+

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

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
    dgu = df["locu"].astype(str).str.extract(r"[ ]{10}(?P<Fu>.{5})", expand=True)
    # Ref [1] : locl
    # --------------
    #     | Br  | J'' | Sym''| F'' |
    # 5X  | A1  | I3  | A1   | A5  |
    dgl = (
        df["locl"]
        .astype(str)
        .str.extract(
            r"[ ]{5}(?P<branch>[\S]{1})(?P<jl>[\d ]{3})(?P<syml>.)(?P<Fl>.{5})",
            expand=True,
        )
    )

    # 2. Convert to numeric

    # dgl['jl'] = dgl.jl.apply(pd.to_numeric)
    cast_to_int64_with_missing_values(dgl, ["jl"])

    # 3. Clean
    del df["locu"]
    del df["locl"]

    return pd.concat([df, dgu, dgl], axis=1)


#                # 10X included in Fu
#               'Fu',     ('a15',  str,  'upper state total angular momentum including nuclear spin'  ,''         )),(
#               #'locl',   ('a15',  str,   'electronic and vibrational local lower quanta'  ,''                      )),(
#               'branch', ('a6',   str,     'O, P, Q, R, S branch symbol'                  ,''                      )),(
#               'jl',     ('a3',   int,    'lower state rotational quantum number'         ,''                      )),(
#               'sym',    ('a1',   str,     'symmetry'                                     ,''                      )),(
#               'Fl',     ('a5',   str,   'lower state total angular momentum including nuclear spin', ''         )),(


def _parse_HITRAN_group3(df, verbose=True):
    r"""Parse Spherical rotors (:py:attr:`~radis.db.classes.HITRAN_GROUP3` ) :
    SF6, CH4

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

    HITRAN syntax [1]_:


    References
    ----------

    .. [1] `Table 4 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """
    if verbose > 2:
        print(
            f"parse_local_quanta not implemented for molecules of HITRAN group 3 ({HITRAN_GROUP3}). Non-LTE calculations will not be possible."
        )
    return df


def _parse_HITRAN_group4(df, verbose=True):
    r"""Parse symmetric rotors (:py:attr:`~radis.db.classes.HITRAN_GROUP4` ):
    CH3D, CH3Cl, C2H6, NH3, PH3, CH3OH

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

    HITRAN syntax [1]_:


    References
    ----------

    .. [1] `Table 4 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """
    if verbose > 2:
        print(
            f"parse_local_quanta not implemented for molecules of HITRAN group 4 ({HITRAN_GROUP4}). Non-LTE calculations will not be possible."
        )
    return df


def _parse_HITRAN_group5(df, verbose=True):
    r"""Parse Triplet-Sigma ground electronic states (:py:attr:`~radis.db.classes.HITRAN_GROUP5` ) :
    O2

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

    HITRAN syntax [1]_:


    References
    ----------

    .. [1] `Table 4 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """
    if verbose > 2:
        print(
            "parse_local_quanta not implemented for molecules of HITRAN group 5 ({HITRAN_GROUP5}). Non-LTE calculations will not be possible."
        )
    return df


def _parse_HITRAN_group6(df, verbose=True):
    r"""Parse Doublet-Pi ground electronic states (:py:attr:`~radis.db.classes.HITRAN_GROUP6` ) :
    NO, OH, ClO

    Parameters
    ----------

    df: pandas Dataframe
        lines read from a HITRAN-like database


    Notes
    -----

    HITRAN syntax [1]_:


    References
    ----------

    .. [1] `Table 4 of Rothman et al. HITRAN 2004 <https://www.cfa.harvard.edu/hitran/Download/HITRAN04paper.pdf>`__


    """
    if verbose > 2:
        print(
            "parse_local_quanta not implemented for molecules of HITRAN group 6 ({HITRAN_GROUP6}). Non-LTE calculations will not be possible."
        )
    return df


# %% Reading function


def parse_local_quanta(df, mol, verbose=True):
    r"""
    Parameters
    ----------

    df: pandas Dataframe

    mol: str
        molecule name
    """

    if mol in HITRAN_GROUP1:
        df = _parse_HITRAN_group1(df, verbose=verbose)
    elif mol in HITRAN_GROUP2:
        df = _parse_HITRAN_group2(df, verbose=verbose)
    elif mol in HITRAN_GROUP3:
        df = _parse_HITRAN_group3(df, verbose=verbose)
    elif mol in HITRAN_GROUP4:
        df = _parse_HITRAN_group4(df, verbose=verbose)
    elif mol in HITRAN_GROUP5:
        df = _parse_HITRAN_group5(df, verbose=verbose)
    elif mol in HITRAN_GROUP6:
        df = _parse_HITRAN_group6(df, verbose=verbose)
    else:
        raise ValueError(
            "Unknown group for molecule {0}. Cant parse local quanta".format(mol)
        )

    return df


def parse_global_quanta(df, mol, verbose=True):
    r"""

    Parameters
    ----------

    df: pandas Dataframe

    mol: str
        molecule name
    """

    if mol in HITRAN_CLASS1:
        df = _parse_HITRAN_class1(df, verbose=verbose)
    elif mol in HITRAN_CLASS2:
        df = _parse_HITRAN_class2(df, verbose=verbose)
    elif mol in HITRAN_CLASS3:
        df = _parse_HITRAN_class3(df, verbose=verbose)
    elif mol in HITRAN_CLASS4:
        df = _parse_HITRAN_class4(df, verbose=verbose)
    elif mol in HITRAN_CLASS5:
        df = _parse_HITRAN_class5(df, verbose=verbose)
    elif mol in HITRAN_CLASS6:
        df = _parse_HITRAN_class6(df, verbose=verbose)
    elif mol in HITRAN_CLASS7:
        df = _parse_HITRAN_class7(df, verbose=verbose)
    elif mol in HITRAN_CLASS8:
        df = _parse_HITRAN_class8(df, verbose=verbose)
    elif mol in HITRAN_CLASS9:
        df = _parse_HITRAN_class9(df, verbose=verbose)
    elif mol in HITRAN_CLASS10:
        df = _parse_HITRAN_class10(df, verbose=verbose)
    else:
        raise ValueError(
            "Unknown class for molecule {0}. Cant parse global quanta".format(mol)
        )

    return df


#%%


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

    def get_filenames(self):
        if self.engine == "vaex":
            return [join(self.local_databases, f"{self.molecule}.hdf5")]
        elif self.engine == "pytables":
            return [join(self.local_databases, f"{self.molecule}.h5")]
        else:
            raise NotImplementedError()

    def download_and_parse(self, local_file, cache=True, parse_quanta=True):
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

        wmin = 1
        wmax = 40000

        def download_all_hitran_isotopes(molecule, directory, extra_params):
            """Blindly try to download all isotpes 1 - 9 for the given molecule

            .. warning::
                this won't be able to download higher isotopes (ex : isotope 10-11-12 for CO2)
                Neglected for the moment, they're irrelevant for most calculations anyway
            .. Isotope Missing:
                When an isotope is missing for a particular molecule then a key error `(molecule_id, isotope_id)
                is raised.

            """
            directory = abspath(expanduser(directory))

            # create temp folder :
            from radis.misc.basics import make_folders

            make_folders(*split(directory))

            db_begin(directory)
            isotope_list = []
            data_file_list = []
            header_file_list = []
            for iso in range(1, 10):
                file = f"{molecule}_{iso}"
                if exists(join(directory, file + ".data")):
                    if cache == "regen":
                        # remove without printing message
                        os.remove(join(directory, file + ".data"))
                    else:
                        from radis.misc.printer import printr

                        printr(
                            "File already exist: {0}. Deleting it.`".format(
                                join(directory, file + ".data")
                            )
                        )
                        os.remove(join(directory, file + ".data"))
                try:
                    if extra_params == "all":
                        fetch(
                            file,
                            get_molecule_identifier(molecule),
                            iso,
                            wmin,
                            wmax,
                            ParameterGroups=[*PARAMETER_GROUPS_HITRAN],
                        )
                    elif extra_params is None:
                        fetch(file, get_molecule_identifier(molecule), iso, wmin, wmax)
                    else:
                        raise ValueError("extra_params can only be 'all' or None ")
                except KeyError as err:
                    list_pattern = ["(", ",", ")"]
                    import re

                    if (
                        set(list_pattern).issubset(set(str(err)))
                        and len(re.findall("\d", str(err))) >= 2
                        and get_molecule_identifier(molecule)
                        == int(
                            re.findall(r"[\w']+", str(err))[0]
                        )  # The regex are cryptic
                    ):
                        # Isotope not defined, go to next isotope
                        continue
                    else:
                        raise KeyError("Error: {0}".format(str(err)))
                else:
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
        isotope_list, data_file_list, header_file_list = download_all_hitran_isotopes(
            molecule, tempdir, extra_params
        )

        writer = self.get_datafile_manager()

        # Create HDF5 cache file for all isotopes
        Nlines = 0
        for iso, data_file in zip(isotope_list, data_file_list):
            df = pd.DataFrame(LOCAL_TABLE_CACHE[data_file.split(".")[0]]["data"])
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

    def register(self):
        """register in ~/radis.json"""

        from radis.db import MOLECULES_LIST_NONEQUILIBRIUM

        local_files = self.get_filenames()

        if self.wmin is None or self.wmax is None:
            print(
                "Somehow wmin and wmax was not given for this database. Reading from the files"
            )
            ##  fix:
            # (can happen if database was downloaded & parsed, but registration failed a first time)
            df_full = self.load(
                local_files,
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

        dict_entries = {
            "info": info,
            "path": local_files,
            "format": "hdf5-radisdb",
            "parfuncfmt": "hapi",
            "wavenumber_min": self.wmin,
            "wavenumber_max": self.wmax,
            "download_date": self.get_today(),
        }

        # Add energy level calculation
        if self.molecule in MOLECULES_LIST_NONEQUILIBRIUM:
            dict_entries[
                "info"
            ] += " and RADIS spectroscopic constants for rovibrational energies (nonequilibrium)"
            dict_entries["levelsfmt"] = "radis"

        super().register(dict_entries)


# ======================================================
# %% Test


if __name__ == "__main__":

    from radis.test.io.test_hitran_cdsd import _run_testcases

    print("Testing HITRAN parsing: ", _run_testcases())
    from radis.test.io.test_query import _run_testcases

    print("Testing HITRAN fetch: ", _run_testcases())
