# -*- coding: utf-8 -*-
"""
Summary
-------

HITRAN database parser


Routine Listing
---------------

- :func:`~radis.io.hitran.hit2df`
- :func:`~radis.io.hitran.parse_local_quanta`
- :func:`~radis.io.hitran.parse_global_quanta`


-------------------------------------------------------------------------------


"""


import sys

# from radis.test.utils import getTestFile
import time
from collections import OrderedDict
from os.path import exists, getmtime

import pandas as pd
from numpy import int64

import radis
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
from radis.io.cache_files import cache_file_name, load_h5_cache_file, save_to_hdf
from radis.io.tools import (
    drop_object_format_columns,
    parse_hitran_file,
    replace_PQR_with_m101,
)

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


def cast_to_int64_with_missing_values(dg, keys):
    """replace missing values of int64 columns with -1"""
    for c in keys:
        if dg.dtypes[c] != int64:
            dg[c] = dg[c].fillna(-1).astype(int64)


def hit2df(
    fname,
    cache=True,
    verbose=True,
    drop_non_numeric=True,
    load_wavenum_min=None,
    load_wavenum_max=None,
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
        wavenumbers above/below the specified value. See :py:func`~radis.io.cache_files.load_h5_cache_file`.
        Default ``'None'``.

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
    fcache = cache_file_name(fname)
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
            valid_if_metadata_is=metadata,
            relevant_if_metadata_above=relevant_if_metadata_above,
            relevant_if_metadata_below=relevant_if_metadata_below,
            current_version=radis.__version__,
            last_compatible_version=radis.config["OLDEST_COMPATIBLE_VERSION"],
            verbose=verbose,
        )
        if df is not None:
            return df

    # Detect the molecule by reading the start of the file
    try:
        with open(fname) as f:
            mol = get_molecule(int(f.read(2)))
    except UnicodeDecodeError as err:
        raise ValueError(
            "You're trying to read a binary file {0} ".format(fname)
            + "instead of an HITRAN file"
        ) from err

    # %% Start reading the full file

    df = parse_hitran_file(fname, columns)

    # %% Post processing

    # assert one molecule per database only. Else the groupbase data reading
    # above doesnt make sense
    nmol = len(set(df["id"]))
    if nmol == 0:
        raise ValueError("Databank looks empty")
    elif nmol != 1:
        # Crash, give explicity error messages
        try:
            secondline = df.iloc[1]
        except IndexError:
            secondline = ""
        raise ValueError(
            "Multiple molecules in database ({0}). Current ".format(nmol)
            + "spectral code only computes 1 species at the time. Use MergeSlabs. "
            + "Verify the parsing was correct by looking at the first row below: "
            + "\n{0}".format(df.iloc[0])
            + "\n----------------\nand the second row "
            + "below: \n{0}".format(secondline)
        )

    # Add local quanta attributes, based on the HITRAN group
    try:
        df = parse_local_quanta(df, mol)
    except ValueError as err:
        # Empty strings (unlabelled lines) have been reported for HITEMP2010-H2O.
        # In this case, do not parse (makes non-equilibrium calculations impossible).
        # see https://github.com/radis/radis/issues/211
        if verbose:
            print(str(err))
            print("-" * 10)
            print(
                f"Impossible to parse local quanta in {fname}, probably an unlabelled line. Ignoring, but nonequilibrium calculations will not be possible. See details above."
            )

    # Add global quanta attributes, based on the HITRAN class
    try:
        df = parse_global_quanta(df, mol)
    except ValueError as err:
        # Empty strings (unlabelled lines) have been reported for HITEMP2010-H2O.
        # In this case, do not parse (makes non-equilibrium calculations impossible).
        # see https://github.com/radis/radis/issues/211
        if verbose:
            print(str(err))
            print("-" * 10)
            print(
                f"Impossible to parse global quanta in {fname}, probably an unlabelled line. Ignoring, but nonequilibrium calculations will not be possible. See details above."
            )

    # Remove non numerical attributes
    if drop_non_numeric:
        if "branch" in df:
            replace_PQR_with_m101(df)
        df = drop_object_format_columns(df, verbose=verbose)

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
        try:
            save_to_hdf(
                df,
                fcache,
                metadata=new_metadata,
                version=radis.__version__,
                key="df",
                overwrite=True,
                verbose=verbose,
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


# %% Hitran global quanta classes


def _parse_HITRAN_class1(df):
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


def _parse_HITRAN_class2(df):
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
    print("parse_global_quanta not implemented for molecules of HITRAN class 2")
    return df


def _parse_HITRAN_class3(df):
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
    print("parse_global_quanta not implemented for molecules of HITRAN class 3")
    return df


def _parse_HITRAN_class4(df):
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


def _parse_HITRAN_class5(df):
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


def _parse_HITRAN_class6(df):
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


def _parse_HITRAN_class7(df):
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
    print("parse_global_quanta not implemented for molecules of HITRAN class 7")
    return df


def _parse_HITRAN_class8(df):
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
    print("parse_global_quanta not implemented for molecules of HITRAN class 8")
    return df


def _parse_HITRAN_class9(df):
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
    print("parse_global_quanta not implemented for molecules of HITRAN class 9")
    return df


def _parse_HITRAN_class10(df):
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
    print("parse_global_quanta not implemented for molecules of HITRAN class 10")
    return df


# %% HITRAN Local quanta


def _parse_HITRAN_group1(df):
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


def _parse_HITRAN_group2(df):
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


def _parse_HITRAN_group3(df):
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
    print("parse_local_quanta not implemented for molecules of HITRAN group 3")
    return df


def _parse_HITRAN_group4(df):
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
    print("parse_local_quanta not implemented for molecules of HITRAN group 4")
    return df


def _parse_HITRAN_group5(df):
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
    print("parse_local_quanta not implemented for molecules of HITRAN group 5")
    return df


def _parse_HITRAN_group6(df):
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
    print("parse_local_quanta not implemented for molecules of HITRAN group 6")
    return df


# %% Reading function


def parse_local_quanta(df, mol):
    r"""
    Parameters
    ----------

    df: pandas Dataframe

    mol: str
        molecule name
    """

    if mol in HITRAN_GROUP1:
        df = _parse_HITRAN_group1(df)
    elif mol in HITRAN_GROUP2:
        df = _parse_HITRAN_group2(df)
    elif mol in HITRAN_GROUP3:
        df = _parse_HITRAN_group3(df)
    elif mol in HITRAN_GROUP4:
        df = _parse_HITRAN_group4(df)
    elif mol in HITRAN_GROUP5:
        df = _parse_HITRAN_group5(df)
    elif mol in HITRAN_GROUP6:
        df = _parse_HITRAN_group6(df)
    else:
        raise ValueError(
            "Unknown group for molecule {0}. Cant parse local quanta".format(mol)
        )

    return df


def parse_global_quanta(df, mol):
    r"""

    Parameters
    ----------

    df: pandas Dataframe

    mol: str
        molecule name
    """

    if mol in HITRAN_CLASS1:
        df = _parse_HITRAN_class1(df)
    elif mol in HITRAN_CLASS2:
        df = _parse_HITRAN_class2(df)
    elif mol in HITRAN_CLASS3:
        df = _parse_HITRAN_class3(df)
    elif mol in HITRAN_CLASS4:
        df = _parse_HITRAN_class4(df)
    elif mol in HITRAN_CLASS5:
        df = _parse_HITRAN_class5(df)
    elif mol in HITRAN_CLASS6:
        df = _parse_HITRAN_class6(df)
    elif mol in HITRAN_CLASS7:
        df = _parse_HITRAN_class7(df)
    elif mol in HITRAN_CLASS8:
        df = _parse_HITRAN_class8(df)
    elif mol in HITRAN_CLASS9:
        df = _parse_HITRAN_class9(df)
    elif mol in HITRAN_CLASS10:
        df = _parse_HITRAN_class10(df)
    else:
        raise ValueError(
            "Unknown class for molecule {0}. Cant parse global quanta".format(mol)
        )

    return df


# ======================================================
# %% Test


if __name__ == "__main__":

    from radis.test.test_io import _run_testcases

    print("Testing HITRAN parsing: ", _run_testcases())
