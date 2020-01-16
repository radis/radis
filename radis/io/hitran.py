# -*- coding: utf-8 -*-
"""
Summary
-------

HITRAN database parser 


Routine Listing
---------------

- :func:`~radis.io.hitran.hit2df`
- :func:`~radis.io.hitran.get_molecule`
- :func:`~radis.io.hitran.get_molecule_identifier`
- :func:`~radis.io.hitran.parse_local_quanta`
- :func:`~radis.io.hitran.parse_global_quanta`


-------------------------------------------------------------------------------


"""

from __future__ import print_function, absolute_import, division, unicode_literals

import sys
import pandas as pd
from collections import OrderedDict
from os.path import exists, splitext
import radis
from radis.io.tools import (
    parse_binary_file,
    drop_object_format_columns,
    replace_PQR_with_m101,
)
from radis.misc.cache_files import save_to_hdf, check_cache_file, get_cache_file


# %% Hitran groups and classes
# As defined in Rothman et al, "The HITRAN 2004 molecular spectroscopic database"
# Tables 3 and 4
# Groups and classes are needed to compute energy in nonequilibrium mode.

# Groups (define local quanta)

HITRAN_GROUP1 = [
    "H2O",
    "O3",
    "SO2",
    "NO2",
    "HNO3",
    "H2CO",
    "HOCl",
    "H2O2",
    "COF2",
    "H2S",
    "HO2",
    "HCOOH",
    "ClONO2",
    "HOBr",
    "C2H4",
]
"""str: asymmetric rotors"""

HITRAN_GROUP2 = [
    "CO2",
    "N2O",
    "CO",
    "HF",
    "HCl",
    "HBr",
    "HI",
    "OCS",
    "N2",
    "HCN",
    "C2H2",
    "NO+",
]
"""str: diatomic and linear molecules"""

HITRAN_GROUP3 = ["SF6", "CH4"]
"""str: Spherical rotors"""

HITRAN_GROUP4 = ["CH3D", "CH3Cl", "C2H6", "NH3", "PH3", "CH3OH"]
"""str: symmetric rotors"""

HITRAN_GROUP5 = ["O2"]
"""str: Triplet-Sigma ground electronic states"""

HITRAN_GROUP6 = ["NO", "OH", "ClO"]
"""str: Doublet-Pi ground electronic states"""

# Classes (define global quanta)

HITRAN_CLASS1 = ["CO", "HF", "HCl", "HBr", "HI", "N2", "NO+"]
"""str: Diatomic molecules with ? """

HITRAN_CLASS2 = ["O2"]
"""str: Diatomic molecules with different electronic levels"""

HITRAN_CLASS3 = ["NO", "OH", "ClO"]
"""str: Diatomic molecules with doublet-Pi electronic state"""

HITRAN_CLASS4 = ["N2O", "OCS", "HCN"]
"""str: Linear triatomic"""

HITRAN_CLASS5 = ["CO2"]
"""str: Linear triatomic with large Fermi resonance"""

HITRAN_CLASS6 = ["H2O", "O3", "SO2", "NO2", "HOCl", "H2S", "HO2", "HOBr"]
"""str: Non-linear triatomic"""

HITRAN_CLASS7 = ["C2H2"]
"""str: Linear tetratomic"""

HITRAN_CLASS8 = ["NH3", "PH3"]
"""str: Pyramidal tetratomic"""

HITRAN_CLASS9 = ["H2CO", "H2O2", "COF2"]
"""str: Non-linear tetratomic"""

HITRAN_CLASS10 = [
    "CH4",
    "CH3Cl",
    "C2H6",
    "HNO3",
    "SF6",
    "HCOOH",
    "ClONO2",
    "C2H4",
    "CH3OH",
]
"""str: Pentatomic or greater polyatomic"""
# 0.9.22: Removed CH3D has no HITRAN identifier

# %% HITRAN ids

trans = {
    "1": "H2O",
    "2": "CO2",
    "3": "O3",
    "4": "N2O",
    "5": "CO",
    "6": "CH4",
    "7": "O2",
    "8": "NO",
    "9": "SO2",
    "10": "NO2",
    "11": "NH3",
    "12": "HNO3",
    "13": "OH",
    "14": "HF",
    "15": "HCl",
    "16": "HBr",
    "17": "HI",
    "18": "ClO",
    "19": "OCS",
    "20": "H2CO",
    "21": "HOCl",
    "22": "N2",
    "23": "HCN",
    "24": "CH3Cl",
    "25": "H2O2",
    "26": "C2H2",
    "27": "C2H6",
    "28": "PH3",
    "29": "COF2",
    "30": "SF6",
    "31": "H2S",
    "32": "HCOOH",
    "33": "HO2",
    "34": "O",
    "35": "ClONO2",
    "36": "NO+",
    "37": "HOBr",
    "38": "C2H4",
    "39": "CH3OH",
    "40": "CH3Br",
    "41": "CH3CN",
    "42": "CF4",
    "43": "C4H2",
    "44": "HC3N",
    "45": "H2",
    "46": "CS",
    "47": "SO3",
    "48": "C2N2",
    "49": "COCl2",
}
HITRAN_MOLECULES = list(trans.values())
""" str: list of [HITRAN-2016]_ molecules. """

# %% Parsing functions

# General case
columns_2004 = OrderedDict(
    [
        (
            # name    # format # type  # description                                 # unit
            "id",
            ("a2", int, "Molecular number", ""),
        ),
        ("iso", ("a1", int, "isotope number", "")),
        ("wav", ("a12", float, "vacuum wavenumber", "cm-1")),
        ("int", ("a10", float, "intensity at 296K", "cm-1/(molecule/cm-2)",)),
        ("A", ("a10", float, "Einstein A coefficient", "s-1")),
        ("airbrd", ("a5", float, "air-broadened half-width at 296K", "cm-1.atm-1")),
        ("selbrd", ("a5", float, "self-broadened half-width at 296K", "cm-1.atm-1")),
        ("El", ("a10", float, "lower-state energy", "cm-1")),
        ("Tdpair", ("a4", float, "temperature-dependance exponent for Gamma air", "")),
        (
            "Pshft",
            ("a8", float, "air pressure-induced line shift at 296K", "cm-1.atm-1"),
        ),
        ("globu", ("a15", str, "electronic and vibrational global upper quanta", "")),
        ("globl", ("a15", str, "electronic and vibrational global lower quanta", "")),
        ("locu", ("a15", str, "electronic and vibrational local upper quanta", "")),
        ("locl", ("a15", str, "electronic and vibrational local lower quanta", "")),
        (
            "ierr",
            (
                "a6",
                str,
                "ordered list of indices corresponding to uncertainty estimates of transition parameters",
                "",
            ),
        ),
        (
            "iref",
            (
                "a12",
                str,
                "ordered list of reference identifiers for transition parameters",
                "",
            ),
        ),
        (
            "lmix",
            (
                "a1",
                str,
                "flag indicating the presence of additional data and code relating to line-mixing",
                "",
            ),
        ),
        ("gp", ("a7", float, "upper state degeneracy", "")),
        ("gpp", ("a7", float, "lower state degeneracy", "")),
    ]
)
""" OrderedDict: parsing order of HITRAN """


def hit2df(fname, count=-1, cache=False, verbose=True, drop_non_numeric=True):
    """ Convert a HITRAN/HITEMP [1]_ file to a Pandas dataframe 

    Parameters    
    ----------

    fname: str
        HITRAN-HITEMP file name 

    count: int
        number of items to read (-1 means all file)

    cache: boolean, or ``'regen'`` or ``'force'``
        if ``True``, a pandas-readable HDF5 file is generated on first access, 
        and later used. This saves on the datatype cast and conversion and
        improves performances a lot (but changes in the database are not 
        taken into account). If False, no database is used. If ``'regen'``, temp
        file are reconstructed. Default ``False``. 

    Other Parameters
    ----------------
    
    drop_non_numeric: boolean
        if ``True``, non numeric columns are dropped. This improves performances, 
        but make sure all the columns you need are converted to numeric formats 
        before hand. Default ``True``. Note that if a cache file is loaded it 
        will be left untouched.

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

    if verbose >= 2:
        print("Opening file {0} (cache={1})".format(fname, cache))

    columns = columns_2004

    # Use cache file if possible
    fcache = splitext(fname)[0] + ".h5"
    check_cache_file(fcache=fcache, use_cached=cache, verbose=verbose)
    if cache and exists(fcache):
        return get_cache_file(fcache, verbose=verbose)

    # Detect the molecule by reading the start of the file
    with open(fname) as f:
        mol = get_molecule(int(f.read(2)))

    # %% Start reading the full file

    df = parse_binary_file(fname, columns, count)

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

    # dd local quanta attributes, based on the HITRAN group
    df = parse_local_quanta(df, mol)

    # Add global quanta attributes, based on the HITRAN class
    df = parse_global_quanta(df, mol)

    # Remove non numerical attributes
    if drop_non_numeric:
        if "branch" in df:
            replace_PQR_with_m101(df)
        df = drop_object_format_columns(df, verbose=verbose)

    # cached file mode but cached file doesn't exist yet (else we had returned)
    if cache:
        if verbose:
            print("Generating cached file: {0}".format(fcache))
        try:
            save_to_hdf(
                df,
                fcache,
                metadata={},
                version=radis.__version__,
                key="df",
                overwrite=True,
                verbose=verbose,
            )
        except:
            if verbose:
                print(sys.exc_info())
                print("An error occured in cache file generation. Lookup access rights")
            pass

    return df


# %% Hitran global quanta classes


def _parse_HITRAN_class1(df):
    r""" Diatomic molecules: CO, HF, HCl, HBr, HI, N2, NO+


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
    dgu = df["globu"].str.extract(r"[ ]{13}(?P<vu>[\d ]{2})", expand=True)
    dgl = df["globl"].str.extract(r"[ ]{13}(?P<vl>[\d ]{2})", expand=True)

    # 2. Convert to numeric
    dgu = dgu.apply(pd.to_numeric)
    dgl = dgl.apply(pd.to_numeric)

    # 3. Clean
    del df["globu"]
    del df["globl"]

    return pd.concat([df, dgu, dgl], axis=1)


def _parse_HITRAN_class2(df):
    r""" Diatomic molecules with different electronic levels: O2


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
    r""" Diatomic molecules with doublet-Pi electronic state: NO, OH, ClO


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
    r""" Parse linear triatomic class in HITRAN [1]_: N2O, OCS, HCN

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
    dgu = df["globu"].str.extract(
        r"[ ]{7}(?P<v1u>[\d ]{2})(?P<v2u>[\d ]{2})(?P<l2u>[\d ]{2})(?P<v3u>[\d ]{2})",
        expand=True,
    )
    dgl = df["globl"].str.extract(
        r"[ ]{7}(?P<v1l>[\d ]{2})(?P<v2l>[\d ]{2})(?P<l2l>[\d ]{2})(?P<v3l>[\d ]{2})",
        expand=True,
    )

    # 2. Convert to numeric
    dgu = dgu.apply(pd.to_numeric)
    dgl = dgl.apply(pd.to_numeric)

    # 3. Clean
    del df["globu"]
    del df["globl"]

    return pd.concat([df, dgu, dgl], axis=1)


def _parse_HITRAN_class5(df):
    r""" Parse linear triatomic with large Fermi resonance in HITRAN [1]_: CO2

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
    dgu = df["globu"].str.extract(
        r"[ ]{6}(?P<v1u>[\d ]{2})(?P<v2u>[\d ]{2})(?P<l2u>[\d ]{2})(?P<v3u>[\d ]{2})(?P<ru>\d)",
        expand=True,
    )
    dgl = df["globl"].str.extract(
        r"[ ]{6}(?P<v1l>[\d ]{2})(?P<v2l>[\d ]{2})(?P<l2l>[\d ]{2})(?P<v3l>[\d ]{2})(?P<rl>\d)",
        expand=True,
    )

    # 2. Convert to numeric
    dgu = dgu.apply(pd.to_numeric)
    dgl = dgl.apply(pd.to_numeric)

    # 3. Clean
    del df["globu"]
    del df["globl"]

    return pd.concat([df, dgu, dgl], axis=1)


def _parse_HITRAN_class6(df):
    r""" Parse non-linear triatomic in HITRAN [1]_: H2O, O3, SO2, NO2, HOCl, H2S, HO2, HOBr

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
    dgu = df["globu"].str.extract(
        #        '[ ]{9}(?P<v1u>[\d ]{2})(?P<v2u>[\d ]{2})(?P<v3u>[\d ]{2})',
        r"[ ]{9}(?P<v1u>[\-\d ]{2})(?P<v2u>[\-\d ]{2})(?P<v3u>[\-\d ]{2})",
        expand=True,
    )
    dgl = df["globl"].str.extract(
        #        '[ ]{9}(?P<v1l>[\d ]{2})(?P<v2l>[\d ]{2})(?P<v3l>[\d ]{2})',
        r"[ ]{9}(?P<v1l>[\-\d ]{2})(?P<v2l>[\-\d ]{2})(?P<v3l>[\-\d ]{2})",
        expand=True,
    )
    # ... note @EP: in HITRAN H2O files, for iso=2, vibrational levels are
    # ... somehow negative. The regex above is adapted to catch negation signs with \-

    # 2. Convert to numeric
    dgu = dgu.apply(pd.to_numeric)
    dgl = dgl.apply(pd.to_numeric)

    # 3. Clean
    del df["globu"]
    del df["globl"]

    return pd.concat([df, dgu, dgl], axis=1)


def _parse_HITRAN_class7(df):
    r""" Parse linear tetratomic in HITRAN [1]_: C2H2

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
    r""" Pyramidal tetratomic in HITRAN [1]_: NH3, PH3


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
    r""" Non-linear tetratomic in HITRAN [1]_: H2CO, H2O2, COF2


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
    r""" Pentatomic or greater polyatomic in HITRAN [1]_


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
    r""" 

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
    dgu = df["locu"].str.extract(
        r"(?P<ju>[\d ]{3})(?P<Kau>[\-\d ]{3})(?P<Kcu>[\-\d ]{3})(?P<Fu>.{5})(?P<symu>.)",
        expand=True,
    )
    # Ref [1] : locl
    # --------------
    # J'' | Ka''| Kc''| F'' | Sym''
    # I3  | I3  | I3  | A5  | A1
    dgl = df["locl"].str.extract(
        r"(?P<jl>[\d ]{3})(?P<Kal>[\-\d ]{3})(?P<Kcl>[\-\d ]{3})(?P<Fl>.{5})(?P<syml>.)",
        expand=True,
    )
    # ... note @EP: in HITRAN H2O files, for iso=2, the Kau, Kcu can somehow
    # ... be negative. The regex above is adapted to catch negation signs with \-

    # 2. Convert to numeric
    for k in ["ju", "Kau", "Kcu"]:
        dgu[k] = pd.to_numeric(dgu[k])
    for k in ["jl", "Kal", "Kcl"]:
        dgl[k] = pd.to_numeric(dgl[k])

    # 3. Clean
    del df["locu"]
    del df["locl"]

    return pd.concat([df, dgu, dgl], axis=1)


def _parse_HITRAN_group2(df):
    r""" 

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
    dgu = df["locu"].str.extract(r"[ ]{10}(?P<Fu>.{5})", expand=True)
    # Ref [1] : locl
    # --------------
    #     | Br  | J'' | Sym''| F'' |
    # 5X  | A1  | I3  | A1   | A5  |
    dgl = df["locl"].str.extract(
        r"[ ]{5}(?P<branch>[\S]{1})(?P<jl>[\d ]{3})(?P<syml>.)(?P<Fl>.{5})", expand=True
    )

    # 2. Convert to numeric

    # dgl['jl'] = dgl.jl.apply(pd.to_numeric)
    dgl["jl"] = pd.to_numeric(dgl.jl)

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
    r""" 

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
    r""" 

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
    r""" 

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
    r""" 

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

    # had some problems with bytes types
    df["locu"] = df.locu.astype(str)
    df["locl"] = df.locl.astype(str)

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

    # had some problems with bytes types
    df["globu"] = df.globu.astype(str)
    df["globl"] = df.globl.astype(str)

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


def get_molecule_identifier(molecule_name):
    r"""
    For a given input molecular formula, return the corresponding HITRAN molecule 
    identifier number [1]_.


    Parameters
    ----------
    molecular_formula : str
        The string describing the molecule.


    Returns
    -------
    M: int
        The HITRAN molecular identified number.


    References
    ----------

    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__

    Function is from from https://github.com/nzhagen/hitran/blob/master/hitran.py

    """

    # Invert the dictionary.
    trans2 = {v: k for k, v in trans.items()}

    try:
        return int(trans2[molecule_name])
    except KeyError:
        raise NotImplementedError(
            "Molecule '{0}' not supported. Choose one of {1}".format(
                molecule_name, list(trans2.keys())
            )
        )


def get_molecule(molecule_id):
    r"""
    For a given input molecular identifier, return the corresponding HITRAN 
    molecule name [1]_.


    Parameters    
    ----------

    molecular_id: str
        Hitran identifier of the molecule.


    References
    ----------

    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__

    """

    # assert str
    id = "{:d}".format(int(molecule_id))

    try:
        return trans[id]
    except KeyError:
        raise NotImplementedError(
            "Molecule ID '{0}' unknown. Choose one of {1}".format(molecule_id, trans)
        )


# ======================================================
# %% Test


if __name__ == "__main__":

    from radis.test.test_io import _run_testcases

    print("Testing HITRAN parsing: ", _run_testcases())
