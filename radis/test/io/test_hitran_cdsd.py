# -*- coding: utf-8 -*-
"""Test parsers.

Notes
-----

Runs tests for radis/io so that they can be accessed by pytest (and hopefully
the CI test suite)

Examples
--------

Run all tests::

        pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)::

        pytest -m fast


"""

import os
from os.path import getmtime
from warnings import warn

import numpy as np
import pandas as pd
import pytest

from radis.api.cdsdapi import cdsd2df
from radis.api.hitranapi import hit2df, parse_global_quanta
from radis.api.tools import _format_dtype
from radis.misc.warning import IrrelevantFileWarning
from radis.test.utils import getTestFile, setup_test_line_databases


@pytest.mark.fast
def test_format_dtype_with_legacy_string_alias():
    assert _format_dtype([("id", "a2")]) == np.dtype([("id", "S2")])


@pytest.mark.fast
@pytest.mark.parametrize("dataframe_type", ["pandas", "vaex"])
@pytest.mark.parametrize("fast_parsing", [False, True])
@pytest.mark.parametrize(
    "molecule,upper,lower,expected",
    [
        ("CO", " " * 13 + " 4", " " * 13 + " 2", {"vu": 4, "vl": 2}),
        (
            "NO",
            " " * 6 + " X3/2" + "  " + " 4",
            " " * 6 + " A1/2" + "  " + " 2",
            {"electronic_stateu": "X", "omegau": "3/2", "vu": 4, "vl": 2},
        ),
        (
            "CO2",
            " " * 6 + " 1 2-2 34",
            " " * 6 + " 0 0 0 01",
            {"v1u": 1, "l2u": -2, "v3u": 3, "ru": 4, "rl": 1},
        ),
        (
            "HCN",
            " " * 7 + " 1-2 3 4",
            " " * 7 + " 0 1 2 3",
            {"v1u": 1, "v2u": -2, "l2u": 3, "v3l": 3},
        ),
        (
            "H2O",
            " " * 9 + " 1-2 3",
            " " * 9 + " 0 1 2",
            {"v1u": 1, "v2u": -2, "v3l": 2},
        ),
        (
            "PH3",
            " " * 5 + " 1 2 3 4A1",
            " " * 5 + " 0 0 0 0A2",
            {"v4u": 4, "vib_symmetryu": "A1", "vib_symmetryl": "A2"},
        ),
        (
            "NH3",
            " 12*4 56 7 A1  ",
            " 0000 00 0 A2  ",
            {"v3u": "*", "l4u": "6", "vib_symmetryl": "A2"},
        ),
        (
            "C2H2",
            " 123 4 5-1 2+ g",
            " 000 0 0 0 0- u",
            {"v4u": 4, "l4u": -1, "vib_parityu": "+", "vib_symmetryl": "u"},
        ),
        (
            "C4H2",
            " 123456789 e A1",
            " 000000000 f A2",
            {"v9u": 9, "vib_parityu": "e", "vib_symmetryl": "A2"},
        ),
        ("HC3N", "  1234567-1 2 3", "  0000000 0 0 0", {"v7u": 7, "l5u": -1, "l7u": 3}),
        (
            "C2N2",
            " 1 2 3 4 5 6+7g",
            " 0 0 0 0 0 0-1u",
            {"v5u": 5, "lu": 6, "ru": 7, "vib_symmetryl": "u"},
        ),
        ("H2CO", "    1 2 3 4 5 6", "    0 0 0 0 0 0", {"v4u": 4, "v6u": 6, "v1l": 0}),
        ("H2O2", "    1 2 334 5 6", "    0 0 000 0 0", {"nu": 3, "tauu": 4, "v5u": 5}),
        (
            "SO3",
            " 1 2 3-1 4 2A1 ",
            " 0 0 0 0 0 0A2 ",
            {"v3u": 3, "l3u": -1, "vib_symmetryl": "A2"},
        ),
        (
            "H3+",
            "      1 2-1 A1 ",
            "      0 0 0 A2 ",
            {"v1u": 1, "l2u": -1, "total_symmetryl": "A2"},
        ),
        (
            "CH4",
            "    1 2 3 401A1",
            "    0 0 0 002A2",
            {"v4u": 4, "multiplicityu": "01", "vib_symmetryl": "A2"},
        ),
        (
            "SF6",
            "          V1+V4",
            "         GROUND",
            {"vib_bandu": "V1+V4", "vib_bandl": "GROUND"},
        ),
    ],
)
def test_hitran_global_class_formats(
    molecule, upper, lower, expected, dataframe_type, fast_parsing
):
    assert len(upper) == len(lower) == 15
    data = pd.DataFrame({"globu": [upper], "globl": [lower]})
    if dataframe_type == "vaex":
        import vaex

        data = vaex.from_pandas(data)
    result = parse_global_quanta(
        data,
        molecule,
        verbose=False,
        dataframe_type=dataframe_type,
        fast_parsing=fast_parsing,
    )
    if dataframe_type == "vaex":
        result = result.to_pandas_df()
    for column, value in expected.items():
        assert result[column][0] == value


@pytest.mark.fast
@pytest.mark.parametrize("dataframe_type", ["pandas", "vaex"])
@pytest.mark.parametrize(
    "molecule,columns", [("H2CO", ("v1u", "v6l")), ("CO2", ("v1u", "rl"))]
)
def test_hitran_new_global_formats_blank_assignments(dataframe_type, molecule, columns):
    data = pd.DataFrame({"globu": [" " * 15], "globl": [" " * 15]})
    if dataframe_type == "vaex":
        import vaex

        data = vaex.from_pandas(data)
    result = parse_global_quanta(data, molecule, dataframe_type=dataframe_type)
    if dataframe_type == "vaex":
        result = result.to_pandas_df()
    for column in columns:
        assert result[column][0] == -1


@pytest.mark.fast
def test_hitran_names_match(verbose=True, warnings=True, *args, **kwargs):
    """Compare that HITRAN species defined in :mod:`radis.api.hitranapi` match
    the nomenclature dictionary : :py:data:`radis.api.hitranapi.trans`.

    This should be ensured by developers when adding new species.
    """
    from radis.db.classes import (
        HITRAN_CLASS1,
        HITRAN_CLASS2,
        HITRAN_CLASS3,
        HITRAN_CLASS4,
        HITRAN_CLASS5,
        HITRAN_CLASS6,
        HITRAN_CLASS7,
        HITRAN_CLASS8,
        HITRAN_CLASS9,
        HITRAN_MOLECULES,
    )
    from radis.misc.basics import compare_lists

    all_hitran = (
        HITRAN_CLASS1
        + HITRAN_CLASS2
        + HITRAN_CLASS3
        + HITRAN_CLASS4
        + HITRAN_CLASS5
        + HITRAN_CLASS6
        + HITRAN_CLASS7
        + HITRAN_CLASS8
        + HITRAN_CLASS9
    )
    all_hitran = list(set(all_hitran))

    # All species in HITRAN groups should be valid HITRAN_MOLECULES names
    for m in all_hitran:
        if m == "CH3D" and "CH4" in HITRAN_MOLECULES:
            pass
        elif not m in HITRAN_MOLECULES:
            raise ValueError(f"{m} is defined in HITRAN groups but has no HITRAN id")

    # Species in 'HITRAN_MOLECULES' should be classified in groups, else nonequilibrium
    # calculations are not possible.
    if warnings and all_hitran != HITRAN_MOLECULES:
        warn(
            "Difference between HITRAN groups (left) and HITRAN id "
            + "dictionary (right). Some HITRAN species are not classified in "
            + "groups. Nonequilibrium calculations wont be possible for these!:\n"
            + f"{compare_lists(all_hitran, HITRAN_MOLECULES, verbose=False, return_string=True)[1]}"
        )


@pytest.mark.fast
def test_local_hitran_co(verbose=True, warnings=True, **kwargs):
    """Analyse some default files to make sure everything still works."""

    # 1. Load
    df = hit2df(getTestFile("hitran_CO_fragment.par"), cache="regen")
    if verbose:
        print("Read hitran_CO_fragment.par")
        print("---------------------------")
        print(df.head())

    # 2. Test
    assert list(df.loc[0, ["vu", "vl"]]) == [4, 4]
    assert df.dtypes["vu"] == np.int64
    assert df.dtypes["vl"] == np.int64


def test_local_hitran_co2(verbose=True, warnings=True, **kwargs):

    # 1. Load
    df = hit2df(getTestFile("hitran_CO2_fragment.par"), cache="regen")
    if verbose:
        print("Read hitran_CO2_fragment.par")
        print("----------------------------")
        print(df.head())

    # 2. Test
    assert list(
        df.loc[0, ["v1u", "v2u", "l2u", "v3u", "v1l", "v2l", "l2l", "v3l"]]
    ) == [4, 0, 0, 0, 0, 0, 0, 1]
    assert df.dtypes["v1l"] == np.int64
    assert df.dtypes["v3u"] == np.int64


def test_local_hitran_h2o(verbose=True, warnings=True, **kwargs):

    # 1. Load
    df = hit2df(getTestFile("hitran_2016_H2O_2iso_2000_2100cm.par"), cache="regen")
    if verbose:
        print("Read hitran_2016_H2O_2iso_2000_2100cm.par")
        print("-----------------------------------------")
        print(df.head())

    # 2. Test
    assert list(df.loc[0, ["v1u", "v2u", "v3u", "v1l", "v2l", "v3l"]]) == [
        0,
        2,
        0,
        0,
        1,
        0,
    ]

    assert df.loc[26, "ju"] == 5  # in .par : line 27, column 99-100
    assert df.loc[27, "ju"] == 18  # in .par : line 28, column 99-100

    assert df.dtypes["v1l"] == np.int64
    assert df.dtypes["v3u"] == np.int64

    assert df.dtypes["ju"] == np.int64
    assert df.dtypes["Kau"] == np.int64
    assert df.dtypes["Kcu"] == np.int64
    assert df.dtypes["jl"] == np.int64
    assert df.dtypes["Kal"] == np.int64
    assert df.dtypes["Kcl"] == np.int64


def test_local_hitemp_file(verbose=True, warnings=True, **kwargs):
    """Analyse some default files to make sure everything still works."""

    # 1. Load
    df = cdsd2df(
        getTestFile("cdsd_hitemp_09_header.txt"), cache="regen", drop_non_numeric=True
    )
    if verbose:
        print(df.head())

    # 2. Tests
    assert df.wav[3] == 2250.00096
    # make sure P Q R is correctly replaced by drop_non_numeric:
    assert "branch" in df
    assert df["branch"].iloc[0] == 0  # Q
    assert df["branch"].iloc[1] == 1  # R


def test_irrelevant_file_loading(*args, **kwargs):
    """check that irrelevant files (irrelevant wavenumber) are not loaded"""

    # For cdsd-hitemp files :

    # Ensures file is not empty
    df = cdsd2df(getTestFile("cdsd_hitemp_09_header.txt"), cache="regen")
    assert len(df) > 0

    # Now test nothing is loaded if asking for outside the wavenumber range
    with pytest.raises(IrrelevantFileWarning):
        df = cdsd2df(getTestFile("cdsd_hitemp_09_header.txt"), load_wavenum_min=100000)
    with pytest.raises(IrrelevantFileWarning):
        df = cdsd2df(getTestFile("cdsd_hitemp_09_header.txt"), load_wavenum_max=0.5)

    # For HITRAN files :

    # Ensures file is not empty
    df = hit2df(getTestFile("hitran_2016_H2O_2iso_2000_2100cm.par"), cache="regen")
    assert len(df) > 0

    # Now test nothing is loaded if asking for outside the wavenumber range
    with pytest.raises(IrrelevantFileWarning):
        df = hit2df(
            getTestFile("hitran_2016_H2O_2iso_2000_2100cm.par"), load_wavenum_min=100000
        )
    with pytest.raises(IrrelevantFileWarning):
        df = hit2df(
            getTestFile("hitran_2016_H2O_2iso_2000_2100cm.par"), load_wavenum_max=0.5
        )


def _run_example(verbose=False):
    from radis import SpectrumFactory

    setup_test_line_databases(
        verbose=verbose
    )  # add HITEMP-CO2-TEST in ~/radis.json if not there
    sf = SpectrumFactory(
        wavelength_min=4165,
        wavelength_max=4200,
        path_length=0.1,
        pressure=20,
        molecule="CO2",
        isotope="1",
        cutoff=1e-25,  # cm/molecule
        truncation=5,  # cm-1
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.load_databank("HITRAN-CO2-TEST")  # this database must be defined in ~/radis.json


def test_cache_regeneration(verbose=True, warnings=True, **kwargs):
    """Checks that if a line database is manually updated (last edited time changes),
    its cache file will be regenerated automatically
    """

    # Initialisation
    # --------------
    # Run case : we generate a first cache file
    _run_example(verbose=verbose)
    # get time when line database file was last modified
    file_last_modification = getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.par")
    )
    # get time when cache file was last modified
    cache_last_modification = getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.h5")
    )

    # To be sure, re run-example and make sure the cache file was not regenerated.
    _run_example(verbose=verbose)
    assert cache_last_modification == getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.h5")
    )

    # Test
    # ----
    # Now we fake the manual editing of the database.
    # change the time when line database was last modified :
    stinfo = os.stat(getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.par"))
    access_time = stinfo.st_atime
    os.utime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.par"),
        (file_last_modification + 1, access_time + 1),
    )

    file_last_modification_again = getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.par")
    )

    assert file_last_modification_again > file_last_modification

    # Run case : this should re-generated the cache file
    _run_example(verbose=verbose)

    cache_last_modification_again = getmtime(
        getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.h5")
    )

    assert cache_last_modification_again > cache_last_modification


def test_hitran_cross_sections(*args, **kwargs):

    datafile = "CH3COCH3_233.4_375.2_700.0-1780.0_13.xsc"

    from radis.api.hitranapi import hitranxsc

    data = hitranxsc(getTestFile(datafile))

    assert data["T"] == 233.4  # K
    assert data["P"] == 375.2  # Pa
    assert data["name"] == "Acetone"
    assert data["wavenumber"].min() == 700.004
    assert data["wavenumber"].max() == 1780.0066


@pytest.mark.fast
@pytest.mark.parametrize(
    "molecule, globu, globl, expected",
    [
        ("CO", "              1", "              0", {"vu": 1, "vl": 0}),
        ("OH", "       X3/2   1", "       X3/2   0", {"vu": 1, "vl": 0}),
        (
            "CO2",
            "       0 0 0 11",
            "       0 0 0 01",
            {"v1u": 0, "v2u": 0, "l2u": 0, "v3u": 1, "ru": 1, "v3l": 0, "rl": 1},
        ),
        (
            "N2O",
            "        0 0 0 1",
            "        0 0 0 0",
            {"v1u": 0, "v2u": 0, "l2u": 0, "v3u": 1, "v3l": 0},
        ),
        ("H2O", "          0 1 0", "          0 0 0", {"v1u": 0, "v2u": 1, "v3u": 0}),
    ],
)
def test_parse_global_quanta_by_class(molecule, globu, globl, expected):
    """Each HITRAN class is read with the parser written for its format"""
    from radis.api.hitranapi import parse_global_quanta

    df = pd.DataFrame({"globu": [globu], "globl": [globl]})
    df = parse_global_quanta(df, molecule, verbose=False)
    for column, value in expected.items():
        assert df[column].iloc[0] == value, column


@pytest.mark.fast
def test_parse_local_quanta_linear():
    from radis.api.hitranapi import parse_local_quanta

    df = pd.DataFrame({"locu": [" " * 15], "locl": ["     R  5      "]})
    df = parse_local_quanta(df, "CO", verbose=False)
    assert df["branch"].iloc[0] == "R"
    assert df["jl"].iloc[0] == 5


@pytest.mark.fast
@pytest.mark.parametrize("dataframe_type", ["pandas", "vaex"])
@pytest.mark.parametrize(
    "molecule,upper,lower,expected",
    [
        (
            "C4H2",
            "-1 2 3 4" + " " * 7,
            " 0-2 1 3" + " R  5e ",
            {"l6u": "-1", "l7l": "-2", "branch": "R", "jl": 5, "syml": "e"},
        ),
        (
            "CH4",
            "  " + " 12A1  7" + "   F1",
            "  " + "  9E   4" + "   F2",
            {"ju": 12, "rot_symmetryu": "A1", "alphal": 4, "Fl": "F2"},
        ),
        (
            "PH3",
            " 12  3-2A1e  F1",
            "  9  2 1A2f  F2",
            {"Ku": 3, "rot_lu": -2, "rot_symmetryl": "A2", "syml": "f"},
        ),
        (
            "NH3",
            f"{12:2}{3:3}{-2:2} {'A1':3}{'E':3} ",
            f"{9:2}{2:3}{1:2} {'A2':3}{'A1':3} ",
            {"ju": 12, "rot_lu": -2, "rot_symmetryu": "A1", "total_symmetryl": "A1"},
        ),
        (
            "C2H6",
            " 12  3-2A1   F1",
            "  9  2 1E    F2",
            {"Ku": 3, "rot_lu": -2, "symu": "A1", "Fl": "F2"},
        ),
        (
            "SO3",
            f"{12:3}{3:3}  {'A1':3}{'E':3} ",
            f"{9:3}{2:3}  {'A2':3}{'A1':3} ",
            {"ju": 12, "rot_symmetryu": "A1", "total_symmetryl": "A1"},
        ),
        (
            "CH3",
            f"{12:3}{3:3}  {'A1':3}{'E':3} ",
            f"{9:3}{2:3}  {'A2':3}{'A1':3} ",
            {"Ku": 3, "rot_symmetryl": "A2"},
        ),
        (
            "H3+",
            f"{12:3}{'A1':2}{3:3}{'u':2}{'K1':3}  ",
            f"{9:3}{'E':2}{2:3}{'l':2}{'K2':3}  ",
            {"ju": 12, "Gu": 3, "Uu": "u", "Kl": "K2", "syml": "E"},
        ),
        (
            "O2",
            " " * 10 + "   F1",
            " P  3R  4   F2d",
            {
                "Fu": "F1",
                "branch_N": "P",
                "Nl": 3,
                "branch_J": "R",
                "jl": 4,
                "transition_type": "d",
            },
        ),
        (
            "NO",
            "m" + " " * 9 + "   F1",
            "  PQ  3.5e   F2",
            {"magnetic": "m", "branch_NJ": "PQ", "jl": 3.5, "syml": "e", "Fl": "F2"},
        ),
        (
            "OH",
            " " * 10 + "   F1",
            " QR  3.5ef   F2",
            {"branch_NJ": "QR", "jl": 3.5, "symmetry_pair": "ef", "Fl": "F2"},
        ),
    ],
)
def test_hitran_local_group_formats(molecule, upper, lower, expected, dataframe_type):
    from radis.api.hitranapi import parse_local_quanta

    assert len(upper) == len(lower) == 15
    data = pd.DataFrame({"locu": [upper, " " * 15], "locl": [lower, " " * 15]})
    if dataframe_type == "vaex":
        import vaex

        data = vaex.from_pandas(data)
    result = parse_local_quanta(data, molecule, dataframe_type=dataframe_type)
    if dataframe_type == "vaex":
        result = result.to_pandas_df()
    for column, value in expected.items():
        assert result[column][0] == value, column
        if isinstance(value, (int, float)):
            assert result[column][1] == -1, column
    assert result["locu"][0] == upper
    assert result["locl"][0] == lower
    if molecule in ("O2", "NO", "OH"):
        assert "branch" not in result


@pytest.mark.fast
@pytest.mark.parametrize("dataframe_type", ["pandas", "vaex"])
@pytest.mark.parametrize("fast_parsing", [False, True])
def test_hitran_group2_magnetic_field(dataframe_type, fast_parsing):
    from radis.api.hitranapi import parse_local_quanta

    data = pd.DataFrame(
        {"locu": ["m" + " " * 9 + "   F1"], "locl": ["     R  5e   F2"]}
    )
    if dataframe_type == "vaex":
        import vaex

        data = vaex.from_pandas(data)
    result = parse_local_quanta(
        data, "CO2", dataframe_type=dataframe_type, fast_parsing=fast_parsing
    )
    if dataframe_type == "vaex":
        result = result.to_pandas_df()
    assert result["magnetic"][0] == "m"
    assert result["branch"][0] == "R"
    assert result["jl"][0] == 5


@pytest.mark.fast
def test_hitran_group1_vaex_lower_labels():
    import vaex

    from radis.api.hitranapi import parse_local_quanta

    data = vaex.from_pandas(
        pd.DataFrame({"locu": ["  1  2  3   F1e"], "locl": ["  4  5  6   F2f"]})
    )
    result = parse_local_quanta(data, "H2O", dataframe_type="vaex").to_pandas_df()
    assert result["Fl"][0] == "F2"
    assert result["syml"][0] == "f"


@pytest.mark.fast
@pytest.mark.parametrize("dataframe_type", ["pandas", "vaex"])
@pytest.mark.parametrize("fast_parsing", [False, True])
@pytest.mark.parametrize(
    "parser,molecule,upper,lower,expected",
    [
        ("global", "CO", " " * 13 + " 4", " " * 13 + " 2", {"vu": 4, "vl": 2}),
        (
            "global",
            "H2O",
            " " * 9 + " 1-2 3",
            " " * 9 + " 0 1 2",
            {"v1u": 1, "v2u": -2, "v3l": 2},
        ),
        (
            "global",
            "HCN",
            " " * 7 + " 1-2 3 4",
            " " * 7 + " 0 1 2 3",
            {"v1u": 1, "v2u": -2, "l2l": 2},
        ),
        (
            "local",
            "H2O",
            "  1 -2  3   F1e",
            "  4  5 -1   F2f",
            {"ju": 1, "Kau": -2, "Kcl": -1},
        ),
        (
            "local",
            "CO2",
            "m" + " " * 9 + "   F1",
            "     R  5e   F2",
            {"jl": 5},
        ),
    ],
)
def test_hitran_legacy_parsers_preserve_unassigned(
    parser, molecule, upper, lower, expected, dataframe_type, fast_parsing
):
    from radis.api.hitranapi import parse_local_quanta

    fields = ("globu", "globl") if parser == "global" else ("locu", "locl")
    data = pd.DataFrame({fields[0]: [upper, " " * 15], fields[1]: [lower, " " * 15]})
    if dataframe_type == "vaex":
        import vaex

        data = vaex.from_pandas(data)
    parse = parse_global_quanta if parser == "global" else parse_local_quanta
    result = parse(
        data, molecule, dataframe_type=dataframe_type, fast_parsing=fast_parsing
    )
    if dataframe_type == "vaex":
        result = result.to_pandas_df()
    for column, value in expected.items():
        assert result[column][0] == value
        assert result[column][1] == -1
    if molecule == "CO2":
        assert result["branch"][0] == "R"
        assert result["branch"][1] == ""


@pytest.mark.fast
@pytest.mark.parametrize("dataframe_type", ["pandas", "vaex"])
def test_hitran_group2_postprocess_unassigned_branch(dataframe_type):
    from radis.api.hitranapi import post_process_hitran_data

    data = pd.DataFrame(
        {
            "id": [2, 2],
            "locu": ["m" + " " * 9 + "   F1", " " * 15],
            "locl": ["     R  5e   F2", " " * 15],
            "globu": ["       0 0 0 11", " " * 15],
            "globl": ["       0 0 0 01", " " * 15],
        }
    )
    if dataframe_type == "vaex":
        import vaex

        data = vaex.from_pandas(data)
    with pytest.warns(UserWarning, match="Invalid values"):
        result = post_process_hitran_data(
            data, "CO2", verbose=False, dataframe_type=dataframe_type
        )
    if dataframe_type == "vaex":
        result = result.to_pandas_df()
    assert result["branch"][0] == 1
    assert np.isnan(result["branch"][1])
    assert result["jl"][1] == -1


@pytest.mark.fast
@pytest.mark.parametrize("dataframe_type", ["pandas", "vaex"])
def test_hitran_c4h2_postprocess_unassigned_branch(dataframe_type):
    from radis.api.hitranapi import post_process_hitran_data

    data = pd.DataFrame(
        {
            "id": [43, 43],
            "locu": ["-1 2 3 4" + " " * 7, " " * 15],
            "locl": [" 0-2 1 3" + " R  5e ", " " * 15],
            "globu": [" 123456789 e A1", " " * 15],
            "globl": [" 000000000 f A2", " " * 15],
        }
    )
    if dataframe_type == "vaex":
        import vaex

        data = vaex.from_pandas(data)
    with pytest.warns(UserWarning, match="Invalid values"):
        result = post_process_hitran_data(
            data, "C4H2", verbose=False, dataframe_type=dataframe_type
        )
    if dataframe_type == "vaex":
        result = result.to_pandas_df()
    assert result["branch"][0] == 1
    assert np.isnan(result["branch"][1])
    assert result["v9u"][0] == 9
    assert result["jl"][1] == -1


@pytest.mark.fast
def test_hitran_local_and_global_symmetry_are_distinct():
    from radis.api.hitranapi import post_process_hitran_data

    data = pd.DataFrame(
        {
            "id": [56],
            "locu": [f"{12:3}{'E':2}{3:3}{'u':2}{'K1':3}  "],
            "locl": [f"{9:3}{'E':2}{2:3}{'l':2}{'K2':3}  "],
            "globu": ["      1 2-1 A1 "],
            "globl": ["      0 0 0 A2 "],
        }
    )
    result = post_process_hitran_data(
        data, "H3+", verbose=False, drop_non_numeric=False
    )
    assert result["symu"][0] == "E"
    assert result["total_symmetryu"][0] == "A1"
    assert result["Kl"][0] == "K2"


@pytest.mark.fast
def test_COFCl_is_an_asymmetric_rotor():
    from radis.db.classes import HITRAN_GROUP1, HITRAN_GROUP2

    assert "COFCl" in HITRAN_GROUP1
    assert "COFCl" not in HITRAN_GROUP2


def _run_testcases(verbose=True, *args, **kwargs):

    test_hitran_names_match(verbose=verbose, *args, **kwargs)
    # test_local_hitran_co(verbose=verbose, *args, **kwargs)
    # test_local_hitran_co2(verbose=verbose, *args, **kwargs)
    # test_local_hitran_h2o(verbose=verbose, *args, **kwargs)
    # test_local_hitemp_file(verbose=verbose, *args, **kwargs)
    # test_irrelevant_file_loading()
    # test_cache_regeneration(verbose=verbose, *args, **kwargs)
    # test_hitran_cross_sections(*args, **kwargs)


if __name__ == "__main__":
    print("Testing test_hitran_cdsd.py: ", _run_testcases(verbose=True))
    test_cache_regeneration(verbose=3)
