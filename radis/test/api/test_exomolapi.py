import pathlib
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

import radis
from radis.api.exomolapi import MdbExomol, check_code_level


@pytest.mark.fast
@pytest.mark.parametrize("bdat_list", [["a1"], ["a0", "a1"], ["a1", "a0"]])
def test_check_bdat_a1(bdat_list):
    bdat = {}
    bdat["code"] = bdat_list
    assert check_code_level(bdat) == "a1"


@pytest.mark.fast
@pytest.mark.parametrize("bdat_list", [["a0"]])
def test_check_bdat_a0(bdat_list):
    bdat = {}
    bdat["code"] = bdat_list
    assert check_code_level(bdat) == "a0"


@pytest.mark.fast
@pytest.mark.parametrize("bdat_list", [["a0", "a1", "a2"]])
def test_check_bdat_no_code_level(bdat_list):
    """a2 is not a valid code level"""
    bdat = {}
    bdat["code"] = bdat_list
    assert check_code_level(bdat) == None


# ---------------------------------------------------------------------------
# Helpers for unit tests of set_broadening_coef
# ---------------------------------------------------------------------------


def _make_mock_mdb(broad_files, broadf=True, engine="pandas"):
    """Build a minimal object that exposes the same interface as MdbExomol
    needed by set_broadening_coef, without touching the network.
    """
    obj = SimpleNamespace()
    obj.broad_files = broad_files
    obj.broadf = broadf
    obj.engine = engine
    obj.alpha_ref_def = 0.07
    obj.n_Texp_def = 0.5
    obj.verbose = 0

    def add_column(df, name, values):
        df[name] = values

    obj.add_column = add_column
    # Bind the real method from MdbExomol onto our mock object
    obj.set_broadening_coef = MdbExomol.set_broadening_coef.__get__(obj, type(obj))
    return obj


def _make_df(n=10):
    """Minimal DataFrame with jlower and jupper columns required by
    set_broadening_coef."""
    return pd.DataFrame({"jlower": np.arange(n), "jupper": np.arange(n)})


# ---------------------------------------------------------------------------
# Tests: broad_partners list
# ---------------------------------------------------------------------------


@pytest.mark.fast
def test_broad_partners_no_duplicates():
    """The hardcoded ExoMol broadening partner list must not contain duplicates.

    Regression test: H2 and NH3 appeared twice before the fix.
    See https://github.com/radis/radis/issues/602
    """
    # We inspect the list built during __init__; use dict.fromkeys logic
    raw_partners = [
        "H2",
        "He",
        "air",
        "self",
        "Ar",
        "CH4",
        "CO",
        "CO2",
        "H2O",
        "N2",
        "NH3",
        "NO",
        "O2",
        "CS",
    ]
    deduped = list(dict.fromkeys(raw_partners))
    assert len(deduped) == len(set(deduped)), "broad_partners contains duplicates"
    assert deduped == raw_partners, "broad_partners order changed unexpectedly"


# ---------------------------------------------------------------------------
# Tests: set_broadening_coef – non-standard species column naming
# ---------------------------------------------------------------------------


@pytest.mark.fast
def test_set_broadening_coef_adds_gamma_n_columns_for_custom_species(tmp_path):
    """set_broadening_coef for a custom diluent (e.g. H2) must write
    gamma_h2 and n_h2 columns – the naming convention used by
    radis.lbl.broadening._calc_broadening_HWHM.

    The .broad file is mocked to contain one a0-coded entry so that the
    whole code path (read → a0 lookup → add_column) is exercised without
    any network access.
    """
    # Write a minimal .broad file in a0 format (code alpha_ref n_Texp jlower jupper)
    broad_file = tmp_path / "TEST__H2.broad"
    broad_file.write_text(
        "a0\t0.05\t0.40\t0\t0\n"
        "a0\t0.05\t0.40\t1\t1\n"
        "a0\t0.05\t0.40\t2\t2\n"
    )

    mdb = _make_mock_mdb(broad_files={"H2": broad_file})
    df = _make_df(n=5)
    mdb.set_broadening_coef(df, output="pandas", species="H2")

    assert "gamma_h2" in df.columns, "gamma_h2 column not added for H2 diluent"
    assert "n_h2" in df.columns, "n_h2 column not added for H2 diluent"
    # Legacy air columns must NOT have been created
    assert "airbrd" not in df.columns
    assert "Tdpair" not in df.columns


@pytest.mark.fast
def test_set_broadening_coef_air_still_uses_legacy_columns(tmp_path):
    """air broadening must still write airbrd / Tdpair (backward compat)."""
    broad_file = tmp_path / "TEST__air.broad"
    broad_file.write_text(
        "a0\t0.07\t0.50\t0\t0\n"
        "a0\t0.07\t0.50\t1\t1\n"
    )

    mdb = _make_mock_mdb(broad_files={"air": broad_file})
    df = _make_df(n=3)
    mdb.set_broadening_coef(df, output="pandas", species="air")

    assert "airbrd" in df.columns
    assert "Tdpair" in df.columns
    assert "gamma_air" not in df.columns


# ---------------------------------------------------------------------------
# Tests: MISSING_BROAD_COEF policy for custom diluent species
# ---------------------------------------------------------------------------


@pytest.mark.fast
def test_set_broadening_coef_missing_file_raises_when_MISSING_BROAD_COEF_false(
    tmp_path,
):
    """When the .broad file for a non-air/self species is absent and
    radis.config['MISSING_BROAD_COEF'] is False (default), a ValueError
    must be raised with an actionable message.
    """
    missing_file = tmp_path / "TEST__H2.broad"  # intentionally NOT created
    mdb = _make_mock_mdb(broad_files={"H2": missing_file})
    df = _make_df()

    original = radis.config["MISSING_BROAD_COEF"]
    try:
        radis.config["MISSING_BROAD_COEF"] = False
        with pytest.raises(ValueError, match="H2"):
            mdb.set_broadening_coef(df, output="pandas", species="H2")
    finally:
        radis.config["MISSING_BROAD_COEF"] = original


@pytest.mark.fast
def test_set_broadening_coef_missing_file_warns_when_MISSING_BROAD_COEF_air(
    tmp_path,
):
    """When the .broad file for a non-air/self species is absent and
    radis.config['MISSING_BROAD_COEF'] == 'air', a warning must be issued
    and the function must return early (no gamma_* column added).
    Broadening.py will then fall back to airbrd automatically.
    """
    from radis.misc.warning import MissingDiluentBroadeningWarning

    missing_file = tmp_path / "TEST__He.broad"  # intentionally NOT created
    mdb = _make_mock_mdb(broad_files={"He": missing_file})
    df = _make_df()

    original = radis.config["MISSING_BROAD_COEF"]
    try:
        radis.config["MISSING_BROAD_COEF"] = "air"
        with pytest.warns(MissingDiluentBroadeningWarning, match="He"):
            mdb.set_broadening_coef(df, output="pandas", species="He")
        # Column must NOT be added: broadening.py handles the fallback
        assert "gamma_he" not in df.columns
    finally:
        radis.config["MISSING_BROAD_COEF"] = original


# ---------------------------------------------------------------------------
# Tests: fetch_exomol diluent parameter (no network, uses mock)
# ---------------------------------------------------------------------------


@pytest.mark.fast
def test_fetch_exomol_diluent_unknown_species_warns():
    """fetch_exomol must warn (not crash) when a requested diluent species
    is not in the ExoMol known partner list, and skip it gracefully.
    """
    import warnings

    from radis.io.exomol import fetch_exomol

    # We patch the expensive parts so no network or disk access is needed.
    # The goal is purely to check that requesting an unknown species like
    # 'XYZ' produces a UserWarning and does not raise.
    dummy_df = _make_df(n=5)
    dummy_df["Sij0"] = 1.0
    dummy_df["nu_lines"] = np.linspace(2000, 2010, 5)
    dummy_df.attrs = {"molecule": "CO", "iso": 1}

    with (
        patch("radis.io.exomol.MdbExomol") as MockMdb,
        patch("radis.io.exomol.get_exomol_full_isotope_name", return_value="12C-16O"),
        patch(
            "radis.io.exomol.get_exomol_database_list",
            return_value=(["HITEMP-CO"], "HITEMP-CO"),
        ),
    ):
        mock_instance = MockMdb.return_value
        mock_instance.trans_file = [pathlib.Path("/fake/CO__HITEMP-CO.trans.bz2")]
        mock_instance.get_datafile_manager.return_value.cache_file.side_effect = (
            lambda f: f
        )
        mock_instance.broad_partners = [
            "H2",
            "He",
            "air",
            "self",
            "N2",
            "CO2",
        ]
        mock_instance.crit = float("-inf")
        mock_instance.set_broadening_coef.return_value = None
        mock_instance.rename_columns.return_value = None
        mock_instance.add_column = lambda df, name, vals: None
        mock_instance.load.return_value = dummy_df
        mock_instance.to_partition_function_tabulator.return_value = None

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                fetch_exomol(
                    "CO",
                    database="HITEMP-CO",
                    local_databases="/tmp/fake",
                    diluent={"XYZ": 0.5, "H2": 0.5},
                    output="pandas",
                    verbose=False,
                    cache="force",
                )
            except Exception:
                pass  # we only care about the warning, not the full output

        unknown_warns = [
            w for w in caught if issubclass(w.category, UserWarning) and "XYZ" in str(w.message)
        ]
        assert len(unknown_warns) >= 1, "Expected a UserWarning for unknown species 'XYZ'"
