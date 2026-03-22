#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for :py:func:`radis.io.chianti.load_chianti`.

Uses the sample ``fe_12.wgfa`` file shipped with RADIS under
``radis/db/chianti/fe/fe_12/``.

Note on units
-------------
``wmin`` / ``wmax`` are in **nanometres (nm)**, consistent with the
RADIS-wide convention used by ``fetch_nist``, ``fetch_astroquery``, etc.
The Fe XII sample data covers ≈ 19.3–19.7 nm  (193–197 Å).
"""

import os

import numpy as np
import pandas as pd
import pytest

import radis
from radis.io.chianti import load_chianti

# Path to the sample CHIANTI database bundled with RADIS
SAMPLE_DATAPATH = os.path.join(
    os.path.dirname(radis.__file__), "db", "chianti"
)

# Convenience nm window that covers the full sample file
WMIN_NM = 19.0   # nm  ≈ 190 Å
WMAX_NM = 20.0   # nm  ≈ 200 Å


# --------------------------------------------------------------------------- #
#  Basic loading
# --------------------------------------------------------------------------- #


def test_load_chianti_returns_dataframe():
    """load_chianti should return a pandas DataFrame."""
    df = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=WMIN_NM,
        wmax=WMAX_NM,
    )
    assert isinstance(df, pd.DataFrame)


def test_load_chianti_nonempty():
    """Sample fe_12.wgfa contains transitions in 19–20 nm; loading that
    range must return at least one row."""
    df = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=WMIN_NM,
        wmax=WMAX_NM,
    )
    assert len(df) > 0, (
        "Expected non-empty DataFrame for Fe XII in 19–20 nm range"
    )


def test_load_chianti_required_columns():
    """The returned DataFrame must contain the documented columns."""
    df = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=WMIN_NM,
        wmax=WMAX_NM,
    )
    for col in ("wav", "gf", "A", "upper", "lower"):
        assert col in df.columns, f"Missing expected column: {col}"


def test_load_chianti_column_types():
    """Numeric columns should be float; index columns should be int."""
    df = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=WMIN_NM,
        wmax=WMAX_NM,
    )
    assert np.issubdtype(df["wav"].dtype, np.floating)
    assert np.issubdtype(df["gf"].dtype, np.floating)
    assert np.issubdtype(df["A"].dtype, np.floating)
    assert np.issubdtype(df["upper"].dtype, np.integer)
    assert np.issubdtype(df["lower"].dtype, np.integer)


# --------------------------------------------------------------------------- #
#  Wavelength unit — most critical test
# --------------------------------------------------------------------------- #


def test_load_chianti_wavelength_unit_is_nm():
    """
    ``wav`` column must be in nm, consistent with RADIS convention.

    The Fe XII sample file covers ≈ 19.3–19.7 nm (193–197 Å).
    If wav values are > 100 it almost certainly means Å leaked through
    without the × 0.1 conversion.
    """
    df = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=WMIN_NM,
        wmax=WMAX_NM,
    )
    assert df["wav"].max() < 100.0, (
        f"wav max = {df['wav'].max():.1f} — looks like Å, expected nm"
    )
    assert df["wav"].min() >= 1.0, (
        f"wav min = {df['wav'].min():.3f} — suspiciously small, check units"
    )
    # Tighter bound specific to the Fe XII sample
    assert df["wav"].min() >= 19.0, "Fe XII wav should start ≥ 19 nm"
    assert df["wav"].max() <= 20.0, "Fe XII wav should end ≤ 20 nm"


# --------------------------------------------------------------------------- #
#  Wavelength filtering
# --------------------------------------------------------------------------- #


def test_load_chianti_empty_range():
    """A wavelength window outside the data range should return an empty
    DataFrame — no error, just zero rows."""
    df = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=0.0,
        wmax=1.0,   # well below Fe XII EUV lines (~19 nm)
    )
    assert len(df) == 0, "Expected 0 rows for window 0–1 nm"


def test_load_chianti_wmin_wmax_filtering_in_nm():
    """
    wmin/wmax must be interpreted as nm.
    Narrowing the window should return fewer (or equal) lines.
    """
    df_wide = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=19.0,   # nm
        wmax=20.0,   # nm
    )
    df_narrow = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=19.3,   # nm  — tighter lower bound
        wmax=19.5,   # nm  — tighter upper bound
    )
    assert len(df_narrow) <= len(df_wide), (
        "Narrower window should not return more lines than wider window"
    )
    # Every line in narrow must also satisfy the filter bounds
    if len(df_narrow) > 0:
        assert df_narrow["wav"].min() >= 19.3 - 1e-9
        assert df_narrow["wav"].max() <= 19.5 + 1e-9


# --------------------------------------------------------------------------- #
#  Physical sanity checks
# --------------------------------------------------------------------------- #


def test_load_chianti_positive_A():
    """All Einstein A coefficients must be strictly positive."""
    df = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=WMIN_NM,
        wmax=WMAX_NM,
    )
    assert (df["A"] > 0).all(), "Einstein A must be > 0 for all transitions"


def test_load_chianti_no_nan():
    """The returned DataFrame should not contain NaN in any column."""
    df = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=WMIN_NM,
        wmax=WMAX_NM,
    )
    assert not df.isna().any().any(), "DataFrame contains unexpected NaN values"


def test_load_chianti_no_sentinel_rows():
    """
    CHIANTI .wgfa files end with a sentinel line whose level indices are -1.
    This sentinel must NOT appear in the output DataFrame.
    """
    df = load_chianti(
        datapath=SAMPLE_DATAPATH,
        ion="fe_12",
        wmin=WMIN_NM,
        wmax=WMAX_NM,
    )
    assert (df["lower"] != -1).all(), (
        "Sentinel line (lower == -1) leaked into output DataFrame"
    )


# --------------------------------------------------------------------------- #
#  Error handling
# --------------------------------------------------------------------------- #


def test_load_chianti_missing_file():
    """Requesting a non-existent path should raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_chianti(
            datapath="/nonexistent/path",
            ion="fe_12",
        )


def test_load_chianti_missing_ion():
    """Requesting an ion not present in the sample database should raise
    FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_chianti(
            datapath=SAMPLE_DATAPATH,
            ion="xx_99",   # does not exist
        )


# --------------------------------------------------------------------------- #
#  Default datapath
# --------------------------------------------------------------------------- #


def test_load_chianti_default_datapath():
    """When datapath=None the function should fall back to the bundled
    sample database and still return data."""
    df = load_chianti(
        datapath=None,
        ion="fe_12",
        wmin=WMIN_NM,
        wmax=WMAX_NM,
    )
    assert len(df) > 0, "Default datapath fallback returned empty DataFrame"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
