# -*- coding: utf-8 -*-
"""
Unit tests for radis.api.hitranapi HITRAN class parsers.
Tests pure parsing functions using synthetic DataFrames — no network, no GPU.
"""
import numpy as np
import pandas as pd
import pytest


@pytest.mark.fast
class TestCastToInt64:
    """Tests for cast_to_int64_with_missing_values."""

    def test_basic_conversion(self):
        from radis.api.hitranapi import cast_to_int64_with_missing_values

        df = pd.DataFrame({"v": [" 1", " 2", " 3"]})
        cast_to_int64_with_missing_values(df, ["v"])
        assert df["v"].dtype == np.int64
        assert list(df["v"]) == [1, 2, 3]

    def test_missing_values_replaced_with_minus_one(self):
        from radis.api.hitranapi import cast_to_int64_with_missing_values

        df = pd.DataFrame({"v": [" 1", "  ", " 3"]})
        cast_to_int64_with_missing_values(df, ["v"])
        assert df["v"].iloc[1] == -1


@pytest.mark.fast
class TestCastAllToInt64:
    """Tests for cast_all_to_int64."""

    def test_basic(self):
        from radis.api.hitranapi import cast_all_to_int64

        df = pd.DataFrame({"a": ["1", "2"], "b": ["3", "4"]})
        cast_all_to_int64(df, ["a", "b"])
        assert df["a"].dtype == np.int64
        assert df["b"].dtype == np.int64


@pytest.mark.fast
class TestParseHITRANClass1FastParsing:
    """Tests for _parse_HITRAN_class1_fast_parsing — diatomic molecules."""

    def test_pandas(self):
        from radis.api.hitranapi import _parse_HITRAN_class1_fast_parsing

        # globu/globl format: 13 spaces + 2-digit vibrational number
        df = pd.DataFrame(
            {
                "globu": ["             01"],
                "globl": ["              2"],
            }
        )
        result = _parse_HITRAN_class1_fast_parsing(df, verbose=False)
        assert "vu" in result.columns
        assert "vl" in result.columns
        assert "globu" not in result.columns
        assert result["vu"].iloc[0] == 1
        assert result["vl"].iloc[0] == 2


@pytest.mark.fast
class TestParseHITRANClass1:
    """Tests for _parse_HITRAN_class1 — diatomic molecules (regex version)."""

    def test_pandas(self):
        from radis.api.hitranapi import _parse_HITRAN_class1

        df = pd.DataFrame(
            {
                "globu": ["             01"],
                "globl": ["              2"],
            }
        )
        result = _parse_HITRAN_class1(df, verbose=False)
        assert "vu" in result.columns
        assert "vl" in result.columns
        assert result["vu"].iloc[0] == 1
        assert result["vl"].iloc[0] == 2


@pytest.mark.fast
class TestParseHITRANClass2:
    """Tests for _parse_HITRAN_class2 — returns df unchanged (not implemented)."""

    def test_returns_unchanged(self):
        from radis.api.hitranapi import _parse_HITRAN_class2

        df = pd.DataFrame({"wav": [100.0]})
        result = _parse_HITRAN_class2(df, verbose=False)
        assert "wav" in result.columns


@pytest.mark.fast
class TestParseHITRANClass3:
    """Tests for _parse_HITRAN_class3 — returns df unchanged (not implemented)."""

    def test_returns_unchanged(self):
        from radis.api.hitranapi import _parse_HITRAN_class3

        df = pd.DataFrame({"wav": [100.0]})
        result = _parse_HITRAN_class3(df, verbose=False)
        assert "wav" in result.columns


@pytest.mark.fast
class TestParseHITRANClass4:
    """Tests for _parse_HITRAN_class4 — linear triatomic: N2O, OCS, HCN."""

    def test_pandas(self):
        from radis.api.hitranapi import _parse_HITRAN_class4

        # globu/globl format: 7 spaces + 4 x I2
        df = pd.DataFrame(
            {
                "globu": ["        1 0 0 1"],
                "globl": ["        0 1 0 0"],
            }
        )
        result = _parse_HITRAN_class4(df, verbose=False)
        assert "v1u" in result.columns
        assert "v3l" in result.columns


@pytest.mark.fast
class TestParseHITRANClass5:
    """Tests for _parse_HITRAN_class5 — CO2."""

    def test_pandas(self):
        from radis.api.hitranapi import _parse_HITRAN_class5

        # CO2 format: 6 spaces + v1(I2) + v2(I2) + l2(I2) + v3(I2) + r(I1)
        df = pd.DataFrame(
            {
                "globu": ["       1 0 0 11"],
                "globl": ["       0 1 0 01"],
            }
        )
        result = _parse_HITRAN_class5(df, verbose=False)
        assert "v1u" in result.columns
        assert "ru" in result.columns
        assert "rl" in result.columns


@pytest.mark.fast
class TestParseHITRANClass6FastParsing:
    """Tests for _parse_HITRAN_class6_fast_parsing — non-linear triatomic (H2O etc.)."""

    def test_pandas(self):
        from radis.api.hitranapi import _parse_HITRAN_class6_fast_parsing

        # Format: 9 spaces + v1(I2) + v2(I2) + v3(I2)
        df = pd.DataFrame(
            {
                "globu": ["          0 1 0"],
                "globl": ["          0 0 0"],
            }
        )
        result = _parse_HITRAN_class6_fast_parsing(df, verbose=False)
        assert "v1u" in result.columns
        assert "v3l" in result.columns
        assert "globu" not in result.columns


@pytest.mark.fast
class TestParseHITRANClass6:
    """Tests for _parse_HITRAN_class6 — non-linear triatomic (regex)."""

    def test_pandas(self):
        from radis.api.hitranapi import _parse_HITRAN_class6

        df = pd.DataFrame(
            {
                "globu": ["          0 1 0"],
                "globl": ["          0 0 0"],
            }
        )
        result = _parse_HITRAN_class6(df, verbose=False)
        assert "v1u" in result.columns
        assert "v3l" in result.columns


@pytest.mark.fast
class TestPostProcessHitranData:
    """Tests for post_process_hitran_data basic paths."""

    def test_empty_databank_raises(self):
        from radis.api.hitranapi import post_process_hitran_data

        df = pd.DataFrame({"id": []})
        with pytest.raises(ValueError, match="empty"):
            post_process_hitran_data(df, molecule="CO", parse_quanta=False)

    def test_multiple_molecules_raises(self):
        from radis.api.hitranapi import post_process_hitran_data

        df = pd.DataFrame({"id": [1, 2], "wav": [100, 200]})
        with pytest.raises(ValueError, match="Multiple molecules"):
            post_process_hitran_data(df, molecule="CO", parse_quanta=False)


if __name__ == "__main__":
    pytest.main([__file__])
