# -*- coding: utf-8 -*-
"""
Unit tests for radis.api.hdf5 DataFileManager.
Tests write/load/read/metadata roundtrips using temp files — no network, no GPU.
"""
import os

import numpy as np
import pandas as pd
import pytest


@pytest.mark.fast
class TestVaexsafeColname:
    """Tests for vaexsafe_colname."""

    def test_replaces_slash(self):
        from radis.api.hdf5 import vaexsafe_colname

        assert vaexsafe_colname("a/b") == "a_b"

    def test_no_slash(self):
        from radis.api.hdf5 import vaexsafe_colname

        assert vaexsafe_colname("abc") == "abc"


@pytest.mark.fast
class TestDataFileManagerInit:
    """Tests for DataFileManager initialization."""

    def test_pytables_engine(self):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("pytables")
        assert mgr.engine == "pytables"

    def test_vaex_engine(self):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("vaex")
        assert mgr.engine == "vaex"

    def test_h5py_engine(self):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("h5py")
        assert mgr.engine == "h5py"


@pytest.mark.fast
class TestDataFileManagerCacheFile:
    """Tests for DataFileManager.cache_file."""

    def test_pytables_cache_extension(self):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("pytables")
        result = str(mgr.cache_file("/data/test.par"))
        assert result.endswith(".h5")

    def test_vaex_cache_extension(self):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("vaex")
        result = str(mgr.cache_file("/data/test.par"))
        assert result.endswith(".hdf5")


@pytest.mark.fast
class TestDataFileManagerWriteRead:
    """Tests for DataFileManager write/read roundtrip."""

    def test_pytables_write_read_roundtrip(self, tmp_path):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("pytables")
        df = pd.DataFrame(
            {
                "wav": [100.0, 200.0, 300.0],
                "int": [1e-20, 2e-20, 3e-20],
                "iso": [1, 1, 2],
            }
        )

        filepath = str(tmp_path / "test.h5")
        mgr.write(filepath, df, key="default", append=False)

        assert os.path.exists(filepath)

        loaded = mgr.read(filepath, key="default")
        assert loaded is not None
        assert len(loaded) == 3
        assert "wav" in loaded.columns
        np.testing.assert_allclose(loaded["wav"].values, [100.0, 200.0, 300.0])

    def test_pytables_write_read_with_columns(self, tmp_path):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("pytables")
        df = pd.DataFrame(
            {
                "wav": [100.0, 200.0],
                "int": [1e-20, 2e-20],
                "iso": [1, 1],
            }
        )

        filepath = str(tmp_path / "test_cols.h5")
        mgr.write(filepath, df, key="default", append=False)

        loaded = mgr.read(filepath, columns=["wav"], key="default")
        assert "wav" in loaded.columns


@pytest.mark.fast
class TestDataFileManagerMetadata:
    """Tests for DataFileManager metadata operations."""

    def test_add_and_read_metadata(self, tmp_path):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("pytables")
        df = pd.DataFrame({"wav": [100.0], "iso": [1]})

        filepath = str(tmp_path / "test_meta.h5")
        mgr.write(filepath, df, key="default", append=False)

        # Add metadata
        metadata = {"version": "1.0", "molecule": "CO"}
        mgr.add_metadata(filepath, metadata, key="default")

        # Read metadata
        result = mgr.read_metadata(filepath, key="default")
        assert result["version"] == "1.0"
        assert result["molecule"] == "CO"


@pytest.mark.fast
class TestDataFileManagerGetColumns:
    """Tests for DataFileManager.get_columns."""

    def test_get_columns(self, tmp_path):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("pytables")
        df = pd.DataFrame({"wav": [1.0], "int": [2.0], "iso": [1]})

        filepath = str(tmp_path / "test_getcols.h5")
        mgr.write(filepath, df, key="default", append=False)

        columns = mgr.get_columns(filepath)
        assert "wav" in columns
        assert "int" in columns
        assert "iso" in columns


@pytest.mark.fast
class TestDataFileManagerLoad:
    """Tests for DataFileManager.load with filtering."""

    def test_load_with_bounds(self, tmp_path):
        from radis.api.hdf5 import DataFileManager

        mgr = DataFileManager("pytables")
        df = pd.DataFrame(
            {
                "wav": [100.0, 200.0, 300.0, 400.0],
                "iso": [1, 1, 2, 2],
            }
        )

        filepath = str(tmp_path / "test_filter.h5")
        mgr.write(filepath, df, key="default", append=False)

        # Load with lower bound
        loaded = mgr.load(
            filepath,
            lower_bound=[("wav", 200.0)],
            output="pandas",
        )
        assert loaded is not None
        assert all(loaded["wav"] >= 200.0)


@pytest.mark.fast
class TestHDF5ManagerDeprecation:
    """Tests for HDF5Manager deprecation warning."""

    def test_deprecated(self):
        from radis.api.hdf5 import HDF5Manager

        with pytest.raises(DeprecationWarning):
            HDF5Manager()


if __name__ == "__main__":
    pytest.main([__file__])
