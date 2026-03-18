"""Tests for DataFrameAdapter and PolarsAdapter — ref issue #658."""

import os
import tempfile
import numpy as np
import pandas as pd
import pytest

try:
    import polars as pl

    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False

from radis.io.dataframe_adapter import DataFrameAdapter


# ---- fixtures ----


@pytest.fixture
def sample_data():
    """Generate synthetic HITRAN-like spectroscopic data."""
    rng = np.random.default_rng(42)
    n = 10000
    return {
        "wav": np.sort(rng.uniform(1000, 5000, n)),
        "int": rng.lognormal(-20, 5, n),
        "A": rng.lognormal(1, 2, n),
        "airbrd": rng.uniform(0.01, 0.15, n),
        "selbrd": rng.uniform(0.01, 0.15, n),
        "El": rng.uniform(0, 5000, n),
        "Tdpair": rng.uniform(0.4, 0.8, n),
        "iso": rng.integers(1, 4, n).astype(float),
    }


@pytest.fixture
def sample_parquet(sample_data):
    """Write sample data to a temporary parquet file."""
    df = pd.DataFrame(sample_data)
    with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
        df.to_parquet(f.name, engine="pyarrow", row_group_size=2000)
        yield f.name
    os.unlink(f.name)


@pytest.fixture
def sample_hdf5(sample_data):
    """Write sample data to a temporary HDF5 file."""
    df = pd.DataFrame(sample_data)
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as f:
        df.to_hdf(f.name, key="df", format="table")
        yield f.name
    os.unlink(f.name)


# ---- PolarsAdapter tests ----


@pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
class TestPolarsAdapter:

    def test_is_subclass(self):
        from radis.io.polars_adapter import PolarsAdapter

        assert issubclass(PolarsAdapter, DataFrameAdapter)

    def test_load_parquet(self, sample_parquet):
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_parquet)
        assert len(adapter) == 10000

    def test_load_hdf5_fallback(self, sample_hdf5):
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_hdf5)
        assert len(adapter) == 10000

    def test_load_with_columns(self, sample_parquet):
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_parquet, columns=["wav", "int", "A"])
        result = adapter.compute()
        assert list(result.columns) == ["wav", "int", "A"]

    def test_filter_range(self, sample_parquet):
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_parquet)
        adapter.filter_range("wav", 1900, 2300)
        result = adapter.compute()
        wav = result["wav"].to_numpy()
        assert np.all(wav >= 1900)
        assert np.all(wav <= 2300)
        assert len(result) < 10000
        assert len(result) > 0

    def test_filter_equals(self, sample_parquet):
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_parquet)
        adapter.filter_equals("iso", 1.0)
        result = adapter.compute()
        assert np.all(result["iso"].to_numpy() == 1.0)

    def test_select_columns(self, sample_parquet):
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_parquet)
        adapter.select_columns(["wav", "int", "A"])
        result = adapter.compute()
        assert list(result.columns) == ["wav", "int", "A"]

    def test_chained_operations(self, sample_parquet):
        """Test filter + select together like calc_spectrum would do."""
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_parquet)
        adapter.filter_range("wav", 2000, 3000)
        adapter.select_columns(["wav", "int", "A", "airbrd", "selbrd", "El"])
        result = adapter.compute()
        assert list(result.columns) == [
            "wav", "int", "A", "airbrd", "selbrd", "El"
        ]
        wav = result["wav"].to_numpy()
        assert np.all(wav >= 2000)
        assert np.all(wav <= 3000)

    def test_to_pandas(self, sample_parquet):
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_parquet)
        adapter.filter_range("wav", 2000, 2100)
        pdf = adapter.to_pandas()
        assert isinstance(pdf, pd.DataFrame)
        assert "wav" in pdf.columns
        assert len(pdf) > 0

    def test_to_numpy(self, sample_parquet):
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_parquet)
        adapter.filter_range("wav", 2000, 2100)
        arrays = adapter.to_numpy()
        assert isinstance(arrays, dict)
        assert isinstance(arrays["wav"], np.ndarray)
        assert np.all(arrays["wav"] >= 2000)

    def test_filter_range_iso_then_pandas(self, sample_parquet):
        """Simulate a typical RADIS loading workflow."""
        from radis.io.polars_adapter import PolarsAdapter

        adapter = PolarsAdapter()
        adapter.load(sample_parquet)
        adapter.filter_range("wav", 1900, 2300)
        adapter.filter_equals("iso", 1.0)
        adapter.select_columns(["wav", "int", "A", "airbrd", "selbrd", "El"])
        pdf = adapter.to_pandas()
        assert isinstance(pdf, pd.DataFrame)
        assert all(pdf["wav"] >= 1900)
        assert all(pdf["wav"] <= 2300)


# ---- PandasAdapter tests ----


class TestPandasAdapter:

    def test_is_subclass(self):
        from radis.io.polars_adapter import PandasAdapter

        assert issubclass(PandasAdapter, DataFrameAdapter)

    def test_load_parquet(self, sample_parquet):
        from radis.io.polars_adapter import PandasAdapter

        adapter = PandasAdapter()
        adapter.load(sample_parquet)
        assert len(adapter) == 10000

    def test_load_hdf5(self, sample_hdf5):
        from radis.io.polars_adapter import PandasAdapter

        adapter = PandasAdapter()
        adapter.load(sample_hdf5)
        assert len(adapter) == 10000

    def test_filter_range(self, sample_parquet):
        from radis.io.polars_adapter import PandasAdapter

        adapter = PandasAdapter()
        adapter.load(sample_parquet)
        adapter.filter_range("wav", 1900, 2300)
        pdf = adapter.to_pandas()
        assert np.all(pdf["wav"] >= 1900)
        assert np.all(pdf["wav"] <= 2300)


# ---- get_adapter factory tests ----


@pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
def test_get_adapter_polars():
    from radis.io.polars_adapter import PolarsAdapter, get_adapter

    adapter = get_adapter("polars")
    assert isinstance(adapter, PolarsAdapter)


def test_get_adapter_pandas():
    from radis.io.polars_adapter import PandasAdapter, get_adapter

    adapter = get_adapter("pytables")
    assert isinstance(adapter, PandasAdapter)


def test_get_adapter_invalid():
    from radis.io.polars_adapter import get_adapter

    with pytest.raises(ValueError):
        get_adapter("invalid_engine")


# ---- HDF5 to Parquet migration tests ----


class TestHDF5ToParquet:

    def test_convert_basic(self, sample_hdf5):
        from radis.io.hdf5_to_parquet import convert_hdf5_to_parquet

        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
            out_path = f.name
        try:
            convert_hdf5_to_parquet(sample_hdf5, out_path, verbose=False)
            assert os.path.exists(out_path)
            # check we can read it back
            check = pd.read_parquet(out_path)
            assert len(check) == 10000
        finally:
            if os.path.exists(out_path):
                os.unlink(out_path)

    def test_sorted_by_wav(self, sample_hdf5):
        from radis.io.hdf5_to_parquet import convert_hdf5_to_parquet

        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
            out_path = f.name
        try:
            convert_hdf5_to_parquet(
                sample_hdf5, out_path, sort_by="wav", verbose=False
            )
            check = pd.read_parquet(out_path)
            wav = check["wav"].values
            # must be sorted
            assert np.all(wav[:-1] <= wav[1:])
        finally:
            if os.path.exists(out_path):
                os.unlink(out_path)

    def test_data_integrity(self, sample_hdf5):
        """Original and converted data should match within tolerance."""
        from radis.io.hdf5_to_parquet import convert_hdf5_to_parquet

        original = pd.read_hdf(sample_hdf5)
        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
            out_path = f.name
        try:
            convert_hdf5_to_parquet(
                sample_hdf5, out_path, validate=True, verbose=False
            )
            converted = pd.read_parquet(out_path)
            assert len(original) == len(converted)
            # check a numeric column
            orig_wav = np.sort(original["wav"].values)
            conv_wav = np.sort(converted["wav"].values)
            assert np.allclose(orig_wav, conv_wav, rtol=1e-14)
        finally:
            if os.path.exists(out_path):
                os.unlink(out_path)

    def test_get_parquet_path(self, sample_hdf5):
        from radis.io.hdf5_to_parquet import get_parquet_path

        parquet_path = get_parquet_path(sample_hdf5)
        assert parquet_path.endswith(".parquet")
        assert os.path.exists(parquet_path)
        # clean up
        os.unlink(parquet_path)

    def test_skip_existing(self, sample_hdf5):
        from radis.io.hdf5_to_parquet import convert_hdf5_to_parquet

        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
            out_path = f.name
        try:
            # first conversion
            convert_hdf5_to_parquet(sample_hdf5, out_path, verbose=False)
            mtime1 = os.path.getmtime(out_path)
            # second call should skip
            convert_hdf5_to_parquet(sample_hdf5, out_path, verbose=False)
            mtime2 = os.path.getmtime(out_path)
            assert mtime1 == mtime2  # file not recreated
        finally:
            if os.path.exists(out_path):
                os.unlink(out_path)