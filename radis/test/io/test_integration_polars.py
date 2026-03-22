"""Integration tests for Polars backend — issue #978"""
import os, tempfile, shutil
import numpy as np
import pandas as pd
import pytest

try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False

def make_db(n=50_000):
    rng = np.random.default_rng(42)
    return pd.DataFrame({
        "wav":    np.sort(rng.uniform(1000, 5000, n)),
        "int":    rng.lognormal(-20, 5, n),
        "A":      rng.lognormal(1, 2, n),
        "airbrd": rng.uniform(0.01, 0.15, n),
        "selbrd": rng.uniform(0.01, 0.15, n),
        "El":     rng.uniform(0, 5000, n),
        "iso":    rng.choice([1.0, 2.0, 3.0], n, p=[0.7, 0.2, 0.1]),
        "id":     np.ones(n),
    })

@pytest.fixture
def h5_and_parquet():
    df = make_db()
    td = tempfile.mkdtemp()
    h5 = os.path.join(td, "test.h5")
    pq = os.path.join(td, "test.parquet")
    df.to_hdf(h5, key="df", format="table")
    df.sort_values("wav").to_parquet(pq, engine="pyarrow",
                                     row_group_size=5000, index=False)
    yield h5, pq, df
    shutil.rmtree(td)

@pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
class TestDataFileManagerPolars:

    def test_engine_is_polars(self, h5_and_parquet):
        from radis.api.hdf5 import DataFileManager
        mgr = DataFileManager(engine="polars")
        assert mgr.engine == "polars"

    def test_load_returns_dataframe(self, h5_and_parquet):
        h5, pq, _ = h5_and_parquet
        from radis.api.hdf5 import DataFileManager
        mgr = DataFileManager(engine="polars")
        result = mgr.load(pq, lower_bound=[("wav", 2000)],
                          upper_bound=[("wav", 3000)])
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_predicate_pushdown_correct_rows(self, h5_and_parquet):
        h5, pq, df = h5_and_parquet
        from radis.api.hdf5 import DataFileManager
        wmin, wmax = 2000.0, 3000.0
        mgr_pl = DataFileManager(engine="polars")
        df_pl = mgr_pl.load(pq, lower_bound=[("wav", wmin)],
                            upper_bound=[("wav", wmax)])
        df_pd = df[(df["wav"] >= wmin) & (df["wav"] <= wmax)]
        assert len(df_pl) == len(df_pd)
        assert df_pl["wav"].min() >= wmin
        assert df_pl["wav"].max() <= wmax

    def test_wav_range_correct(self, h5_and_parquet):
        _, pq, _ = h5_and_parquet
        from radis.api.hdf5 import DataFileManager
        mgr = DataFileManager(engine="polars")
        result = mgr.load(pq, lower_bound=[("wav", 1500)],
                          upper_bound=[("wav", 2500)])
        assert result["wav"].min() >= 1500
        assert result["wav"].max() <= 2500

class TestConfigIntegration:

    def test_get_engine_default(self):
        import os
        os.environ.pop("RADIS_DATAFRAME_ENGINE", None)
        from radis.misc.config import get_engine
        engine = get_engine()
        assert engine in ("polars", "pytables", "pandas")

    def test_set_engine(self):
        import os
        from radis.misc.config import set_engine, get_engine
        set_engine("polars")
        assert get_engine() == "polars"
        set_engine("pytables")
        assert get_engine() == "pytables"
        os.environ.pop("RADIS_DATAFRAME_ENGINE", None)
