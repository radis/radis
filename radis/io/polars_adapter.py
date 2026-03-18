import polars as pl
from pathlib import Path
from radis.io.dataframe_adapter import DataFrameAdapter


class PolarsAdapter(DataFrameAdapter):

    def __init__(self):
        self._lf = None
        self._path = None

    def load(self, path, columns=None):
        self._path = Path(path)

        if self._path.suffix == ".parquet":
            self._lf = pl.scan_parquet(path)
        elif self._path.suffix in (".hdf5", ".h5"):
            # fallback for old HDF5 files that haven't been migrated yet
            import pandas as pd

            df = pd.read_hdf(path)
            self._lf = pl.from_pandas(df).lazy()
        else:
            raise ValueError(f"Unsupported file format: {self._path.suffix}")

        if columns:
            self._lf = self._lf.select(columns)
        return self

    def filter_range(self, column, vmin, vmax):
        # predicate pushdown — only reads matching row groups from parquet
        self._lf = self._lf.filter(pl.col(column).is_between(vmin, vmax))
        return self

    def filter_equals(self, column, value):
        self._lf = self._lf.filter(pl.col(column) == value)
        return self

    def select_columns(self, columns):
        self._lf = self._lf.select(columns)
        return self

    def compute(self):
        return self._lf.collect()

    def to_pandas(self):
        return self.compute().to_pandas()

    def to_numpy(self):
        df = self.compute()
        return {col: df[col].to_numpy() for col in df.columns}

    def __len__(self):
        return self._lf.select(pl.len()).collect().item()


class PandasAdapter(DataFrameAdapter):
    """Pandas fallback — wraps existing pytables behavior."""

    def __init__(self):
        self._df = None
        self._columns = None

    def load(self, path, columns=None):
        import pandas as pd

        path = Path(path)
        if path.suffix == ".parquet":
            self._df = pd.read_parquet(path, columns=columns)
        elif path.suffix in (".hdf5", ".h5"):
            self._df = pd.read_hdf(path, columns=columns)
        else:
            raise ValueError(f"Unsupported file format: {path.suffix}")
        return self

    def filter_range(self, column, vmin, vmax):
        self._df = self._df[(self._df[column] >= vmin) & (self._df[column] <= vmax)]
        return self

    def filter_equals(self, column, value):
        self._df = self._df[self._df[column] == value]
        return self

    def select_columns(self, columns):
        self._df = self._df[columns]
        return self

    def compute(self):
        return self._df

    def to_pandas(self):
        return self._df.copy()

    def to_numpy(self):
        return {col: self._df[col].values for col in self._df.columns}

    def __len__(self):
        return len(self._df)


def get_adapter(engine="default"):
    """Return the right adapter based on config or explicit engine name."""
    if engine == "default":
        try:
            from radis.misc.config import config

            engine = config.get("DATAFRAME_ENGINE", "pytables")
        except Exception:
            engine = "pytables"

    if engine == "polars":
        return PolarsAdapter()
    elif engine in ("pytables", "pandas"):
        return PandasAdapter()
    else:
        raise ValueError(
            f"Unknown engine: {engine}. Use 'polars' or 'pytables'."
        )