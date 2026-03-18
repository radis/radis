from abc import ABC, abstractmethod


class DataFrameAdapter(ABC):
    """Base adapter for swapping DataFrame backends in RADIS.
    Configured via config["DATAFRAME_ENGINE"] in ~/radis.json
    """

    @abstractmethod
    def load(self, path, columns=None):
        ...

    @abstractmethod
    def filter_range(self, column, vmin, vmax):
        ...

    @abstractmethod
    def filter_equals(self, column, value):
        ...

    @abstractmethod
    def select_columns(self, columns):
        ...

    @abstractmethod
    def compute(self):
        ...

    @abstractmethod
    def to_pandas(self):
        ...

    @abstractmethod
    def to_numpy(self):
        ...

    @abstractmethod
    def __len__(self):
        ...