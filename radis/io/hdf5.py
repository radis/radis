# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 21:27:15 2021

@author: erwan
"""

import sys
from time import time

import pandas as pd


class HDF5Manager(object):
    def __init__(self, engine="pytables"):
        self.engine = engine

    def open(self, file, mode="w"):
        if self.engine == "pytables":
            return pd.HDFStore(file, mode=mode, complib="blosc", complevel=9)
        else:
            raise NotImplementedError(self.engine)

    def write(self, handler, df, append=True):
        if self.engine == "pytables":
            DATA_COLUMNS = ["iso", "wav"]
            """
            list : only these column names will be searchable directly on disk to
            only load certain lines. See :py:func:`~radis.io.hdf5.hdf2df`
            """
            handler.put(
                key="df",
                value=df,
                append=append,
                format="table",
                data_columns=DATA_COLUMNS,
            )
        else:
            raise NotImplementedError(self.engine)

    def load(self, fname, columns, where=None, **store_kwargs) -> pd.DataFrame:
        """
        Parameters
        ----------
        columns: list of str
            list of columns to load. If ``None``, returns all columns in the file.
        """
        if self.engine == "pytables":
            try:
                df = pd.read_hdf(fname, columns=columns, where=where, **store_kwargs)
            except TypeError as err:
                if "reading from a Fixed format store" in str(err):
                    raise TypeError(
                        f"radis.io.hdf5.hdf2df can only be used to load specific HDF5 files generated in a 'Table' which allows to select only certain columns or rows. Here the file {fname} is in 'Fixed' format. Regenerate it ? If it's a cache file of a .par file, load the .par file directly ?"
                    )
        else:
            raise NotImplementedError(self.engine)

        return df

    def add_metadata(self, local_file: str, metadata: dict):

        if self.engine == "pytables":
            with pd.HDFStore(local_file, mode="a", complib="blosc", complevel=9) as f:

                f.get_storer("df").attrs.metadata = metadata
        else:
            raise NotImplementedError(self.engine)

    def read_metadata(self, local_file: str) -> dict:

        if self.engine == "pytables":
            with pd.HDFStore(local_file, mode="r", complib="blosc", complevel=9) as f:

                metadata = f.get_storer("df").attrs.metadata
        else:
            raise NotImplementedError(self.engine)

        return metadata


def hdf2df(
    fname,
    columns=None,
    isotope=None,
    load_wavenum_min=None,
    load_wavenum_max=None,
    verbose=True,
    store_kwargs={},
    engine="pytables",
):
    """Load a HDF5 line databank into a Pandas DataFrame.

    Adds HDF5 metadata in ``df.attrs``

    Parameters
    ----------
    fname : str
        HDF5 file name
    columns: list of str
        list of columns to load. If ``None``, returns all columns in the file.
    isotope: str
        load only certain isotopes : ``'2'``, ``'1,2'``, etc. If ``None``, loads
        everything. Default ``None``.
    load_wavenum_min, load_wavenum_max: float (cm-1)
        load only specific wavelength.

    Other Parameters
    ----------------
    store_kwargs: dict
        arguments forwarded to :py:meth:`~pandas.io.pytables.read_hdf`


    Returns
    -------
    df: pandas Dataframe
        dataframe containing all lines or energy levels

    Examples
    --------

    ::


        path = getDatabankEntries("HITEMP-OH")['path'][0]
        df = hdf2df(path)

        df = hdf2df(path, columns=['wav', 'int'])

        df = hdf2df(path, isotope='2')
        df = hdf2df(path, isotope='1,2)

        df = hdf2df(path, load_wavenum_min=2300, load_wavenum_max=2500)

    Notes
    -----

    DataFrame metadata in ``df.attrs`` is still experimental in Pandas and can be lost
    during ``groupby, pivot, join or loc`` operations on the Dataframe.
    See https://stackoverflow.com/questions/14688306/adding-meta-information-metadata-to-pandas-dataframe

    Always check for existence !

    """
    where = []
    if load_wavenum_min is not None:
        where.append(f"wav > {load_wavenum_min}")
    if load_wavenum_max is not None:
        where.append(f"wav < {load_wavenum_max}")
    if isotope:
        where.append(f'iso in {isotope.split(",")}')

    # Load :
    t0 = time()

    manager = HDF5Manager(engine)
    df = manager.load(fname, columns=columns, where=where, **store_kwargs)
    metadata = manager.read_metadata(fname)

    # Sanity Checks
    if "total_lines" in metadata:
        assert len(df) == metadata["total_lines"]
    if "wavenumber_min" in metadata:
        assert df["wav"].min() == metadata["wavenumber_min"]
    if "wavenumber_max" in metadata:
        assert df["wav"].max() == metadata["wavenumber_max"]

    df.attrs.update(metadata)

    if verbose >= 3:
        from radis.misc.printer import printg

        printg(
            f"Generated dataframe from {fname} in {time()-t0:.2f}s ({len(df)} rows, {len(df.columns)} columns, {sys.getsizeof(df)*1e-6:.2f} MB)"
        )

    return df


#%%

if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_hdf5.py"]))
