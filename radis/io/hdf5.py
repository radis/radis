# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 21:27:15 2021

@author: erwan
"""

import sys
from time import time

import pandas as pd


def hdf2df(
    fname,
    columns=None,
    isotope=None,
    load_wavenum_min=None,
    load_wavenum_max=None,
    verbose=True,
    store_kwargs={},
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
    try:
        df = pd.read_hdf(fname, columns=columns, where=where, **store_kwargs)
    except TypeError as err:
        if "reading from a Fixed format store" in str(err):
            raise TypeError(
                f"radis.io.hdf5.hdf2df can only be used to load specific HDF5 files generated in a 'Table' which allows to select only certain columns or rows. Here the file {fname} is in 'Fixed' format. Regenerate it ? If it's a cache file of a .par file, load the .par file directly ?"
            )
    with pd.HDFStore(fname, mode="r") as store:
        df.attrs.update(store.get_storer("df").attrs.metadata)

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
