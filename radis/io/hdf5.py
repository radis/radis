# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 21:27:15 2021

@author: erwan
"""

import pandas as pd


def hdf2df(
    fname,
    columns=None,
    isotope=None,
    load_only_wavenum_above=None,
    load_only_wavenum_below=None,
    verbose=True,
    store_kwargs={},
):
    """
    Load a HDF5 line databank into a Pandas DataFrame.
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
    load_only_wavenum_above, load_only_wavenum_below: float (cm-1)
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

        df = hdf2df(path, load_only_wavenum_above=2300,
                    load_only_wavenum_below=2500)

    Notes
    -----

    DataFrame metadata in ``df.attrs`` is still experimental in Pandas and can be lost
    during ``groupby, pivot, join or loc`` operations on the Dataframe.
    See https://stackoverflow.com/questions/14688306/adding-meta-information-metadata-to-pandas-dataframe

    Always check for existence !

    """

    where = []
    if load_only_wavenum_above is not None:
        where.append(f"wav > {load_only_wavenum_above}")
    if load_only_wavenum_below is not None:
        where.append(f"wav < {load_only_wavenum_below}")
    if isotope:
        where.append(f'iso in {isotope.split(",")}')

    # Load :
    try:
        df = pd.read_hdf(fname, columns=columns, where=where, **store_kwargs)
    except TypeError as err:
        if "reading from a Fixed format store" in str(err):
            raise TypeError(
                f"radis.io.hdf5.hdf2df can only be used to load specific HDF5 files generated in a 'Table' which allows to select only certain columns or rows. Here the file {fname} is in 'Fixed' format. Regenerate it ? If it's a cache file of a .par file, load the .par file directly ?"
            )

    with pd.HDFStore(fname, mode="r+") as store:
        df.attrs.update(store.get_storer("df").attrs.metadata)

    return df
