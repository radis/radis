# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 21:27:15 2021

@author: erwan
"""

import pandas as pd


def hdf2df(fname, verbose=True, **kwargs):
    """
    Load a HDF5 line databank into a Pandas DataFrame.
    Adds HDF5 metadata in ``df.attrs``

    Parameters
    ----------
    fname : str
        HDF5 file name

    Other Parameters
    ----------------

    **kwargs: dict
        arguments forwarded to :py:meth:`~pandas.io.pytables.HDFStore`

    Returns
    -------

    df: pandas Dataframe
        dataframe containing all lines or energy levels

    Notes
    -----

    DataFrame metadata in ``df.attrs`` is still experimental in Pandas and can be lost
    during ``groupby, pivot, join or loc`` operations on the Dataframe.
    See https://stackoverflow.com/questions/14688306/adding-meta-information-metadata-to-pandas-dataframe

    Always check for existence !

    """

    # df = pd.read_hdf(fname)

    with pd.HDFStore(fname, mode="r+", **kwargs) as store:
        df = store["df"]
        df.attrs.update(store.get_storer("df").attrs.metadata)

    return df
