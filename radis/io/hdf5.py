# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 21:27:15 2021

@author: erwan
"""

import pandas as pd


def hdf2df(fname, count=-1, cache=False, verbose=True, drop_non_numeric=True):
    """
    Load a HDF5 line databank into a Pandas DataFrame

    Parameters
    ----------
    fname : TYPE
        DESCRIPTION.
    count : TYPE, optional
        DESCRIPTION. The default is -1.
    cache : TYPE, optional
        DESCRIPTION. The default is False.
    verbose : TYPE, optional
        DESCRIPTION. The default is True.
    drop_non_numeric : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------

    df: pandas Dataframe
        dataframe containing all lines and parameters


    """

    return pd.read_hdf(fname)
