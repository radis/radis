# -*- coding: utf-8 -*-
"""
Functions to read/write Spectrum objects under HDF5 format
"""

import pandas as pd

try:
    from .hdf5 import HDF5Manager
except ImportError:
    from radis.io.hdf5 import HDF5Manager


# Convert


def to_pandas(s, copy=True):
    """Convert a Spectrum to a Pandas DataFrame

    Returns
    -------
    pd.DataFrame : pandas DataFrame where columns are spectral arrays, and
        units are stored in attributes ``df.attrs``

    Notes
    -----
    Pandas does not support units yet. pint-pandas is an advanced
    project but not fully working.
    See discussion in https://github.com/pandas-dev/pandas/issues/10349

    For the moment, we store units as metadata"""
    df = pd.DataFrame(s._q, copy=copy)

    df.attrs = s.units

    return df


def spec2hdf(s, file, engine="pytables"):

    df_arrays = to_pandas(s, copy=False)

    # Store
    mgr = HDF5Manager(engine=engine)

    # Store spectral arrays and lines as Datasets

    mgr.write(file, df_arrays, key="arrays", append=False)

    if s.lines:
        mgr.write(file, s.lines, key="lines", append=True)

    # Store the rest as metadata

    metadata = {}
    metadata["units"] = s.units

    if s.populations:
        metadata["populations"] = s.populations

    metadata["conditions"] = s.conditions

    if s.references:
        metadata["references"] = s.references

    mgr.add_metadata(file, metadata, key="arrays")


def hdf2spec(file, engine="pytables"):
    from radis import Spectrum

    mgr = HDF5Manager(engine=engine)

    df2 = mgr.load(file, key="arrays")
    metadata = mgr.read_metadata(file, key="arrays")

    units2 = metadata["units"]
    conditions2 = metadata["conditions"]
    populations2 = metadata.get("populations", None)
    references2 = metadata.get("references", None)

    try:
        lines2 = mgr.load(file, key="lines")
    except KeyError:
        lines2 = None

    s2 = Spectrum(
        {k: v.values for k, v in dict(df2).items()},
        lines=lines2,
        units=units2,
        conditions=conditions2,
        references=references2,
        populations=populations2,
    )

    return s2
