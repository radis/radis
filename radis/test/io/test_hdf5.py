# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 20:36:38 2021

@author: erwan
"""

import pytest

from radis.api.hdf5 import DataFileManager, hdf2df
from radis.io.hitemp import fetch_hitemp
from radis.misc.config import getDatabankEntries


@pytest.mark.fast
def test_hdf5_io_engines(*args, **kwargs):
    """Test different engines implemented in :py:class:`radis.api.hdf5.DataFileManager`"""

    import os
    from os.path import exists

    for file in ["test_pytables.h5", "test_h5py.h5"]:
        if exists(file):
            os.remove(file)

    # Test data
    import h5py
    import numpy as np
    import pandas as pd

    df0 = pd.DataFrame({"a": np.arange(10) ** 2, "b": np.arange(10) ** 3})
    metadata0 = {"some_metadata": True}
    # ... a Pandas HDFStore file
    df0.to_hdf("test_pytables.h5", "df")
    # ... a h5py file with same content
    with h5py.File("test_h5py.h5", "w") as f:
        # group = f.create_group("df")
        for c in df0.columns:
            f.create_dataset(c, data=df0[c])

    # Test Pytables engine : add_metadata ; load ; read_metadata
    manager = DataFileManager(engine="pytables")
    manager.add_metadata("test_pytables.h5", metadata0)

    df = manager.read("test_pytables.h5")
    assert (df == df0).all().all()
    assert manager.read_metadata("test_pytables.h5") == metadata0

    # Test h5py engine : add_metadata ; load ; read_metadata
    manager = DataFileManager(engine="h5py")
    manager.add_metadata("test_h5py.h5", metadata0, key=None)
    df = manager.read("test_h5py.h5", key=None)
    assert (df == df0).all().all()
    assert manager.read_metadata("test_h5py.h5") == metadata0

    # Test vaex engine : add_metadata ; load ; read_
    # .. also able to read h5py files
    manager = DataFileManager(engine="vaex")
    manager.add_metadata("test_h5py.h5", metadata0, key=None)
    df = manager.read("test_h5py.h5", key=None)
    # this time; df is a vaex DataFrame
    # ... it keeps the file open so we should close the file handle
    assert (df.to_pandas_df() == df0).all().all()
    df.close()
    assert manager.read_metadata("test_h5py.h5", key=None) == metadata0

    # Test get_columns function of DataFileManager
    # For vaex
    manager = DataFileManager(engine="vaex")
    assert list(manager.get_columns("test_h5py.h5")) == ["a", "b"]

    # For pytables
    manager = DataFileManager(engine="pytables")
    assert list(manager.get_columns("test_pytables.h5")) == ["a", "b"]


@pytest.mark.needs_connection
def test_local_hdf5_lines_loading(*args, **kwargs):
    """
    We use the OH HITEMP line database to test :py:func:`~radis.io.hitemp.fetch_hitemp`
    and :py:func:`~radis.api.hdf5.hdf2df`

    - Partial loading (only specific wavenumbers)
    - Only certain isotopes
    - Only certain columns

    """

    fetch_hitemp("OH")  # to initialize the database

    path = getDatabankEntries("HITEMP-OH")["path"]

    # Initialize the database
    fetch_hitemp("OH")
    path = getDatabankEntries("HITEMP-OH")["path"][0]
    df = hdf2df(path)
    wmin, wmax = df.wav.min(), df.wav.max()
    assert wmin < 2300  # needed for next test to be valid
    assert wmax > 2500  # needed for next test to be valid
    assert len(df.columns) > 5  # many columns loaded by default
    assert len(df.iso.unique()) > 1

    # Test loading only certain columns
    df = hdf2df(path, columns=["wav", "int"])
    assert len(df.columns) == 2 and "wav" in df.columns and "int" in df.columns

    # Test loading only certain isotopes
    df = hdf2df(path, isotope="2")
    assert df.iso.unique() == 2

    # Test partial loading of wavenumbers
    df = hdf2df(path, load_wavenum_min=2300, load_wavenum_max=2500)
    assert df.wav.min() >= 2300
    assert df.wav.max() <= 2500

    # Test with only one
    assert hdf2df(path, load_wavenum_min=2300).wav.min() >= 2300

    # Test with the other
    assert hdf2df(path, load_wavenum_max=2500).wav.max() <= 2500


if __name__ == "__main__":
    test_hdf5_io_engines()
    test_local_hdf5_lines_loading()
