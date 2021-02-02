# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 22:40:51 2021

@author: erwan


https://stackoverflow.com/questions/55610891/numpy-load-from-io-bytesio-stream
https://stupidpythonideas.blogspot.com/2014/07/three-ways-to-read-files.html

"""

import os
import sys
from datetime import date
from os.path import abspath, exists, expanduser, join
from time import time

import numpy as np
import pandas as pd
from numpy import DataSource

from radis.db import MOLECULES_LIST_NONEQUILIBRIUM
from radis.io.hdf5 import hdf2df
from radis.io.hitran import columns_2004
from radis.io.tools import _create_dtype, _get_linereturnformat, _ndarray2df
from radis.misc.config import addDatabankEntries, getDatabankList
from radis.misc.warning import DatabaseAlreadyExists

BASE_URL = "https://hitran.org/hitemp/data/bzip2format/"
HITEMP_SOURCE_FILES = {
    "CO2": "",  # NotImplemented
    "N2O": "04_HITEMP2019.par.bz2",
    "CO": "05_HITEMP2019.par.bz2",
    "CH4": "06_HITEMP2020.par.bz2",
    "NO": "08_HITEMP2019.par.bz2",
    "NO2": "10_HITEMP2019.par.bz2",
    "OH": "13_HITEMP2020.par.bz2",
}
"""
dict: list of available HITEMP source files

Last update : 02 Feb 2021
Compare with files available on https://hitran.org/hitemp/
"""

DATA_COLUMNS = ["iso", "wav"]
"""
list : only these column names will be searchable directly on disk to
only load certain lines. See :py:func:`~radis.io.hdf5.hdf2df`
"""
# TODO: WIP. Maybe move somewhere else to be also used by HITRAN queries


def fetch_hitemp(
    molecule,
    local_databases="~/.radisdb/",
    databank_name="HITEMP-{molecule}",
    isotope=None,
    load_only_wavenum_above=None,
    load_only_wavenum_below=None,
    cache=True,
    verbose=True,
):
    """Stream HITEMP file from HITRAN website. Unzip and build a HDF5 file directly.
    Returns a Pandas DataFrame containing all lines.

    Parameters
    ----------
    molecule: `"CO2", "N2O", "CO", "CH4", "NO", "NO2", "OH"`
        HITEMP molecule. See :py:attr:`~radis.io.hitemp.HITEMP_SOURCE_FILES`
    local_databases: str
        where to download the files
    databank_name: str
        name of the databank in RADIS :ref:`Configuration file <label_lbl_config_file>`
        Default ``"HITEMP-{molecule}"``
    isotope: str
        load only certain isotopes : ``'2'``, ``'1,2'``, etc. If ``None``, loads
        everything. Default ``None``.
    load_only_wavenum_above, load_only_wavenum_below: float (cm-1)
        load only specific wavelength.

    Other Parameters
    ----------------

    cache: bool, or ``'regen'``
        if ``True``, use existing HDF5 file. If ``False`` or ``'regen'``, rebuild it.
    verbose: bool

    Returns
    -------

    df: pd.DataFrame
        Line list
        A HDF5 file is also created in ``local_databases`` and referenced
        in the :ref:`RADIS config file <label_lbl_config_file>` with name
        ``databank_name``

    Notes
    -----

    if using ``load_only_wavenum_above/below`` or ``isotope``, the whole
    database is anyway downloaded and uncompressed in ``local_databases``
    fast .HDF5 files (which will take a long time on first call). Only
    the expected wavenumber range & isotopes are returned with :py:func:`~radis.io.hdf5.hdf2df`

    See Also
    --------
    :py:func:`~radis.io.hdf5.hdf2df`

    """
    # TODO: unzip only parts of the database
    # see https://github.com/radis/radis/pull/194

    # TODO: clean DataSource downloads after HDF5 is generated

    if databank_name == "HITEMP-{molecule}":
        databank_name = databank_name.format(**{"molecule": molecule})
    local_databases = local_databases.replace("~", expanduser("~"))

    if molecule in ["H2O", "CO2"]:
        raise NotImplementedError(
            "Multiple files. Download manually on https://hitran.org/hitemp/ "
        )

    inputf = HITEMP_SOURCE_FILES[molecule]
    urlname = BASE_URL + inputf

    try:
        os.mkdir(local_databases)
    except OSError:
        pass
    else:
        if verbose:
            print("Created folder :", local_databases)

    output = join(local_databases, molecule + "-" + inputf.replace(".par.bz2", ".h5"))

    if not cache or cache == "regen":
        # Delete existing HDF5 file
        if exists(output):
            if verbose:
                print("Removing existing file ", output)
                # TODO: also clean the getDatabankList? Todo once it is in JSON format. https://github.com/radis/radis/issues/167
            os.remove(output)

    if exists(output):
        # TODO check metadata here?
        try:
            assert databank_name in getDatabankList()
        except AssertionError:
            raise AssertionError(
                f"{databank_name} not in {getDatabankList()} although {output} exists. Maybe try regenerating cache files with `SpectrumFactory.fetch_databank(..., use_cached='regen'`"
            )
        if verbose:
            print(f"Using existing database {databank_name}")
        return hdf2df(
            output,
            isotope=isotope,
            load_only_wavenum_below=load_only_wavenum_below,
            load_only_wavenum_above=load_only_wavenum_above,
            verbose=verbose,
        )

    # Doesnt exist : download
    ds = DataSource(local_databases)

    if verbose:
        print(f"Downloading {inputf} for {molecule} in {output}")
    download_date = date.today().strftime("%d %b %Y")

    columns = columns_2004

    # Get linereturn (depends on OS, but file may also have been generated
    # on a different OS. Here we simply read the file to find out)
    with ds.open(urlname) as gfile:  # locally downloaded file

        dt = _create_dtype(columns, "a2")  # 'a2' allocates space to get \n or \n\r
        b = np.empty(161, dtype=dt)
        gfile.readinto(b)
        linereturnformat = _get_linereturnformat(b, columns)

    with ds.open(urlname) as gfile:  # locally downloaded file

        CHUNK = 1000  # TODO: adjust if needed (RAM dependant?)

        dt = _create_dtype(columns, linereturnformat)
        b = np.empty(
            dt.itemsize * CHUNK, dtype=dt
        )  # receives the HITRAN 160-format data.

        t0 = time()
        wmin = np.inf
        wmax = 0
        if verbose:
            print("Building", output)

        Nlines = 0

        with pd.HDFStore(output, mode="a", complib="blosc", complevel=9) as f:
            for nbytes in iter(lambda: gfile.readinto(b), 0):

                if not b[-1]:
                    # End of file within the chunk (but does not start with end of file)
                    # so nbytes != 0
                    b = get_last(b)

                df = _ndarray2df(b, columns, linereturnformat)

                # df.to_hdf(
                #     output, "df", format="table", append=True, complib="blosc", complevel=9
                # )
                f.put(key="df", value=df, format="table", data_columns=DATA_COLUMNS)

                wmin = np.min((wmin, df.wav.min()))
                wmax = np.max((wmax, df.wav.max()))
                Nlines += len(df)
                if verbose:
                    sys.stdout.write(
                        "\r({0:.0f}s) Built {1} database from {2:.1f} to {3:.1f} cm-1 ({4} lines)".format(
                            time() - t0, molecule, wmin, wmax, Nlines
                        )
                    )

                # Reinitialize for next read
                b = np.empty(
                    dt.itemsize * CHUNK, dtype=dt
                )  # receives the HITRAN 160-format data.

            f.get_storer("df").attrs.metadata = {
                "wmin": wmin,
                "wmax": wmax,
                "download_date": download_date,
                "download_url": urlname,
            }
    if verbose:
        sys.stdout.write("\n")

    # Add to registered Databanks in RADIS config file
    dict_entries = {
        "info": f"HITEMP {molecule} lines ({wmin:.1f}-{wmax:.1f} cm-1) with TIPS-2017 (through HAPI) for partition functions",
        "path": abspath(output),
        "format": "hdf5",
        "parfuncfmt": "hapi",
        "download_date": download_date,
        "download_url": urlname,
    }
    if molecule in MOLECULES_LIST_NONEQUILIBRIUM:
        dict_entries[
            "info"
        ] += " and RADIS spectroscopic constants for rovibrational energies (nonequilibrium)"
        dict_entries["levelsfmt"] = "radis"

    try:
        addDatabankEntries(databank_name, dict_entries)
    except DatabaseAlreadyExists as err:
        raise DatabaseAlreadyExists(
            f"{databank_name} already exists in your ~/.radis config file. "
            + "Choose a different name for the downloaded database with `fetch_hitemp(databank_name=...)`"
        ) from err

    return hdf2df(
        output,
        isotope=isotope,
        load_only_wavenum_below=load_only_wavenum_below,
        load_only_wavenum_above=load_only_wavenum_above,
        verbose=verbose,
    )


#%%
def get_last(b):
    """ Get non-empty lines of a chunk b, parsing the bytes """

    element_length = np.vectorize(lambda x: len(x.__str__()))(b)
    non_zero = element_length > element_length[-1]
    threshold = non_zero.argmin() - 1
    assert (non_zero[: threshold + 1] == 1).all()
    assert (non_zero[threshold + 1 :] == 0).all()
    return b[non_zero]


#%%

if __name__ == "__main__":

    output = fetch_hitemp("OH")
