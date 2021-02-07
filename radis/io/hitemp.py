# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 22:40:51 2021

@author: erwan


https://stackoverflow.com/questions/55610891/numpy-load-from-io-bytesio-stream
https://stupidpythonideas.blogspot.com/2014/07/three-ways-to-read-files.html

"""

import os
from datetime import date
from os.path import abspath, exists, expanduser, join

import numpy as np
import pandas as pd
from numpy import DataSource

import radis
from radis.db import MOLECULES_LIST_NONEQUILIBRIUM
from radis.io.hdf5 import hdf2df
from radis.io.hitran import columns_2004
from radis.io.tools import _create_dtype, _get_linereturnformat, _ndarray2df
from radis.misc.cache_files import check_not_deprecated, check_relevancy
from radis.misc.config import addDatabankEntries, getDatabankList
from radis.misc.progress_bar import ProgressBar
from radis.misc.warning import DatabaseAlreadyExists

BASE_URL = "https://hitran.org/hitemp/data/bzip2format/"
HITEMP_SOURCE_FILES = {
    "H2O": "",  # NotImplemented
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

INFO_HITEMP_LINE_COUNT = {
    "H2O": 114241164,
    "CO2": 11193608,
    "N2O": 3626425,
    "CO": 752976,
    "CH4": 31880412,
    "NO": 1137192,
    "NO2": 1108709,
    "OH": 57019,
}
"""
dict: total line count according to https://hitran.org/hitemp/
Used for tests and to give an idea of the uncompress progress
in :py:func:`~radis.io.hitemp.fetch_hitemp

.. warning::
    may change with new HITEMP updates. Only use as an informative
    param, do not rely on this in the code.
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
    load_wavenum_min=None,
    load_wavenum_max=None,
    cache=True,
    verbose=True,
    chunksize=100000,
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
    load_wavenum_min, load_wavenum_max: float (cm-1)
        load only specific wavenumbers.

    Other Parameters
    ----------------
    cache: bool, or ``'regen'``
        if ``True``, use existing HDF5 file. If ``False`` or ``'regen'``, rebuild it.
    verbose: bool
    chunksize: int
        number of lines to process at a same time. Higher is usually faster
        but can create Memory problems and keep the user uninformed of the progress.

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
    database is anyway downloaded and uncompressed to ``local_databases``
    fast access .HDF5 files (which will take a long time on first call). Only
    the expected wavenumber range & isotopes are returned. The .HFD5 parsing uses
    :py:func:`~radis.io.hdf5.hdf2df`

    See Also
    --------
    :py:func:`~radis.io.hdf5.hdf2df`

    """
    # TODO: unzip only parts of the database
    # see https://github.com/radis/radis/pull/194

    # TODO: clean DataSource downloads after HDF5 is generated

    if databank_name == "HITEMP-{molecule}":
        databank_name = databank_name.format(**{"molecule": molecule})
    local_databases = abspath(local_databases.replace("~", expanduser("~")))

    if molecule in ["H2O", "CO2"]:
        raise NotImplementedError(
            "Automatic HITEMP download not implemented for {0} : multiple files. Download manually on https://hitran.org/hitemp/ ".format(
                molecule
            )
        )

    try:
        inputf = HITEMP_SOURCE_FILES[molecule]
    except KeyError as err:
        raise KeyError(
            f"Please choose one of HITEMP molecules : {list(HITEMP_SOURCE_FILES.keys())}. Got '{molecule}'"
        ) from err
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
        check_not_deprecated(
            output,
            metadata_is={},
            metadata_keys_contain=["wavenumber_min", "wavenumber_max"],
        )
        check_relevancy(
            output,
            relevant_if_metadata_above={"wavenumber_max": load_wavenum_min}
            if load_wavenum_min
            else {},
            relevant_if_metadata_below={"wavenumber_min": load_wavenum_max}
            if load_wavenum_max
            else {},
        )
        try:
            assert databank_name in getDatabankList()
        except AssertionError:
            raise AssertionError(
                f"{databank_name} not declared in your RADIS ~/.config file although {output} exists. If you know this file, add it to ~/.radisdb manually. Else regenerate the database with: \n\t>>> radis.SpectrumFactory().fetch_databank(..., use_cached='regen') \nor\n\t>>> radis.io.hitemp.fetch_hitemp({molecule}, cache='regen')\n\n⚠️ It will re-download & uncompress the whole database from HITEMP.\n\nList of declared databanks: {getDatabankList()} "
            )
        if verbose:
            print(f"Using existing database {databank_name}")
        return hdf2df(
            output,
            isotope=isotope,
            load_wavenum_min=load_wavenum_min,
            load_wavenum_max=load_wavenum_max,
            verbose=verbose,
        )

    # Doesnt exist : download
    ds = DataSource(local_databases)

    if verbose:
        print(f"Downloading {inputf} for {molecule}.")
    download_date = date.today().strftime("%d %b %Y")

    columns = columns_2004

    # Get linereturn (depends on OS, but file may also have been generated
    # on a different OS. Here we simply read the file to find out)
    with ds.open(urlname) as gfile:  # locally downloaded file

        dt = _create_dtype(
            columns, "a2"
        )  # 'a2' allocates space to get \n or \n\r for linereturn character
        b = np.zeros(1, dtype=dt)
        gfile.readinto(b)
        linereturnformat = _get_linereturnformat(b, columns)

    with ds.open(urlname) as gfile:  # locally downloaded file

        dt = _create_dtype(columns, linereturnformat)
        b = np.zeros(chunksize, dtype=dt)  # receives the HITRAN 160-character data.

        wmin = np.inf
        wmax = 0
        if verbose:
            print(f"Download complete. Building {molecule} database to {output}")

        with pd.HDFStore(output, mode="a", complib="blosc", complevel=9) as f:
            Nlines = 0
            Ntotal_lines_expected = INFO_HITEMP_LINE_COUNT[molecule]
            pb = ProgressBar(N=Ntotal_lines_expected, active=verbose)
            for nbytes in iter(lambda: gfile.readinto(b), 0):

                if not b[-1]:
                    # End of file flag within the chunk (but does not start
                    # with End of file flag) so nbytes != 0
                    b = get_last(b)

                df = _ndarray2df(b, columns, linereturnformat)

                # df.to_hdf(
                #     output, "df", format="table", append=True, complib="blosc", complevel=9
                # )
                f.put(
                    key="df",
                    value=df,
                    append=True,
                    format="table",
                    data_columns=DATA_COLUMNS,
                )

                wmin = np.min((wmin, df.wav.min()))
                wmax = np.max((wmax, df.wav.max()))
                Nlines += len(df)
                pb.update(
                    Nlines,
                    message=f"Got {Nlines} / {Ntotal_lines_expected} total lines. {wmin:.3f}-{wmax:.3f} cm-1 wavenumber range is complete.",
                )

                # Reinitialize for next read
                b = np.zeros(
                    chunksize, dtype=dt
                )  # receives the HITRAN 160-character data.

            f.get_storer("df").attrs.metadata = {
                "wavenumber_min": wmin,
                "wavenumber_max": wmax,
                "download_date": download_date,
                "download_url": urlname,
                "version": radis.get_version(add_git_number=False),
            }
            pb.done()

    # Add to registered Databanks in RADIS config file
    dict_entries = {
        "info": f"HITEMP {molecule} lines ({wmin:.1f}-{wmax:.1f} cm-1) with TIPS-2017 (through HAPI) for partition functions",
        "path": [abspath(output)],
        "format": "hdf5",
        "parfuncfmt": "hapi",
        "wavenumber_min": wmin,
        "wavenumber_max": wmax,
        "download_date": download_date,
        "download_url": urlname,
    }
    if molecule in MOLECULES_LIST_NONEQUILIBRIUM:
        dict_entries[
            "info"
        ] += " and RADIS spectroscopic constants for rovibrational energies (nonequilibrium)"
        dict_entries["levelsfmt"] = "radis"

    # Register database in ~/.radis to be able to use it with load_databank()
    try:
        addDatabankEntries(databank_name, dict_entries)
    except DatabaseAlreadyExists as err:
        # Check that the existing database had the same entries
        try:
            from radis.misc.config import getDatabankEntries

            for k, v in getDatabankEntries(
                databank_name,
                get_extra_keys=["download_url", "wavenumber_min", "wavenumber_max"],
            ).items():
                assert dict_entries[k] == v
            # TODO @dev : replace once we have configfile as JSON (https://github.com/radis/radis/issues/167)
        except AssertionError:
            raise DatabaseAlreadyExists(
                f"{databank_name} already exists in your ~/.radis config file. "
                + "Remove it from your config file, or choose a different name "
                + "for the downloaded database with `fetch_hitemp(databank_name=...)`"
            ) from err
        else:  # no other error raised
            if verbose:
                print(
                    f"{databank_name} already registered in ~/.radis config file, with the same parameters."
                )

    return hdf2df(
        output,
        isotope=isotope,
        load_wavenum_min=load_wavenum_min,
        load_wavenum_max=load_wavenum_max,
        verbose=verbose,
    )


#%%
def get_last(b):
    """Get non-empty lines of a chunk b, parsing the bytes."""
    element_length = np.vectorize(lambda x: len(x.__str__()))(b)
    non_zero = element_length > element_length[-1]
    threshold = non_zero.argmin() - 1
    assert (non_zero[: threshold + 1] == 1).all()
    assert (non_zero[threshold + 1 :] == 0).all()
    return b[non_zero]


#%%

if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_hitemp.py"]))
