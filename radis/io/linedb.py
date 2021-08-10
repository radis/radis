# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 19:42:51 2021

@author: erwan
"""
import os
from io import BytesIO
from os.path import abspath, exists, splitext
from zipfile import ZipFile

from radis.db import MOLECULES_LIST_NONEQUILIBRIUM
from radis.misc.config import addDatabankEntries, getDatabankEntries, getDatabankList
from radis.misc.warning import DatabaseAlreadyExists, DeprecatedFileWarning

try:
    from .hdf5 import hdf2df
except ImportError:
    from radis.io.hdf5 import hdf2df

from datetime import date
from os.path import join

import numpy as np
import pandas as pd
from dateutil.parser import parse as parse_date
from numpy import DataSource

LAST_VALID_DATE = (
    "01 Jan 2010"  # set to a later date to force re-download of all HITEMP databases
)

# Add a zip opener to the datasource _file_openers
def open_zip(zipname, mode="r", encoding=None, newline=None):
    output = BytesIO()
    with ZipFile(zipname, mode[0]) as myzip:
        fnames = myzip.namelist()
        for fname in fnames:
            output.write(myzip.read(fname))
    output.seek(0)
    return output


np.lib._datasource._file_openers._file_openers[".zip"] = open_zip


class DatabaseManager(object):
    """Line Database parser

    Parameters
    ----------
    name: str
        database name as registered in ~/radis.json
    molecule: str
    local_databases: str
        path to local database
    """

    # Should be as a close as possible to the content of the corresponding ~/radis.json entry
    # Essentially a FileManager

    def __init__(self, name, molecule, local_databases, verbose=False):

        self.name = name
        self.molecule = molecule
        self.local_databases = local_databases
        self.downloadable = False  # by default

        self.ds = DataSource(join(self.local_databases, "downloads"))

        self.verbose = verbose

    def get_filenames(self):
        verbose = self.verbose
        local_databases = self.local_databases

        # First; checking if the database is registered in radis.json
        if self.is_registered():
            entries = getDatabankEntries(self.name)
            try:
                assert "download_url" in entries
                assert "download_date" in entries
                assert parse_date(entries["download_date"]) > parse_date(
                    LAST_VALID_DATE
                )
            except AssertionError as err:
                raise DeprecatedFileWarning(
                    "Database file {0} not valid anymore"
                ) from err
            else:
                local_files = entries["path"]
            urlnames = None

        elif self.is_downloadable():
            # local_files = self.fetch_filenames()
            urlnames = self.fetch_urlnames()
            local_fnames = [
                (
                    splitext(splitext(url.split("/")[-1])[0])[
                        0
                    ]  # twice to remove .par.bz2
                    + ".h5"
                )
                for url in urlnames
            ]

            try:
                os.mkdir(local_databases)
            except OSError:
                pass
            else:
                if verbose:
                    print("Created folder :", local_databases)

            local_files = [
                abspath(
                    join(
                        local_databases,
                        self.molecule + "-" + local_fname,
                    )
                )
                for local_fname in local_fnames
            ]

        else:
            raise NotImplementedError

        return local_files, urlnames

    def remove_local_files(self, local_files):
        # Delete existing HDF5 file
        for local_file in local_files:
            if exists(local_file):
                if self.verbose:
                    print("Removing existing file ", local_file)
                    # TODO: also clean the getDatabankList? Todo once it is in JSON format. https://github.com/radis/radis/issues/167
                os.remove(local_file)

    def is_downloadable(self):
        return self.downloadable

    def is_registered(self):
        return self.name in getDatabankList()

    def fetch_filenames(self):
        raise NotImplementedError

    def get_today(self):
        return date.today().strftime("%d %b %Y")

    def register(self, local_files, urlname):
        download_date = self.get_today()
        verbose = self.verbose

        # Add database to  ~/radis.json
        return register_database(
            self.name,
            local_files,
            self.molecule,
            self.wmin,
            self.wmax,
            download_date,
            urlname,
            verbose,
        )

    def download_and_parse(self, urlnames, local_files):
        all_local_files, _ = self.get_filenames()

        verbose = self.verbose
        molecule = self.molecule

        # self.parse_to_local_file(ds, urlname, local_file, molecule)

        from time import time

        t0 = time()
        pbar_Ntot_estimate_factor = None
        if len(urlnames) != len(all_local_files):
            # we're only downloading a part of the database
            # expected number of lines is approximately Ntot * N_files/N_files_total
            # (this is just to give the user an expected download & parse time)
            pbar_Ntot_estimate_factor = len(urlnames) / len(all_local_files)
        else:
            pbar_Ntot_estimate_factor = None
        Nlines_total = 0
        Ndownload = 1
        Ntotal_downloads = len(local_files)
        for urlname, local_file in zip(urlnames, local_files):

            if verbose:
                inputf = urlname.split("/")[-1]
                print(
                    f"Downloading {inputf} for {molecule} ({Ndownload}/{Ntotal_downloads})."
                )

            Nlines = self.parse_to_local_file(
                self.ds,
                urlname,
                local_file,
                pbar_t0=time() - t0,
                pbar_Ntot_estimate_factor=pbar_Ntot_estimate_factor,
                pbar_Nlines_already=Nlines_total,
            )
            Ndownload += 1
            Nlines_total += Nlines

    def clean_download_files(self, urlnames):
        # Fully unzipped (and working, as it was reloaded): clean
        if not isinstance(urlnames, list):
            urlnames = [urlnames]
        for urlname in urlnames:
            os.remove(self.ds._findfile(urlname))
            if self.verbose >= 3:
                from radis.misc.printer import printg

                printg(f"... removed downloaded cache file for {urlname}")

    def load(
        self,
        local_files,
        isotope,
        load_wavenum_min,
        load_wavenum_max,
        engine="pytables",
    ):
        if engine == "pytables":
            df_all = []
            for local_file in local_files:
                df_all.append(
                    hdf2df(
                        local_file,
                        isotope=isotope,
                        load_wavenum_min=load_wavenum_min,
                        load_wavenum_max=load_wavenum_max,
                        verbose=self.verbose,
                    )
                )
            return pd.concat(df_all)

        elif engine in ["vaex", "h5py"]:
            raise NotImplementedError
        else:
            raise ValueError(engine)

    def get_nrows(self, local_file, engine="pytables"):
        if engine == "pytables":
            with pd.HDFStore(local_file, "r") as store:
                nrows = store.get_storer("df").nrows

        elif engine in ["vaex", "h5py"]:
            raise NotImplementedError
        else:
            raise ValueError(engine)
        return nrows


def register_database(
    databank_name, path_list, molecule, wmin, wmax, download_date, urlname, verbose
):
    """Add to registered databases in RADIS config file.

    If database exists, assert it has the same entries.

    Parameters
    ----------
    databank_name: str
        name of the database in :ref:`~/radis.json config file <label_lbl_config_file>`

    Other Parameters
    ----------------
    verbose: bool

    Returns
    -------
    None:
        adds to :ref:`~/radis.json <label_lbl_config_file>` with all the input
        parameters. Also adds ::

            format : "hitemp-radisdb"
            parfuncfmt : "hapi"   # TIPS-2017 for equilibrium partition functions

        And if the molecule is in :py:attr:`~radis.db.MOLECULES_LIST_NONEQUILIBRIUM`::

            levelsfmt : "radis"   # use RADIS spectroscopic constants for rovibrational energies (nonequilibrium)

    """
    dict_entries = {
        "info": f"HITEMP {molecule} lines ({wmin:.1f}-{wmax:.1f} cm-1) with TIPS-2017 (through HAPI) for partition functions",
        "path": path_list,
        "format": "hitemp-radisdb",
        "parfuncfmt": "hapi",
        "wavenumber_min": f"{wmin}",
        "wavenumber_max": f"{wmax}",
        "download_date": download_date,
        "download_url": urlname,
    }

    if molecule in MOLECULES_LIST_NONEQUILIBRIUM:
        dict_entries[
            "info"
        ] += " and RADIS spectroscopic constants for rovibrational energies (nonequilibrium)"
        dict_entries["levelsfmt"] = "radis"

    # Register database in ~/radis.json to be able to use it with load_databank()
    try:
        addDatabankEntries(databank_name, dict_entries)
    except DatabaseAlreadyExists as err:
        # Check that the existing database had the same entries
        try:
            from radis.misc.config import getDatabankEntries

            for k, v in getDatabankEntries(databank_name).items():
                if k == "download_date":
                    continue
                assert dict_entries[k] == v
            # TODO @dev : replace once we have configfile as JSON (https://github.com/radis/radis/issues/167 )
        except AssertionError:
            raise DatabaseAlreadyExists(
                f"{databank_name} already exists in your ~/radis.json config file, "
                + f"with different key `{k}` : `{v}` (~/radis.json) ≠ `{dict_entries[k]}` (new). "
                + "If you're sure of what you're doing, fix the registered database in ~/radis.json. "
                + "Else, remove it from your config file, or choose a different name "
                + "for the downloaded database with `fetch_hitemp(databank_name=...)`, "
                + "and restart."
            ) from err
        else:  # no other error raised
            if verbose:
                print(
                    f"{databank_name} already registered in ~/radis.json config file, with the same parameters."
                )
