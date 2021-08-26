# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 19:42:51 2021

@author: erwan
"""
import os
import shutil
from io import BytesIO
from os.path import abspath, exists, splitext
from zipfile import ZipFile

from radis.misc.config import addDatabankEntries, getDatabankEntries, getDatabankList
from radis.misc.printer import printr
from radis.misc.warning import DatabaseAlreadyExists, DeprecatedFileWarning

try:
    from .cache_files import check_not_deprecated
    from .hdf5 import HDF5Manager, hdf2df
except ImportError:
    from radis.io.hdf5 import hdf2df, HDF5Manager
    from radis.io.cache_files import check_not_deprecated

from datetime import date
from os.path import join

import numpy as np
import pandas as pd
from dateutil.parser import parse as parse_date
from numpy import DataSource

LAST_VALID_DATE = (
    "01 Jan 2010"  # set to a later date to force re-download of all databases
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
        self.format = ""

        self.tempdir = join(self.local_databases, "downloads(can_be_deleted)")
        self.ds = DataSource(self.tempdir)

        self.verbose = verbose

    def get_filenames(self, engine):
        """Get names of all files in the database (even if not downloaded yet)

        See Also
        --------
        :py:meth:`~radis.io.linedb.get_files_to_download`"""
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
                    "Database entry {self.name} not valid anymore. See above."
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

        if engine == "vaex":
            local_files = [fname.replace(".h5", ".hdf5") for fname in local_files]

        return local_files, urlnames

    def fetch_urlnames(self) -> list:
        """ "This function should be overwritten by the DatabaseManager subclass

        Returns
        -------
        list: list of urlnames

        See for instance :py:meth:`radis.io.hitemp.HITEMPDatabaseManager"""

        raise NotImplementedError(
            "This function should be overwritten by the DatabaseManager subclass"
        )

    def check_deprecated_files(self, local_files, engine, remove=True):
        """Check metadata of files and remove the deprecated ones

        Unless remove=False: Then raise an error"""
        verbose = self.verbose
        for local_file in local_files:
            try:
                check_not_deprecated(
                    local_file,
                    metadata_is={},
                    metadata_keys_contain=["wavenumber_min", "wavenumber_max"],
                    engine=engine,
                )
            except DeprecatedFileWarning as err:
                if not remove:
                    raise err
                else:  # delete file to regenerate it in the end of the script
                    if verbose:
                        printr(
                            "File {0} deprecated:\n{1}\nDeleting it!".format(
                                local_file, str(err)
                            )
                        )
                    os.remove(local_file)

    def get_existing_files(self, files):
        """Return files that exist among ``files``

        See Also
        --------
        :py:meth:`~radis.io.linedb.get_filenames`"""
        return [k for k in files if exists(k)]

    def get_missing_files(self, files):
        """Return files that do not exist among ``files``

        Note : in 'vaex' mode; if "FILE.hdf5" does not exist
        but "FILE.h5" does (a likely 'pytables' file), does
        not consider it missing so it can be converted
        automatically

        See Also
        --------
        :py:meth:`~radis.io.linedb.get_filenames`"""
        return [k for k in files if not exists(k)]

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
        return date.today().strftime(
            "%Y-%m-%d"
        )  # only numbers avoids locale problems. Old : ("%d %b %Y")

    def register(self, dict_entries):

        # Add database to  ~/radis.json
        return register_database(
            self.name,
            dict_entries,
            self.verbose,
        )

    def get_hdf5_manager(self, engine):
        return HDF5Manager(engine=engine)

    def download_and_parse(self, urlnames, local_files, engine="pytables"):
        all_local_files, _ = self.get_filenames(engine)

        verbose = self.verbose
        molecule = self.molecule

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

            # Check we can open the file, give the path if there is an error
            try:
                self.ds.open(urlname)
            except Exception as err:
                raise OSError(
                    f"Problem opening : {self.ds._findfile(urlname)}. See above. You may want to delete the downloaded file."
                ) from err

            # try:
            Nlines = self.parse_to_local_file(
                self.ds,
                urlname,
                local_file,
                pbar_t0=time() - t0,
                pbar_Ntot_estimate_factor=pbar_Ntot_estimate_factor,
                pbar_Nlines_already=Nlines_total,
                pbar_last=(Ndownload == Ntotal_downloads),
                engine=engine,
            )
            # except Exception as err:
            #     raise IOError("Problem parsing `{0}`. Check the error above. It may arise if the file wasn't properly downloaded. Try to delete it".format(self.ds._findfile(urlname))) from err

            Ndownload += 1
            Nlines_total += Nlines

    def clean_download_files(self):
        """Fully unzipped (and working, as it was reloaded): clean files

        Note : we do not let np.DataSource clean its own tempfiles; as
        they may be downloaded but not fully parsed / registered."""
        if exists(self.tempdir):
            try:
                shutil.rmtree(self.tempdir)
            except PermissionError as err:
                if self.verbose >= 3:
                    from radis.misc.printer import printr

                    printr(
                        f"... couldnt delete downloaded cache files {self.tempdir}: {str(err)}"
                    )
            else:
                if self.verbose >= 3:
                    from radis.misc.printer import printg

                    printg(f"... removed downloaded cache files for {self.tempdir}")

    def load(
        self,
        local_files,
        isotope,
        columns,
        load_wavenum_min,
        load_wavenum_max,
        engine="pytables",
    ):
        """
        Other Parameters
        ----------------
        columns: list of str
            list of columns to load. If ``None``, returns all columns in the file.
        """
        if engine == "pytables":
            df_all = []
            for local_file in local_files:
                df_all.append(
                    hdf2df(
                        local_file,
                        columns=columns,
                        isotope=isotope,
                        load_wavenum_min=load_wavenum_min,
                        load_wavenum_max=load_wavenum_max,
                        verbose=self.verbose,
                        engine="pytables",
                    )
                )
            return pd.concat(df_all)

        elif engine == "vaex":
            # vaex can open several files at the same time:
            return hdf2df(
                local_files,
                columns=columns,
                isotope=isotope,
                load_wavenum_min=load_wavenum_min,
                load_wavenum_max=load_wavenum_max,
                verbose=self.verbose,
                engine="vaex",
            )
        else:
            raise NotImplementedError(engine)

    def plot(self, local_files, isotope, wavenum_min, wavenum_max):
        """Convenience function to plot linestrengths of the database"""

        df = self.load(
            local_files, isotope, wavenum_min, wavenum_max, columns=["wav", "int"]
        )
        df.plot("wav", "int")

    def get_nrows(self, local_file, engine="pytables"):
        if engine == "pytables":
            with pd.HDFStore(local_file, "r") as store:
                nrows = store.get_storer("df").nrows

        elif engine in ["vaex", "h5py"]:
            raise NotImplementedError
        else:
            raise ValueError(engine)
        return nrows


def register_database(databank_name, dict_entries, verbose):
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
