# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 19:42:51 2021

@author: erwan
"""
import os
import shutil
from io import BytesIO
from os.path import abspath, dirname, exists, expanduser, join, split, splitext
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

import numpy as np
import pandas as pd
from dateutil.parser import parse as parse_date
from joblib import Parallel, delayed
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
    engine: 'vaex', 'pytables', 'h5py', or 'default'
        memory-mapping library to use with this database. If 'default' use
        the value from ~/radis.json

    Other Parameters
    ----------------
    *input for :class:`~joblib.parallel.Parallel` loading of database*
    parallel: bool
        if ``True``, use parallel loading.
        Default ``True``.
    nJobs: int
        Number of processors to use to load a database (useful for big
        databases). BE CAREFUL, no check is done on processor use prior
        to the execution ! Default ``-2``: use all but 1 processors.
        Use ``1`` for single processor.
    batch_size: int or ``'auto'``
        The number of atomic tasks to dispatch at once to each
        worker. When individual evaluations are very fast, dispatching
        calls to workers can be slower than sequential computation because
        of the overhead. Batching fast computations together can mitigate
        this. Default: ``'auto'``
    More information in :class:`joblib.parallel.Parallel`

    """

    # Should be as a close as possible to the content of the corresponding ~/radis.json entry
    # Essentially a FileManager

    def __init__(
        self,
        name,
        molecule,
        local_databases,
        engine,
        verbose=False,
        parallel=True,
        nJobs=-2,
        batch_size="auto",
    ):
        from os import environ

        if engine == "default":
            from radis import config

            engine = config["MEMORY_MAPPING_ENGINE"]  # 'pytables', 'vaex', 'feather'
            # Quick fix for #401
            if engine == "auto":
                # "auto" uses "vaex" in most cases unless you're using the Spyder IDE (where it may result in freezes).
                # see https://github.com/spyder-ide/spyder/issues/16183.
                # and https://github.com/radis/radis/issues/401
                if any("SPYDER" in name for name in environ):
                    engine = "pytables"  # for HITRAN and HITEMP databases
                    if verbose >= 3:
                        print(
                            f"Spyder IDE detected. Memory-mapping-engine set to '{engine}' (less powerful than 'vaex' but Spyder user experience freezes). See https://github.com/spyder-ide/spyder/issues/16183. Change this behavior by setting the radis.config['MEMORY_MAPPING_ENGINE'] key"
                        )
                # temp fix for vaex not building on RTD
                # see https://github.com/radis/radis/issues/404
                elif any("READTHEDOCS" in name for name in environ):
                    engine = "pytables"  # for HITRAN and HITEMP databases
                    if verbose >= 3:
                        print(
                            f"ReadTheDocs environment detected. Memory-mapping-engine set to '{engine}'. See https://github.com/radis/radis/issues/404"
                        )
                else:
                    engine = "vaex"

        # vaex processes are stuck if ran from Spyder. See https://github.com/spyder-ide/spyder/issues/16183
        if engine == "vaex" and any("SPYDER" in name for name in environ):
            from radis.misc.log import printwarn

            printwarn(
                "Spyder IDE detected while using memory_mapping_engine='vaex'.\nVaex is the fastest way to read database files in RADIS, but Vaex processes may be stuck if ran from Spyder. See https://github.com/spyder-ide/spyder/issues/16183. Quick fix: starting a new console releases the lock, usually for the rest of your session. You may consider using another IDE, or using a different `memory_mapping_engine` such as 'pytables' or 'feather'. You can change the engine in Spectrum.fetch_databank() calls, or globally by setting the 'MEMORY_MAPPING_ENGINE' key in your ~/radis.json \n"
            )

        self.name = name
        self.molecule = molecule
        self.local_databases = local_databases
        # create folder if needed
        if not exists(local_databases):
            from radis.misc.basics import make_folders

            make_folders(*split(abspath(dirname(local_databases))))
            make_folders(*split(abspath(local_databases)))

        if self.is_registered():
            registered_paths = getDatabankEntries(self.name)["path"]
            for registered_path in registered_paths:
                if (
                    not abspath(expanduser(registered_path))
                    .lower()
                    .startswith(abspath(expanduser(local_databases).lower()))
                ):  # TODO: replace with pathlib
                    raise ValueError(
                        f"Databank `{self.name}` is already registered in radis.json but the declared path ({registered_path}) is not in the expected local databases folder ({local_databases}). Please fix/delete the radis.json entry, or change the default local databases path entry 'DEFAULT_DOWNLOAD_PATH' in `radis.config` or ~/radis.json"
                    )

        self.downloadable = False  # by default
        self.format = ""
        self.engine = engine

        self.tempdir = join(self.local_databases, "downloads__can_be_deleted")
        self.ds = DataSource(self.tempdir)

        self.verbose = verbose

        self.parallel = parallel
        self.nJobs = nJobs
        self.batch_size = batch_size
        self.minimum_nfiles = (
            4  #: type: int. If there are less files, don't use parallel mode.
        )

    def get_filenames(self):
        """Get names of all files in the database (even if not downloaded yet)

        See Also
        --------
        :py:meth:`~radis.io.linedb.get_files_to_download`"""
        verbose = self.verbose
        local_databases = self.local_databases
        engine = self.engine

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
                    f"Database entry {self.name} not valid anymore. See above which keys are missing/wrong. Update or delete the entry in your ~/radis.json"
                ) from err
            else:
                local_files = entries["path"]
            urlnames = None

            # Check that local files are the one we expect :
            for f in local_files:
                if (
                    not abspath(expanduser(f))
                    .lower()
                    .startswith(abspath(expanduser(local_databases)).lower())
                ):
                    raise ValueError(
                        f"Database {self.name} is inconsistent : it should be stored in {local_databases} but files registered in ~/radis.json contains {f}. Please fix or delete the ~/radis.json entry."
                    )

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

    def check_deprecated_files(self, local_files, auto_remove=True):
        """Check metadata of files and remove the deprecated ones

        Unless auto_remove=False: Then raise an error"""
        verbose = self.verbose
        engine = self.engine
        for local_file in local_files:
            try:
                check_not_deprecated(
                    local_file,
                    metadata_is={},
                    metadata_keys_contain=["wavenumber_min", "wavenumber_max"],
                    engine=engine,
                )
            except DeprecatedFileWarning as err:
                if not auto_remove:
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

    def get_hdf5_manager(self):
        return HDF5Manager(engine=self.engine)

    def download_and_parse(self, urlnames, local_files):
        all_local_files, _ = self.get_filenames()

        verbose = self.verbose
        molecule = self.molecule
        parallel = self.parallel

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
        Ntotal_downloads = len(local_files)

        def download_and_parse_one_file(urlname, local_file, Ndownload):
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
                pbar_active=(not parallel),
                pbar_t0=time() - t0,
                pbar_Ntot_estimate_factor=pbar_Ntot_estimate_factor,
                pbar_Nlines_already=Nlines_total,
                pbar_last=(Ndownload == Ntotal_downloads),
            )
            # except Exception as err:
            #     raise IOError("Problem parsing `{0}`. Check the error above. It may arise if the file wasn't properly downloaded. Try to delete it".format(self.ds._findfile(urlname))) from err

            return Nlines

        if parallel and len(local_files) > self.minimum_nfiles:
            nJobs = self.nJobs
            batch_size = self.batch_size
            if self.verbose:
                print(
                    f"Downloading and parsing {urlnames} to {local_files} "
                    + f"({len(local_files)}) files), in parallel ({nJobs} jobs)"
                )
            Nlines_total = sum(
                Parallel(n_jobs=nJobs, batch_size=batch_size, verbose=self.verbose)(
                    delayed(download_and_parse_one_file)(urlname, local_file, Ndownload)
                    for urlname, local_file, Ndownload in zip(
                        urlnames, local_files, range(1, len(local_files) + 1)
                    )
                )
            )
        else:
            for urlname, local_file, Ndownload in zip(
                urlnames, local_files, range(1, len(local_files) + 1)
            ):
                download_and_parse_one_file(urlname, local_file, Ndownload)

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
    ):
        """
        Other Parameters
        ----------------
        columns: list of str
            list of columns to load. If ``None``, returns all columns in the file.
        """
        engine = self.engine
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
                        engine=engine,
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
                engine=engine,
            )
        else:
            raise NotImplementedError(engine)

    def plot(self, local_files, isotope, wavenum_min, wavenum_max):
        """Convenience function to plot linestrengths of the database"""

        df = self.load(
            local_files, isotope, wavenum_min, wavenum_max, columns=["wav", "int"]
        )
        df.plot("wav", "int")

    def get_nrows(self, local_file):
        """ Get number of rows (without loading all DataFrame)"""
        engine = self.engine
        local_file = expanduser(local_file)
        if engine == "pytables":
            with pd.HDFStore(local_file, "r") as store:
                nrows = store.get_storer("df").nrows

        elif engine == "vaex":
            import vaex

            # by default vaex does not load everything
            df = vaex.open(local_file)
            nrows = len(df)
            df.close()

        elif engine in ["h5py"]:
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
