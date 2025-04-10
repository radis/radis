# -*- coding: utf-8 -*-
"""

"""
import os
import shutil
from os.path import abspath, dirname, exists, expanduser, join, split, splitext

try:
    from ..misc.config import addDatabankEntries, getDatabankEntries, getDatabankList
    from ..misc.printer import printr
    from ..misc.utils import NotInstalled, not_installed_vaex_args
    from ..misc.warning import DatabaseAlreadyExists, DeprecatedFileWarning
    from .cache_files import check_not_deprecated
    from .hdf5 import DataFileManager
except ImportError:
    if __name__ == "__main__":  # running from this file, as a script
        from radis.api.cache_files import check_not_deprecated
        from radis.api.hdf5 import DataFileManager
        from radis.misc.config import (
            addDatabankEntries,
            getDatabankEntries,
            getDatabankList,
        )
        from radis.misc.printer import printr
        from radis.misc.utils import NotInstalled, not_installed_vaex_args
        from radis.misc.warning import DatabaseAlreadyExists, DeprecatedFileWarning
    else:
        raise

try:
    import vaex
except ImportError:
    vaex = NotInstalled(*not_installed_vaex_args)

from datetime import date

import pandas as pd
import requests
from dateutil.parser import parse as parse_date
from joblib import Parallel, delayed

LAST_VALID_DATE = (
    "01 Jan 2010"  # set to a later date to force re-download of all databases
)


def get_auto_MEMORY_MAPPING_ENGINE():
    """see https://github.com/radis/radis/issues/653

    Use Vaex by default if it exists (only Python <= 3.11 as of June 2024) ,
    else use PyTables"""
    try:
        import vaex

        vaex
    except ImportError:
        return "pytables"
    else:
        return "vaex"


# Class to mimic numpy.DataSource behavior
class RequestsFileOpener:
    """Simple file opener class that mimics the behavior of numpy.DataSource
    using the requests pip package.
    """

    def __init__(self, file_path):
        self.file_path = file_path

    def open(self, url_or_path=None):
        """Open the file for reading"""
        if self.file_path.endswith(".zip"):
            from io import BytesIO
            from zipfile import ZipFile

            output = BytesIO()
            with ZipFile(self.file_path, "r") as myzip:
                fnames = myzip.namelist()
                for fname in fnames:
                    output.write(myzip.read(fname))
            output.seek(0)
            return output
        elif self.file_path.endswith(".bz2"):
            import bz2
            from io import BytesIO

            # Read compressed file
            with open(self.file_path, "rb") as f:
                compressed_data = f.read()

            # Decompress the data
            decompressed_data = bz2.decompress(compressed_data)

            # Create BytesIO object with decompressed data
            output = BytesIO(decompressed_data)
            return output
        else:
            return open(self.file_path, "rb")

    def abspath(self, url_or_path=None):
        """Return the absolute path of the file"""
        return self.file_path


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
        extra_params=None,
        verbose=False,
        parallel=True,
        nJobs=-2,
        batch_size="auto",
    ):
        if engine == "default":
            from radis import config

            engine = config["MEMORY_MAPPING_ENGINE"]  # 'pytables', 'vaex', 'feather'
            if engine == "auto":
                engine = get_auto_MEMORY_MAPPING_ENGINE()

        self.name = name
        self.molecule = molecule
        self.local_databases = local_databases
        self.extra_params = extra_params
        # create folder if needed
        if not exists(local_databases):
            from radis.misc.basics import make_folders

            make_folders(*split(abspath(dirname(expanduser(local_databases)))))
            make_folders(*split(abspath(expanduser(local_databases))))

        if self.is_registered():
            registered_paths = getDatabankEntries(self.name)["path"]
            for registered_path in registered_paths:
                registered_path_abspath = abspath(expanduser(registered_path)).lower()
                local_databases_abspath = abspath(expanduser(local_databases).lower())
                if not registered_path_abspath.startswith(
                    local_databases_abspath
                ):  # TODO: replace with pathlib
                    raise ValueError(
                        f"Databank `{self.name}` is already registered in radis.json but the declared path ({registered_path_abspath}) is not in the expected local databases folder ({local_databases_abspath}). Please fix/delete the radis.json entry, change the `databank_name`, or change the default local databases path entry 'DEFAULT_DOWNLOAD_PATH' in `radis.config` or ~/radis.json"
                    )

        self.downloadable = False  # by default
        self.format = ""
        self.engine = engine

        self.tempdir = join(self.local_databases, "downloads__can_be_deleted")
        # Create the temp directory if it doesn't exist
        if not exists(self.tempdir):
            os.makedirs(self.tempdir, exist_ok=True)

        self.verbose = verbose

        self.parallel = parallel
        self.nJobs = nJobs
        self.batch_size = batch_size
        self.minimum_nfiles = (
            4  #: type: int. If there are less files, don't use parallel mode.
        )

    def get_filenames(self, return_reg_urls=False):
        """Get names of all files in the database (even if not downloaded yet)

        Parameters
        ----------
        return_reg_urls: (boolean)
            When the database is registered, whether to return the registered urls (``True``) or ``None`` (``False``)

        See Also
        --------
        :py:meth:`~radis.api.dbmanager.DatabaseManager.get_files_to_download`"""
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
                    f"Database entry {self.name} not valid anymore. See above which keys are missing/wrong. Update or delete the entry in your ~/radis.json"
                ) from err
            else:
                local_files = entries["path"]
            if return_reg_urls:
                urlnames = entries["download_url"]
            else:
                urlnames = []

            # Check that local files are the one we expect :
            for f in local_files:
                local_file_abspath = abspath(expanduser(f)).lower()
                local_databases_abspath = abspath(expanduser(local_databases)).lower()
                if not local_file_abspath.startswith(local_databases_abspath):
                    raise ValueError(
                        f"Database {self.name} is inconsistent : it should be stored in {local_databases_abspath} but files registered in ~/radis.json contains {local_file_abspath}. Please fix or delete the ~/radis.json entry."
                    )

        elif self.is_downloadable():
            urlnames = self.fetch_urlnames()
            local_files = self.fetch_filenames(urlnames, local_databases)
        else:
            raise NotImplementedError

        local_files = [expanduser(f) for f in local_files]

        return local_files, urlnames

    def fetch_urlnames(self) -> list:
        """ "This function should be overwritten by the DatabaseManager subclass

        Returns
        -------
        list: list of urlnames

        See for instance :py:class:`radis.api.hitempapi.HITEMPDatabaseManager"""

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
        :py:meth:`~radis.api.dbmanager.DatabaseManager.get_filenames`"""
        return [k for k in files if exists(k)]

    def get_missing_files(self, files):
        """Return files that do not exist among ``files``

        Note : in 'vaex' mode; if "FILE.hdf5" does not exist
        but "FILE.h5" does (a likely 'pytables' file), does
        not consider it missing so it can be converted
        automatically

        See Also
        --------
        :py:meth:`~radis.api.dbmanager.DatabaseManager.get_filenames`"""
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

    def fetch_filenames(self, urlnames, local_databases):
        engine = self.engine
        verbose = self.verbose

        local_fnames = [
            (
                splitext(splitext(url.split("/")[-1])[0])[0]  # twice to remove .par.bz2
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

        if engine == "vaex":
            local_files = [fname.replace(".h5", ".hdf5") for fname in local_files]

        local_files = [expanduser(f) for f in local_files]

        return local_files

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

    def get_datafile_manager(self, engine=None):
        if engine is None:
            engine = self.engine
        return DataFileManager(engine=engine)

    def download_and_parse(self, urlnames, local_files, N_files_total=None):
        if N_files_total is None:
            all_local_files, _ = self.get_filenames()
            N_files_total = len(all_local_files)

        verbose = self.verbose
        molecule = self.molecule
        parallel = self.parallel

        from time import time

        t0 = time()
        pbar_Ntot_estimate_factor = None
        if len(urlnames) != N_files_total:
            # we're only downloading a part of the database
            # expected number of lines is approximately Ntot * N_files/N_files_total
            # (this is just to give the user an expected download & parse time)
            pbar_Ntot_estimate_factor = len(urlnames) / N_files_total
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

            # Download file with requests
            try:
                # Get session from HITEMP API
                from radis.api.hitempapi import login_to_hitran

                session = login_to_hitran(verbose=verbose)

                # Set headers to indicate we want the actual file
                headers = {
                    "Accept": "application/zip, application/octet-stream",
                    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
                }

                # First check if we can access the file
                head_response = session.head(
                    urlname, headers=headers, allow_redirects=True
                )
                if head_response.status_code != 200:
                    raise OSError(f"Failed to access file: {head_response.status_code}")

                # Check if we got redirected to login page
                if "text/html" in head_response.headers.get("content-type", "").lower():
                    raise OSError(
                        "Got HTML response instead of file. Please ensure you're logged in and have access to the file."
                    )

                # Now download the file
                response = session.get(
                    urlname, headers=headers, stream=True, allow_redirects=True
                )
                response.raise_for_status()  # Raise an error if request fails

                # Create a temporary file to store the downloaded content
                temp_file_path = join(self.tempdir, urlname.split("/")[-1])
                with open(temp_file_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:  # filter out keep-alive new chunks
                            f.write(chunk)

                # Create an opener object using the global RequestsFileOpener class
                opener = RequestsFileOpener(temp_file_path)

                # Pass the opener to the parse function (maintaining backward compatibility)
                Nlines = self.parse_to_local_file(
                    opener,
                    urlname,
                    local_file,
                    pbar_active=(not parallel),
                    pbar_t0=time() - t0,
                    pbar_Ntot_estimate_factor=pbar_Ntot_estimate_factor,
                    pbar_Nlines_already=Nlines_total,
                    pbar_last=(Ndownload == Ntotal_downloads),
                )
            except requests.RequestException as err:
                raise OSError(
                    f"Problem downloading: {urlname}. Error: {str(err)}"
                ) from err
            except Exception as err:
                raise IOError(
                    f"Problem parsing downloaded file from {urlname}. Check the error above. It may arise if the file wasn't properly downloaded. Try to delete it."
                ) from err

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

    def parse_to_local_file(
        self,
        opener,
        urlname,
        local_file,
        pbar_active=True,
        pbar_t0=0,
        pbar_Ntot_estimate_factor=None,
        pbar_Nlines_already=0,
        pbar_last=True,
    ) -> list:
        """This function should be overwritten by the DatabaseManager subclass

        Uncompress ``urlname`` into ``local_file``.
        Also add metadata

        Returns
        -------
        list: list of urlnames

        See for instance :py:class:`radis.api.hitempapi.HITEMPDatabaseManager"""

        raise NotImplementedError(
            "This function should be overwritten by the DatabaseManager subclass"
        )

    def clean_download_files(self):
        """Fully unzipped (and working, as it was reloaded): clean files"""
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
        columns=None,
        lower_bound=[],
        upper_bound=[],
        within=[],
        output="pandas",
    ):
        """
        Other Parameters
        ----------------
        columns: list of str
            list of columns to load. If ``None``, returns all columns in the file.
        output: 'pandas', 'vaex', 'jax'
            format of the output DataFrame. If ``'jax'``, returns a dictionary of
            jax arrays.
        lower_bound: list of tuples [(column, lower_bound), etc.]
            ::

                lower_bound =[("wav", load_wavenum_min)]
        upper_bound_bound: list of tuples [(column, upper_bound), etc.]
            ::

                upper_bound=[("wav", load_wavenum_max)]
        within: list of tuples [(column, within_list), etc.]
            ::

                within=[("iso", isotope.split(","))]
        """
        engine = self.engine
        mgr = self.get_datafile_manager()
        if engine in ["pytables", "feather"]:
            df_all = []
            for local_file in local_files:
                df_all.append(
                    mgr.load(
                        local_file,
                        columns=columns,
                        lower_bound=lower_bound,
                        upper_bound=upper_bound,
                        within=within,
                        output=output,
                    )
                )
            return pd.concat(df_all)

        elif engine == "vaex":
            # vaex can open several files at the same time:
            return mgr.load(
                local_files,
                columns=columns,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                within=within,
                output=output,
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
        """Get number of rows (without loading all DataFrame)"""
        engine = self.engine
        local_file = expanduser(local_file)
        if engine == "pytables":
            with pd.HDFStore(local_file, "r") as store:
                nrows = store.get_storer("df").nrows

        elif engine == "vaex":
            # by default vaex does not load everything
            df = vaex.open(local_file)
            nrows = len(df)
            df.close()

        elif engine in ["h5py"]:
            raise NotImplementedError
        else:
            raise ValueError(engine)
        return nrows

    def get_columns(self, local_file):
        """Get all columns using DataFileManager class get_columns function"""
        return self.get_datafile_manager().get_columns(local_file)

    def add_column(self, df, key, value):
        """Create column ``key`` in DataFrame or dictionary ``df`` with value ``value``"""
        from radis.misc.basics import is_number

        if is_number(self.alpha_ref):
            if isinstance(df, vaex.dataframe.DataFrameLocal):
                # see https://github.com/vaexio/vaex/pull/1570
                df[key] = vaex.vconstant(float(value), length=df.length_unfiltered())
            else:
                df[key] = value
        else:
            df[key] = value

    def rename_columns(self, df, rename_dict):
        """Example::

        mdb.rename_columns(df, {"nu_lines":"wav"})
        """
        if not isinstance(vaex, NotInstalled) and isinstance(
            df, vaex.dataframe.DataFrameLocal
        ):
            for k, v in rename_dict.items():
                df.rename(k, v)
        elif isinstance(df, pd.DataFrame):
            df.rename(columns=rename_dict, inplace=True)
        elif isinstance(df, dict):
            for k, v in rename_dict.items():
                df[v] = df.pop(k)


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
