# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 22:40:51 2021

@author: erwan


https://stackoverflow.com/questions/55610891/numpy-load-from-io-bytesio-stream
https://stupidpythonideas.blogspot.com/2014/07/three-ways-to-read-files.html

"""

import os
import re
import urllib.request
from os.path import abspath, exists, expanduser

import numpy as np
import pandas as pd

import radis
from radis.misc.printer import printr

try:
    from .cache_files import check_not_deprecated
    from .hitran import columns_2004, parse_global_quanta, parse_local_quanta
    from .linedb import DatabaseManager
    from .tools import (
        _create_dtype,
        _get_linereturnformat,
        _ndarray2df,
        replace_PQR_with_m101,
    )
except ImportError:  # ran from here
    from radis.io.linedb import DatabaseManager
    from radis.io.cache_files import check_not_deprecated
    from radis.io.hitran import columns_2004, parse_global_quanta, parse_local_quanta
    from radis.io.tools import (
        _create_dtype,
        _get_linereturnformat,
        _ndarray2df,
        replace_PQR_with_m101,
    )

from radis.misc.progress_bar import ProgressBar
from radis.misc.warning import DeprecatedFileWarning

HITEMP_MOLECULES = ["H2O", "CO2", "N2O", "CO", "CH4", "NO", "NO2", "OH"]
DATA_COLUMNS = ["iso", "wav"]
"""
list : only these column names will be searchable directly on disk to
only load certain lines. See :py:func:`~radis.io.hdf5.hdf2df`
"""
# TODO: WIP. Maybe move somewhere else to be also used by HITRAN queries


def keep_only_relevant(
    inputfiles,
    wavenum_min=None,
    wavenum_max=None,
):
    """Parser file names for ``wavenum_format`` (min and max) and only keep
    relevant files if the requested range is ``[wavenum_min, wavenum_max]``

    Returns
    -------
    relevant: list of relevant files
    files_wmin, files_wmax: (float, float) : wavenum min & max of relevant range
    """
    wavenum_format = r"\d{5}"
    relevantfiles = []
    files_wmin = np.inf
    files_wmax = 0
    for file in inputfiles:
        fname_wmin, fname_wmax = re.findall(wavenum_format, file)
        relevant = False
        if wavenum_min is not None and wavenum_max is not None:
            if (float(fname_wmax) >= wavenum_min) and (
                float(fname_wmin) <= wavenum_max
            ):
                relevant = True
        elif wavenum_min is not None:
            if float(fname_wmax) >= wavenum_min:
                relevant = True
        elif wavenum_max is not None:
            if float(fname_wmin) <= wavenum_max:
                relevant = True
        else:
            relevant = True
        if relevant:
            relevantfiles.append(file)
            files_wmin = min(float(fname_wmin), files_wmin)
            files_wmax = max(float(fname_wmax), files_wmax)
    # small checks
    if relevantfiles:
        if wavenum_min:
            assert files_wmin <= wavenum_min
        if wavenum_max:
            assert files_wmax >= wavenum_max
    return relevantfiles, files_wmin, files_wmax


#%%
def get_last(b):
    """Get non-empty lines of a chunk b, parsing the bytes."""
    element_length = np.vectorize(lambda x: len(x.__str__()))(b)
    non_zero = element_length > element_length[-1]
    threshold = non_zero.argmin() - 1
    assert (non_zero[: threshold + 1] == 1).all()
    assert (non_zero[threshold + 1 :] == 0).all()
    return b[non_zero]


class HITEMPDatabaseManager(DatabaseManager):
    def __init__(
        self,
        name,
        molecule,
        local_databases="~/.radisdb/",
        verbose=True,
        chunksize=100000,
    ):
        super().__init__(name, molecule, local_databases, verbose=verbose)
        self.chunksize = chunksize
        self.downloadable = True
        self.base_url = None
        self.Nlines = None

    def fetch_url_and_Nlines(self, hitemp_url="https://hitran.org/hitemp/"):
        """requires connexion"""

        molecule = self.molecule

        if self.base_url is not None and self.Nlines is not None:
            return self.base_url, self.Nlines

        else:

            response = urllib.request.urlopen(hitemp_url)
            text = response.read().decode()
            text = text[
                text.find(
                    '<table id="hitemp-molecules-table" class="selectable-table list-table">'
                ) : text.find("</table>")
            ]
            text = re.sub(r"<!--.+?-->\s*\n", "", text)  # remove commented lines
            html_molecule = re.sub(r"(\d{1})", r"(<sub>\1</sub>)", molecule)
            text = text[
                re.search(
                    "<td>(?:<strong>)?" + html_molecule + "(?:</strong>)?</td>", text
                ).start() :
            ]
            lines = text.splitlines()

            Nlines = int(re.findall(r"(\d+)", lines[3].replace("&nbsp;", ""))[0])
            url = "https://hitran.org" + re.findall(r'href="(.+?)"', lines[7])[0]

            self.base_url, self.Nlines = url, Nlines

        return url, Nlines

    def fetch_urlnames(self):
        """requires connexion"""

        molecule = self.molecule

        if molecule in ["H2O", "CO2"]:

            base_url, Ntotal_lines_expected = self.fetch_url_and_Nlines()
            response = urllib.request.urlopen(base_url)
            response_string = response.read().decode()
            inputfiles = re.findall('href="(\S+.zip)"', response_string)

            # inputfiles = keep_only_relevant(
            #     inputfiles, load_wavenum_min, load_wavenum_max, wavenum_format=r"\d{5}"
            # )
            # if verbose:
            #     print("relevant files:", inputfiles)

            urlnames = [base_url + f for f in inputfiles]

        elif molecule in HITEMP_MOLECULES:
            url, Ntotal_lines_expected = self.fetch_url_and_Nlines()
            urlnames = [url]
        else:
            raise KeyError(
                f"Please choose one of HITEMP molecules : {HITEMP_MOLECULES}. Got '{molecule}'"
            )

        return urlnames

    def parse_to_local_file(
        self,
        opener,
        urlname,
        local_file,
        pbar_t0=0,
        pbar_Ntot_estimate_factor=None,
        pbar_Nlines_already=0,
    ):
        """
        Parameters
        ----------
        opener: an opener with an .open() command
        gfile : file handler. Filename: for info"""

        # Get linereturn (depends on OS, but file may also have been generated
        # on a different OS. Here we simply read the file to find out)
        columns = columns_2004
        chunksize = self.chunksize
        verbose = self.verbose
        molecule = self.molecule

        with opener.open(urlname) as gfile:  # locally downloaded file
            dt = _create_dtype(
                columns, "a2"
            )  # 'a2' allocates space to get \n or \n\r for linereturn character
            b = np.zeros(1, dtype=dt)
            try:
                gfile.readinto(b)
            except EOFError as err:
                raise ValueError(
                    f"End of file while parsing file {opener.abspath(urlname)}. May be due to download error. Delete file ?"
                ) from err
            linereturnformat = _get_linereturnformat(b, columns)

        Nlines = pbar_Nlines_already
        if verbose:
            _, Ntotal_lines_expected = self.fetch_url_and_Nlines()
            if pbar_Ntot_estimate_factor:
                # multiply Ntotal_lines_expected by pbar_Ntot_estimate_factor
                # (accounts for total lines divided in number of files, and
                # not all files downloaded)
                Ntotal_lines_expected = int(
                    Ntotal_lines_expected * pbar_Ntot_estimate_factor
                )
        pb = ProgressBar(N=Ntotal_lines_expected, active=verbose, t0=pbar_t0)
        wmin = np.inf
        wmax = 0

        with opener.open(urlname) as gfile:  # locally downloaded file

            dt = _create_dtype(columns, linereturnformat)
            b = np.zeros(chunksize, dtype=dt)  # receives the HITRAN 160-character data.

            if verbose:
                print(f"Download complete. Parsing {molecule} database to {local_file}")

            with pd.HDFStore(local_file, mode="a", complib="blosc", complevel=9) as f:
                # TODO : implement with engines = 'h5py' too

                for nbytes in iter(lambda: gfile.readinto(b), 0):

                    if not b[-1]:
                        # End of file flag within the chunk (but does not start
                        # with End of file flag) so nbytes != 0
                        b = get_last(b)

                    df = _ndarray2df(b, columns, linereturnformat)

                    # Post-processing :
                    # ... Add local quanta attributes, based on the HITRAN group
                    df = parse_local_quanta(df, molecule)

                    # ... Add global quanta attributes, based on the HITRAN class
                    df = parse_global_quanta(df, molecule)

                    # Switch 'P', 'Q', 'R' to -1, 0, 1
                    if "branch" in df:
                        replace_PQR_with_m101(df)

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
                    if pbar_Ntot_estimate_factor is None:
                        pbar_Ntot_message = f"{Ntotal_lines_expected:,} lines"
                    else:
                        pbar_Ntot_message = (
                            f"~{Ntotal_lines_expected:,} lines (estimate)"
                        )
                    pb.update(
                        Nlines,
                        message=f"Parsed {Nlines:,} / {pbar_Ntot_message}. Wavenumber range {wmin:.2f}-{wmax:.2f} cm-1 is complete.",
                    )
                    # Reinitialize for next read
                    b = np.zeros(
                        chunksize, dtype=dt
                    )  # receives the HITRAN 160-character data.

        # Add metadata
        with pd.HDFStore(local_file, mode="a", complib="blosc", complevel=9) as f:

            f.get_storer("df").attrs.metadata = {
                "wavenumber_min": wmin,
                "wavenumber_max": wmax,
                "download_date": self.get_today(),
                "download_url": urlname,
                "version": radis.__version__,
            }

        return Nlines


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
    clean_cache_files=True,
    return_local_path=False,
):
    """Stream HITEMP file from HITRAN website. Unzip and build a HDF5 file directly.

    Returns a Pandas DataFrame containing all lines.

    Parameters
    ----------
    molecule: `"H2O", "CO2", "N2O", "CO", "CH4", "NO", "NO2", "OH"`
        HITEMP molecule. See https://hitran.org/hitemp/
    local_databases: str
        where to create the RADIS HDF5 files. Default ``"~/.radisdb/"``
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
        if ``True``, use existing HDF5 file. If ``False`` or ``'regen'``, rebuild it..
        Default ``True``.
    verbose: bool
    chunksize: int
        number of lines to process at a same time. Higher is usually faster
        but can create Memory problems and keep the user uninformed of the progress.
    clean_cache_files: bool
        if ``True`` clean downloaded cache files after HDF5 are created.
    return_local_path: bool
        if ``True``, also returns the path of the local database file.

    Returns
    -------
    df: pd.DataFrame
        Line list
        A HDF5 file is also created in ``local_databases`` and referenced
        in the :ref:`RADIS config file <label_lbl_config_file>` with name
        ``databank_name``
    local_path: str
        path of local database file if ``return_local_path``

    Examples
    --------
    ::

        from radis import fetch_hitemp
        df = fetch_hitemp("CO")
        print(df.columns)
        >>> Index(['id', 'iso', 'wav', 'int', 'A', 'airbrd', 'selbrd', 'El', 'Tdpair',
            'Pshft', 'ierr', 'iref', 'lmix', 'gp', 'gpp', 'Fu', 'branch', 'jl',
            'syml', 'Fl', 'vu', 'vl'],
            dtype='object')

    .. minigallery:: radis.io.hitemp.fetch_hitemp

    .. minigallery:: radis.fetch_hitemp

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

    if databank_name == "HITEMP-{molecule}":
        databank_name = databank_name.format(**{"molecule": molecule})

    local_databases = abspath(local_databases.replace("~", expanduser("~")))

    ldb = HITEMPDatabaseManager(
        databank_name,
        molecule=molecule,
        local_databases=local_databases,
        verbose=verbose,
        chunksize=chunksize,
    )

    # Get list of all expected local files for this database:
    local_files, urlnames = ldb.get_filenames()

    # Delete files if needed:
    if cache == "regen":
        ldb.remove_local_files(local_files)

    # Get number of lines to download
    download_files = []
    for local_file in local_files:
        if exists(local_file):
            try:
                check_not_deprecated(
                    local_file,
                    metadata_is={},
                    metadata_keys_contain=["wavenumber_min", "wavenumber_max"],
                )
            except DeprecatedFileWarning as err:
                if cache == "force":
                    raise err
                else:  # delete file to regenerate it in the end of the script
                    if verbose:
                        printr(
                            "File {0} deprecated:\n{1}\nDeleting it!".format(
                                local_file, str(err)
                            )
                        )
                    os.remove(local_file)
                    download_files.append(local_file)
        else:
            download_files.append(local_file)

    # Download files
    # url_to_download =
    if len(download_files) > 0:
        _, Ntotal_lines_expected = ldb.fetch_url_and_Nlines()
        if urlnames is None:
            urlnames = ldb.fetch_urlnames()
        filesmap = dict(zip(local_files, urlnames))
        download_urls = [filesmap[k] for k in download_files]
        ldb.download_and_parse(download_urls, download_files)

        # Done: add final checks
        if molecule not in ["CO2", "H2O"]:
            # ... check on the created file that all lines are there :
            nrows = ldb.get_nrows(local_files[0])
            # assert nrows == Nlines
            if nrows != Ntotal_lines_expected:
                raise AssertionError(
                    f"Number of lines in local database ({nrows:,}) "
                    + "differ from the expected number of lines for "
                    + f"HITEMP {molecule}: {Ntotal_lines_expected}"
                )

    # # Open and run

    # if local_files_exist:

    # check that the database is correctly registered in ~/radis.json
    # if not databank_name in getDatabankList():
    #     # if not, register it.
    #     if len(local_files) == 1:
    #         local_file = local_files[0]
    #         # ... First check number of lines is correct :
    #         if Ntotal_lines_expected is None:
    #             _, Ntotal_lines_expected = get_url_and_Nlines(molecule)
    #         error_msg = ""
    #         with pd.HDFStore(local_file, "r") as store:
    #             nrows = store.get_storer("df").nrows
    #             # TODO: replace with Database.get_rows()  # which would work for any backend (pytables / h5py)
    #             if nrows != Ntotal_lines_expected:
    #                 error_msg += (
    #                     f"\nNumber of lines in local database ({nrows:,}) "
    #                     + "differ from the expected number of lines for "
    #                     + f"HITEMP {molecule}: {Ntotal_lines_expected}"
    #                 )
    #             file_metadata = store.get_storer("df").attrs.metadata
    #             # TODO: replace with Database.get_metadata()  # which would work for any backend (pytables / h5py)
    #             for k in [
    #                 "wavenumber_min",
    #                 "wavenumber_max",
    #                 "download_url",
    #                 "download_date",
    #             ]:
    #                 if k not in file_metadata:
    #                     error_msg += (
    #                         "\nMissing key in file metadata to register the database "
    #                         + f"automatically : {k}"
    #                     )

    #         if error_msg:
    #             raise ValueError(
    #                 f"{databank_name} not declared in your RADIS ~/.config file although "
    #                 + f"{local_file} exists. {error_msg}\n"
    #                 + "If you know this file, add it to ~/radis.json manually. "
    #                 + "Else regenerate the database with:\n\t"
    #                 + ">>> radis.SpectrumFactory().fetch_databank(..., use_cached='regen')"
    #                 + "\nor\n\t"
    #                 + ">>> radis.io.hitemp.fetch_hitemp({molecule}, cache='regen')"
    #                 + "\n\n⚠️ It will re-download & uncompress the whole database "
    #                 + "from HITEMP.\n\nList of declared databanks: {getDatabankList()}.\n"
    #                 + f"{local_file} metadata: {file_metadata}"
    #             )

    #         # ... Database looks ok : register it
    #         if verbose:
    #             print(
    #                 f"{databank_name} not declared in your RADIS ~/.config file although "
    #                 + f"{local_file} exists. Registering the database automatically."
    #             )

    #         register_database(
    #             databank_name,
    #             [local_file],
    #             molecule=molecule,
    #             wmin=file_metadata["wavenumber_min"],
    #             wmax=file_metadata["wavenumber_max"],
    #             download_date=file_metadata["download_date"],
    #             urlname=file_metadata["download_url"],
    #             verbose=verbose,
    #         )
    #     else:    # case of CO2, H2O

    #         _, files_wmin, files_wmax = keep_only_relevant(inputfiles)

    #         # just register the full database
    #         register_database(
    #             databank_name,
    #             local_files,  # write all filenames in the database
    #             molecule=molecule,
    #             wmin=files_wmin,
    #             wmax=files_wmax,
    #             download_date=file_metadata["download_date"],
    #             urlname=file_metadata["download_url"],
    #             verbose=verbose,
    #         )

    # Database exists, and is registered : we can return it directly
    # ... unless it's CO2 / H2O : there are many files and we don't need all of them

    # if molecule in ["CO2", "H2O"]:
    #     # ... Check database is complete
    #     entries = getDatabankEntries(databank_name)
    #     if load_wavenum_min is not None and load_wavenum_max is not None:
    #         if (float(entries["wavenum_max"]) < load_wavenum_min) or (
    #             float(entries["wavenum_min"]) > load_wavenum_max
    #         ):
    #             # downloaded range is enough
    #             full_range_downloaded = False
    #             # TODO / Unless we're asking for a range beyond the maximum range given...
    #             # That's relevant for all databases... What to do ? BeyondTheRangeWarning (checking if also complete)
    #     # if only one of the two extrema is given > NotImplemented > we
    #     # ... require to download the full database.
    #     # elif load_wavenum_min is not None:
    #     #     if float(fname_wmax) >= wavenum_min:
    #     #         relevant.append(file)
    #     # elif load_wavenum_max is not None:
    #     #     if float(fname_wmin) <= wavenum_max:
    #     #         relevant.append(file)
    #     else:
    #         # check number of lines is complete
    #         with pd.HDFStore(local_file, "r") as store:
    #             nrows = store.get_storer("df").nrows
    #             # TODO: replace with Database.get_rows()  # which would work for any backend (pytables / h5py)
    #             if nrows < Ntotal_lines_expected:
    #                 raise DeprecatedFileWarning(
    #                     f"\nNumber of lines in local database ({nrows:,}) "
    #                     + "differ from the expected number of lines for "
    #                     + f"HITEMP {molecule}: {Ntotal_lines_expected}"
    #                 )
    # if full_range_downloaded:
    #     if verbose:
    #         print(f"Using existing database {databank_name}")
    #     df = hdf2df(
    #         local_file,
    #         isotope=isotope,
    #         load_wavenum_min=load_wavenum_min,
    #         load_wavenum_max=load_wavenum_max,
    #         verbose=verbose,
    #     )
    #     return (df, local_file) if return_local_path else df

    # else:
    #     # keep a continuous database (i.e. : all lines from [wmin] to [wmax] even if we only need two separate sections)

    #     # Note : with Vaex it's very easy to open many different files with
    #     # only a minimal overhead. Why put everything in a large file then ?
    #     raise NotImplementedError
    #     # ... TODO

    #     # use is_relevant function of Corentin ?

    #     # therefore, only download relevant_files...

    # TODO here: for CO2; H2O : if exists : re-check if valid ; compare with relevant range ;
    # download only missing!

    ###################################################
    # Doesnt exist : download    missing files

    if not ldb.is_registered():
        ldb.register(local_files, urlnames)

    if len(download_files) > 0 and clean_cache_files:
        ldb.clean_download_files(urlnames)

    # Load and return
    df = ldb.load(
        local_files,
        isotope=isotope,
        load_wavenum_min=load_wavenum_min,
        load_wavenum_max=load_wavenum_max,
    )

    return (df, local_files) if return_local_path else df


#%%

if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_hitemp.py"]))
