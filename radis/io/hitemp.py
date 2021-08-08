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
from datetime import date
from io import BytesIO
from os.path import abspath, exists, expanduser, join, splitext
from zipfile import ZipFile

import numpy as np
import pandas as pd
from dateutil.parser import parse as parse_date
from numpy import DataSource

import radis
from radis.db import MOLECULES_LIST_NONEQUILIBRIUM
from radis.io.cache_files import check_not_deprecated
from radis.io.hdf5 import hdf2df
from radis.io.hitran import columns_2004, parse_global_quanta, parse_local_quanta
from radis.io.tools import (
    _create_dtype,
    _get_linereturnformat,
    _ndarray2df,
    replace_PQR_with_m101,
)
from radis.misc.config import addDatabankEntries, getDatabankEntries, getDatabankList
from radis.misc.progress_bar import ProgressBar
from radis.misc.warning import DatabaseAlreadyExists, DeprecatedFileWarning

LAST_VALID_DATE = (
    "01 Jan 2010"  # set to a later date to force re-download of all HITEMP databases
)
HITEMP_MOLECULES = ["H2O", "CO2", "N2O", "CO", "CH4", "NO", "NO2", "OH"]
DATA_COLUMNS = ["iso", "wav"]
"""
list : only these column names will be searchable directly on disk to
only load certain lines. See :py:func:`~radis.io.hdf5.hdf2df`
"""
# TODO: WIP. Maybe move somewhere else to be also used by HITRAN queries


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


def get_url_and_Nlines(molecule, hitemp_url="https://hitran.org/hitemp/"):
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

    return url, Nlines


def keep_only_relevant(
    inputfiles, wavenum_min=None, wavenum_max=None, wavenum_format=r"\d{5}"
):
    """Parser file names for ``wavenum_format`` (min and max) and only keep
    relevant files if the requested range is ``[wavenum_min, wavenum_max]``"""
    relevant = []
    for file in inputfiles:
        fname_wmin, fname_wmax = re.findall(wavenum_format, file)
        if wavenum_min is not None and wavenum_max is not None:
            if (float(fname_wmax) >= wavenum_min) and (
                float(fname_wmin) <= wavenum_max
            ):
                relevant.append(file)
        elif wavenum_min is not None:
            if float(fname_wmax) >= wavenum_min:
                relevant.append(file)
        elif wavenum_max is not None:
            if float(fname_wmin) <= wavenum_max:
                relevant.append(file)
        else:
            relevant.append(file)
    return relevant


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
    # TODO ? : unzip only parts of the database
    # see https://github.com/radis/radis/pull/194

    if databank_name == "HITEMP-{molecule}":
        databank_name = databank_name.format(**{"molecule": molecule})

    Ntotal_lines_expected = None

    # First; checking if the database is registered in radis.json
    def fetch_urlnames(molecule):

        if molecule in ["H2O", "CO2"]:

            base_url, Ntotal_lines_expected = get_url_and_Nlines(molecule)
            response = urllib.request.urlopen(base_url)
            response_string = response.read().decode()
            inputfiles = re.findall('href="(\S+.zip)"', response_string)

            inputfiles = keep_only_relevant(
                inputfiles, load_wavenum_min, load_wavenum_max, wavenum_format=r"\d{5}"
            )
            if verbose:
                print("relevant files:", inputfiles)

            urlnames = [base_url + f for f in inputfiles]

            # harcoded local name; regroups all 0X_[wmin]_[wmax]_HITEMP2010* files :
            assert "HITEMP2010" in inputfiles[0]
            local_fname = inputfiles[0][:3] + "HITEMP2010.h5"

        elif molecule in HITEMP_MOLECULES:
            url, Ntotal_lines_expected = get_url_and_Nlines(molecule)
            urlnames = [url]
            local_fname = (
                splitext(splitext(url.split("/")[-1])[0])[0]  # twice to remove .par.bz2
                + ".h5"
            )
        else:
            raise KeyError(
                f"Please choose one of HITEMP molecules : {HITEMP_MOLECULES}. Got '{molecule}'"
            )

        try:
            os.mkdir(local_databases)
        except OSError:
            pass
        else:
            if verbose:
                print("Created folder :", local_databases)

        local_file = abspath(
            join(
                local_databases,
                molecule + "-" + local_fname,
            )
        )
        return local_file, urlnames

    if databank_name in getDatabankList():
        entries = getDatabankEntries(databank_name)
        try:
            assert "download_url" in entries
            assert "download_date" in entries
            assert parse_date(entries["download_date"]) > parse_date(LAST_VALID_DATE)
        except AssertionError as err:
            raise DeprecatedFileWarning("Database file {0} not valid anymore") from err
        else:
            local_file = entries["path"]
            assert len(local_file) == 1
            local_file = local_file[0]
    else:
        local_file, urlnames = fetch_urlnames(molecule)

    # Now, check if the local file (as registered in radis.json, or fetched from the website)
    # exists
    local_databases = abspath(local_databases.replace("~", expanduser("~")))

    if not cache or cache == "regen":
        # Delete existing HDF5 file
        if exists(local_file):
            if verbose:
                print("Removing existing file ", local_file)
                # TODO: also clean the getDatabankList? Todo once it is in JSON format. https://github.com/radis/radis/issues/167
            os.remove(local_file)

    if exists(local_file):
        # Read and return from local file

        # check metadata :
        check_not_deprecated(
            local_file,
            metadata_is={},
            metadata_keys_contain=["wavenumber_min", "wavenumber_max"],
        )
        # check database is registered in ~/radis.json
        if not databank_name in getDatabankList():
            # if not, register it.
            # ... First check number of lines is correct :
            if Ntotal_lines_expected is None:
                _, Ntotal_lines_expected = get_url_and_Nlines(molecule)
            error_msg = ""
            with pd.HDFStore(local_file, "r") as store:
                nrows = store.get_storer("df").nrows
                # TODO: replace with Database.get_rows()  # which would work for any backend (pytables / h5py)
                if nrows != Ntotal_lines_expected:
                    error_msg += (
                        f"\nNumber of lines in local database ({nrows:,}) "
                        + "differ from the expected number of lines for "
                        + f"HITEMP {molecule}: {Ntotal_lines_expected}"
                    )
                file_metadata = store.get_storer("df").attrs.metadata
                # TODO: replace with Database.get_metadata()  # which would work for any backend (pytables / h5py)
                for k in [
                    "wavenumber_min",
                    "wavenumber_max",
                    "download_url",
                    "download_date",
                ]:
                    if k not in file_metadata:
                        error_msg += (
                            "\nMissing key in file metadata to register the database "
                            + f"automatically : {k}"
                        )

            if error_msg:
                raise ValueError(
                    f"{databank_name} not declared in your RADIS ~/.config file although "
                    + f"{local_file} exists. {error_msg}\n"
                    + "If you know this file, add it to ~/radis.json manually. "
                    + "Else regenerate the database with:\n\t"
                    + ">>> radis.SpectrumFactory().fetch_databank(..., use_cached='regen')"
                    + "\nor\n\t"
                    + ">>> radis.io.hitemp.fetch_hitemp({molecule}, cache='regen')"
                    + "\n\n⚠️ It will re-download & uncompress the whole database "
                    + "from HITEMP.\n\nList of declared databanks: {getDatabankList()}.\n"
                    + f"{local_file} metadata: {file_metadata}"
                )

            # ... Database looks ok : register it
            if verbose:
                print(
                    f"{databank_name} not declared in your RADIS ~/.config file although "
                    + f"{local_file} exists. Registering the database automatically."
                )

            register_database(
                databank_name,
                [local_file],
                molecule=molecule,
                wmin=file_metadata["wavenumber_min"],
                wmax=file_metadata["wavenumber_max"],
                download_date=file_metadata["download_date"],
                urlname=file_metadata["download_url"],
                verbose=verbose,
            )

        # Database exists, and is registered : we can return it directly
        # ... unless it's CO2 / H2O : there are many files and maybe it's not complete.

        full_range_downloaded = True

        if molecule in ["CO2", "H2O"]:
            # ... Check database is complete
            entries = getDatabankEntries(databank_name)
            if load_wavenum_min is not None and load_wavenum_max is not None:
                if (float(entries["wavenum_max"]) < load_wavenum_min) or (
                    float(entries["wavenum_min"]) > load_wavenum_max
                ):
                    # downloaded range is enough
                    full_range_downloaded = False
                    # TODO / Unless we're asking for a range beyond the maximum range given...
                    # That's relevant for all databases... What to do ? BeyondTheRangeWarning (checking if also complete)
            # if only one of the two extrema is given > NotImplemented > we
            # ... require to download the full database.
            # elif load_wavenum_min is not None:
            #     if float(fname_wmax) >= wavenum_min:
            #         relevant.append(file)
            # elif load_wavenum_max is not None:
            #     if float(fname_wmin) <= wavenum_max:
            #         relevant.append(file)
            else:
                # check number of lines is complete
                with pd.HDFStore(local_file, "r") as store:
                    nrows = store.get_storer("df").nrows
                    # TODO: replace with Database.get_rows()  # which would work for any backend (pytables / h5py)
                    if nrows < Ntotal_lines_expected:
                        raise DeprecatedFileWarning(
                            f"\nNumber of lines in local database ({nrows:,}) "
                            + "differ from the expected number of lines for "
                            + f"HITEMP {molecule}: {Ntotal_lines_expected}"
                        )
        if full_range_downloaded:
            if verbose:
                print(f"Using existing database {databank_name}")
            df = hdf2df(
                local_file,
                isotope=isotope,
                load_wavenum_min=load_wavenum_min,
                load_wavenum_max=load_wavenum_max,
                verbose=verbose,
            )
            return (df, local_file) if return_local_path else df

        else:
            # keep a continuous database (i.e. : all lines from [wmin] to [wmax] even if we only need two separate sections)

            # Note : with Vaex it's very easy to open many different files with
            # only a minimal overhead. Why put everything in a large file then ?
            raise NotImplementedError
            # ... TODO

            # use is_relevant function of Corentin ?

            # therefore, only download relevant_files...

    # TODO here: for CO2; H2O : if exists : re-check if valid ; compare with relevant range ;
    # download only missing!

    ###################################################
    # Doesnt exist : download
    if Ntotal_lines_expected is None:
        _, Ntotal_lines_expected = get_url_and_Nlines(molecule)
    ds = DataSource(join(local_databases, "downloads"))
    Ndownload = 1
    Ntotal_downloads = len(urlnames)

    Nlines = 0
    pb = ProgressBar(N=Ntotal_lines_expected, active=verbose)

    wmin = np.inf
    wmax = 0

    for urlname in urlnames:

        if verbose:
            inputf = urlname.split("/")[-1]
            print(
                f"Downloading {inputf} for {molecule} ({Ndownload}/{Ntotal_downloads})."
            )
        download_date = date.today().strftime("%d %b %Y")

        columns = columns_2004

        # Get linereturn (depends on OS, but file may also have been generated
        # on a different OS. Here we simply read the file to find out)
        with ds.open(urlname) as gfile:  # locally downloaded file
            dt = _create_dtype(
                columns, "a2"
            )  # 'a2' allocates space to get \n or \n\r for linereturn character
            b = np.zeros(1, dtype=dt)
            try:
                gfile.readinto(b)
            except EOFError as err:
                raise ValueError(
                    f"End of file while parsing file {ds.abspath(urlname)}. May be due to download error. Delete file ?"
                ) from err
            linereturnformat = _get_linereturnformat(b, columns)

        with ds.open(urlname) as gfile:  # locally downloaded file

            dt = _create_dtype(columns, linereturnformat)
            b = np.zeros(chunksize, dtype=dt)  # receives the HITRAN 160-character data.

            if verbose:
                print(
                    f"Download complete. Building {molecule} database to {local_file}"
                )

            with pd.HDFStore(local_file, mode="a", complib="blosc", complevel=9) as f:

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
                    pb.update(
                        Nlines,
                        message=f"Parsed {Nlines:,} / {Ntotal_lines_expected:,} lines. Wavenumber range {wmin:.2f}-{wmax:.2f} cm-1 is complete.",
                    )

                    # Reinitialize for next read
                    b = np.zeros(
                        chunksize, dtype=dt
                    )  # receives the HITRAN 160-character data.

        Ndownload += 1

    pb.done()

    url_store = urlnames[0] if len(urlnames) == 1 else urlnames
    with pd.HDFStore(local_file, mode="a", complib="blosc", complevel=9) as f:

        f.get_storer("df").attrs.metadata = {
            "wavenumber_min": wmin,
            "wavenumber_max": wmax,
            "download_date": download_date,
            "download_url": url_store,
            "version": radis.__version__,
        }

    # Done: add final checks
    # ... check on the created file that all lines are there :
    with pd.HDFStore(local_file, "r") as store:
        nrows = store.get_storer("df").nrows
        assert nrows == Nlines
        if nrows != Ntotal_lines_expected:
            raise AssertionError(
                f"Number of lines in local database ({nrows:,}) "
                + "differ from the expected number of lines for "
                + f"HITEMP {molecule}: {Ntotal_lines_expected}"
            )

    # Add database to  ~/radis.json
    register_database(
        databank_name,
        [local_file],
        molecule,
        wmin,
        wmax,
        download_date,
        url_store,
        verbose,
    )

    df = hdf2df(
        local_file,
        isotope=isotope,
        load_wavenum_min=load_wavenum_min,
        load_wavenum_max=load_wavenum_max,
        verbose=verbose,
    )

    # Fully unzipped (and working, as it was reloaded): clean
    if clean_cache_files:
        os.remove(ds._findfile(urlname))
        if verbose >= 3:
            from radis.misc.printer import printg

            printg("... removed downloaded cache file")

    return (df, local_file) if return_local_path else df


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
            # TODO @dev : replace once we have configfile as JSON (https://github.com/radis/radis/issues/167)
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
