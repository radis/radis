# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 22:40:51 2021

@author: erwan


https://stackoverflow.com/questions/55610891/numpy-load-from-io-bytesio-stream
https://stupidpythonideas.blogspot.com/2014/07/three-ways-to-read-files.html

"""

import re
import urllib.request
from os.path import abspath, basename, exists, expanduser, join
from typing import Union

import numpy as np

try:
    from .dbmanager import DatabaseManager
    from .hdf5 import update_pytables_to_vaex
    from .hitran import columns_2004, parse_global_quanta, parse_local_quanta
    from .tools import (
        _create_dtype,
        _get_linereturnformat,
        _ndarray2df,
        replace_PQR_with_m101,
    )
except ImportError:  # ran from here
    from radis.io.dbmanager import DatabaseManager
    from radis.io.hitran import columns_2004, parse_global_quanta, parse_local_quanta
    from radis.io.hdf5 import update_pytables_to_vaex
    from radis.io.tools import (
        _create_dtype,
        _get_linereturnformat,
        _ndarray2df,
        replace_PQR_with_m101,
    )

from radis.db import MOLECULES_LIST_NONEQUILIBRIUM
from radis.misc.progress_bar import ProgressBar

HITEMP_MOLECULES = ["H2O", "CO2", "N2O", "CO", "CH4", "NO", "NO2", "OH"]


def keep_only_relevant(
    inputfiles,
    wavenum_min=None,
    wavenum_max=None,
) -> Union[list, float, float]:
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
    for filepath in inputfiles:
        file = basename(filepath)
        fname_wmin, fname_wmax = re.findall(wavenum_format, file)
        relevant = False
        if wavenum_min is not None and wavenum_max is not None:
            if (float(fname_wmax) > wavenum_min) and (float(fname_wmin) < wavenum_max):
                # strict '>' :  we exclude "CO2-02_02250-02500_HITEMP2010.h5'" if calculating 2500 - 3000 cm-1
                # strict '<' :  we exclude "CO2-02_03000-03250_HITEMP2010.h5" if calculating 2500 - 3000 cm-1
                relevant = True
        elif wavenum_min is not None:
            if float(fname_wmax) > wavenum_min:
                # strict '>' :  we exclude "CO2-02_02250-02500_HITEMP2010.h5'" if calculating 2500 - 3000 cm-1
                relevant = True
        elif wavenum_max is not None:
            if float(fname_wmin) < wavenum_max:
                # strict '<' :  we exclude "CO2-02_03000-03250_HITEMP2010.h5" if calculating 2500 - 3000 cm-1
                relevant = True
        else:
            relevant = True
        if relevant:
            relevantfiles.append(filepath)
            files_wmin = min(float(fname_wmin), files_wmin)
            files_wmax = max(float(fname_wmax), files_wmax)

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
        local_databases,
        engine="default",
        verbose=True,
        chunksize=100000,
        parallel=True,
    ):
        super().__init__(
            name,
            molecule,
            local_databases,
            engine,
            verbose=verbose,
            parallel=parallel,
        )
        self.chunksize = chunksize
        self.downloadable = True
        self.base_url = None
        self.Nlines = None
        self.wmin = None  # available on HITEMP website. See HITEMPDatabaseManager.fetch_url_Nlines_wmin_wmax
        self.wmax = None  # available on HITEMP website. See HITEMPDatabaseManager.fetch_url_Nlines_wmin_wmax
        self.urlnames = None

    def fetch_url_Nlines_wmin_wmax(self, hitemp_url="https://hitran.org/hitemp/"):
        """requires connexion"""

        molecule = self.molecule

        if (
            self.base_url is not None
            and self.Nlines is not None
            and self.wmin is not None
            and self.wmax is not None
        ):
            return self.base_url, self.Nlines, self.wmin, self.wmax

        else:

            response = urllib.request.urlopen(hitemp_url)

            # Alternative to return a Pandas Dataframe :
            # ... Doesnt work because missing <tr> in HITEMP website table for N2O

            # soup = BeautifulSoup(
            #         response, features="lxml"
            #     )
            # table = soup.find(lambda tag: tag.name=='table' and tag.has_attr('id') and tag['id']=="hitemp-molecules-table")

            # def tableDataText(table):
            #     """Parses a html segment started with tag <table> followed
            #     by multiple <tr> (table rows) and inner <td> (table data) tags.
            #     It returns a list of rows with inner columns.
            #     Accepts only one <th> (table header/data) in the first row.

            #     From https://stackoverflow.com/a/58274853/5622825
            #     """

            #

            #     def rowgetDataText(tr, coltag='td'): # td (data) or th (header)
            #         return [td.get_text(strip=True) for td in tr.find_all(coltag)]
            #     rows = []
            #     trs = table.find_all('tr')
            #     headerow = rowgetDataText(trs[0], 'th')
            #     if headerow: # if there is a header row include first
            #         rows.append(headerow)
            #         trs = trs[1:]
            #     for tr in trs: # for every table row
            #         rows.append(rowgetDataText(tr, 'td') ) # data row

            #     df = pd.DataFrame(rows[1:], columns=rows[0])
            #     df.index = df.Formula

            #     return df

            # df = tableDataText(table)

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
            wmin = int(re.findall(r"(\d+)", lines[4].replace("&nbsp;", ""))[0])
            wmax = int(re.findall(r"(\d+)", lines[5].replace("&nbsp;", ""))[0])
            url = "https://hitran.org" + re.findall(r'href="(.+?)"', lines[7])[0]

            self.base_url, self.Nlines, self.wmin, self.wmax = url, Nlines, wmin, wmax

        return url, Nlines, wmin, wmax

    def fetch_urlnames(self):
        """requires connexion"""

        if self.urlnames is not None:
            return self.urlnames

        molecule = self.molecule

        if molecule in ["H2O", "CO2"]:

            base_url, Ntotal_lines_expected, _, _ = self.fetch_url_Nlines_wmin_wmax()
            response = urllib.request.urlopen(base_url)
            response_string = response.read().decode()
            inputfiles = re.findall('href="(\S+.zip)"', response_string)

            urlnames = [join(base_url, f) for f in inputfiles]

        elif molecule in HITEMP_MOLECULES:
            url, Ntotal_lines_expected, _, _ = self.fetch_url_Nlines_wmin_wmax()
            urlnames = [url]
        else:
            raise KeyError(
                f"Please choose one of HITEMP molecules : {HITEMP_MOLECULES}. Got '{molecule}'"
            )

        self.urlnames = urlnames

        return urlnames

    def keep_only_relevant(
        self,
        inputfiles,
        wavenum_min=None,
        wavenum_max=None,
    ) -> list:
        """For CO2 and H2O, return only relevant files for given wavenumber range.

        If other molecule, return the file anyway.
        see :py:func:`radis.io.hitemp.keep_only_relevant`"""
        if self.molecule in ["CO2", "H2O"]:
            inputfiles, _, _ = keep_only_relevant(inputfiles, wavenum_min, wavenum_max)
        return inputfiles

    def get_linereturn_format(self, opener, urlname, columns):

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
        return linereturnformat

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
    ):
        """Uncompress ``urlname`` into ``local_file``.
        Also add metadata

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

        if not verbose:
            pbar_active = False

        linereturnformat = self.get_linereturn_format(opener, urlname, columns)

        Nlines = 0
        Nlines_raw = 0
        Nlines_tot = Nlines + pbar_Nlines_already
        _, Ntotal_lines_expected, _, _ = self.fetch_url_Nlines_wmin_wmax()
        if pbar_Ntot_estimate_factor:
            # multiply Ntotal_lines_expected by pbar_Ntot_estimate_factor
            # (accounts for total lines divided in number of files, and
            # not all files downloaded)
            Ntotal_lines_expected = int(
                Ntotal_lines_expected * pbar_Ntot_estimate_factor
            )
        pb = ProgressBar(N=Ntotal_lines_expected, active=pbar_active, t0=pbar_t0)
        wmin = np.inf
        wmax = 0

        writer = self.get_hdf5_manager()

        with opener.open(urlname) as gfile:  # locally downloaded file

            dt = _create_dtype(columns, linereturnformat)

            if verbose:
                print(f"Download complete. Parsing {molecule} database to {local_file}")

            # assert not(exists(local_file))

            b = np.zeros(chunksize, dtype=dt)  # receives the HITRAN 160-character data.

            for nbytes in iter(lambda: gfile.readinto(b), 0):

                if not b[-1]:
                    # End of file flag within the chunk (but does not start
                    # with End of file flag) so nbytes != 0
                    b = get_last(b)

                df = _ndarray2df(b, columns, linereturnformat)

                # Post-processing :
                # ... Add local quanta attributes, based on the HITRAN group
                df = parse_local_quanta(df, molecule, verbose=verbose)

                # ... Add global quanta attributes, based on the HITRAN class
                df = parse_global_quanta(df, molecule, verbose=verbose)

                # Switch 'P', 'Q', 'R' to -1, 0, 1
                if "branch" in df:
                    replace_PQR_with_m101(df)

                writer.write(local_file, df, append=True)

                wmin = np.min((wmin, df.wav.min()))
                wmax = np.max((wmax, df.wav.max()))

                Nlines += len(df)
                Nlines_tot += len(df)
                Nlines_raw += len(b)
                if pbar_Ntot_estimate_factor is None:
                    pbar_Ntot_message = f"{Ntotal_lines_expected:,} lines"
                else:
                    pbar_Ntot_message = f"~{Ntotal_lines_expected:,} lines (estimate)"
                pb.update(
                    Nlines_tot,
                    message=f"  Parsed {Nlines_tot:,} / {pbar_Ntot_message}. Wavenumber range {wmin:.2f}-{wmax:.2f} cm-1 is complete.",
                )
                # Reinitialize for next read
                b = np.zeros(
                    chunksize, dtype=dt
                )  # receives the HITRAN 160-character data.
        writer.combine_temp_batch_files(local_file)  # used for vaex mode only
        if pbar_last:
            pb.update(
                Nlines_tot,
                message=f"  Parsed {Nlines_tot:,} / {Nlines_tot:,} lines. Wavenumber range {wmin:.2f}-{wmax:.2f} cm-1 is complete.",
            )
            pb.done()
        else:
            print("")

        # Check number of lines is consistent
        assert Nlines == Nlines_raw

        # Add metadata
        from radis import __version__

        writer.add_metadata(
            local_file,
            {
                "wavenumber_min": wmin,
                "wavenumber_max": wmax,
                "download_date": self.get_today(),
                "download_url": urlname,
                "total_lines": Nlines_raw,
                "version": __version__,
            },
        )

        return Nlines

    def register(self):
        """register in ~/radis.json"""

        local_files, urlnames = self.get_filenames()
        info = f"HITEMP {self.molecule} lines ({self.wmin:.1f}-{self.wmax:.1f} cm-1) with TIPS-2017 (through HAPI) for partition functions"

        if self.molecule in ["CO2", "H2O"]:
            info = "(registered files will be downloaded only when required) " + info

        dict_entries = {
            "info": info,
            "path": local_files,
            "format": "hitemp-radisdb",
            "parfuncfmt": "hapi",
            "wavenumber_min": self.wmin,
            "wavenumber_max": self.wmax,
            "download_date": self.get_today(),
            "download_url": urlnames,
        }

        # Add energy level calculation
        if self.molecule in MOLECULES_LIST_NONEQUILIBRIUM:
            dict_entries[
                "info"
            ] += " and RADIS spectroscopic constants for rovibrational energies (nonequilibrium)"
            dict_entries["levelsfmt"] = "radis"

        super().register(dict_entries)


def fetch_hitemp(
    molecule,
    local_databases=None,
    databank_name="HITEMP-{molecule}",
    isotope=None,
    load_wavenum_min=None,
    load_wavenum_max=None,
    columns=None,
    cache=True,
    verbose=True,
    chunksize=100000,
    clean_cache_files=True,
    return_local_path=False,
    engine="default",
    parallel=True,
):
    """Stream HITEMP file from HITRAN website. Unzip and build a HDF5 file directly.

    Returns a Pandas DataFrame containing all lines.

    Parameters
    ----------
    molecule: `"H2O", "CO2", "N2O", "CO", "CH4", "NO", "NO2", "OH"`
        HITEMP molecule. See https://hitran.org/hitemp/
    local_databases: str
        where to create the RADIS HDF5 files. Default ``"~/.radisdb/hitemp"``.
        Can be changed in ``radis.config["DEFAULT_DOWNLOAD_PATH"]`` or in ~/radis.json config file
    databank_name: str
        name of the databank in RADIS :ref:`Configuration file <label_lbl_config_file>`
        Default ``"HITEMP-{molecule}"``
    isotope: str, int or None
        load only certain isotopes : ``'2'``, ``'1,2'``, etc. If ``None``, loads
        everything. Default ``None``.
    load_wavenum_min, load_wavenum_max: float (cm-1)
        load only specific wavenumbers.
    columns: list of str
        list of columns to load. If ``None``, returns all columns in the file.

    Other Parameters
    ----------------
    cache: ``True``, ``False``, ``'regen'`` or ``'force'``
        if ``True``, use existing HDF5 file. If ``False`` or ``'regen'``, rebuild it.
        If ``'force'``, raise an error if cache file cannot be used (useful for debugging).
        Default ``True``.
    verbose: bool
    chunksize: int
        number of lines to process at a same time. Higher is usually faster
        but can create Memory problems and keep the user uninformed of the progress.
    clean_cache_files: bool
        if ``True`` clean downloaded cache files after HDF5 are created.
    return_local_path: bool
        if ``True``, also returns the path of the local database file.
    engine: 'pytables', 'vaex', 'default'
        which HDF5 library to use. If 'default' use the value from ~/radis.json
    parallel: bool
        if ``True``, uses joblib.parallel to load database with multiple processes

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

    Notes
    -----
    if using ``load_only_wavenum_above/below`` or ``isotope``, the whole
    database is anyway downloaded and uncompressed to ``local_databases``
    fast access .HDF5 files (which will take a long time on first call). Only
    the expected wavenumber range & isotopes are returned. The .HFD5 parsing uses
    :py:func:`~radis.io.hdf5.hdf2df`

    See Also
    --------
    :py:func:`~radis.io.hdf5.hdf2df`, :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank`

    """

    if r"{molecule}" in databank_name:
        databank_name = databank_name.format(**{"molecule": molecule})

    if local_databases is None:
        import radis

        local_databases = join(radis.config["DEFAULT_DOWNLOAD_PATH"], "hitemp")
    local_databases = abspath(expanduser(local_databases))

    ldb = HITEMPDatabaseManager(
        databank_name,
        molecule=molecule,
        local_databases=local_databases,
        verbose=verbose,
        chunksize=chunksize,
        parallel=parallel,
        engine=engine,
    )

    # Get list of all expected local files for this database:
    local_files, urlnames = ldb.get_filenames()

    # Delete files if needed:
    relevant_files = ldb.keep_only_relevant(
        local_files, load_wavenum_min, load_wavenum_max
    )
    if cache == "regen":
        ldb.remove_local_files(relevant_files)
    ldb.check_deprecated_files(
        ldb.get_existing_files(relevant_files),
        auto_remove=True if cache != "force" else False,
    )

    # Get missing files
    download_files = ldb.get_missing_files(local_files)
    download_files = ldb.keep_only_relevant(
        download_files, load_wavenum_min, load_wavenum_max
    )
    # do not re-download files if they exist in another format :
    if engine in ["vaex", "auto", "default"]:
        # ... convert files if asked:
        from radis import config

        if config["AUTO_UPDATE_DATABASE"]:
            converted = []
            for f in download_files:
                if exists(f.replace(".hdf5", ".h5")):
                    update_pytables_to_vaex(f.replace(".hdf5", ".h5"))
                    converted.append(f)
            download_files = [f for f in download_files if f not in converted]
        # do not re-download remaining files that exist. Let user decide what to do.
        # (download & re-parsing is a long solution!)
        download_files = [
            f for f in download_files if not exists(f.replace(".hdf5", ".h5"))
        ]

    # Download files
    if len(download_files) > 0:
        if urlnames is None:
            urlnames = ldb.fetch_urlnames()
        filesmap = dict(zip(local_files, urlnames))
        download_urls = [filesmap[k] for k in download_files]
        ldb.download_and_parse(download_urls, download_files)

        # Done: add final checks
        if molecule not in ["CO2", "H2O"]:
            # check on the created file that all lines are there
            # ... (for CO2, H2O, database is split in many files so it's not possible)
            nrows = ldb.get_nrows(local_files[0])
            # assert nrows == Nlines
            _, Ntotal_lines_expected, _, _ = ldb.fetch_url_Nlines_wmin_wmax()
            if nrows != Ntotal_lines_expected:
                raise AssertionError(
                    f"Number of lines in local database ({nrows:,}) "
                    + "differ from the expected number of lines for "
                    + f"HITEMP {molecule}: {Ntotal_lines_expected}"
                )

    # Register
    if not ldb.is_registered():
        ldb.register()

    if len(download_files) > 0 and clean_cache_files:
        ldb.clean_download_files()

    # Load and return
    files_loaded = ldb.keep_only_relevant(
        local_files, load_wavenum_min, load_wavenum_max
    )

    if isotope and type(isotope) == int:
        isotope = str(isotope)

    df = ldb.load(
        files_loaded,  # filter other files,
        columns=columns,
        isotope=isotope,
        load_wavenum_min=load_wavenum_min,  # for relevant files, get only the right range
        load_wavenum_max=load_wavenum_max,
    )

    return (df, files_loaded) if return_local_path else df


#%%

if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_hitemp.py"]))
