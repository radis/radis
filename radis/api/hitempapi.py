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
from os.path import basename, commonpath, join
from typing import Union

import numpy as np
import requests
from bs4 import BeautifulSoup
from cryptography.fernet import Fernet
from dotenv import load_dotenv
from getpass4 import getpass
from tqdm import tqdm

try:
    from .dbmanager import DatabaseManager
    from .hitranapi import columns_2004, parse_global_quanta, parse_local_quanta
    from .tools import (
        _create_dtype,
        _get_linereturnformat,
        _ndarray2df,
        replace_PQR_with_m101,
    )
except ImportError:  # ran from here
    if __name__ == "__main__":  # running from this file, as a script
        from radis.api.dbmanager import DatabaseManager
        from radis.api.hitranapi import (
            columns_2004,
            parse_global_quanta,
            parse_local_quanta,
        )
        from radis.io.tools import (
            _create_dtype,
            _get_linereturnformat,
            _ndarray2df,
            replace_PQR_with_m101,
        )
    else:
        raise

from radis.db import MOLECULES_LIST_NONEQUILIBRIUM
from radis.misc.progress_bar import ProgressBar

HITEMP_MOLECULES = ["H2O", "CO2", "N2O", "CO", "CH4", "NO", "NO2", "OH"]


def keep_only_relevant(
    inputfiles, wavenum_min=None, wavenum_max=None, verbose=True
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

    if verbose and relevantfiles != []:
        if len(relevantfiles) > 1:
            folder = commonpath(relevantfiles)
            # file_list = [x.replace(folder+'\\', '') for x in relevantfiles]
            # print(f"In {folder} keep only relevant input files: {file_list}")
            print(f"In ``{folder}`` keep only relevant input files:")
            for file in relevantfiles:
                print(file.replace(folder + "\\", ""))
        else:
            print(f"Keep only relevant input file: {relevantfiles}")

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


def setup_credentials():
    """Set up HITRAN credentials and store them in .env file"""
    # Get credentials from user
    username = input("Enter HITRAN username: ")
    password = getpass("Enter HITRAN password: ")
    return username, password


def get_encryption_key():
    """Get or create encryption key for HITRAN credentials"""
    from radis.misc.utils import getProjectRoot

    key_path = os.path.join(getProjectRoot(), "db", ".radis.key")
    if os.path.exists(key_path):
        with open(key_path, "rb") as f:
            return f.read()
    else:
        # Generate a new key
        key = Fernet.generate_key()
        # Ensure directory exists
        os.makedirs(os.path.dirname(key_path), exist_ok=True)
        # Save the key
        with open(key_path, "wb") as f:
            f.write(key)
        # Set restrictive permissions on key file
        os.chmod(key_path, 0o600)
        return key


def encrypt_password(password):
    """Encrypt password using Fernet symmetric encryption"""
    key = get_encryption_key()
    f = Fernet(key)
    return f.encrypt(password.encode()).decode()


def decrypt_password(encrypted_password):
    """Decrypt password using Fernet symmetric encryption"""
    key = get_encryption_key()
    f = Fernet(key)
    return f.decrypt(encrypted_password.encode()).decode()


def store_credentials(username, password):
    """Store HITRAN credentials in radis.env file in db folder with encrypted username and password"""
    from radis.misc.utils import getProjectRoot

    env_path = os.path.join(getProjectRoot(), "db", "radis.env")

    # Create db directory if it doesn't exist
    os.makedirs(os.path.dirname(env_path), exist_ok=True)

    # Encrypt both username and password before storing
    encrypted_username = encrypt_password(username)  # reuse same encryption function
    encrypted_password = encrypt_password(password)

    print(
        f"Your HITRAN credentials will be saved securely in {env_path}. You can delete this file if you wish but you will have to prompt your credentials at next download."
    )
    with open(env_path, "w") as f:
        f.write(f"HITRAN_USERNAME={encrypted_username}\n")
        f.write(f"HITRAN_PASSWORD={encrypted_password}\n")

    # Set restrictive permissions on env file
    os.chmod(env_path, 0o600)


def login_to_hitran(verbose=False):
    """Login to HITRAN using stored credentials or prompt if not available"""
    login_url = "https://hitran.org/login/"
    session = requests.Session()

    def attempt_login(username, password):
        """Attempt to login with provided credentials"""
        # Get CSRF token
        response = session.get(login_url)
        soup = BeautifulSoup(response.text, "html.parser")
        csrf = soup.find("input", {"name": "csrfmiddlewaretoken"})["value"]

        login_data = {
            "csrfmiddlewaretoken": csrf,
            "email": username,
            "password": password,
        }

        headers = {
            "Referer": login_url,
            "Origin": "https://hitran.org",
            "Cookie": f"csrftoken={csrf}",
        }

        login_response = session.post(
            login_url, data=login_data, headers=headers, allow_redirects=False
        )

        return login_response, session

    def is_login_successful(response):
        """Check if login was successful by looking for specific elements"""
        return response.status_code == 302 or "Logout" in response.text

    # Check if radis.env exists in db folder
    from radis.misc.utils import getProjectRoot

    env_path = os.path.join(getProjectRoot(), "db", "radis.env")

    if os.path.exists(env_path):
        load_dotenv(env_path)
        encrypted_username = os.getenv("HITRAN_USERNAME")
        encrypted_password = os.getenv("HITRAN_PASSWORD")

        if encrypted_username and encrypted_password:
            try:
                # Decrypt both username and password
                username = decrypt_password(encrypted_username)
                password = decrypt_password(encrypted_password)
                login_response, session = attempt_login(username, password)
                if is_login_successful(login_response):
                    if verbose:
                        print("Login successful.")
                    return session
            except Exception as e:
                if verbose:
                    print(f"Error decrypting credentials: {str(e)}")
                # Delete invalid credentials from radis.env
                if os.path.exists(env_path):
                    os.remove(env_path)
                print(
                    "Invalid stored credentials. Please enter your HITRAN credentials again."
                )
                # Continue to new user flow

    # First time use or no stored credentials
    username, password = setup_credentials()
    login_response, session = attempt_login(username, password)

    if is_login_successful(login_response):
        if verbose:
            print("Login successful.")
        store_credentials(username, password)
        return session
    else:
        if verbose:
            print(f"Login failed: {login_response.status_code}")
        raise OSError(
            "HITRAN login failed. Please ensure you entered correct credentials from https://hitran.org/login/"
        )


def download_hitemp_file(session, file_url, output_filename, verbose=False):
    print(f"Starting download from {file_url} to {output_filename}")
    file_response = session.get(file_url, stream=True)
    if file_response.status_code == 200:
        total_size = int(file_response.headers.get("content-length", 0))
        print(f"Total size to download: {total_size} bytes")

        with open(output_filename, "wb") as f, tqdm(
            total=total_size, unit="B", unit_scale=True, desc=output_filename
        ) as pbar:
            for chunk in file_response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))

        print("\nDownload complete!")
    else:
        print(f"Download failed: {file_response.status_code}")
        print("Response:", file_response.text[:500])
        raise Warning(
            f"Failed to download {file_url}. Please download manually and place it in the following location:"
        )
        temp_folder = os.path.join(
            os.path.dirname(output_filename),
            "downloads__can_be_deleted",
            "hitran.org",
            "files",
            "HITEMP",
            "HITEMP-2024",
            "CO2_line_list",
        )
        print(f"{file_url} ==> {temp_folder} \n")


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
        r"""
        See Also
        --------
        HITEMPDatabaseManager is compatible with Exojax :py:class:`exojax.spec.api.MdbHitemp`

        """
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
        r"""requires connexion"""

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
            #     headerrow = rowgetDataText(trs[0], 'th')
            #     if headerrow: # if there is a header row include first
            #         rows.append(headerrow)
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
        r"""requires connection"""

        if self.urlnames is not None:
            return self.urlnames

        molecule = self.molecule

        if molecule in ["H2O"]:  # CO2 is a single file since 01/2025

            base_url, Ntotal_lines_expected, _, _ = self.fetch_url_Nlines_wmin_wmax()

            # response = urllib.request.urlopen(base_url)
            # response_string = response.read().decode()
            # inputfiles = re.findall(r'href="(\S+.zip)"', response_string)
            # urlnames = [join(base_url, f) for f in inputfiles]

            from radis.misc.utils import getProjectRoot

            with open(
                join(getProjectRoot(), "db", "H2O", "HITRANpage_january2025.htm")
            ) as file:
                response_string = file.read()

            inputfiles = re.findall(r'href="(\S+.zip)"', response_string)
            base_url = "https://hitran.org"
            urlnames = [f"{base_url}{f}" for f in inputfiles]

        elif molecule in HITEMP_MOLECULES:
            session = login_to_hitran(verbose=self.verbose)
            if session:
                url, Ntotal_lines_expected, _, _ = self.fetch_url_Nlines_wmin_wmax()
                download_hitemp_file(session, url, basename(url))
                urlnames = [basename(url)]
            else:
                return []  # Exit if login failed
        else:
            raise KeyError(
                f"Please choose one of HITEMP molecules : {HITEMP_MOLECULES}. Got '{molecule}'"
            )

        self.urlnames = urlnames

        return urlnames

    def keep_only_relevant(
        self, inputfiles, wavenum_min=None, wavenum_max=None, verbose=True
    ) -> list:
        r"""For CO2 and H2O, return only relevant files for given wavenumber range.

        If other molecule, return the file anyway.
        see :py:func:`radis.api.hitempapi.keep_only_relevant`"""
        if self.molecule in ["H2O"]:  # CO2 is a single file since 01/2025
            inputfiles, _, _ = keep_only_relevant(
                inputfiles, wavenum_min, wavenum_max, verbose
            )
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
        r"""Uncompress ``urlname`` into ``local_file``.
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

        writer = self.get_datafile_manager()

        if molecule == "CO2":
            session = login_to_hitran()
            download_hitemp_file(
                session,
                "https://hitran.org/files/HITEMP/bzip2format/02_HITEMP2024.par.bz2",
                "02_HITEMP2024.par.bz2",
            )
            urlname = "02_HITEMP2024.par.bz2"

        with opener.open(urlname) as gfile:  # locally downloaded file

            dt = _create_dtype(columns, linereturnformat)

            if verbose:
                print(f"Download complete. Parsing {molecule} database to {local_file}")
                print(
                    "This step is executed only ONCE and will considerably accelerate the computation of spectra. It will also dramatically reduce the memory usage. The parsing/conversion can be very fast (e.g. HITEMP OH takes a few seconds) or extremely long (e.g. HITEMP CO2 takes approximately 1 hour)."
                )
            if molecule == "CO2":
                from warnings import warn

                warn(
                    "Parsing will take approximately 1 hour for HITEMP CO2 (compressed = 6 GB",
                    UserWarning,
                )

            # assert not(exists(local_file))

            b = np.zeros(chunksize, dtype=dt)  # receives the HITRAN 160-character data.

            for nbytes in iter(lambda: gfile.readinto(b), 0):

                if not b[-1]:
                    # End of file flag within the chunk (but does not start
                    # with End of file flag) so nbytes != 0
                    b = get_last(b)

                df = _ndarray2df(b, columns, linereturnformat, molecule=self.molecule)

                # 19722
                if molecule == "CO2":
                    df["iso"] = (
                        df["iso"].replace({"A": 10, "B": 11, "C": 12}).astype(int)
                    )  # in HITEMP2024, isotopologue 10, 11, 12 are A, B, C.

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
        r"""register in ~/radis.json"""

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


#%%

if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_hitemp.py"]))
