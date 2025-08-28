# -*- coding: utf-8 -*-
"""
Summary
-----------------------

Kurucz database parser

Based largely on the `Exojax <https://github.com/HajimeKawahara/exojax>`__ code

-----------------------


"""

import io
import os
import urllib
from os.path import abspath, expanduser, join, splitext

import pandas as pd

import radis
from radis.api.dbmanager import DatabaseManager
from radis.api.kuruczapi import get_atomic_number, get_ionization_state
from radis.misc.utils import getProjectRoot


class KuruczPFManager(DatabaseManager):
    def __init__(
        self,
        name,
        molecule,
        local_databases,
        engine="default",
        verbose=True,
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

        self.atomic_number = f"{get_atomic_number(molecule):02}"
        self.ionization_state = f"{get_ionization_state(molecule):02}"

    def get_pf_path(self):
        """returns the url at which the dedicated partition function file is expected to be located, and the derived file path at which it would be saved after downloading and parsing"""
        code = f"{self.atomic_number}{self.ionization_state}"
        url = "http://kurucz.harvard.edu/atoms/" + code + "/partfn" + code + ".dat"
        fname = splitext(url.split("/")[-1])[0] + ".h5"
        path = expanduser(
            abspath(join(self.local_databases, self.molecule + "-" + fname))
        )
        return path, url

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
        """overwrites parent class method as required"""

        writer = self.get_datafile_manager()
        df = read_kurucz_pfs(opener.abspath(urlname))
        writer.write(local_file, df, append=False)

        writer.combine_temp_batch_files(local_file)

        # Add metadata
        from radis import __version__

        writer.add_metadata(
            local_file,
            {
                "download_date": self.get_today(),
                "download_url": urlname,
                "version": __version__,
            },
        )
        Nlines = len(df)  # irrelevant
        return Nlines


def read_kurucz_pfs(file):
    """Convert a Kurucz partfnxxyy.dat file containing a table of partition functions by temperature and potential lowering to a Pandas DataFrame

    Parameters
    ----------
    fname: str
        file name

    Returns
    ----------
    df: Pandas DataFrame
    """
    df = pd.read_csv(file, sep=r"\s+", header=2)
    df.set_index("LOG10(T)", inplace=True)
    return df


def fetch_kurucz_pfs(
    molecule,
    local_databases=None,
    databank_name="Kurucz-PF-{molecule}",
    cache=True,
    verbose=True,
    clean_cache_files=True,
    engine="default",
    output="pandas",
    parallel=True,
):
    """
    See e.g. :func:`~radis.io.hitemp.fetch_hitemp` for an explanation of the parameters largely applicable to `fetch_kurucz_pfs`
    """

    # largely based on :py:func:`~radis.io.fetch_geisa`
    if r"{molecule}" in databank_name:
        databank_name = databank_name.format(**{"molecule": molecule})

    if local_databases is None:

        local_databases = join(
            radis.config["DEFAULT_DOWNLOAD_PATH"],
            "atomic_partition_functions",
            "kurucz",
        )
    local_databases = abspath(expanduser(local_databases))

    ldb = KuruczPFManager(
        databank_name,
        molecule=molecule,
        local_databases=local_databases,
        verbose=verbose,
        parallel=parallel,
        engine=engine,
    )

    path, url = ldb.get_pf_path()

    if cache == "regen":
        ldb.remove_local_files([path])

    if ldb.get_missing_files([path]):
        try:
            ldb.download_and_parse([url], [path], 1)
        except OSError:
            print("a partition function file specific to this species was not found")
            return None, None

    if clean_cache_files:
        ldb.clean_download_files()

    df = ldb.load(
        [path],
        output=output,
    )
    return df, path


def load_pf_Barklem2016():
    """Load a table of the partition functions for 284 atomic species.

    Returns
    -------
    pfTdat:  pd.DataFrame
        Steps of temperature (K)
    pfdat:  pd.DataFrame
        Partition functions for 284 atomic species

    References
    ----------
    `Barklem & Collet (2016), Table 8 <https://doi.org/10.1051/0004-6361/201526961>`_

    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, "./../db/pfTKurucz_values.txt")
    with open(file_path, "r") as file:
        pfT_str = file.read()
    pfTdat = pd.read_csv(io.StringIO(pfT_str), sep=r"\s+")
    pfTdat = pd.Series(pfTdat.columns[1:]).astype(
        "float64"
    )  # Converts the values to float64, skipping the first value

    with open(os.path.join(getProjectRoot(), "db", "kuruczpartfn.txt"), "r") as f:
        pfdat = pd.read_csv(f, sep=r"\s+", comment="#", names=pfTdat.index)

    return pfTdat, pfdat


class NISTPFManager(DatabaseManager):
    def __init__(
        self,
        name,
        molecule,
        local_databases,
        engine="default",
        verbose=True,
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
        self.downloadable = True

    def fetch_urlnames(self):

        payload = {
            "spectrum": self.molecule,  # should accept spectroscopic notation with an underscore as is the current standard for atoms
            "units": 0,  # the unit for energy, 1 for cm-1
            "format": 3,  # output format, 3 for tab-delimited
            "output": 0,  # return output in its entirety rather than in pages
            "level_out": "on",  # output energies of each level
            "g_out": "on",  # output statistical weights each level
        }

        query = urllib.parse.urlencode(payload, doseq=True)

        url = urllib.parse.urlunsplit(
            ("https", "physics.nist.gov", "/cgi-bin/ASD/energy1.pl", query, "")
        )

        return [url]

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
        """overwrites parent class method as required"""

        writer = self.get_datafile_manager()
        from io import StringIO

        # Use the opener's open method which should be available in all implementations
        with opener.open() as file:
            file_content = file.read().decode("utf-8")
        file = StringIO(file_content)

        df = read_NIST_pfs(file)
        writer.write(local_file, df, append=False)

        writer.combine_temp_batch_files(local_file)

        # Add metadata
        from radis import __version__

        writer.add_metadata(
            local_file,
            {
                "download_date": self.get_today(),
                "download_url": urlname,
                "version": __version__,
            },
        )
        Nlines = len(df)  # irrelevant
        return Nlines


def read_NIST_pfs(file):
    df = pd.read_csv(file, sep="\t", index_col=False)
    df.rename(columns={"Level (cm-1)": "E"}, inplace=True)
    df["E"] = df["E"].astype("float")
    return df


def fetch_NIST_pfs(
    molecule,
    local_databases=None,
    databank_name="NIST-PF-{molecule}",
    cache=True,
    verbose=True,
    clean_cache_files=True,
    engine="default",
    output="pandas",
    parallel=True,
):
    """
    See e.g. :func:`~radis.io.hitemp.fetch_hitemp` for an explanation of the parameters largely applicable to `fetch_NIST_pfs`
    """

    # largely based on :py:func:`~radis.io.fetch_geisa`
    if r"{molecule}" in databank_name:
        databank_name = databank_name.format(**{"molecule": molecule})

    if local_databases is None:

        local_databases = join(
            radis.config["DEFAULT_DOWNLOAD_PATH"], "atomic_partition_functions", "NIST"
        )
    local_databases = abspath(expanduser(local_databases))

    ldb = NISTPFManager(
        databank_name,
        molecule=molecule,
        local_databases=local_databases,
        verbose=verbose,
        parallel=parallel,
        engine=engine,
    )

    local_files, urlnames = ldb.get_filenames()

    if cache == "regen":
        ldb.remove_local_files(local_files)

    if ldb.get_missing_files(local_files):
        ldb.download_and_parse(urlnames, local_files, 1)

    if clean_cache_files:
        ldb.clean_download_files()

    df = ldb.load(
        local_files,
        output=output,
    )

    return df, local_files[0]
