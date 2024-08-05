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
from os.path import abspath, expanduser, join, splitext

import numpy as np
import pandas as pd
import periodictable

from radis.api.dbmanager import DatabaseManager
from radis.db.classes import get_ielem_charge, roman_to_int
from radis.misc.config import getDatabankEntries
from radis.misc.utils import getProjectRoot
from radis.misc.warning import DatabaseAlreadyExists
from radis.phys.air import air2vacuum
from radis.phys.constants import c_CGS, e_CGS, m_e_CGS


def pick_ionE(ielem, charge):
    # This method was extracted from exojax/src/exojax/spec/atomllapi.py
    # (https://github.com/HajimeKawahara/exojax.git)

    """Pick out the ionization energy of a specific atomic species from NIST.

    Parameters
    ----------
    ielem: int
        atomic number (e.g., Fe=26)
    charge: int
        charge of the species

    Returns
    -------
    ionE: float
        ionization energy

    """

    fn_IonE = os.path.join(
        getProjectRoot(), "db", "NIST_Atomic_Ionization_Energies.txt"
    )  # pkgutil.get_data(
    #     "exojax", "data/atom/NIST_Atomic_Ionization_Energies.txt"
    # )
    df_ionE = pd.read_csv(fn_IonE, sep="|", skiprows=6, header=0)

    def f_droppare(x):
        return (
            x.str.replace("\(", "", regex=True)
            .str.replace("\)", "", regex=True)
            .str.replace("\[", "", regex=True)
            .str.replace("\]", "", regex=True)
            .str.replace("                                      ", "0", regex=True)
        )

    ionE = f_droppare(
        df_ionE[(df_ionE["At. num "] == ielem) & (df_ionE[" Ion Charge "] == charge)][
            "      Ionization Energy (a) (eV)      "
        ]
    )

    assert ionE.size == 1

    return float(ionE.iloc[0])


def get_atomic_number(species):
    """
    Extracts the atomic number from the species given in spectroscopic notation
    """
    atomic_symbol = species.split("_")[0]
    el = getattr(periodictable, atomic_symbol)

    atomic_number = el.number
    return atomic_number


def get_ionization_state(species):
    """
    Extracts the ionization state from the species given in spectroscopic notation
    """
    ionization_str = species.split("_")[1]
    ionization_int = roman_to_int(ionization_str) - 1
    # formatted_str = f"{ionization_int:02}"

    return ionization_int  # formatted_str


def get_element_symbol(species):
    """
    Extracts the element symbol from the species given in spectroscopic notation
    """

    atomic_symbol = species.split("_")[0]
    el = getattr(periodictable, atomic_symbol)
    return el


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
    df = pd.read_csv(file, sep="\s+", header=2)
    df.set_index("LOG10(T)", inplace=True)
    return df


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
    pfTdat = pd.read_csv(io.StringIO(pfT_str), sep="\s+")
    pfTdat = pd.Series(pfTdat.columns[1:]).astype(
        "float64"
    )  # Converts the values to float64, skipping the first value

    with open(os.path.join(getProjectRoot(), "db", "kuruczpartfn.txt"), "r") as f:
        pfdat = pd.read_csv(f, sep="\s+", comment="#", names=pfTdat.index)

    return pfTdat, pfdat


class KuruczDatabaseManager(DatabaseManager):
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
        self.wmin = None
        self.wmax = None

        self.atomic_number = f"{get_atomic_number(molecule):02}"
        self.ionization_state = f"{get_ionization_state(molecule):02}"

        self.actual_file = None
        self.actual_url = None
        self.pf_path = None

    def fetch_urlnames(self):
        """returns the possible urls of the desired linelist file, pending confirmation upon attempting download as to which is correct"""

        code = f"{self.atomic_number}{self.ionization_state}"

        urlnames = [
            "http://kurucz.harvard.edu/atoms/"
            + code
            + "/gf"
            + code
            + ".all",  # 1. new linelist including lab lines
            "http://kurucz.harvard.edu/atoms/"
            + code
            + "/gf"
            + code
            + ".pos",  # 2. new linelist without lab lines
            "http://kurucz.harvard.edu/linelists/gfall/gf"
            + code
            + ".all",  # 3. old linelist including lab lines
        ]

        return urlnames

    def get_pf_path(self):
        """returns the url at which the dedicated partition function file is expected to be located, and the derived file path at which it would be saved after downloading and parsing"""
        code = f"{self.atomic_number}{self.ionization_state}"
        url = "http://kurucz.harvard.edu/atoms/" + code + "/partfn" + code + ".dat"
        fname = splitext(url.split("/")[-1])[0] + ".h5"
        path = expanduser(
            abspath(join(self.local_databases, self.molecule + "-" + fname))
        )
        return path, url

    def get_possible_files(self, urlnames=None):
        """returns the urls from fretch_urlnames and the derived file paths at which each would be saved after downloading and parsing"""
        verbose = self.verbose
        local_databases = self.local_databases
        engine = self.engine

        # copied from DatabaseManager.get_filenames:
        if urlnames is None:
            urlnames = self.fetch_urlnames()
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

        return local_files, urlnames

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
        if local_file == self.pf_path:
            df = read_kurucz_pfs(opener.abspath(urlname))
            with pd.HDFStore(local_file, mode="w", complib="blosc", complevel=9) as f:
                f.put(key="df", value=df)
            Nlines = len(df)  # irrelevant
        else:
            df = read_kurucz(opener.abspath(urlname))

            writer.write(local_file, df, append=False)

            self.wmin = df.wav.min()
            self.wmax = df.wav.max()

            Nlines = len(df)

            writer.combine_temp_batch_files(local_file)

            # Add metadata
            from radis import __version__

            writer.add_metadata(
                local_file,
                {
                    "wavenumber_min": self.wmin,
                    "wavenumber_max": self.wmax,
                    "download_date": self.get_today(),
                    "download_url": urlname,
                    "total_lines": Nlines,
                    "version": __version__,
                },
            )

        return Nlines

    def register(self, get_main_files, get_pf_files):

        if self.is_registered():
            dict_entries = getDatabankEntries(
                self.name
            )  # just update previous register details
        else:
            dict_entries = {}

        if get_main_files:
            files = [self.actual_file]
            urls = [self.actual_url]
            dict_entries.update(
                {
                    "path": files,
                    "format": "hdf5-radisdb",
                    "download_date": self.get_today(),
                    "download_url": urls,
                }
            )
            if self.wmin and self.wmin:
                info = f"Kurucz {self.molecule} lines ({self.wmin:.1f}-{self.wmax:.1f} cm-1)."
                dict_entries.update(
                    {
                        "info": info,
                        "wavenumber_min": self.wmin,
                        "wavenumber_max": self.wmax,
                    }
                )
        # else:
        #     raise Exception("Kurucz database can't be registered until the correct url to use is known")

        if get_pf_files:
            # files.append(self.pf_path)
            # urls.append(self.pf_url)
            dict_entries.update(
                {
                    "parfunc": self.pf_path,
                    "parfuncfmt": "kurucz",
                }
            )

        dict_entries.update(
            {
                "parfuncfmt": "kurucz",
            }
        )

        try:
            super().register(dict_entries)
        except DatabaseAlreadyExists as e:
            raise Exception(
                'If you want RADIS to overwrite the existing entry for a registered databank, set the config option "ALLOW_OVERWRITE" to True.'
            ) from e


def read_kurucz(kuruczf, preserve_orig_levels=False):
    """
    Parse a Kurucz linelist, process its columns as required, and return as a Pandas DataFrame

    Inspired by: https://github.com/rasmus98/NLTE-Helium/blob/fc6161a30ebecfcf59f36386b1dc7f02ff749905/Grotrian_diagram_Helium.ipynb cell 2. Also see:

    - https://github.com/DBerke/varconlib/blob/e57250ca359026ae8b8059dae179fb0ad9625aa2/varconlib/scripts/select_line_pairs.py#L711

    - https://github.com/followthesheep/specutils/blob/e2873719f3820e83ab29c79062d7a59a7664fa2f/specutils/read_kurucz_linelist.py

    - https://github.com/ajwheeler/Korg.jl/blob/53420d38c23e21a7fe202b8680bdb201c9a62a2a/src/linelist.jl#L153

    Parameters
    ----------
    kuruczf: str
        file path
    preserve_orig_levels: boolean
        whether to preserve the original columns pertaining to the levels prior to transforming them, with '_orig' appended to their name; note the columns are still reversed so that the database index increases in wavenumber rather than wavelength

    Returns
    ----------
        Pandas DataFrame containing the required columns for the Kurucz linelist database
    """
    colspecs = (
        (0, 11),
        (11, 18),
        (18, 24),
        (24, 36),
        (36, 41),
        (42, 52),
        (52, 64),
        (64, 69),
        (70, 80),
        (80, 86),
        (86, 92),
        (92, 98),
        # (98,102),(102,104),(104,106),
        (106, 109),
    )
    # (109,115),(115,118),(118,124),(124,129),(129,134),(135,136),(136,137),(138,139),(139,140),(140,141),(141,144),(144,149),(149,154),(154,160))
    dtypes = (
        np.float64,
        np.float64,
        "object",  # make species object rather than float for easier parsing as a string
        np.float64,
        np.float64,
        "object",
        np.float64,
        np.float64,
        "object",
        np.float64,
        np.float64,
        np.float64,
        # str,np.int64,np.int64,
        np.int64,
    )
    # np.float64,np.int64,np.float64,np.int64,np.int64,str,str,str,str,np.int64,str,np.int64,np.int64,np.int64)
    colnames = (
        "orig_wavelen",
        "loggf",
        "species",
        "elower_orig",
        "jlower_orig",
        "labellower_orig",
        "eupper_orig",
        "jupper_orig",
        "labelupper_orig",
        "gamRad",
        "gamSta",
        "gamvdW",
        #'ref','NLTElo','NLTEhi',
        "iso",
    )
    #'hyperstrength','isotope2','isoabund','hyperlo','hyperhi','hyperFlo','hypernotelo','hyperFhi','hypernotehi','strengthclassnote','trans_code','lande_g_even','lande_g_odd','iso_shift')
    df = pd.read_fwf(
        kuruczf,
        colspecs=colspecs,
        dtype=dict(zip(colnames, dtypes)),
        names=colnames,
        skip_blank_lines=True,
    )  # , converters={'labellower_orig': str.strip, 'labelupper_orig': str.strip}) # appears to strip automatically
    species_unique = df["species"].unique()
    assert len(species_unique) == 1
    ielem, charge = get_ielem_charge(species_unique[0])
    df["ionE"] = pick_ionE(ielem, charge)

    # simple method to preserve dtypes of each column:
    cond = (df["eupper_orig"] - df["elower_orig"]) > 0
    oldlower = ["jlower_orig", "labellower_orig", "elower_orig"]
    oldupper = ["jupper_orig", "labelupper_orig", "eupper_orig"]
    newlower = ["jl", "labellower", "El"]
    newupper = ["ju", "labelupper", "Eu"]
    for oldcol1, oldcol2, newcol1, newcol2 in zip(
        oldlower, oldupper, newlower, newupper
    ):
        df[newcol1] = df[oldcol1].where(cond, df[oldcol2])
        df[newcol2] = df[oldcol2].where(cond, df[oldcol1])
    # df1[['jlower_orig', 'jupper_orig', 'labellower_orig', 'labelupper_orig', 'elower_orig', 'eupper_orig']] = df[['jupper_orig', 'jlower_orig', 'labelupper_orig', 'labellower_orig', 'eupper_orig', 'elower_orig']]
    # df.where(cond, df1, inplace=True) # returns TypeError: Cannot do inplace boolean setting on mixed-types with a non np.nan value
    # np.where(pd.concat([cond]*6, axis=1), df[['jlower_orig', 'jupper_orig', 'labellower_orig', 'labelupper_orig', 'elower_orig', 'eupper_orig']], df[['jupper_orig', 'jlower_orig', 'labelupper_orig', 'labellower_orig', 'eupper_orig', 'elower_orig']]) # doesn't preserve dtype

    df["wav"] = 1e7 / df["orig_wavelen"].where(
        df["orig_wavelen"] < 200, air2vacuum(df["orig_wavelen"])
    )
    df["gu"] = df["ju"] * 2 + 1
    df["A"] = (
        10 ** df["loggf"]
        / df["gu"]
        * (c_CGS * df["wav"]) ** 2
        * (8 * np.pi**2 * e_CGS**2)
        / (m_e_CGS * c_CGS**3)
    )

    df[:] = df[::-1]

    if not preserve_orig_levels:
        df.drop(
            [
                "elower_orig",
                "eupper_orig",
                "jlower_orig",
                "jupper_orig",
                "labellower_orig",
                "labelupper_orig",
            ],
            axis=1,
            inplace=True,
        )

    return df
