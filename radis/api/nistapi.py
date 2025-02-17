import urllib

import pandas as pd

from radis.api.dbmanager import DatabaseManager
from radis.misc.config import getDatabankEntries
from radis.misc.warning import DatabaseAlreadyExists


class NISTDatabaseManager(DatabaseManager):
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
        self.wmin = None
        self.wmax = None
        self.downloadable = True

    def fetch_urlnames(self, low_w="", upp_w=""):

        payload = {
            "spectra": self.molecule,  # should accept spectroscopic notation with an underscore as is the current standard for atoms
            "low_w": low_w,  # the lower wavelength bound (not wavenumber, even if using 'show_wn=1' instead of 'show_calc_wl=1')
            "upp_w": upp_w,  # the upper wavelength bound
            "unit": 1,  # the unit for 'low_w' and 'upp_w', 1 for nm
            "format": 3,  # output format, 3 for tab-delimited
            "line_out": 1,  # output only lines with transition probabilities, otherwise it includes lines that don't have any of Aki, fik, Sik, loggf
            "en_unit": 0,  # the unit for energy returned, 0 for cm-1
            "output": 0,  # return output in its entirety rather than in pages
            "show_calc_wl": 1,  # return the Ritz wavelength
            "order_out": 0,  # order output by wavelength
            "show_av": 3,  # use vacuum wavelengths throughout
            "A_out": 0,  # return Einstein A coefficients rather than gA
            "allowed_out": 1,  # include 'Allowed (E1)' transitions
            "forbid_out": 1,  # include 'Forbidden (M1,E2,...)' transitions
            "enrg_out": "on",  # output energies of upper and lower level
            "g_out": "on",  # output statistical weights of upper and lower levels
        }

        query = urllib.parse.urlencode(payload, doseq=True)

        url = urllib.parse.urlunsplit(
            ("https", "physics.nist.gov", "/cgi-bin/ASD/lines1.pl", query, "")
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
        df = nist2df(opener.abspath(urlname))

        df["species"] = self.molecule
        df["iso"] = 0

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

    def register(self, local_files, urlnames):

        if self.is_registered():
            dict_entries = getDatabankEntries(
                self.name
            )  # just update previous register details
        else:
            dict_entries = {}

        files = local_files
        urls = urlnames
        dict_entries.update(
            {
                "path": files,
                "format": "hdf5-radisdb",
                "download_date": self.get_today(),
                "download_url": urls,
            }
        )
        if self.wmin and self.wmax:
            info = f"NIST {self.molecule} lines ({self.wmin:.1f}-{self.wmax:.1f} cm-1)."
            dict_entries.update(
                {
                    "info": info,
                    "wavenumber_min": self.wmin,
                    "wavenumber_max": self.wmax,
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


def nist2df(file):
    df = pd.read_csv(file, sep="\t")

    # based on Kurucz method:
    cond = (df["Ek(cm-1)"] - df["Ei(cm-1)"]) > 0
    icols = ["g_i", "Ei(cm-1)"]
    kcols = ["g_k", "Ek(cm-1)"]
    lcols = ["gl", "El"]
    ucols = ["gu", "Eu"]
    for icol, kcol, lcol, ucol in zip(icols, kcols, lcols, ucols):
        df[lcol] = df[icol].where(cond, df[kcol])
        df[ucol] = df[kcol].where(cond, df[icol])

    df["wav"] = 1e7 / df["ritz_wl_vac(nm)"]
    df[["jl", "ju"]] = (
        df[["gl", "gu"]] - 1
    ) / 2  # calculate J from g rather than requesting from NIST to avoid having to deal with half-integer '1/2' etc notation that NIST returns
    df["A"] = df["Aki(s^-1)"]

    df[:] = df[::-1]

    df = df[["ritz_wl_vac(nm)", "wav", "A", "gl", "El", "gu", "Eu", "jl", "ju"]]

    return df
