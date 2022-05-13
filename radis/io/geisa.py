# -*- coding: utf-8 -*-
"""
Summary
-----------------------

GEISA database parser

-----------------------


"""


import re
import time
from collections import OrderedDict
from os.path import abspath, basename, exists, expanduser, getmtime, join
from typing import Union

import numpy as np

import radis
from radis.io.cache_files import cache_file_name, load_h5_cache_file, save_to_hdf
from radis.io.dbmanager import DatabaseManager
from radis.io.tools import (
    _create_dtype,
    _get_linereturnformat,
    _ndarray2df,
    drop_object_format_columns,
    parse_hitran_file,
)
from radis.misc.progress_bar import ProgressBar

from .hdf5 import update_pytables_to_vaex

# %% Parsing functions

# General case : GEISA 2020
# Provided by Thibault Delahaye

# BE ADVISED: the columns "int" and "ierrB" should be in format FLOAT instead of str as below.
# However, because of Fortran's double-precision floating-point values used in these 2 columns,
# which numpy cannot typecast to float. Thus, we set them as string here and will typecast them
# to float again later in geisa_utils.

# fmt: off
columns_GEISA = OrderedDict(
    [
        # name  # format  # type  # description  # unit
        ("wav", ("a12", float, "vacuum wavenumber", "cm-1")),
        ("int", ("a11", str, "intensity at 296K", "cm-1/(molecule/cm-2)")),
        ("airbrd", ("a6", float, "air-broadened half-width at 296K", "cm-1.atm-1")),
        ("El", ("a10", float, "lower-state energy", "cm-1")),
        ("globu", ("a25", str, "electronic and vibrational global upper quanta", "")),
        ("globl", ("a25", str, "electronic and vibrational global lower quanta", "")),
        ("locu", ("a15", str, "electronic and vibrational local upper quanta", "")),
        ("locl", ("a15", str, "electronic and vibrational local lower quanta", "")),
        ("Tdpgair", ("a4", float, "temperature-dependance exponent for Gamma air", "")),
        ("isoG", ("a3", int, "GEISA isotope number", "")),
        ("mol", ("a3", int, "GEISA molecular number", "")),
        ("idG", ("a3", str, "Internal GEISA code for the data identification", "")),
        ("id", ("a2", int, "Hitran molecular number", "")),
        ("iso", ("a1", int, "Hitran isotope number", "")),
        ("A", ("a10", float, "Einstein A coefficient", "s-1")),
        ("selbrd", ("a7", float, "self-broadened half-width at 296K", "cm-1.atm-1")),
        ("Pshft", ("a9", float, "air pressure-induced line shift at 296K", "cm-1.atm-1")),
        ("Tdppair", ("a6", float, "temperature-dependance exponent for air pressure-induced line shift", "")),
        ("ierrA", ("a10", float, "estimated accuracy on the line position", "cm-1")),
        ("ierrB", ("a11", str, "estimated accuracy on the intensity of the line", "cm-1/(molecule/cm-2)")),
        ("ierrC", ("a6", float, "estimated accuracy on the air collision halfwidth", "cm-1.atm-1")),
        ("ierrF", ("a4", float, "estimated accuracy on the temperature dependence coefficient of the air-broadening halfwidth", "")),
        ("ierrO", ("a9", float, "estimated accuracy on the air pressure shift of the line transition at 296K", "cm-1.atm-1")),
        ("ierrR", ("a6", float, "estimated accuracy on the temperature dependence coefficient of the air pressure shift", "")),
        ("ierrN", ("a7", float, "estimated accuracy on the self-broadened at 296K", "cm-1.atm-1")),
        ("Tdpgself", ("a4", float, "temperature-dependance exponent for self-broadening halfwidth", "")),
        ("ierrS", ("a4", float, "estimated accuracy on the temperature dependence coefficient of the self-broadening halfwidth", "")),
        ("Pshfts", ("a8", float, "self pressure-induced line shift at 296K", "cm-1.atm-1")),
        ("ierrT", ("a8", float, "estimated accuracy on the self-pressure shift of the line transition at 296K", "cm-1.atm-1")),
        ("Tdppself", ("a4", float, "temperature-dependance exponent for self pressure-induced line shift", "")),
        ("ierrU", ("a4", float, "estimated accuracy on the temperature dependence coefficient of the self pressure shift", "")),
    ]
)
""" OrderedDict: parsing order of GEISA2020 format """

# The columns required for equilibrium calculations only
equilibrium_columns = OrderedDict(
    [
        # name  # format  # type  # description  # unit
        ("wav", ("a12", float, "vacuum wavenumber", "cm-1")),
        ("int", ("a11", float, "intensity at 296K", "cm-1/(molecule/cm-2)")),
        ("airbrd", ("a6", float, "air-broadened half-width at 296K", "cm-1.atm-1")),
        ("iso", ("a3", int, "GEISA isotope number", "")),
        ("mol", ("a3", int, "GEISA molecular number", "")),
        ("id", ("a3", str, "Internal GEISA code for the data identification", "")),
        ("selbrd", ("a7", float, "self-broadened half-width at 296K", "cm-1.atm-1")),
        ("Pshft", ("a9", float, "air pressure-induced line shift at 296K", "cm-1.atm-1")),
        ("Pshfts", ("a8", float, "self pressure-induced line shift at 296K", "cm-1.atm-1")),
    ]
)
# fmt: on

GEISA_MOLECULES_Nlines = {
    "H2O": 75699,
    "CO2": 532533,
    "O3": 389378,
    "N2O": 50633,
    "CO": 13515,
    "CH4": 240858,
    "O2": 6428,
    "NO": 105079,
    "SO2": 68728,
    "NO2": 113880,
    "NH3": 29082,
    "PH3": 34542,
    "HNO3": 691161,
    "OH": 42866,
    "HF": 20010,
    "HCL": 35985,
    "HBR": 8980,
    "HI": 4751,
    "CLO": 7230,
    "OCS": 33809,
    "H2CO": 37050,
    "C2H6": 28439,
    "CH3D": 49237,
    "C2H2": 11340,
    "C2H4": 18378,
    "GEH4": 32372,
    "HCN": 81889,
    "C3H8": 8983,
    "C2N2": 2577,
    "C4H2": 119480,
    "HC3N": 179347,
    "HOCL": 17862,
    "N2": 120,
    "CH3CL": 18344,
    "H2O2": 126983,
    "H2S": 20788,
    "HCOOH": 62684,
    "COF2": 70904,
    "SF6": 92398,
    "C3H4": 19001,
    "HO2": 38804,
    "CLONO2": 356899,
    "CH3BR": 36911,
    "CH3OH": 19897,
    "NO+": 1206,
    "HNC": 5619,
    "C6H6": 9797,
    "C2HD": 15512,
    "CF4": 60033,
    "CH3CN": 17172,
    "HDO": 63641,
    "SO3": 10881,
    "HONO": 26041,
    "COFCL": 215639,
    "CH3I": 70291,
    "CH3F": 1499,
    "RUO4": 30205,
    "H2C3H2": 31686,
}

GEISA_MOLECULES = list(GEISA_MOLECULES_Nlines.keys())


def gei2df(
    fname,
    cache=True,
    load_columns=None,
    verbose=True,
    drop_non_numeric=True,
    load_wavenum_min=None,
    load_wavenum_max=None,
    engine="pytables",
):
    """Convert a GEISA [1]_ file to a Pandas dataframe.
    Parameters
    ----------
    fname: str
        GEISA file name.
    cache: boolean, or 'regen'
        if ``True``, a pandas-readable HDF5 file is generated on first access,
        and later used. This saves on the datatype cast and conversion and
        improves performances a lot (but changes in the database are not
        taken into account). If ``False``, no database is used. If 'regen', temp
        file are reconstructed. Default ``True``.
    load_columns: list
        columns to load. If ``None``, loads everything.
        .. note::
            this is only relevant when loading from a cache file. To generate
            the cache file, all columns are loaded anyway.
    Other Parameters
    ----------------
    drop_non_numeric: boolean
        if ``True``, non numeric columns are dropped. This improves performances,
        but make sure all the columns you need are converted to numeric formats
        before hand. Default ``True``. Note that if a cache file is loaded it
        will be left untouched.
    load_wavenum_min, load_wavenum_max: float
        if not ``'None'``, only load the cached file if it contains data for
        wavenumbers above/below the specified value. See :py:func`~radis.io.cache_files.load_h5_cache_file`.
        Default ``'None'``.
    engine: 'pytables', 'vaex'
        format for Hdf5 cache file, `pytables` by default.
    Returns
    -------
    df: pandas Dataframe
        dataframe containing all lines and parameters.
    Notes
    -----
    GEISA Database 2020 release can be downloaded from [2]_
    References
    ----------
    .. [1] `The 2020 edition of the GEISA spectroscopic database, Thibault Delahaye et al., 2021 <https://www.sciencedirect.com/science/article/abs/pii/S0022285221000928>`_
    .. [2] `GEISA Database 2020 release <https://geisa.aeris-data.fr/interactive-access/?db=2020&info=ftp>`_
    See Also
    --------
    :func:`~radis.io.hitran.hit2df`
    :func:`~radis.io.cdsd.cdsd2df`
    """

    # Try to access last modification time of original file
    metadata = {}
    metadata["last_modification"] = time.ctime(getmtime(fname))
    # Make sure the wavenums are legit
    if load_wavenum_min and load_wavenum_max:
        assert load_wavenum_min < load_wavenum_max
    # Verbose
    if verbose >= 2:
        print("Opening file {0}, cache={1})".format(fname, cache))
        print("Last Modification time: {0}".format(metadata["last_modification"]))

    parse_columns = columns_GEISA

    # Attempt to use cache file
    fcache = cache_file_name(fname, engine=engine)
    if cache and exists(fcache):
        relevant_if_metadata_above = (
            {"wavenum_max": load_wavenum_min} if load_wavenum_min else {}
        )  # not relevant if wavenum_max of file is < wavenum min required
        relevant_if_metadata_below = (
            {"wavenum_min": load_wavenum_max} if load_wavenum_max else {}
        )  # not relevant if wavenum_min of file is > wavenum max required
        df = load_h5_cache_file(
            fcache,
            cache,
            columns=load_columns,
            valid_if_metadata_is=metadata,
            relevant_if_metadata_above=relevant_if_metadata_above,
            relevant_if_metadata_below=relevant_if_metadata_below,
            current_version=radis.__version__,
            last_compatible_version=radis.config["OLDEST_COMPATIBLE_VERSION"],
            verbose=verbose,
            engine=engine,
        )
        if df is not None:
            return df

    # If cache files are not found, commence reading of full file
    df = parse_hitran_file(fname, parse_columns)

    # Commence "D to E" conversion on 2nd and 20th columns, by using
    # df.columns.str.replace() which is much faster than Python loop.
    # Finally, typecast them to float.
    df["int"] = df["int"].astype(str).str.replace("D", "E").astype(float)
    df["ierrB"] = df["ierrB"].astype(str).str.replace("D", "E").astype(float)

    # Remove non numerical attributes
    if drop_non_numeric:
        # The function replace_PQR_with_m101 below will be developed
        # in later updates for GEISA's non-equilibrium calculations
        # replace_PQR_with_m101(df)

        # Basically ripping out 4 quanta columns
        df = drop_object_format_columns(df, verbose=verbose)

    # Generate cache file for later use
    if cache:
        new_metadata = {
            # Last modification time of the original file :
            "last_modification": time.ctime(getmtime(fname)),
            "wavenum_min": df.wav.min(),
            "wavenum_max": df.wav.max(),
        }
        if verbose:
            print(
                "Generating cache file {0} with metadata :\n{1}".format(
                    fcache, new_metadata
                )
            )
        try:
            save_to_hdf(
                df,
                fcache,
                metadata=new_metadata,
                version=radis.__version__,
                key="df",
                overwrite=True,
                verbose=verbose,
                engine=engine,
            )
        except PermissionError:
            if verbose:
                print(
                    "An error occured in cache file generation. Lookup access rights."
                )
            pass

    return df


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


class GEISADatabaseManager(DatabaseManager):
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
        self.wmin = None
        self.wmax = None
        self.urlnames = None

    def fetch_url_Nlines(
        self,
        geisa_url="https://aeris-geisa.ipsl.fr/geisa_files/2020/Lines/line_GEISA2020_asc_gs08_v1.0_",
    ):
        """requires connexion"""

        molecule = self.molecule

        if self.base_url is not None and self.Nlines is not None:
            return self.base_url, self.Nlines

        else:

            Nlines = GEISA_MOLECULES_Nlines[molecule]
            url = geisa_url + molecule.lower()

            self.base_url, self.Nlines = url, Nlines

        print(url)
        print(Nlines)

        return url, Nlines

    def fetch_urlnames(self):
        """requires connexion"""

        if self.urlnames is not None:
            return self.urlnames

        molecule = self.molecule

        print(f"Molecule: {molecule}")

        if molecule.upper() in GEISA_MOLECULES:
            url, Ntotal_lines_expected = self.fetch_url_Nlines
            urlnames = [url]
        else:
            raise KeyError(
                f"Please choose one of GEISA molecules : {GEISA_MOLECULES}. Got '{molecule}'"
            )

        self.urlnames = urlnames

        return urlnames

    def keep_only_relevant(
        self,
        inputfiles,
        wavenum_min=None,
        wavenum_max=None,
    ) -> list:

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
        columns = columns_GEISA
        chunksize = self.chunksize
        verbose = self.verbose
        molecule = self.molecule

        if not verbose:
            pbar_active = False

        linereturnformat = self.get_linereturn_format(opener, urlname, columns)

        Nlines = 0
        Nlines_raw = 0
        Nlines_tot = Nlines + pbar_Nlines_already
        _, Ntotal_lines_expected = self.fetch_url_Nlines
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
        info = f"GEISA {self.molecule} lines ({self.wmin:.1f}-{self.wmax:.1f} cm-1)"

        dict_entries = {
            "info": info,
            "path": local_files,
            "format": "geisa-radisdb",
            "parfuncfmt": "hapi",
            "wavenumber_min": self.wmin,
            "wavenumber_max": self.wmax,
            "download_date": self.get_today(),
            "download_url": urlnames,
        }

        super().register(dict_entries)


def fetch_geisa(
    molecule,
    local_databases=None,
    databank_name="GEISA-{molecule}",
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
    output="pandas",
    parallel=True,
):
    """Stream GEISA file from GEISA website. Unzip and build a HDF5 file directly.

    Returns a Pandas DataFrame containing all lines.

    Parameters
    ----------
    molecule: all 58 GEISA 2020 molecules. See here https://geisa.aeris-data.fr/interactive-access/?db=2020&info=ftp
    local_databases: str
        where to create the RADIS HDF5 files. Default ``"~/.radisdb/hitemp"``.
        Can be changed in ``radis.config["DEFAULT_DOWNLOAD_PATH"]`` or in ~/radis.json config file
    databank_name: str
        name of the databank in RADIS :ref:`Configuration file <label_lbl_config_file>`
        Default ``"GEISA-{molecule}"``
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
        which HDF5 library to use to parse local files. If 'default' use the value from ~/radis.json
    output: 'pandas', 'vaex', 'jax'
        format of the output DataFrame. If ``'jax'``, returns a dictionary of
        jax arrays.
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
    :py:func:`~radis.io.hitran.fetch_hitran`, :py:func:`~radis.io.exomol.fetch_exomol`
    :py:func:`~radis.io.hitemp.fetch_hitemp`, :py:func:`~radis.io.hdf5.hdf2df`
    :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank`

    """

    if r"{molecule}" in databank_name:
        databank_name = databank_name.format(**{"molecule": molecule})

    if local_databases is None:
        import radis

        local_databases = join(radis.config["DEFAULT_DOWNLOAD_PATH"], "geisa")
    local_databases = abspath(expanduser(local_databases))

    ldb = GEISADatabaseManager(
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
        output=output,
    )

    return (df, files_loaded) if return_local_path else df


#%%

if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_geisa.py"]))
