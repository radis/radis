# -*- coding: utf-8 -*-
"""
Summary
-----------------------

GEISA database parser

-----------------------


"""


# import re
import time
from collections import OrderedDict
from os.path import abspath, exists, expanduser, getmtime, join

import numpy as np

import radis

try:
    from .cache_files import load_h5_cache_file, save_to_hdf
    from .dbmanager import DatabaseManager
    from .hdf5 import DataFileManager
    from .tools import drop_object_format_columns, parse_hitran_file
except ImportError:
    from radis.io.cache_files import load_h5_cache_file, save_to_hdf
    from radis.io.dbmanager import DatabaseManager
    from radis.io.hdf5 import DataFileManager
    from radis.io.tools import drop_object_format_columns, parse_hitran_file
# from typing import Union


# %% Parsing functions

# General case : GEISA 2020
# Provided by Thibault Delahaye

# BE ADVISED #1: the columns "int" and "ierrB" should be in format FLOAT instead of str as below.
# However, because of Fortran's double-precision floating-point values used in these 2 columns,
# which numpy cannot typecast to float. Thus, we set them as string here and will typecast them
# to float again later in geisa_utils.

# BE ADVISED #2: the columns "idH" for HITRAN molecular number and "isoH" for HITRAN isotope number
# have been changed to "id" and "iso" respectively. Meanwhile, original columns "id" and "iso" have
# also been changed to "idG" and "isoG" to avoid duplication. This is to make sure the parsing process
# and spectrum calculation using HAPI partition function can operate normally.

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
        ("Tdpair", ("a4", float, "temperature-dependance exponent for Gamma air", "")),
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
        ("Tdpself", ("a4", float, "temperature-dependance exponent for self-broadening halfwidth", "")),
        ("ierrS", ("a4", float, "estimated accuracy on the temperature dependence coefficient of the self-broadening halfwidth", "")),
        ("Pshfts", ("a8", float, "self pressure-induced line shift at 296K", "cm-1.atm-1")),
        ("ierrT", ("a8", float, "estimated accuracy on the self-pressure shift of the line transition at 296K", "cm-1.atm-1")),
        ("Tdpnself", ("a4", float, "temperature-dependance exponent for self pressure-induced line shift", "")),
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
    "H2O": 362222,
    "CO2": 532533,
    "O3": 389378,
    "N2O": 50633,
    "CO": 14985,
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

    # Attempt to use cache file
    fcache = DataFileManager(engine).cache_file(fname)
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
    df = parse_hitran_file(fname, columns_GEISA)

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

    def fetch_urlnames(self):
        """requires connexion"""

        if self.urlnames is not None:
            return self.urlnames

        molecule = self.molecule
        geisa_url = "https://aeris-geisa.ipsl.fr/geisa_files/2020/Lines/line_GEISA2020_asc_gs08_v1.0_"

        print(f"Molecule: {molecule}")

        if molecule.upper() in GEISA_MOLECULES:
            url = geisa_url + molecule.lower()
            urlnames = [url]
        else:
            raise KeyError(
                f"Please choose one of GEISA molecules : {GEISA_MOLECULES}. Got '{molecule}'"
            )

        self.urlnames = urlnames

        return urlnames

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

        # molecule = self.molecule
        # Ntotal_lines_expected = GEISA_MOLECULES_Nlines[molecule.upper()]

        writer = self.get_datafile_manager()

        with opener.open(urlname) as gfile:  # locally downloaded file

            gfile  #  so the linter doesn't annoy us. We're not using this file anyway, just unzipping the cache file directly :
            df = gei2df(opener.abspath(urlname), drop_non_numeric=False, cache=False)

            writer.write(local_file, df, append=False)

            self.wmin = df.wav.min()
            self.wmax = df.wav.max()

            Nlines = len(df)

        writer.combine_temp_batch_files(local_file)  # used for vaex mode only

        # Check number of lines is consistent
        # assert Nlines == Ntotal_lines_expected

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
        where to create the RADIS HDF5 files. Default ``"~/.radisdb/geisa"``.
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
        jax arrays. If ``'vaex'``, output is a :py:class:`vaex.dataframe.DataFrameLocal`

        .. note::
            Vaex DataFrames are memory-mapped. They do not take any space in RAM
            and are extremelly useful to deal with the largest databases.

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

        from radis import fetch_geisa
        df = fetch_geisa("CO")
        print(df.columns)
        >>> Index(['wav', 'int', 'airbrd', 'El', 'globu', 'globl', 'locu', 'locl',
            'Tdpgair', 'isoG', 'mol', 'idG', 'id', 'iso', 'A', 'selbrd', 'Pshft',
            'Tdpair', 'ierrA', 'ierrB', 'ierrC', 'ierrF', 'ierrO', 'ierrR', 'ierrN',
            'Tdpgself', 'ierrS', 'Pshfts', 'ierrT', 'Tdppself', 'ierrU'],
            dtype='object')

    .. minigallery:: radis.fetch_geisa

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

    if cache == "regen":
        ldb.remove_local_files(local_files)
    ldb.check_deprecated_files(
        ldb.get_existing_files(local_files),
        auto_remove=True if cache != "force" else False,
    )

    # Get missing files
    download_files = ldb.get_missing_files(local_files)

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

    if isotope and type(isotope) == int:
        isotope = str(isotope)

    df = ldb.load(
        local_files,  # filter other files,
        columns=columns,
        within=[("iso", isotope)] if isotope is not None else [],
        # for relevant files, get only the right range :
        lower_bound=[("wav", load_wavenum_min)] if load_wavenum_min is not None else [],
        upper_bound=[("wav", load_wavenum_max)] if load_wavenum_max is not None else [],
        output=output,
    )

    return (df, local_files) if return_local_path else df


#%%

if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_geisa.py"]))
