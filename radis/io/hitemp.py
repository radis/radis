# -*- coding: utf-8 -*-
"""
Created on Sun May 22 17:35:05 2022

@author: erwan
"""

from os.path import abspath, exists, expanduser, join

from radis import config
from radis.api.hdf5 import update_pytables_to_vaex
from radis.api.hitempapi import HITEMPDatabaseManager, get_recent_hitemp_database_year
from radis.misc.config import getDatabankEntries


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
    output="pandas",
    parallel=True,
    database="most_recent",
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
        which HDF5 library to use to parse local files. If 'default' use the value from ~/radis.json
    output: 'pandas', 'vaex', 'jax'
        format of the output DataFrame. If ``'jax'``, returns a dictionary of
        jax arrays. If ``'vaex'``, output is a :py:class:`vaex.dataframe.DataFrameLocal`

        .. note::
            Vaex DataFrames are memory-mapped. They do not take any space in RAM
            and are extremely useful to deal with the largest databases.

    parallel: bool
        if ``True``, uses joblib.parallel to load database with multiple processes
    database: ``str``
        The database version to retrieve. Options include:
        - `"most_recent"`: Fetches the latest available database version.
        - A four-digit year (e.g., `"2010"`): Selects a specific version, such as the 2010 or 2019 database for CO.

        If not provided, the default is `"most_recent"`.
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
    :py:func:`~radis.api.hdf5.hdf2df`

    if a registered entry already exists and `radis.config["ALLOW_OVERWRITE"]` is `True`:
    - if any situation arises where the databank needs to be re-downloaded, the possible urls are attempted in their usual order of preference, as if the databank hadn't been registered, rather than directly re-downloading from the same url that was previously registered, in case e.g. a new linelist has been uploaded since the databank was previously registered
    - If no partition function file is registered, e.g because one wasn't available server-side when the databank was last registered, an attempt is still made again to download it, to account for e.g. the case where one has since been uploaded

    See Also
    --------
    :py:func:`~radis.io.hitran.fetch_hitran`, :py:func:`~radis.io.exomol.fetch_exomol`
    :py:func:`~radis.io.geisa.fetch_geisa`, :py:func:`~radis.io.kurucz.fetch_kurucz`
    :py:func:`~radis.io.hdf5.hdf2df`
    :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank`

    """

    recent_database = get_recent_hitemp_database_year(molecule)
    if database == "most_recent":
        database = recent_database
    available_years = {"2010", str(recent_database)}
    if database not in available_years:
        raise KeyError(
            f"The database '{database}' is not recognized as an available HITEMP database for '{molecule}'. Please choose from the following available years: {available_years}."
        )

    if r"{molecule}" in databank_name:
        databank_name = databank_name.format(**{"molecule": molecule})
        databank_name += "-2010" if database == "2010" else ""

    if local_databases is None:

        local_databases = join(config["DEFAULT_DOWNLOAD_PATH"], "hitemp")
    local_databases = abspath(expanduser(local_databases))

    ldb = HITEMPDatabaseManager(
        databank_name,
        molecule=molecule,
        local_databases=local_databases,
        verbose=verbose,
        chunksize=chunksize,
        parallel=parallel,
        engine=engine,
        database=database,
    )
    # Check if the database is registered in radis.json
    local_files, urlnames = [], []
    if ldb.is_registered():
        entries = getDatabankEntries(ldb.name)
        local_files, urlnames = entries["path"], entries["download_url"]

    if ldb.is_registered() and not config["ALLOW_OVERWRITE"]:
        print(ldb.is_registered())
        error = False
        if cache == "regen" or not local_files:
            error = True
        files_to_check = local_files
        if ldb.get_missing_files(files_to_check):
            error = True
        if error:
            raise Exception(
                'Changes are required to the local database, and hence updating the registered entry, but "ALLOW_OVERWRITE" is False. Set `radis.config["ALLOW_OVERWRITE"]=True` to allow the changes to be made and config file to be automatically updated accordingly.'
            )

    # Delete files if needed:
    relevant_files = ldb.keep_only_relevant(
        local_files, load_wavenum_min, load_wavenum_max, verbose=(verbose > 1)
    )
    if cache == "regen":
        ldb.remove_local_files(relevant_files)
    ldb.check_deprecated_files(
        ldb.get_existing_files(relevant_files),
        auto_remove=True if cache != "force" else False,
    )

    # Get list of all expected local files for this database if there is no registry
    downloaded = True
    if len(local_files) < 1:
        local_files, urlnames = ldb.get_filenames()

    # Get missing files
    download_files = ldb.get_missing_files(local_files)
    download_files = ldb.keep_only_relevant(
        download_files, load_wavenum_min, load_wavenum_max, verbose=(verbose > 1)
    )
    # do not re-download files if they exist in another format :
    if engine in ["vaex", "auto", "default"]:
        # ... convert files if asked:

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

    if len(download_files) < 1:
        downloaded = False

    # Download files
    if downloaded:
        if len(urlnames) == 0:
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
            # Do not raise an exception for the 2010 version since it is not available on the HITEMP website.
            if nrows != Ntotal_lines_expected and Ntotal_lines_expected != 0:
                raise AssertionError(
                    f"Number of lines in local database ({nrows:,}) "
                    + "differ from the expected number of lines for "
                    + f"HITEMP {molecule}: {Ntotal_lines_expected}"
                )

    # Register
    if downloaded or not ldb.is_registered():
        ldb.register(downloaded)

    if len(download_files) > 0 and clean_cache_files:
        ldb.clean_download_files()

    # Load and return
    files_loaded = ldb.keep_only_relevant(
        local_files, load_wavenum_min, load_wavenum_max, verbose=verbose
    )

    if isotope and type(isotope) == int:
        isotope = str(isotope)

    df = ldb.load(
        files_loaded,  # filter other files,
        columns=columns,
        within=[("iso", isotope)] if isotope is not None else [],
        # for relevant files, get only the right range :
        lower_bound=[("wav", load_wavenum_min)] if load_wavenum_min is not None else [],
        upper_bound=[("wav", load_wavenum_max)] if load_wavenum_max is not None else [],
        output=output,
    )

    return (df, files_loaded) if return_local_path else df


if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_hitemp.py"]))
