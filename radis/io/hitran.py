# -*- coding: utf-8 -*-
"""
RADIS Hitran Functions; based on Common API hitranapi.py
"""

from os.path import abspath, expanduser, join

from radis.api.hitranapi import HITRANDatabaseManager
from radis.misc.warning import AccuracyWarning


# TODO: implement parallel=True for all isotopes ?
def fetch_hitran(
    molecule,
    extra_params=None,
    local_databases=None,
    databank_name="HITRAN-{molecule}",
    isotope=None,
    load_wavenum_min=None,
    load_wavenum_max=None,
    columns=None,
    cache=True,
    verbose=True,
    clean_cache_files=True,
    return_local_path=False,
    engine="default",
    output="pandas",
    parallel=True,
    parse_quanta=True,
):
    """Download all HITRAN lines from HITRAN website. Unzip and build a HDF5 file directly.

    Returns a Pandas DataFrame containing all lines.

    Parameters
    ----------
    molecule: str
        one specific molecule name, listed in HITRAN molecule metadata.
        See https://hitran.org/docs/molec-meta/
        Example: "H2O", "CO2", etc.
    local_databases: str
        where to create the RADIS HDF5 files. Default ``"~/.radisdb/hitran"``.
        Can be changed in ``radis.config["DEFAULT_DOWNLOAD_PATH"]`` or in ~/radis.json config file
    databank_name: str
        name of the databank in RADIS :ref:`Configuration file <label_lbl_config_file>`
        Default ``"HITRAN-{molecule}"``
    isotope: str
        load only certain isotopes : ``'2'``, ``'1,2'``, etc. If ``None``, loads
        everything. Default ``None``.
    load_wavenum_min, load_wavenum_max: float (cm-1)
        load only specific wavenumbers.
    columns: list of str
        list of columns to load. If ``None``, returns all columns in the file.
    extra_params: 'all' or None
        Downloads all additional columns available in the HAPI database for the molecule including
        parameters like `gamma_co2`, `n_co2` that are required to calculate spectrum in co2 diluent.
        For eg:
        ::

            from radis.io.hitran import fetch_hitran
            df = fetch_hitran('CO', extra_params='all', cache='regen') # cache='regen' to regenerate new database with additional columns

    Other Parameters
    ----------------
    cache: ``True``, ``False``, ``'regen'`` or ``'force'``
        if ``True``, use existing HDF5 file. If ``False`` or ``'regen'``, rebuild it.
        If ``'force'``, raise an error if cache file cannot be used (useful for debugging).
        Default ``True``.
    verbose: bool
    clean_cache_files: bool
        if ``True`` clean downloaded cache files after HDF5 are created.
    return_local_path: bool
        if ``True``, also returns the path of the local database file.
    engine: 'pytables', 'vaex', 'default'
        which HDF5 library to use. If 'default' use the value from ~/radis.json
    output: 'pandas', 'vaex', 'jax'
        format of the output DataFrame. If ``'jax'``, returns a dictionary of
        jax arrays. If ``'vaex'``, output is a :py:class:`vaex.dataframe.DataFrameLocal`

        .. note::
            Vaex DataFrames are memory-mapped. They do not take any space in RAM
            and are extremelly useful to deal with the largest databases.

    parallel: bool
        if ``True``, uses joblib.parallel to load database with multiple processes
    parse_quanta: bool
        if ``True``, parse local & global quanta (required to identify lines
        for non-LTE calculations ; but sometimes lines are not labelled.)


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

        from radis.io.hitran import fetch_hitran
        df = fetch_hitran("CO")
        print(df.columns)
        >>> Index(['id', 'iso', 'wav', 'int', 'A', 'airbrd', 'selbrd', 'El', 'Tdpair',
            'Pshft', 'gp', 'gpp', 'branch', 'jl', 'vu', 'vl'],
            dtype='object')

    .. minigallery:: radis.fetch_hitran

    Notes
    -----
    if using ``load_only_wavenum_above/below`` or ``isotope``, the whole
    database is anyway downloaded and uncompressed to ``local_databases``
    fast access .HDF5 files (which will take a long time on first call). Only
    the expected wavenumber range & isotopes are returned. The .HFD5 parsing uses
    :py:func:`~radis.io.hdf5.hdf2df`

    See Also
    --------
    :py:func:`~radis.io.hitemp.fetch_hitemp`, :py:func:`~radis.io.exomol.fetch_exomol`,
    :py:func:`~radis.io.geisa.fetch_geisa`,
    :py:func:`~radis.api.hdf5.hdf2df`, :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank`

    """

    if r"{molecule}" in databank_name:
        databank_name = databank_name.format(**{"molecule": molecule})

    if local_databases is None:
        import radis

        local_databases = join(radis.config["DEFAULT_DOWNLOAD_PATH"], "hitran")
    local_databases = abspath(local_databases.replace("~", expanduser("~")))

    ldb = HITRANDatabaseManager(
        databank_name,
        molecule=molecule,
        local_databases=local_databases,
        engine=engine,
        extra_params=extra_params,
        verbose=verbose,
        parallel=parallel,
    )

    # Get expected local files for this database:
    local_file = ldb.get_filenames()

    # Delete files if needed:
    if cache == "regen":
        ldb.remove_local_files(local_file)
    else:
        # Raising AccuracyWarning if local_file exists and doesn't have extra columns in it
        if ldb.get_existing_files(local_file) and extra_params == "all":
            columns = ldb.get_columns(local_file[0])
            extra_columns = ["y_", "gamma_", "n_"]
            found = False
            for key in extra_columns:
                for column_name in columns:
                    if key in column_name:
                        found = True
                        break

            if not found:
                import warnings

                warnings.warn(
                    AccuracyWarning(
                        "All columns are not downloaded currently, please use cache = 'regen' and extra_params='all' to download all columns."
                    )
                )

    ldb.check_deprecated_files(
        ldb.get_existing_files(local_file),
        auto_remove=True if cache != "force" else False,
    )

    # Download files
    download_files = ldb.get_missing_files(local_file)
    if download_files:
        ldb.download_and_parse(download_files, cache=cache, parse_quanta=parse_quanta)

    # Register
    if not ldb.is_registered():
        ldb.register()

    if len(download_files) > 0 and clean_cache_files:
        ldb.clean_download_files()

    # Load and return
    df = ldb.load(
        local_file,
        columns=columns,
        within=[("iso", isotope)] if isotope is not None else [],
        # for relevant files, get only the right range :
        lower_bound=[("wav", load_wavenum_min)] if load_wavenum_min is not None else [],
        upper_bound=[("wav", load_wavenum_max)] if load_wavenum_max is not None else [],
        output=output,
    )

    return (df, local_file) if return_local_path else df


# ======================================================
# %% Test


if __name__ == "__main__":

    from radis.test.io.test_hitran_cdsd import _run_testcases

    print("Testing HITRAN parsing: ", _run_testcases())
    from radis.test.io.test_query import _run_testcases

    print("Testing HITRAN fetch: ", _run_testcases())
