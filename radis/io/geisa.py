# -*- coding: utf-8 -*-
"""

Summary
-----------------------

GEISA database parser

-----------------------

Defines :py:func:`~radis.io.fetch_geisa` based on :py:class:`~radis.api.geisaapi.GEISADatabaseManager`


"""

from os.path import abspath, expanduser, join

from radis.api.geisaapi import GEISADatabaseManager


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
    :py:func:`~radis.api.hdf5.hdf2df`

    See Also
    --------
    :py:func:`~radis.io.hitran.fetch_hitran`, :py:func:`~radis.io.exomol.fetch_exomol`
    :py:func:`~radis.io.hitemp.fetch_hitemp`, :py:func:`~radis.api.hdf5.hdf2df`
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
        isotope=isotope,
        load_wavenum_min=load_wavenum_min,  # for relevant files, get only the right range
        load_wavenum_max=load_wavenum_max,
        output=output,
    )

    return (df, local_files) if return_local_path else df


#%%

if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_geisa.py"]))
