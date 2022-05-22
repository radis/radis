# -*- coding: utf-8 -*-
"""
Created on Sun May 22 17:35:05 2022

@author: erwan
"""

from os.path import abspath, exists, expanduser, join

from radis.api.hdf5 import update_pytables_to_vaex
from radis.api.hitempapi import HITEMPDatabaseManager


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

    See Also
    --------
    :py:func:`~radis.io.hitran.fetch_hitran`, :py:func:`~radis.io.exomol.fetch_exomol`
    :py:func:`~radis.io.geisa.fetch_geisa`, :py:func:`~radis.io.hdf5.hdf2df`
    :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank`

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
