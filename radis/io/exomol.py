# -*- coding: utf-8 -*-
"""
Created on Sun May 22 18:11:23 2022

@author: erwan
"""

import pathlib

import numpy as np

from radis.api.exomolapi import (
    MdbExomol,
    get_exomol_database_list,
    get_exomol_full_isotope_name,
)
from radis.db.classes import get_molecule_identifier


def fetch_exomol(
    molecule,
    database=None,
    local_databases=None,
    databank_name="EXOMOL-{molecule}",
    isotope="1",
    load_wavenum_min=None,
    load_wavenum_max=None,
    columns=None,
    cache=True,
    verbose=True,
    clean_cache_files=True,
    return_local_path=False,
    return_partition_function=False,
    engine="default",
    output="pandas",
    skip_optional_data=True,
):
    """Stream ExoMol file from EXOMOL website. Unzip and build a HDF5 file directly.

    Returns a Pandas DataFrame containing all lines.

    Parameters
    ----------
    molecule: ``str``
        ExoMol molecule
    database: ``str``
        database name. Ex:: ``POKAZATEL`` or ``BT2`` for ``H2O``. See
        :py:data:`~radis.api.exomolapi.KNOWN_EXOMOL_DATABASE_NAMES`. If ``None`` and
        there is only one database available, use it.
    local_databases: ``str``
        where to create the RADIS HDF5 files. Default ``"~/.radisdb/exomol"``.
        Can be changed in ``radis.config["DEFAULT_DOWNLOAD_PATH"]`` or in ~/radis.json config file
    databank_name: ``str``
        name of the databank in RADIS :ref:`Configuration file <label_lbl_config_file>`
        Default ``"EXOMOL-{molecule}"``
    isotope: ``str`` or ``int``
        load only certain isotopes, sorted by terrestrial abundances : ``'1'``, ``'2'``,
        etc. Default ``1``.

        .. note::

            In RADIS, isotope abundance is included in the line intensity
            calculation. However, the terrestrial abundances used may not be
            relevant to non-terrestrial applications.
            By default, the abundance is given reading HITRAN data. If the molecule
            does not exist in the HITRAN database, the abundance is read from
            the ``radis/radis_default.json`` configuration file, which can be
            modified by editing ``radis.config`` after import or directly
            by editing the user ``~/radis.json`` user configuration file
            (overwrites ``radis_default.json``). In the ``radis/radis_default.json``
            file, values were calculated with a simple model based on the
            terrestrial isotopic abundance of each element.

    load_wavenum_min, load_wavenum_max: float (cm-1)
        load only specific wavenumbers.
    columns: list of str
        list of columns to load. If ``None``, returns all columns in the file.

    Other Parameters
    ----------------
    cache: bool, or ``'regen'`` or ``'force'``
        if ``True``, use existing HDF5 file. If ``False`` or ``'regen'``, rebuild it.
        If ``'force'``, crash if not cache file found. Default ``True``.
    verbose: bool
    clean_cache_files: bool
        if ``True`` clean downloaded cache files after HDF5 are created.
    return_local_path: bool
        if ``True``, also returns the path of the local database file.
    return_partition_function: bool
        if ``True``, also returns a :py:class:`~radis.levels.partfunc.PartFuncExoMol` object.
    engine: 'vaex', 'feather'
        which memory-mapping library to use. If 'default' use the value from ~/radis.json
    output: 'pandas', 'vaex', 'jax'
        format of the output DataFrame. If ``'jax'``, returns a dictionary of
        jax arrays. If ``'vaex'``, output is a :py:class:`vaex.dataframe.DataFrameLocal`

        .. note::
            Vaex DataFrames are memory-mapped. They do not take any space in RAM
            and are extremelly useful to deal with the largest databases.

    skip_optional_data : bool
        If False, fetch all fields which are marked as available in the ExoMol definition
        file. If True, load only the first 4 columns of the states file
        ("i", "E", "g", "J"). The structure of the columns above 5 depend on the
        the definitions file (*.def) and the Exomol version.
        If ``skip_optional_data=False``, two errors may occur:

            - a field is marked as present/absent in the *.def field but is
              absent/present in the *.states file (ie both files are inconsistent).
            - in the updated version of Exomol, new fields have been added in the
              states file of some species. But it has not been done for all species,
              so both structures exist. For instance, the states file of
              https://exomol.com/data/molecules/HCl/1H-35Cl/HITRAN-HCl/ follows the
              structure described in [1]_, unlike the states file of
              https://exomol.com/data/molecules/NO/14N-16O/XABC/ which follows the
              structure described in [2]_.

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

    .. minigallery:: radis.fetch_exomol

    Notes
    -----
    if using ``load_only_wavenum_above/below`` or ``isotope``, the whole
    database is anyway downloaded and uncompressed to ``local_databases``
    fast access .HDF5 files (which will take a long time on first call). Only
    the expected wavenumber range & isotopes are returned. The .HFD5 parsing uses
    :py:func:`~radis.api.hdf5.hdf2df`

    References
    ----------

    .. [1] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Barton, E. J., Chubb, K. L., Coles, P. A., … Zak, E. (2016). The ExoMol database: molecular line lists for exoplanet and other hot atmospheres. https://doi.org/10.1016/j.jms.2016.05.002
    .. [2] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Clark, V. H. J., Chubb, K. L., Conway, E. K., … Yurchenko, O. P. (2020). The 2020 release of the ExoMol database: Molecular line lists for exoplanet and other hot atmospheres. Journal of Quantitative Spectroscopy and Radiative Transfer, 255, 107228. https://doi.org/10.1016/j.jqsrt.2020.107228

    See Also
    --------
    :py:func:`~radis.io.hitran.fetch_hitran`, :py:func:`~radis.io.hitemp.fetch_hitemp`
    :py:func:`~radis.api.hdf5.hdf2df`

    """
    # TODO: implement columns= ... to load only specific columns.
    # refactor with "self._quantumNumbers" (which serves the same purpose)

    # Ensure isotope format:
    try:
        isotope = int(isotope)
    except:
        raise ValueError(
            f"In fetch_exomol, ``isotope`` must be an integer. Got `{isotope}` "
            + "Only one isotope can be queried at a time. "
        )

    full_molecule_name = get_exomol_full_isotope_name(molecule, isotope)
    known_exomol_databases, recommended_database = get_exomol_database_list(
        molecule, full_molecule_name
    )

    _exomol_use_hint = "Select one of them with `radis.fetch_exomol(DATABASE_NAME)`, `SpectrumFactory.fetch_databank('exomol', exomol_database=DATABASE_NAME')`, or `calc_spectrum(..., databank=('exomol', DATABASE_NAME))` \n"
    if database is None or database == "default":
        if len(known_exomol_databases) == 1:
            database = known_exomol_databases[0]
            if verbose:
                print(f"Using ExoMol database {database} for {full_molecule_name}")
        elif recommended_database:
            database = recommended_database
            if verbose:
                print(
                    f"Using ExoMol database {database} for {full_molecule_name} (recommended by the ExoMol team). All available databases are {known_exomol_databases}. {_exomol_use_hint}"
                )
        else:
            raise KeyError(
                f"Choose one of the several databases available for {full_molecule_name} in ExoMol: {known_exomol_databases}. ({recommended_database} is recommended by the ExoMol team). {_exomol_use_hint}"
            )
    else:
        if database not in known_exomol_databases:
            raise KeyError(
                f"{database} is not of the known available ExoMol databases for {full_molecule_name}. Choose one of : {known_exomol_databases}. ({recommended_database} is recommended by the ExoMol team). {_exomol_use_hint}"
            )

    if local_databases is None:
        import radis

        local_databases = pathlib.Path(radis.config["DEFAULT_DOWNLOAD_PATH"]) / "exomol"

    local_path = (
        pathlib.Path(local_databases).expanduser()
        / molecule
        / full_molecule_name
        / database
    )

    # TODO: add deprecation if missing columns in cache file

    # Init database, download files if needed.
    mdb = MdbExomol(
        local_path,
        molecule=molecule,
        name=databank_name,
        local_databases=local_databases,
        nurange=[
            load_wavenum_min if load_wavenum_min is not None else 0.0,
            load_wavenum_max if load_wavenum_max is not None else np.inf,
        ],
        engine=engine,
        cache=cache,
        skip_optional_data=skip_optional_data,
    )

    # Get local files
    local_files = mdb.trans_file
    if not isinstance(local_files, list):
        local_files = [local_files]
    mgr = mdb.get_datafile_manager()
    local_files = [mgr.cache_file(f) for f in local_files]

    # Specific for RADIS : rename columns
    radis2exomol_columns = {
        "wav": "nu_lines",
        "airbrd": "alpha_ref",
        "El": "elower",
        "ju": "jupper",
        "jl": "jlower",
        "gp": "gupper",
        "gpp": "glower",
        "Tdpair": "n_Texp",
    }
    # get column name converting to exomol/exojax format if possible, else use the same
    if columns is not None:
        columns_exomol = [radis2exomol_columns.get(c, c) for c in columns] + [
            "Sij0",
            "jlower",
            "jupper",
        ]  # needed for broadening
    else:
        columns_exomol = None

    df = mdb.load(
        local_files,
        columns=columns_exomol,
        lower_bound=([("nu_lines", load_wavenum_min)] if load_wavenum_min else [])
        + ([("Sij0", mdb.crit)] if not np.isneginf(mdb.crit) else []),
        upper_bound=([("nu_lines", load_wavenum_max)] if load_wavenum_max else []),
        output=output,
    )

    if "jlower" not in df:
        raise KeyError(
            f"jlower not found. Maybe try to delete cache file {local_files} and restart?"
        )

    # Add broadening
    mdb.set_broadening(df, output=output)

    # Specific for RADIS :
    # ... Get RADIS column names:
    exomol2radis_columns = {v: k for k, v in radis2exomol_columns.items()}
    mdb.rename_columns(df, {k: v for k, v in exomol2radis_columns.items() if k in df})

    assert "wav" in df

    # ... include isotopic abundance in linestrength :
    # Note : ExoMol treats isotopes as independant molecules ; linestrength is not
    # corrected by isotopic abundance.
    # Below, replace Linestrength with Line Intensity taking into account
    # Terrestrial isotopic abundance (to be compatible with HITRAN/HITEMP/etc. )
    from radis.db.molparam import MOLPARAMS_EXTRA_PATH, MolParams

    Ia = MolParams(extra_file_json=MOLPARAMS_EXTRA_PATH).get(
        molecule, isotope, "abundance"
    )

    if output == "jax":
        try:
            import jax.numpy as jnp
        except:
            import numpy as jnp
        df["logsij0"] += jnp.log(Ia)
    else:
        df["Sij0"] *= Ia
        mdb.rename_columns(df, {"Sij0": "int"})

    # Add Attributes of the DataFrame
    if output == "pandas":  # no attribtes in "Jax" or "Vaex" mode
        from radis.db.classes import HITRAN_MOLECULES

        attrs = {}
        if molecule in HITRAN_MOLECULES:
            attrs["id"] = get_molecule_identifier(
                molecule
            )  # HITRAN id-number (if available)
        attrs["molecule"] = molecule
        attrs["iso"] = isotope

        for k, v in attrs.items():
            df.attrs[k] = v

    # Return:
    out = df
    if return_local_path or return_partition_function:
        out = [out]
    if return_local_path:
        out.append(str(mdb.path))
    if return_partition_function:
        assert return_local_path
        out.append(mdb.to_partition_function_tabulator())

    return out
