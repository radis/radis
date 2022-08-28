"""Molecular database (MDB) class

   * MdbExomol is the MDB for ExoMol

Initial code borrowed from the `Exojax <https://github.com/HajimeKawahara/exojax>`__
code (which you should also have a look at !), by @HajimeKawahara, under MIT License.

"""
import os
import pathlib
import warnings

import numpy as np

try:
    from . import exomolapi
    from .dbmanager import DatabaseManager
    from .exomol_utils import e2s
except ImportError:  # if local import
    from radis.io import exomolapi
    from radis.io.dbmanager import DatabaseManager
    from radis.io.exomol_utils import e2s

from radis.db.classes import (
    EXOMOL_MOLECULES,
    EXOMOL_ONLY_ISOTOPES_NAMES,
    get_molecule_identifier,
)

EXOMOL_URL = "http://www.exomol.com/db/"

from urllib.request import HTTPError, urlopen

from bs4 import BeautifulSoup


def get_exomol_full_isotope_name(molecule, isotope):
    """Get full isotope name for ``molecule`` and ``isotope`` number

    Parameters
    ----------
    molecule: str
    isotope: int
        terrestrial abundance

    Examples
    --------
    ::

        get_exomol_full_isotope_name("CH4", 1)
        >>> '12C-1H4'

    See Also
    --------
    :py:func:`~radis.io.exomol.get_exomol_database_list`"""

    if molecule not in EXOMOL_MOLECULES:
        raise ValueError(
            f"Molecule {molecule} not available in ExoMol. Change database, or choose one of ExoMol available molecules: {sorted(EXOMOL_MOLECULES)}"
        )

    if (molecule, isotope) in EXOMOL_ONLY_ISOTOPES_NAMES:
        return EXOMOL_ONLY_ISOTOPES_NAMES[(molecule, isotope)]
    else:
        # Read and convert from HITRAN molecules
        from radis.db.molparam import MolParams

        mp = MolParams()
        return mp.get(molecule, isotope, "isotope_name_exomol")


def get_list_of_known_isotopes(molecule):
    """find all isotopes until error ensues"""

    i = 1
    isotope_list = []
    while True:
        try:
            iso_name = get_exomol_full_isotope_name(molecule, i)
        except:
            break
        else:
            isotope_list.append(iso_name)
        finally:
            i += 1
    return isotope_list


def get_exomol_database_list(molecule, isotope_full_name=None):
    """Parse ExoMol website and return list of available databases, and recommended database

    Parameters
    ----------
    molecule: str
    isotope_full_name: str
        isotope full name (ex. ``12C-1H4`` for CH4,1). Get it from
        :py:func:`radis.io.exomol.get_exomol_full_isotope_name`

    Returns
    -------
    list of databases, database recommended by ExoMol

    Examples
    --------
    Get CH4 from ExoMol :
    ::

        databases, recommended = get_exomol_database_list("CH4", "12C-1H4")
        >>> ['xsec-YT10to10', 'YT10to10', 'YT34to10'], 'YT34to10'

    Or combine with :py:func:`~radis.io.exomol.get_exomol_full_isotope_name` to
    get the isopologue (sorted by terrestrial abundance) ::

        from radis.io.exomol import get_exomol_database_list, get_exomol_full_isotope_name
        databases, recommended = get_exomol_database_list("CH4", get_exomol_full_isotope_name("CH4", 1))
        >>> ['xsec-YT10to10', 'YT10to10', 'YT34to10'], 'YT34to10'


    .. minigallery:: radis.io.exomol.get_exomol_database_list


    See Also
    --------
    :py:func:`~radis.io.exomol.get_exomol_full_isotope_name`
    """

    if isotope_full_name is None:
        raise ValueError(
            f"Give isotope name. List of known isotopes for {molecule} : {get_list_of_known_isotopes(molecule)}"
        )

    url = f"https://exomol.com/data/molecules/{molecule}/{isotope_full_name}"
    try:
        response = urlopen(url).read()
    except HTTPError as err:
        if isotope_full_name not in get_list_of_known_isotopes(molecule):
            extra = f". Isotope name {isotope_full_name} is not in list of known isotopes : {get_list_of_known_isotopes(molecule)}"
        else:
            extra = ""
        raise ValueError(f"HTTPError opening url={url}" + extra) from err

    soup = BeautifulSoup(
        response, features="lxml"
    )  # make soup that is parse-able by bs

    # Recommended database
    rows = soup.find_all(
        "a", {"class": "list-group-item link-list-group-item recommended"}
    )
    databases_recommended = [r.get_attribute_list("title")[0] for r in rows]

    # All others
    rows = soup.find_all("a", {"class": "list-group-item link-list-group-item"})
    databases = [r.get_attribute_list("title")[0] for r in rows]

    if len(databases_recommended) > 1:
        # Known exceptions :
        if (
            isotope_full_name == "28Si-16O"
            and databases_recommended[0] == "xsec-SiOUVenIR"
        ):
            # this is a cross-section dataset, shouldn't be used. Reverse and use the other one:
            databases_recommended = databases_recommended[::-1]
        else:
            print(
                f"Multiple recommended databases found for {molecule} in ExoMol : {databases_recommended}. This is unexpected. Using the first"
            )

    databases = databases + databases_recommended

    if len(databases_recommended) > 0:
        recommended_database = databases_recommended[0]
    else:
        recommended_database = False

    return databases, recommended_database


# def fetch_exomol_molecule_list():
#     """Parse ExoMol website and return list of available databases, and recommended database

#     Parameters
#     ----------
#     molecule: str
#     isotope_full_name: str
#         isotope full name (ex. ``12C-1H4`` for CH4,1). Get it from
#         :py:func:`radis.io.exomol.get_exomol_full_isotope_name`

#     Returns
#     -------

#     Examples
#     --------
#     Get CH4 from ExoMol :
#     ::
#         databases, recommended = get_exomol_database_list("CH4", "12C-1H4")
#         >>> ['xsec-YT10to10', 'YT10to10', 'YT34to10'], 'YT34to10'

#     Or combine with :py:func:`~radis.io.exomol.get_exomol_full_isotope_name` to
#     get the isopologue (sorted by terrestrial abundance) ::

#         from radis.io.exomol import get_exomol_database_list, get_exomol_full_isotope_name
#         databases, recommended = get_exomol_database_list("CH4", get_exomol_full_isotope_name("CH4", 1))
#         >>> ['xsec-YT10to10', 'YT10to10', 'YT34to10'], 'YT34to10'


#     .. minigallery:: radis.io.exomol.get_exomol_database_list


#     See Also
#     --------
#     :py:func:`~radis.io.exomol.get_exomol_full_isotope_name`
#     """

# url = f"https://exomol.com/data/molecules/"
# try:
#     response = urlopen(url).read()
# except HTTPError as err:
#     raise ValueError(f"HTTPError opening url={url}") from err

# soup = BeautifulSoup(
#     response, features="lxml"
# )  # make soup that is parse-able by bs

# # Recommended database
# rows = soup.find_all(
#     "a", {"class": "list-group-item link-list-group-item recommended"}
# )
# databases_recommended = [r.get_attribute_list("title")[0] for r in rows]

# # All others
# rows = soup.find_all("a", {"class": "list-group-item link-list-group-item"})
# databases = [r.get_attribute_list("title")[0] for r in rows]

# if len(databases_recommended) > 1:
#     # Known exceptions :
#     if (
#         isotope_full_name == "28Si-16O"
#         and databases_recommended[0] == "xsec-SiOUVenIR"
#     ):
#         # this is a cross-section dataset, shouldn't be used. Reverse and use the other one:
#         databases_recommended = databases_recommended[::-1]
#     else:
#         print(
#             f"Multiple recommended databases found for {molecule} in ExoMol : {databases_recommended}. This is unexpected. Using the first"
#         )

# databases = databases + databases_recommended

# return databases, databases_recommended[0]


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
        :py:data:`~radis.io.exomol.KNOWN_EXOMOL_DATABASE_NAMES`. If ``None`` and
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
    :py:func:`~radis.io.hdf5.hdf2df`

    References
    ----------

    .. [1] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Barton, E. J., Chubb, K. L., Coles, P. A., … Zak, E. (2016). The ExoMol database: molecular line lists for exoplanet and other hot atmospheres. https://doi.org/10.1016/j.jms.2016.05.002
    .. [2] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Clark, V. H. J., Chubb, K. L., Conway, E. K., … Yurchenko, O. P. (2020). The 2020 release of the ExoMol database: Molecular line lists for exoplanet and other hot atmospheres. Journal of Quantitative Spectroscopy and Radiative Transfer, 255, 107228. https://doi.org/10.1016/j.jqsrt.2020.107228

    See Also
    --------
    :py:func:`~radis.io.hitran.fetch_hitran`, :py:func:`~radis.io.hitemp.fetch_hitemp`
    :py:func:`~radis.io.hdf5.hdf2df`

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


class MdbExomol(DatabaseManager):

    """molecular database of ExoMol

    MdbExomol is a class for ExoMol.

    Parameters
    ----------
    path: str
        path for Exomol data directory/tag. For instance, "/home/CO/12C-16O/Li2015"
    nurange: array
        wavenumber range list (cm-1) or wavenumber array
    margin: float
        margin for nurange (cm-1)
    crit: float
        line strength lower limit for extraction
    bkgdatm: str
        background atmosphere for broadening. e.g. H2, He,
    broadf: bool
        if False, the default broadening parameters in .def file is used

    Other Parameters
    ----------------
    engine : str
        which memory mapping engine to use : 'vaex', 'pytables' (HDF5), 'feather'
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

    Notes
    -----

    The trans/states files can be very large. For the first time to read it,
    we convert it to the feather or hdf5-format. After the second-time,
    we use the feather/hdf5 format instead.

    Examples
    --------
    ::

        # Init database, download files if needed.
        mdb = MdbExomol(
            local_path,
            molecule=molecule,
            name=databank_name,
            local_databases=local_databases,
            # nurange=[load_wavenum_min, load_wavenum_max],
            engine="vaex",
        )

        # Get cache files to load :
        mgr = mdb.get_datafile_manager()
        local_files = [mgr.cache_file(f) for f in mdb.trans_file]

        # Load files
        df = mdb.load(
            local_files,
            columns=columns_exomol,
            lower_bound=([('nu_lines', load_wavenum_min)] if load_wavenum_min else []) + ([("Sij0", mdb.crit)] if not np.isneginf(mdb.crit) else []),
            upper_bound=([('nu_lines', load_wavenum_max)] if load_wavenum_max else []),
            output="jax", # or "pytables", "vaex"
        )

    .. minigallery:: radis.fetch_exomol


    DataFrame columns
    -----------------

    nu_lines (nd array): line center (cm-1)
    Sij0 (nd array): line strength at T=Tref (cm)
    dev_nu_lines (np array): line center in device (cm-1)
    logsij0 (np array): log line strength at T=Tref
    A (np array): Einstein A coefficient
    elower (np array): the lower state energy (cm-1)
    gpp (np array): statistical weight
    jlower (np array): J_lower
    jupper (np array): J_upper
    n_Tref (np array): temperature exponent
    alpha_ref (np array): alpha_ref (gamma0)
    n_Tref_def: default temperature exponent in .def file, used for jlower not given in .broad
    alpha_ref_def: default alpha_ref (gamma0) in .def file, used for jlower not given in .broad

    References
    ----------

    .. [1] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Barton, E. J., Chubb, K. L., Coles, P. A., … Zak, E. (2016). The ExoMol database: molecular line lists for exoplanet and other hot atmospheres. https://doi.org/10.1016/j.jms.2016.05.002
    .. [2] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Clark, V. H. J., Chubb, K. L., Conway, E. K., … Yurchenko, O. P. (2020). The 2020 release of the ExoMol database: Molecular line lists for exoplanet and other hot atmospheres. Journal of Quantitative Spectroscopy and Radiative Transfer, 255, 107228. https://doi.org/10.1016/j.jqsrt.2020.107228

    """

    # TODO : inherit from DatabaseManager or similar

    # @dev: In exojax this class is defined in exojax/spec/moldb.py
    # see https://github.com/HajimeKawahara/exojax/blob/develop/src/exojax/spec/moldb.py
    # It uses Jax arrays (jnp). Here RADIS uses Numpy arrays.

    def __init__(
        self,
        path,
        molecule,
        database=None,
        local_databases=None,
        name="EXOMOL-{molecule}",
        nurange=[0.0, np.inf],
        margin=1.0,
        crit=-np.inf,
        bkgdatm="Air",  # TODO: use Air whenever possible, to be consistent with HITRAN/HITEMP
        broadf=True,
        engine="vaex",
        verbose=True,
        cache=True,
        skip_optional_data=True,
    ):
        super().__init__(
            name,
            molecule,
            local_databases,
            engine,
            verbose=verbose,
        )

        assert cache  # cache only used for cache='regen' or cache='force' modes, cache=False is not expected

        if engine == "default":
            from radis import config

            engine = config["MEMORY_MAPPING_ENGINE"]
            if engine == "auto":
                engine = "vaex"

        self.path = pathlib.Path(path)

        if local_databases is not None:
            self.path = pathlib.Path(local_databases).expanduser() / self.path

        t0 = self.path.parents[0].stem
        molec = t0 + "__" + str(self.path.stem)
        self.bkgdatm = bkgdatm
        print("Background atmosphere: ", self.bkgdatm)
        molecbroad = t0 + "__" + self.bkgdatm

        self.crit = crit
        self.margin = margin
        self.nurange = [np.min(nurange), np.max(nurange)]
        self.wmin, self.wmax = self.nurange
        self.broadf = broadf
        # Where exomol files are
        self.states_file = self.path / pathlib.Path(molec + ".states.bz2")
        self.pf_file = self.path / pathlib.Path(molec + ".pf")
        self.def_file = self.path / pathlib.Path(molec + ".def")
        self.broad_file = self.path / pathlib.Path(molecbroad + ".broad")

        if not self.def_file.exists():
            self.download(molec, extension=[".def"])
        if not self.pf_file.exists():
            self.download(molec, extension=[".pf"])
        if not self.states_file.exists():
            self.download(molec, extension=[".states.bz2"])
        if not self.broad_file.exists():
            self.download(molec, extension=[".broad"])

        # Add molecule name
        tag = molec.split("__")
        self.isotope_fullname = tag[0]
        self.molecule = e2s(tag[0])
        # self.isotope = 1  # Placeholder. TODO : impement parsing of other isotopes.

        # load def
        dic_def = exomolapi.read_def(self.def_file)  # approx. 3 ms
        self.n_Texp_def = dic_def["n_Texp"]
        self.alpha_ref_def = dic_def["alpha_ref"]

        #  default n_Texp value if not given
        if self.n_Texp_def is None:
            self.n_Texp_def = 0.5
            warnings.warn(
                Warning(
                    f"""
                    No default broadening exponent in def file. Assigned n = {self.n_Texp_def}
                    """
                )
            )
        #  default alpha_ref value if not given
        if self.alpha_ref_def is None:
            self.alpha_ref_def = 0.07
            warnings.warn(
                Warning(
                    f"""
                    No default broadening in def file. Assigned alpha_ref = {self.alpha_ref_def}
                    """
                )
            )

        # load states
        mgr = self.get_datafile_manager()
        if cache == "regen" and mgr.cache_file(self.states_file).exists():
            if verbose:
                print("Removing existing file ", mgr.cache_file(self.states_file))
            os.remove(mgr.cache_file(self.states_file))
        if mgr.cache_file(self.states_file).exists():
            states = mgr.read(mgr.cache_file(self.states_file))
        else:
            if cache == "force":
                raise ValueError(
                    f"Cache file {str(mgr.cache_file(self.states_file))} does not exist"
                )
            print(
                f"Note: Caching states data to the {engine} format. After the second time, it will become much faster."
            )
            states = exomolapi.read_states(
                self.states_file,
                dic_def,
                engine="vaex" if engine == "vaex" else "csv",
                skip_optional_data=skip_optional_data,
            )
            mgr.write(mgr.cache_file(self.states_file), states)

        # load pf
        pf = exomolapi.read_pf(self.pf_file)
        self.gQT = pf["QT"].to_numpy()  # grid QT
        self.T_gQT = pf["T"].to_numpy()  # T forgrid QT

        # trans file(s)
        print("Reading transition file")

        # Compute linestrengths or retrieve them from cache
        self.Tref = 296.0
        self.QTref = np.array(self.QT_interp(self.Tref))

        # Download files

        # Generate list of files
        # ---------------------
        # Case1 : transitions are stored in a single file
        if dic_def["numinf"] is None:
            self.trans_file = [self.path / pathlib.Path(molec + ".trans.bz2")]
            self.num_tag = [None]

        # Case2 : Transitions are stored in multiple files:
        else:  # dic_def["numinf"] is not None
            # dic_def["numinf"] contains the limit points ``[w(0), w(1), ..., w(n)]``
            # (n+1 elements) defining the spectral ranges appearing in the list
            # dic_def["numtag"] which looks like ``["w(0)-w(1)", "w(1)-w(2)", ..., w(n-1)-w(n)]``
            # (n elements)
            # imin :
            # index i of the "w(i)-w(i+1)" element in dic_def["numtag"] such
            # that nurange[0]<=w(i)
            imin = np.searchsorted(dic_def["numinf"], nurange[0], side="right") - 1
            # imax :
            # index i of the "w(i)-w(i+1)" element in dic_def["numtag"] such
            # that w(i+1)<=nurange[1]
            imax = np.searchsorted(dic_def["numinf"], nurange[1], side="right") - 2
            self.trans_file = []
            self.num_tag = []
            for k, i in enumerate(range(imin, imax + 1)):
                trans_file = self.path / pathlib.Path(
                    molec + "__" + dic_def["numtag"][i] + ".trans.bz2"
                )
                self.trans_file.append(trans_file)
                self.num_tag.append(dic_def["numtag"][i])

        # Look-up missing parameters and write file
        # -----------------------------------------
        for trans_file, num_tag in zip(self.trans_file, self.num_tag):

            if cache == "regen" and mgr.cache_file(trans_file).exists():
                if verbose:
                    print("Removing existing file ", mgr.cache_file(trans_file))
                os.remove(mgr.cache_file(trans_file))

            if not mgr.cache_file(trans_file).exists():

                if cache == "force":
                    raise ValueError(
                        f"Cache file {str(mgr.cache_file(trans_file))} does not exist"
                    )

                if not trans_file.exists():
                    self.download(molec, extension=[".trans.bz2"], numtag=num_tag)
                    # TODO: add option to delete file at the end

                print(
                    f"Note: Caching line transition data to the {engine} format. After the second time, it will become much faster."
                )
                trans = exomolapi.read_trans(
                    trans_file, engine="vaex" if engine == "vaex" else "csv"
                )

                # Complete transition data with lookup on upper & lower state :
                # In particular, compute gup and elower

                exomolapi.pickup_gE(
                    states,
                    trans,
                    trans_file,
                    dic_def,
                    skip_optional_data=skip_optional_data,
                )

                ##Recompute Line strength:
                from radis.lbl.base import (  # TODO: move elsewhere
                    linestrength_from_Einstein,
                )

                self.Sij0 = linestrength_from_Einstein(
                    A=trans["A"],
                    gu=trans["gup"],
                    El=trans["elower"],
                    Ia=1,  #  Sij0 is a linestrength calculated without taking into account isotopic abundance (unlike line intensity parameter of HITRAN. In RADIS this is corrected for in fetch_exomol()  )
                    nu=trans["nu_lines"],
                    Q=self.QTref,
                    T=self.Tref,
                )

                trans["Sij0"] = self.Sij0

                mgr.write(mgr.cache_file(trans_file), trans)

        # Database ready to be loaded.
        # Proceed with mdb.load()

    def set_broadening(self, df, alpha_ref_def=None, n_Texp_def=None, output=None):
        """setting broadening parameters

        Parameters
        ----------
        alpha_ref: set default alpha_ref and apply it. None=use self.alpha_ref_def
        n_Texp_def: set default n_Texp and apply it. None=use self.n_Texp_def

        Returns
        -------
        None. Updates ``df`` with ``alpha_ref``, ``n_Texp_def``
        """
        if alpha_ref_def:
            self.alpha_ref_def = alpha_ref_def
        if n_Texp_def:
            self.n_Texp_def = n_Texp_def

        if self.broadf:
            try:
                print(".broad is used.")
                bdat = exomolapi.read_broad(self.broad_file)
                codelv = exomolapi.check_bdat(bdat)
                print("Broadening code level=", codelv)
                if codelv == "a0":
                    j2alpha_ref, j2n_Texp = exomolapi.make_j2b(
                        bdat,
                        alpha_ref_default=self.alpha_ref_def,
                        n_Texp_default=self.n_Texp_def,
                        jlower_max=np.max(df["jlower"]),
                    )
                    self.alpha_ref = np.array(j2alpha_ref[df["jlower"]])
                    self.n_Texp = np.array(j2n_Texp[df["jlower"]])
                elif codelv == "a1":
                    j2alpha_ref, j2n_Texp = exomolapi.make_j2b(
                        bdat,
                        alpha_ref_default=self.alpha_ref_def,
                        n_Texp_default=self.n_Texp_def,
                        jlower_max=np.max(df["jlower"]),
                    )
                    jj2alpha_ref, jj2n_Texp = exomolapi.make_jj2b(
                        bdat,
                        j2alpha_ref_def=j2alpha_ref,
                        j2n_Texp_def=j2n_Texp,
                        jupper_max=np.max(df["jupper"]),
                    )
                    self.alpha_ref = np.array(jj2alpha_ref[df["jlower"], df["jupper"]])
                    self.n_Texp = np.array(jj2n_Texp[df["jlower"], df["jupper"]])
            except FileNotFoundError:
                print(
                    "Warning: Cannot load .broad. The default broadening parameters are used."
                )
                self.alpha_ref = np.array(
                    self.alpha_ref_def * np.ones_like(df["jlower"])
                )
                self.n_Texp = np.array(self.n_Texp_def * np.ones_like(df["jlower"]))

        else:
            print("The default broadening parameters are used.")
            self.alpha_ref = np.array(self.alpha_ref_def * np.ones_like(df["jlower"]))
            self.n_Texp = np.array(self.n_Texp_def * np.ones_like(df["jlower"]))

        # Add values
        self.add_column(df, "alpha_ref", self.alpha_ref)
        self.add_column(df, "n_Texp", self.n_Texp)

    def QT_interp(self, T):
        """interpolated partition function

        Parameters
        ----------
        T: temperature

        Returns
        -------
        Q(T) interpolated in jnp.array

        """
        return np.interp(T, self.T_gQT, self.gQT)

    def qr_interp(self, T):
        """interpolated partition function ratio

        Parameters
        ----------
        T: temperature

        Returns
        -------
        qr(T)=Q(T)/Q(Tref) interpolated in jnp.array

        """
        return self.QT_interp(T) / self.QT_interp(self.Tref)

    def download(self, molec, extension, numtag=None):
        """Downloading Exomol files

        Parameters
        ----------
        molec: like "12C-16O__Li2015"
        extension: extension list e.g. [".pf",".def",".trans.bz2",".states.bz2",".broad"]
        numtag: number tag of transition file if exists. e.g. "11100-11200"

        Notes
        -----
        The download URL is written in exojax.utils.url.

        """
        import os
        import urllib.request

        tag = molec.split("__")
        molname_simple = e2s(tag[0])

        for ext in extension:
            if ext == ".trans.bz2" and numtag is not None:
                ext = "__" + numtag + ext

            if ext == ".broad":
                pfname_arr = [
                    tag[0] + "__H2" + ext,
                    tag[0] + "__He" + ext,
                    tag[0] + "__air" + ext,
                ]
                url = EXOMOL_URL + molname_simple + "/" + tag[0] + "/"
            else:
                pfname_arr = [molec + ext]
                url = EXOMOL_URL + molname_simple + "/" + tag[0] + "/" + tag[1] + "/"

            for pfname in pfname_arr:
                pfpath = url + pfname
                os.makedirs(str(self.path), exist_ok=True)
                print("Downloading " + pfpath)
                try:
                    urllib.request.urlretrieve(pfpath, str(self.path / pfname))
                except HTTPError:
                    print(f"Error: Couldn't download {ext} file at {pfpath} and save.")

    def to_partition_function_tabulator(self):
        """Generate a :py:class:`~radis.levels.partfunc.PartFuncExoMol` object"""
        from radis.levels.partfunc import PartFuncExoMol

        return PartFuncExoMol(self.isotope_fullname, self.T_gQT, self.gQT)


if __name__ == "__main__":

    from radis.test.io.test_exomol import (
        test_calc_exomol_spectrum,
        test_calc_exomol_vs_hitemp,
        test_exomol_parsing_functions,
    )

    test_exomol_parsing_functions()
    test_calc_exomol_spectrum()
    test_calc_exomol_vs_hitemp()
