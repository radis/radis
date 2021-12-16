"""Molecular database (MDB) class

   * MdbExomol is the MDB for ExoMol

Initial code borrowed from the `Exojax <https://github.com/HajimeKawahara/exojax>`__
code (which you should also have a look at !), by @HajimeKawahara, under MIT License.

"""
import pathlib
import warnings

import numpy as np

# from exojax.spec import hapi, exomolapi, exomol
# from exojax.spec.hitran import gamma_natural as gn
import pandas as pd

try:
    from . import exomolapi
    from .exomol_utils import e2s
except ImportError:  # if local import
    from radis.io import exomolapi
    from radis.io.exomol_utils import e2s

from radis.db.classes import get_molecule_identifier

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

    from radis.db.classes import EXOMOL_ONLY_ISOTOPES_NAMES

    if (molecule, isotope) in EXOMOL_ONLY_ISOTOPES_NAMES:
        return EXOMOL_ONLY_ISOTOPES_NAMES[(molecule, isotope)]
    else:
        # Read and convert from HITRAN molecules
        from radis.db.molparam import MolParams

        mp = MolParams()
        return mp.get(molecule, isotope, "isotope_name_exomol")


def get_exomol_database_list(molecule, isotope_full_name):
    """Parse ExoMol website and return list of available databases, and recommended database

    Parameters
    ----------
    molecule: str
    isotope_full_name: str
        isotope full name (ex. ``12C-1H4`` for CH4,1). Get it from
        :py:func:`radis.io.exomol.get_exomol_full_isotope_name`

    Returns
    -------

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

    url = f"https://exomol.com/data/molecules/{molecule}/{isotope_full_name}"
    try:
        response = urlopen(url).read()
    except HTTPError as err:
        raise ValueError(f"HTTPError opening url={url}") from err

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

    return databases, databases_recommended[0]


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
    cache: bool, or ``'regen'``
        if ``True``, use existing HDF5 file. If ``False`` or ``'regen'``, rebuild it.
        Default ``True``.
    verbose: bool
    clean_cache_files: bool
        if ``True`` clean downloaded cache files after HDF5 are created.
    return_local_path: bool
        if ``True``, also returns the path of the local database file.
    return_partition_function: bool
        if ``True``, also returns a :py:class:`~radis.levels.partfunc.PartFuncExoMol` object.
    engine: 'vaex', 'feather'
        which memory-mapping library to use. If 'default' use the value from ~/radis.json

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

    .. minigallery:: radis.io.exomol.fetch_exomol

    Notes
    -----
    if using ``load_only_wavenum_above/below`` or ``isotope``, the whole
    database is anyway downloaded and uncompressed to ``local_databases``
    fast access .HDF5 files (which will take a long time on first call). Only
    the expected wavenumber range & isotopes are returned. The .HFD5 parsing uses
    :py:func:`~radis.io.hdf5.hdf2df`

    See Also
    --------
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

    if load_wavenum_min is None:
        load_wavenum_min = -np.inf
    if load_wavenum_max is None:
        load_wavenum_max = np.inf
    mdb = MdbExomol(local_path, [load_wavenum_min, load_wavenum_max], engine=engine)

    # Attributes of the DataFrame
    from radis.db.classes import HITRAN_MOLECULES

    attrs = {}
    if molecule in HITRAN_MOLECULES:
        attrs["id"] = get_molecule_identifier(
            molecule
        )  # HITRAN id-number (if available)
    attrs["molecule"] = molecule
    attrs["iso"] = isotope

    df = mdb.to_df(attrs=attrs)

    # Replace Linestrength with Line Intensity taking into account Terrestrial isotopic abundance
    from radis.db.molparam import MolParams

    try:
        Ia = MolParams().get(molecule, isotope, "abundance")
    except NotImplementedError:
        from radis import config

        Ia = config["molparams"]["abundances"][molecule][str(isotope)]

    df["Sij0"] *= Ia
    df.rename(columns={"Sij0": "int"}, inplace=True)

    # Return:
    out = [df]

    if return_local_path:
        out.append(str(mdb.path))
    if return_partition_function:
        assert return_local_path
        out.append(mdb.to_partition_function_tabulator())

    return out


class MdbExomol(object):
    """molecular database of ExoMol

    MdbExomol is a class for ExoMol.

    Attributes:
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

    """

    # TODO : inherit from DatabaseManager or similar

    # @dev: In exojax this class is defined in exojax/spec/moldb.py
    # see https://github.com/HajimeKawahara/exojax/blob/develop/src/exojax/spec/moldb.py
    # It uses Jax arrays (jnp). Here RADIS uses Numpy arrays.

    def __init__(
        self,
        path,
        nurange=[-np.inf, np.inf],
        margin=1.0,
        crit=-np.inf,
        bkgdatm="H2",
        broadf=True,
        engine="vaex",
        verbose=True,
    ):
        """Molecular database for Exomol form

        Args:
           path: path for Exomol data directory/tag. For instance, "/home/CO/12C-16O/Li2015"
           nurange: wavenumber range list (cm-1) or wavenumber array
           margin: margin for nurange (cm-1)
           crit: line strength lower limit for extraction
           bkgdatm: background atmosphere for broadening. e.g. H2, He,
           broadf: if False, the default broadening parameters in .def file is used

        Note:
           The trans/states files can be very large. For the first time to read it, we convert it to the feather-format. After the second-time, we use the feather format instead.

        """
        from os import environ

        if engine == "default":
            from radis import config

            engine = config["MEMORY_MAPPING_ENGINE"]
            # Quick fix for #401
            if engine == "auto":
                # "auto" uses "vaex" in most cases unless you're using the Spyder IDE (where it may result in freezes).
                # see https://github.com/spyder-ide/spyder/issues/16183.
                # and https://github.com/radis/radis/issues/401
                if any("SPYDER" in name for name in environ):
                    if verbose >= 3:
                        print(
                            "Spyder IDE detected. Memory-mapping-engine set to 'feather' (less powerful than 'vaex' but Spyder user experience freezes). See https://github.com/spyder-ide/spyder/issues/16183. Change this behavior by setting the radis.config['MEMORY_MAPPING_ENGINE'] key"
                        )
                    engine = "feather"  # for ExoMol database
                # temp fix for vaex not building on RTD
                # see https://github.com/radis/radis/issues/404
                elif any("READTHEDOCS" in name for name in environ):
                    engine = "feather"
                    if verbose >= 3:
                        print(
                            f"ReadTheDocs environment detected. Memory-mapping-engine set to '{engine}'. See https://github.com/radis/radis/issues/404"
                        )
                else:
                    engine = "vaex"

        if engine == "vaex":
            import vaex
        elif engine == "feather":
            pass
            # vaex will be the future default, but it still fails on some systems
            # Ex : on Spyder : See https://github.com/spyder-ide/spyder/issues/16183
            # We keep "feather" as backup)
        else:
            raise NotImplementedError(f"{engine} is not implemented")
        self.engine = engine  # which memory mapping engine to use : 'vaex', 'pytables' (HDF5), 'feather'

        self.path = pathlib.Path(path)
        t0 = self.path.parents[0].stem
        molec = t0 + "__" + str(self.path.stem)
        self.bkgdatm = bkgdatm
        print("Background atmosphere: ", self.bkgdatm)
        molecbroad = t0 + "__" + self.bkgdatm

        self.crit = crit
        self.margin = margin
        self.nurange = [np.min(nurange), np.max(nurange)]
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
        self.isotope = 1  # Placeholder. TODO : impement parsing of other isotopes.

        # load def
        # @minou: I feel this next line should just update the self.stuff
        dic_def = exomolapi.read_def(self.def_file)  # approx. 3 ms
        self.n_Texp_def = dic_def["n_Texp"]
        self.alpha_ref_def = dic_def["alpha_ref"]
        self.molmass = dic_def["molmass"]

        #  default n_Texp value if not given
        if self.n_Texp_def is None:
            warnings.warn(
                Warning(
                    """
                    No default broadening exponent in def file. Assigned n = 0.5
                    """
                )
            )
            self.n_Texp_def = 0.5
        #  default alpha_ref value if not given
        if self.alpha_ref_def is None:
            warnings.warn(
                Warning(
                    """
                    No default broadening in def file. Assigned alpha_ref = 0.07
                    """
                )
            )
            self.alpha_ref_def = 0.07

        # load states
        if engine == "feather":
            if self.states_file.with_suffix(".feather").exists():
                states = pd.read_feather(self.states_file.with_suffix(".feather"))
            else:
                print(
                    "Note: Caching states data to the feather format. After the second time, it will become much faster."
                )
                states = exomolapi.read_states(self.states_file, dic_def, engine="csv")
                states.to_feather(self.states_file.with_suffix(".feather"))
            ndstates = states.to_numpy()[
                :, :4
            ]  # the i, E, g, J are in the 4 first columns
        elif engine == "vaex":
            if self.states_file.with_suffix(".bz2.hdf5").exists():
                states = vaex.open(self.states_file.with_suffix(".bz2.hdf5"))
                ndstates = vaex.array_types.to_numpy(states)
            else:
                print(
                    "Note: Caching states data to the hdf5 format with vaex. After the second time, it will become much faster."
                )
                states = exomolapi.read_states(self.states_file, dic_def, engine="vaex")
                ndstates = vaex.array_types.to_numpy(states)

        # load pf
        pf = exomolapi.read_pf(self.pf_file)
        self.gQT = pf["QT"].to_numpy()  # grid QT
        self.T_gQT = pf["T"].to_numpy()  # T forgrid QT

        # trans file(s)
        print("Reading transition file")

        # Compute linestrengths or retrieve them from cache
        self.Tref = 296.0
        self.QTref = np.array(self.QT_interp(self.Tref))

        if dic_def["numinf"] is None:
            self.trans_file = self.path / pathlib.Path(molec + ".trans.bz2")
            if not self.trans_file.exists():
                self.download(molec, [".trans.bz2"])

            if engine == "vaex":
                if self.trans_file.with_suffix(".hdf5").exists():
                    trans = vaex.open(self.trans_file.with_suffix(".hdf5"))
                    cdt = 1
                    if not np.isneginf(self.nurange[0]):
                        cdt *= trans.nu_lines > self.nurange[0] - self.margin
                    if not np.isinf(self.nurange[1]):
                        cdt *= trans.nu_lines < self.nurange[1] - self.margin
                    if not np.isneginf(self.crit):
                        cdt *= trans.Sij0 > self.crit
                    if cdt != 1:
                        trans = trans[cdt]
                    ndtrans = vaex.array_types.to_numpy(trans)

                    # mask has been alraedy applied
                    mask_needed = False
                else:
                    print(
                        "Note: Caching line transition data to the HDF5 format with vaex. After the second time, it will become much faster."
                    )
                    trans = exomolapi.read_trans(self.trans_file, engine="vaex")
                    ndtrans = vaex.array_types.to_numpy(trans)

                    # mask needs to be applied
                    mask_needed = True
            elif engine == "feather":
                if self.trans_file.with_suffix(".feather").exists():
                    trans = pd.read_feather(self.trans_file.with_suffix(".feather"))
                else:
                    print(
                        "Note: Caching line transition data to the feather format. After the second time, it will become much faster."
                    )
                    trans = exomolapi.read_trans(self.trans_file, engine="csv")
                    trans.to_feather(self.trans_file.with_suffix(".feather"))
                ndtrans = trans.to_numpy()
                # mask needs to be applied   (in feather mode we don't sleect wavneumbers)
                mask_needed = True

            # compute gup and elower
            (
                self._A,
                self.nu_lines,
                self._elower,
                self._gpp,
                self._jlower,
                self._jupper,
                mask_zeronu,
                self._quantumNumbers,
            ) = exomolapi.pickup_gE(ndstates, ndtrans, self.trans_file, dic_def)

            if engine == "vaex" and self.trans_file.with_suffix(".hdf5").exists():
                # if engine == 'feather' we recompute all the time
                # Todo : we should get the column name instead ?
                self.Sij0 = ndtrans[:, 4]
            else:
                ##Line strength:
                from radis.lbl.base import (  # TODO: move elsewhere
                    linestrength_from_Einstein,
                )

                self.Sij0 = linestrength_from_Einstein(
                    A=self._A,
                    gu=self._gpp,
                    El=self._elower,
                    Ia=1,  #  Sij0 is a linestrength calculated without taking into account isotopic abundance (unlike line intensity parameter of HITRAN. In RADIS this is corrected for in fetch_exomol()  )
                    nu=self.nu_lines,
                    Q=self.QTref,
                    T=self.Tref,
                )

                trans["nu_positive"] = mask_zeronu
                if engine == "vaex":
                    # exclude the lines whose nu_lines evaluated inside exomolapi.pickup_gE (thus sometimes different from the "nu_lines" column in trans) is not positive

                    trans = trans[trans.nu_positive].extract()
                    trans.drop("nu_positive", inplace=True)
                else:
                    if False in mask_zeronu:
                        raise NotImplementedError(
                            "some wavenumber is not defined;  masking not impleemtend so far in 'feather' engine"
                        )

                trans["nu_lines"] = self.nu_lines
                trans["Sij0"] = self.Sij0

                if engine == "vaex":
                    trans.export(self.trans_file.with_suffix(".hdf5"))
                #  TODO : implement masking in 'feather' mode

        else:  # dic_def["numinf"] is not None
            imin = (
                np.searchsorted(dic_def["numinf"], nurange[0], side="right") - 1
            )  # left side
            imax = (
                np.searchsorted(dic_def["numinf"], nurange[1], side="right") - 1
            )  # left side
            self.trans_file = []
            for k, i in enumerate(range(imin, imax + 1)):
                trans_file = self.path / pathlib.Path(
                    molec + "__" + dic_def["numtag"][i] + ".trans.bz2"
                )
                if not trans_file.exists():
                    self.download(
                        molec, extension=[".trans.bz2"], numtag=dic_def["numtag"][i]
                    )
                if engine == "feather":
                    trans_feather_path = trans_file.with_suffix(".bz2.feather")
                    if trans_feather_path.exists():
                        trans = pd.read_feather(trans_feather_path)
                    else:
                        print(
                            "Note: Caching line transition data to the feather format. After the second time, it will become much faster."
                        )
                        trans = exomolapi.read_trans(trans_file, engine="csv")
                        trans.to_feather(trans_feather_path)
                        #!!TODO:restrict NOW the trans size to avoid useless overload of memory and CPU
                        # trans = trans[(trans['nu'] > self.nurange[0] - self.margin) & (trans['nu'] < self.nurange[1] + self.margin)]
                    ndtrans = trans.to_numpy()
                    # mask needs to be applied   (in feather mode we don't sleect wavneumbers)
                    mask_needed = True
                elif engine == "vaex":
                    if trans_file.with_suffix(".hdf5").exists():
                        trans = vaex.open(trans_file.with_suffix(".hdf5"))
                        cdt = 1
                        if not np.isneginf(self.nurange[0]):
                            cdt *= trans.nu_lines > self.nurange[0] - self.margin
                        if not np.isinf(self.nurange[1]):
                            cdt *= trans.nu_lines < self.nurange[1] - self.margin
                        if not np.isneginf(self.crit):
                            cdt *= trans.Sij0 > self.crit
                        if cdt != 1:
                            trans = trans[cdt]
                        ndtrans = vaex.array_types.to_numpy(trans)

                        # mask has been already applied
                        mask_needed = False
                    else:
                        print(
                            "Note: Caching line transition data to the HDF5 format with vaex. After the second time, it will become much faster."
                        )
                        trans = exomolapi.read_trans(trans_file, engine="vaex")
                        ndtrans = vaex.array_types.to_numpy(trans)

                        # mask needs to be applied
                        mask_needed = True
                self.trans_file.append(trans_file)

                # compute gup and elower
                if k == 0:
                    (
                        self._A,
                        self.nu_lines,
                        self._elower,
                        self._gpp,
                        self._jlower,
                        self._jupper,
                        mask_zeronu,
                        self._quantumNumbers,
                    ) = exomolapi.pickup_gE(ndstates, ndtrans, trans_file, dic_def)

                    if engine == "vaex" and trans_file.with_suffix(".hdf5").exists():
                        # Sij0x already stored in cache file
                        self.Sij0 = ndtrans[:, 4]
                    else:
                        ##Recompute Line strength:
                        from radis.lbl.base import (  # TODO: move elsewhere
                            linestrength_from_Einstein,
                        )

                        self.Sij0 = linestrength_from_Einstein(
                            A=self._A,
                            gu=self._gpp,
                            El=self._elower,
                            Ia=1,  #  Sij0 is a linestrength calculated without taking into account isotopic abundance (unlike line intensity parameter of HITRAN. In RADIS this is corrected for in fetch_exomol()  )
                            nu=self.nu_lines,
                            Q=self.QTref,
                            T=self.Tref,
                        )

                        # exclude the lines whose nu_lines evaluated inside exomolapi.pickup_gE (thus sometimes different from the "nu_lines" column in trans) is not positive
                        if engine == "vaex":
                            trans["nu_positive"] = mask_zeronu
                            trans = trans[trans.nu_positive].extract()
                            trans.drop("nu_positive", inplace=True)
                        else:
                            if False in mask_zeronu:
                                raise NotImplementedError(
                                    "some wavenumber is not defined;  masking not impleemtend so far in 'feather' engine"
                                )

                        trans["nu_lines"] = self.nu_lines
                        trans["Sij0"] = self.Sij0
                else:  # k!=0

                    (
                        Ax,
                        nulx,
                        elowerx,
                        gppx,
                        jlowerx,
                        jupperx,
                        mask_zeronu,
                        quantumNumbersx,
                    ) = exomolapi.pickup_gE(ndstates, ndtrans, trans_file, dic_def)
                    if engine == "vaex" and trans_file.with_suffix(".hdf5").exists():
                        # Sij0x already stored in cache file
                        Sij0x = ndtrans[:, 4]
                    else:
                        ##Recompute Line strength:
                        from radis.lbl.base import (  # TODO: move elsewhere
                            linestrength_from_Einstein,
                        )

                        Sij0x = linestrength_from_Einstein(
                            A=Ax,
                            gu=gppx,
                            El=elowerx,
                            Ia=1,  #  Sij0 is a linestrength calculated without taking into account isotopic abundance (unlike line intensity parameter of HITRAN. In RADIS this is corrected for in fetch_exomol()  )
                            nu=nulx,
                            Q=self.QTref,
                            T=self.Tref,
                        )

                        # exclude the lines whose nu_lines evaluated inside exomolapi.pickup_gE (thus sometimes different from the "nu_lines" column in trans) is not positive
                        if engine == "vaex":
                            trans["nu_positive"] = mask_zeronu
                            trans = trans[trans.nu_positive].extract()
                            trans.drop("nu_positive", inplace=True)
                        else:
                            if False in mask_zeronu:
                                raise NotImplementedError(
                                    "some wavenumber is not defined;  masking not impleemtend so far in 'feather' engine"
                                )

                        trans["nu_lines"] = nulx
                        trans["Sij0"] = Sij0x

                    self._A = np.hstack([self._A, Ax])
                    self.nu_lines = np.hstack([self.nu_lines, nulx])
                    self.Sij0 = np.hstack([self.Sij0, Sij0x])
                    self._elower = np.hstack([self._elower, elowerx])
                    self._gpp = np.hstack([self._gpp, gppx])
                    self._jlower = np.hstack([self._jlower, jlowerx])
                    self._jupper = np.hstack([self._jupper, jupperx])
                    for key in self._quantumNumbers.keys():
                        self._quantumNumbers[key] = np.hstack(
                            [self._quantumNumbers[key], quantumNumbersx[key]]
                        )

                    if (
                        engine == "vaex"
                        and not trans_file.with_suffix(".hdf5").exists()
                    ):
                        trans.export(trans_file.with_suffix(".hdf5"))

        ### MASKING ###
        mask = (
            (self.nu_lines > self.nurange[0] - self.margin)
            * (self.nu_lines < self.nurange[1] + self.margin)
            * (self.Sij0 > self.crit)
        )

        self.masking(mask, mask_needed)

    def masking(self, mask, mask_needed=True):
        """applying mask and (re)generate jnp.arrays

        Args:
           mask: mask to be applied. self.mask is updated.
           mask_needed: whether mask needs to be applied or not

        """
        # TODO : with Vaex-HDF5, selection should happen on disk.

        # numpy float 64 Do not convert them jnp array
        if mask_needed:
            self.nu_lines = self.nu_lines[mask]
            self.Sij0 = self.Sij0[mask]
            self._A = self._A[mask]
            self._elower = self._elower[mask]
            self._gpp = self._gpp[mask]
            self._jlower = self._jlower[mask]
            self._jupper = self._jupper[mask]
            for key, array in self._quantumNumbers.items():
                self._quantumNumbers[key] = array[mask]

        # jnp arrays
        self.dev_nu_lines = np.array(self.nu_lines)
        self.logsij0 = np.array(np.log(self.Sij0))
        self.A = np.array(self._A)
        # self.gamma_natural = gn(self.A)   #  Natural Broadening neglected in RADIS
        self.elower = np.array(self._elower)
        self.gpp = np.array(self._gpp)
        self.jlower = np.array(self._jlower, dtype=int)
        self.jupper = np.array(self._jupper, dtype=int)
        ##Broadening parameters
        self.set_broadening()

    def set_broadening(self, alpha_ref_def=None, n_Texp_def=None):
        """setting broadening parameters

        Args:
           alpha_ref: set default alpha_ref and apply it. None=use self.alpha_ref_def
           n_Texp_def: set default n_Texp and apply it. None=use self.n_Texp_def
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
                        jlower_max=np.max(self._jlower),
                    )
                    self.alpha_ref = np.array(j2alpha_ref[self._jlower])
                    self.n_Texp = np.array(j2n_Texp[self._jlower])
                elif codelv == "a1":
                    j2alpha_ref, j2n_Texp = exomolapi.make_j2b(
                        bdat,
                        alpha_ref_default=self.alpha_ref_def,
                        n_Texp_default=self.n_Texp_def,
                        jlower_max=np.max(self._jlower),
                    )
                    jj2alpha_ref, jj2n_Texp = exomolapi.make_jj2b(
                        bdat,
                        j2alpha_ref_def=j2alpha_ref,
                        j2n_Texp_def=j2n_Texp,
                        jupper_max=np.max(self._jupper),
                    )
                    self.alpha_ref = np.array(jj2alpha_ref[self._jlower, self._jupper])
                    self.n_Texp = np.array(jj2n_Texp[self._jlower, self._jupper])
            except:
                print(
                    "Warning: Cannot load .broad. The default broadening parameters are used."
                )
                self.alpha_ref = np.array(
                    self.alpha_ref_def * np.ones_like(self._jlower)
                )
                self.n_Texp = np.array(self.n_Texp_def * np.ones_like(self._jlower))

        else:
            print("The default broadening parameters are used.")
            self.alpha_ref = np.array(self.alpha_ref_def * np.ones_like(self._jlower))
            self.n_Texp = np.array(self.n_Texp_def * np.ones_like(self._jlower))

    def QT_interp(self, T):
        """interpolated partition function

        Args:
           T: temperature

        Returns:
           Q(T) interpolated in jnp.array

        """
        return np.interp(T, self.T_gQT, self.gQT)

    def qr_interp(self, T):
        """interpolated partition function ratio

        Args:
           T: temperature

        Returns:
           qr(T)=Q(T)/Q(Tref) interpolated in jnp.array

        """
        return self.QT_interp(T) / self.QT_interp(self.Tref)

    def download(self, molec, extension, numtag=None):
        """Downloading Exomol files

        Args:
           molec: like "12C-16O__Li2015"
           extension: extension list e.g. [".pf",".def",".trans.bz2",".states.bz2",".broad"]
           numtag: number tag of transition file if exists. e.g. "11100-11200"

        Note:
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
                except:
                    print("Error: Couldn't download " + ext + " file and save.")

    def to_partition_function_tabulator(self):
        from radis.levels.partfunc import PartFuncExoMol

        return PartFuncExoMol(self.isotope_fullname, self.T_gQT, self.gQT)

    def to_df(self, attrs={}):
        """Export Line Database to a RADIS-friendly Pandas DataFrame

        Parameters
        ----------
        attrs: dict
            add these as attributes of the DataFrame
        """

        ##Broadening parameters
        self.set_broadening()

        # TODO : define RADIS database format somewhere else; with description of the column names.
        df = pd.DataFrame(
            {
                "wav": self.nu_lines,
                "Sij0": self.Sij0,  #  linestrength (not corrected for isotopic abundance)
                "A": self._A,
                "airbrd": self.alpha_ref,  # temperature dependance exponent for air broadening
                # "selbrd": None,                  # self-temperature dependant exponent. Not implementedi in ExoMol
                "El": self._elower,
                # "Tdpair": None,
                # "Pshft": None,
                "ju": self._jupper,
                "jl": self._jlower,
                "gpp": self._gpp,
                "Tdpair": self.n_Texp,  # temperature dependance exponent. Here we use Tdair; no Tdpsel. TODO
            }
        )

        # Check format :
        for key, array in self._quantumNumbers.items():
            try:
                df[key] = pd.to_numeric(array, errors="raise")
            except:
                df[key] = array
                df[key] = df[key].astype("str")

        for k, v in attrs.items():
            df.attrs[k] = v

        assert self.Tref == 296  # default in RADIS

        return df


if __name__ == "__main__":
    # mdb=MdbExomol("/home/kawahara/exojax/data/CO/12C-16O/Li2015/")
    # mdb=MdbExomol("/home/kawahara/exojax/data/CH4/12C-1H4/YT34to10/",nurange=[6050.0,6150.0])

    # %% Test  by overriding Spectrumfactory's DataFrame df0
    # mdb = MdbExomol(".database/H2O/1H2-16O/POKAZATEL", [4310.0, 4320.0], crit=1.0e-45)
    # df = mdb.to_df()
    # from radis import SpectrumFactory
    # sf = SpectrumFactory(
    #     4310,
    #     4320,
    #     molecule="H2O",
    #     isotope="1",
    # )
    # sf.fetch_databank(
    #     "hitran"
    # )  # placeholder. Load lines (will be replaced), load partition function.
    # sf.df0 = df  # override.
    # s = sf.eq_spectrum(500, name="ExoMol")
    # # sf.fetch_databank('hitran')  # placeholder. Load lines (will be replaced), load partition function.
    # # s_hit = sf.eq_spectrum(500, name='HITRAN')

    #%% Test by direct caclulation
    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_exomol.py"]))
