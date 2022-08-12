# -*- coding: utf-8 -*-
"""
Summary
-------

Functions to parse the ``~/radis.json`` :ref:`Configuration file <label_lbl_config_file>`

Notes
-----

Create a ~/radis.json file in your HOME that contains all machine-specific information
(e.g: path to databanks). See :data:`~radis.misc.config.DBFORMAT` for expected
format

Routine Listing
---------------

- :func:`~radis.misc.config.get_user_config_configformat`
- :func:`~radis.misc.config.getDatabankEntries_configformat`
- :func:`~radis.misc.config.getDatabankList_configformat`
- :func:`~radis.misc.config.addDatabankEntries_configformat`
- :func:`~radis.misc.config.diffDatabankEntries`
- :func:`~radis.misc.config.printDatabankEntries_configformat`
- :func:`~radis.misc.config.printDatabankList_configformat`
- :func:`~radis.misc.config.get_user_config`
- :func:`~radis.misc.config.convertRadisToJSON`
- :func:`~radis.misc.config.getDatabankEntries`
- :func:`~radis.misc.config.getDatabankList`
- :func:`~radis.misc.config.addDatabankEntries`
- :func:`~radis.misc.config.printDatabankEntries`
- :func:`~radis.misc.config.printDatabankList`

-------------------------------------------------------------------------------


"""


import configparser
import json
import os
import warnings
from json import JSONDecodeError
from os.path import dirname, exists, expanduser, join

import hjson

from radis.misc.basics import compare_dict, compare_lists, stdpath
from radis.misc.utils import DatabankNotFound, getProjectRoot
from radis.misc.warning import DatabaseAlreadyExists

# %% Functions to parse radis/config.json

CONFIG_PATH_DEFAULT = join(getProjectRoot(), "default_radis.json")
CONFIG_PATH_JSON = join(expanduser("~"), "radis.json")
CONFIG_PATH_OLD = join(expanduser("~"), ".radis")

assert exists(CONFIG_PATH_DEFAULT)


def get_config(configpath=CONFIG_PATH_JSON):
    """Read the default RADIS config.json file `configpath` (default
    :py:attr:`~radis.misc.config.CONFIG_PATH_JSON` and override it with the
    entries of the user config file ``~/.radis``
    (:py:attr:`~radis.misc.config.CONFIG_PATH_DEFAULT`"""

    jsonfile = CONFIG_PATH_DEFAULT
    with open(jsonfile) as f:
        try:
            config = hjson.load(f)
            # read with Hjson to allow comments in CONFIG_PATH_DEFAULT
            # we do not allow comments in (user) CONFIG_PATH_JSON
            # because it is read/written automatically. So far there is no
            # two-way JSON editor with comments. See https://github.com/radis/radis/issues/328
        except JSONDecodeError as err:
            raise JSONDecodeError(
                "Error reading '{0}' (line {2} col {3}): \n{1}".format(
                    jsonfile, err.msg, err.lineno, err.colno
                ),
                err.doc,
                err.pos,
            ) from err

    user_config = get_user_config(configpath)

    config.update(user_config)

    return config


# %% Functions to parse ~/.radis file

DBFORMAT = r"""
--------------------------

[MY-HITEMP-CO2]                  #  your databank name: use this in calc_spectrum()
                                 #  or SpectrumFactory.load_databank()
info = HITEMP 2010 databank      #  whatever you want
path =                           #  no "", multipath allowed
       D:\Databases\HITEMP-CO2\hitemp_07
       D:\Databases\HITEMP-CO2\hitemp_08
       D:\Databases\HITEMP-CO2\hitemp_09
format = hitran                  #  'hitran' (HITRAN/HITEMP), 'cdsd-hitemp', 'cdsd-4000'
                                 # databank text file format. More info in
                                 # SpectrumFactory.load_databank function.
parfuncfmt:                      #  'cdsd', 'hapi', etc.
                                 # format to read tabulated partition function
                                 # file. If `hapi`, then HAPI (HITRAN Python
                                 # interface) is used to retrieve them (valid if
                                 # your databank is HITRAN data). HAPI is embedded
                                 # into RADIS. Check the version.
# Optional
# ----------
parfunc:                         #  path to tabulated partition function to use.
                                 # If `parfuncfmt` is `hapi` then `parfunc`
                                 # should be the link to the hapi.py file. If
                                 # not given, then the hapi.py embedded in RADIS
                                 # is used (check version)
levels_iso1                      #  path to energy levels (needed for non-eq
                                 # calculations). Default None
levels_iso2                      # etc
levels_iso4                      # etc
levelsfmt:                       #  'cdsd', etc.
                                 # how to read the previous file. Default None.
levelsZPE:                       #  zero-point-energy (cm-1): offset for all level
                                 # energies. Default 0 (if not given)

--------------------------"""
r"""str: Typical expected format of a ~/.radis entry::

    --------------------------

    [MY-HITEMP-CO2]                  #  your databank name: use this in calc_spectrum()
                                     #  or SpectrumFactory.load_databank()
    info = HITEMP 2010 databank      #  whatever you want
    path =                           #  no "", multipath allowed
           D:\Databases\HITEMP-CO2\hitemp_07
           D:\Databases\HITEMP-CO2\hitemp_08
           D:\Databases\HITEMP-CO2\hitemp_09
    format = hitran                  #  'hitran' (HITRAN/HITEMP), 'cdsd-hitemp', 'cdsd-4000'
                                     # databank text file format. List of all
                                     # formats in :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`
                                     # More info in
                                     # :py:meth:`~radis.lbl.loader.DatabankLoader.load_databank` function.
    parfuncfmt:                      #  'cdsd', 'hapi', etc.
                                     # format to read tabulated partition function
                                     # file. If `hapi`, then HAPI (HITRAN Python
                                     # interface) is used to retrieve them (valid if
                                     # your databank is HITRAN data). HAPI is embedded
                                     # into RADIS. Check the version.
                                     # List of all formats in :py:data:`~radis.lbl.loader.KNOWN_LVLFORMAT`
    # Optional
    # ----------
    parfunc:                         #  path to tabulated partition function to use.
                                     # If `parfuncfmt` is `hapi` then `parfunc`
                                     # should be the link to the hapi.py file. If
                                     # not given, then the hapi.py embedded in RADIS
                                     # is used (check version)
    levels_iso1                      #  path to energy levels (needed for non-eq
                                     # calculations). Default None
    levels_iso2                      # etc
    levels_iso4                      # etc
    levelsfmt:                       #  'cdsd', etc.
                                     # how to read the previous file. Default None.
    levelsZPE:                       #  zero-point-energy (cm-1): offset for all level
                                     # energies. Default 0 (if not given)

    --------------------------

For more information refer to the documentation: :ref:`Configuration file <label_lbl_config_file>` :

- Setup test databases with :py:func:`~radis.test.utils.setup_test_line_databases`
- `format` :  format: :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`


See Also
--------

:ref:`Configuration file <label_lbl_config_file>`,
:py:func:`~radis.misc.config.get_user_config`,
:py:meth:`~radis.lbl.loader.DatabankLoader.load_databank`
"""

DBFORMATJSON = r"""
--------------------------

{
"database": {                       # database key, stores all databases
    "MY-HITEMP-CO2": {              # your databank name: use this in calc_spectrum()
                                    # or SpectrumFactory.load_databank()
    "info": "HITEMP 2010 databank", # whatever you want
    "path":                         # no "", multipath allowed
       ["D:\Databases\HITEMP-CO2\hitemp_07",
        "D:\Databases\HITEMP-CO2\hitemp_08",
        "D:\Databases\HITEMP-CO2\hitemp_09"],
    "format": "hitran",             # 'hitran' (HITRAN/HITEMP), 'cdsd-hitemp', 'cdsd-4000'
                                    # databank text file format. More info in
                                    # SpectrumFactory.load_databank function.
    "parfuncfmt":                   # 'cdsd', 'hapi', etc.
                                    # format to read tabulated partition function
                                    # file. If `hapi`, then HAPI (HITRAN Python
                                    # interface) is used to retrieve them (valid if
                                    # your databank is HITRAN data). HAPI is embedded
                                    # into RADIS. Check the version.
    # Optional
    # ----------
    "parfunc":                      # path to tabulated partition function to use.
                                    # If `parfuncfmt` is `hapi` then `parfunc`
                                    # should be the link to the hapi.py file. If
                                    # not given, then the hapi.py embedded in RADIS
                                    # is used (check version)
    "levels_iso1":                  # path to energy levels (needed for non-eq
                                    # calculations). Default None
    "levels_iso2":                  # etc
    "levels_iso4":                  # etc
    "levelsfmt":                    # 'cdsd', etc.
                                    # how to read the previous file. Default None.
    "levelsZPE":                    # zero-point-energy (cm-1): offset for all level
    }                               # energies. Default 0 (if not given)
  }
}

--------------------------"""
r"""str: Typical expected format of a ~/radis.json entry in "database" key::

    --------------------------

    "MY-HITEMP-CO2": {              # your databank name: use this in calc_spectrum()
                                    # or SpectrumFactory.load_databank()
    "info": "HITEMP 2010 databank", # whatever you want
    "path":[                        # no "", multipath allowed
        "D:\Databases\HITEMP-CO2\hitemp_07"
        "D:\Databases\HITEMP-CO2\hitemp_08"
        "D:\Databases\HITEMP-CO2\hitemp_09"
    ],
    "format": "hitran",             # 'hitran' (HITRAN/HITEMP), 'cdsd-hitemp', 'cdsd-4000'
                                    # databank text file format. List of all
                                    # formats in :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`
                                    # More info in
                                    # :py:meth:`~radis.lbl.loader.DatabankLoader.load_databank` function.
    "parfuncfmt":                   # 'cdsd', 'hapi', etc.
                                    # format to read tabulated partition function
                                    # file. If `hapi`, then HAPI (HITRAN Python
                                    # interface) is used to retrieve them (valid if
                                    # your databank is HITRAN data). HAPI is embedded
                                    # into RADIS. Check the version.
                                    # List of all formats in :py:data:`~radis.lbl.loader.KNOWN_LVLFORMAT`
    # Optional
    # ----------
    "parfunc":                      # path to tabulated partition function to use.
                                    # If `parfuncfmt` is `hapi` then `parfunc`
                                    # should be the link to the hapi.py file. If
                                    # not given, then the hapi.py embedded in RADIS
                                    # is used (check version)
    "levels_iso1":                  # path to energy levels (needed for non-eq
                                    # calculations). Default None
    "levels_iso2":                  # etc
    "levels_iso4":                  # etc
    "levelsfmt":                    # 'cdsd', etc.
                                    # how to read the previous file. Default None.
    "levelsZPE":                    # zero-point-energy (cm-1): offset for all level
    }                               # energies. Default 0 (if not given)

    --------------------------

For more information refer to the documentation: :ref:`Configuration file <label_lbl_config_file>` :

- Setup test databases with :py:func:`~radis.test.utils.setup_test_line_databases`
- `format` :  format: :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`


See Also
--------

:ref:`Configuration file <label_lbl_config_file>`,
:py:func:`~radis.misc.config.get_user_config`,
:py:meth:`~radis.lbl.loader.DatabankLoader.load_databank`
"""


def get_user_config_configformat():
    r"""Read config file and returns it.

    Config file name is harcoded: :ref:`~/.radis <label_lbl_config_file>`
    """

    config = configparser.ConfigParser()
    configpath = CONFIG_PATH_OLD

    # Test ~/.radis exists
    if not exists(configpath):

        raise FileNotFoundError(
            "Create a `.radis` file in {0} to store links to ".format(
                dirname(configpath)
            )
            + "your local databanks. Format must be:\n {0}".format(DBFORMAT)
        )
    config.read(configpath)

    return config


def convertRadisToJSON(config_path_json, config_path_old=CONFIG_PATH_OLD):
    r"""Converts the ~/.radis file into json formatted file ~/radis.json

    Example
    -------
    original ~/.radis file format::

        [HITRAN-CO2-TEST]
        info = HITRAN 2016 database, CO2, 1 main isotope (CO2-626), bandhead: 2380-2398 cm-1 (4165-4200 nm)
        path = PATH_TO\radis\radis\test\files\hitran_co2_626_bandhead_4165_4200nm.par
        format = hitran
        parfuncfmt = hapi
        levelsfmt = radis

    -----------
    Converted ~/radis.json file format::

        {
        "database": {
            "HITRAN-CO2-TEST": {
                "info": "HITRAN 2016 database, CO2, 1 main isotope (CO2-626), bandhead: 2380-2398 cm-1 (4165-4200 nm)",
                "path": "PATH_TO\\radis\\radis\\test\\files\\hitran_co2_626_bandhead_4165_4200nm.par",
                "format": "hitran",
                "parfuncfmt": "hapi",
                "levelsfmt": "radis"
            }
          }
        }
    """

    # Loads configuration file ~/.radis
    config = get_user_config_configformat(config_path_old)

    # Variable to store data in JSON format
    config_json = {}

    # Converting configuration into JSON format and storing in config_json variable
    for i in config.sections():
        temp = {}
        for j in config[i]:
            # Checking if the `path` is multiple then creating list otherwise string
            if j == "path":
                if "\n" in config[i][j]:
                    store_list = config[i][j].split("\n")
                    # using remove() to remove empty list elements
                    while "" in store_list:
                        store_list.remove("")

                    if len(store_list) == 1:
                        temp[j] = store_list[0]
                    else:
                        temp[j] = store_list
                else:
                    temp[j] = config[i][j]
            else:
                # Adding to `temp` dictionaru
                temp[j] = config[i][j]

        config_json[i] = temp

    # Adding all config element to a `database` key
    config_final = {}
    config_final["database"] = config_json

    # Creating json file
    config_json_dir = config_path_json
    with open(config_json_dir, "w") as outfile:
        json.dump(config_final, outfile, indent=3)
    outfile.close()

    return


def init_radis_json(config_path_json, config_path_old=CONFIG_PATH_OLD):
    r"""Checks whether ``config_path_json`` (usually `~/radis.json`) exists.

    If not then we use
    :func:`~radis.misc.config.convertRadisToJSON`
    to convert `~/.radis` into `~/radis.json`
    """
    # If `~/.radis` and `~/radis.json` both found, raises DeprecationWarning
    if exists(config_path_json) and exists(config_path_old):
        raise (
            DeprecationWarning(
                f"{config_path_old} and {config_path_json} both found. "
                + "`~/.radis` config file was replaced by `~/radis.json` in version '0.9.29'."
                + f" Make sure all file content is in {config_path_json} then remove {config_path_old}.",
            )
        )
    elif exists(config_path_json):
        pass
    elif exists(config_path_old):
        convertRadisToJSON(config_path_json, config_path_old)
    else:
        # Create it from the default config path

        with open(CONFIG_PATH_DEFAULT) as f:
            try:
                default_config = hjson.load(f)
                # read with Hjson to allow comments in CONFIG_PATH_DEFAULT
                # we do not allow comments in (user) CONFIG_PATH_JSON
                # because it is read/written automatically. So far there is no
                # two-way JSON editor with comments. See https://github.com/radis/radis/issues/328
            except JSONDecodeError as err:
                raise JSONDecodeError(
                    "Error reading '{0}' (line {2} col {3}): \n{1}".format(
                        CONFIG_PATH_DEFAULT, err.msg, err.lineno, err.colno
                    ),
                    err.doc,
                    err.pos,
                ) from err

        # Creating json file
        config = {
            "_HELLO!": "This is your RADIS configuration file. Any key here will overwrite the values of the radis/default_radis.json file on your local computer (see online version of the latest release: https://github.com/radis/radis/blob/master/radis/default_radis.json). More information on the online docs : https://radis.readthedocs.io/en/latest/lbl/lbl.html#label-lbl-config-file"
        }
        config.update(default_config)
        with open(config_path_json, "w") as outfile:
            json.dump(config, outfile, indent=3)
        # outfile.close()

        # try:
        #     shutil.copy2(CONFIG_PATH_DEFAULT, config_path_json)
        # except Exception as err:
        #     print(
        #         f"Couldn't create RADIS configuration in {config_path_json} ({str(err)})"
        #     )
        # else:
        #     print(f"RADIS configuration file created in {config_path_json}")

    # Check it user file exists
    if not exists(config_path_json):
        raise FileNotFoundError(
            "Create a `radis.json file in {0} to store links to ".format(
                dirname(config_path_json)
            )
            + "your local databanks. Format must be:\n {0}".format(DBFORMATJSON)
            + "\n(it can be empty too)"
        )

    return


def get_user_config(configpath=CONFIG_PATH_JSON):
    r"""Read config file and returns it."""

    # First checking `~radis.json` exist or not, if not create it
    init_radis_json(configpath)

    # Checking file is empty
    if os.stat(configpath).st_size == 0:
        return {}

    # Load the JSON file
    with open(configpath) as f:
        try:
            config = json.load(f)
        except JSONDecodeError as err:
            raise JSONDecodeError(
                "Error reading '{0}' (line {2} col {3}): \n{1}".format(
                    configpath, err.msg, err.lineno, err.colno
                ),
                err.doc,
                err.pos,
            ) from err

    return config


def createConfigFile(config_path, verbose=True):

    # Safety check : file shouldnt exist
    if exists(config_path):
        raise FileExistsError(
            f"Config file already exist: {config_path}. You shouldnt be re-creating it. If sure, check it and delete it manually"
        )

    # generate ~/radis.json:
    open(config_path, "w").close()
    if verbose:
        print("Created ~/radis.json in {0}".format(dirname(config_path)))


def getDatabankEntries_configformat(dbname, get_extra_keys=[]):
    r"""Read :ref:`~/.radis <label_lbl_config_file>` config file and returns a dictionary of entries.

    Parameters
    ----------
    dbname: str
        database name in :ref:`~/.radis <label_lbl_config_file>` config file
    get_extra_keys: list
        read additional parameters on top of the usual Databank format keys :

    Notes
    -----
    Databank format::

        [MY-HITEMP-CO2]                  #  your databank name: use this in calc_spectrum()
                                         #  or SpectrumFactory.load_databank()
        info = HITEMP 2010 databank      # whatever you want
        path =                           # no "", multipath allowed
               D:\Databases\HITEMP-CO2\hitemp_07
               D:\Databases\HITEMP-CO2\hitemp_08
               D:\Databases\HITEMP-CO2\hitemp_09
        format = hitemp                  # 'hitran' (HITRAN / HITEMP), 'cdsd-hitemp', 'cdsd-4000'
                                         # Databank text file format. More info in
                                         # SpectrumFactory.load_databank function.

        # Optional:

        parfunc                          # path or 'USE_HAPI'
                                         # path to tabulated partition functions. If
                                         # `USE_HAPI`, then HAPI (HITRAN Python
                                         interface) [1]_ is used to retrieve them (valid
                                         if your databank is HITRAN data). HAPI
                                         is embedded into RADIS. Check the version.

        parfuncfmt:                      # 'cdsd'
                                         # format to read tabulated partition function
                                         # file. If `USE_HAPI` is given as `parfunc`
                                         # parameter then this line should not be used.

        levels_iso1                      # path to energy levels (needed for non-eq)
                                         # calculations.
        levels_iso2                      # etc
        levels_iso4                      # etc

        levelsfmt                        # 'cdsd'
                                         # how to read the previous file.

    References
    ----------
    .. [1] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_


    """
    warnings.warn(DeprecationWarning("config format changed to JSON in radis 0.9.29"))
    config = get_user_config_configformat()

    # Make sure it looks like a databank (path and format are given)
    try:
        config.get(dbname, "path")
        config.get(dbname, "format")
    except configparser.NoSectionError:
        msg = (
            "{1}\nDBFORMAT\n{0}\n".format(DBFORMAT, dbname)
            + "No databank named {0} in `{1}`. ".format(dbname, CONFIG_PATH_OLD)
            + "Available databanks: {0}. ".format(getDatabankList_configformat())
            + "See databank format above. More information in "
            + "https://radis.readthedocs.io/en/latest/lbl/lbl.html#configuration-file"
        )
        raise DatabankNotFound(msg)

    entries = dict(config.items(dbname))

    # Parse paths correctly
    entries["path"] = entries["path"].strip("\n").split("\n")

    # Merge all isotope-dependant levels into one dict
    iso_list = [k for k in entries if k.startswith("levels_iso")]
    if len(iso_list) > 0:
        levels = {}
        for k in iso_list:
            iso = float(k[10:])
            levels[iso] = entries.pop(k)
        entries["levels"] = levels
    else:
        if "levels" in entries:
            raise SyntaxError(
                "in {0}: `levels` replaced with ".format(CONFIG_PATH_OLD)
                + "`levels_iso#` where # is the isotope number"
            )

    return entries


def getDatabankList_configformat():
    r"""Get all databanks available in :ref:`~/.radis <label_lbl_config_file>`."""
    warnings.warn(DeprecationWarning("config format changed to JSON in radis 0.9.29"))

    config = get_user_config_configformat()

    # Get databank path and format
    validdb = []
    for dbname in config.sections():
        try:
            config.get(dbname, "path")
            config.get(dbname, "format")
        except configparser.NoSectionError:
            # not a db
            continue
        except configparser.NoOptionError:
            # not a db
            continue
        # looks like a db. add to db
        validdb.append(dbname)
    return validdb


def addDatabankEntries_configformat(dbname, dict_entries, verbose=True):
    r"""Add database dbname with entries from dict_entries.

    If database already exists in :ref:`~/.radis <label_lbl_config_file>`, raises an error

    Examples
    --------

    ::

        addDatabankEntries("HITEMP2010-CO2",
            {
                "info": "HITEMP2020 CO2 lines with TIPS-2017 for partition functions (equilibrium) and RADIS for rovibrational energies (nonequilibrium) ",
                "path": ["PATH/TO/HITEMP/CO2/*.par"],
                "format": "hitran",
                "parfuncfmt": "hapi",
                "levelsfmt": "radis",
            })
    """
    warnings.warn(DeprecationWarning("config format changed to JSON in radis 0.9.29"))

    # Get ~/.radis if exists, else create it
    try:
        dbnames = getDatabankList_configformat()
    except FileNotFoundError:
        # generate ~/.radis:
        dbnames = []
        open(CONFIG_PATH_OLD, "a").close()
        if verbose:
            print("Created ~/.radis in {0}".format(dirname(CONFIG_PATH_OLD)))

    # Check database doesnt exist
    if dbname in dbnames:
        raise DatabaseAlreadyExists(
            "Database already exists: {0}".format(dbname) + ". Cant add it"
        )

    # Add entries to parser
    config = configparser.ConfigParser()
    config[dbname] = {}

    if "info" in dict_entries:  # optional
        config[dbname]["info"] = dict_entries.pop("info")

    # ... Parse paths correctly
    if isinstance(dict_entries["path"], str):
        config[dbname]["path"] = dict_entries.pop("path")
    else:  # list
        config[dbname]["path"] = "\n       ".join(dict_entries.pop("path"))

    config[dbname]["format"] = dict_entries.pop("format")
    config[dbname]["parfuncfmt"] = dict_entries.pop("parfuncfmt")

    # Optional:
    # ... partition functions:
    if "parfunc" in dict_entries:
        config[dbname]["parfunc"] = dict_entries.pop("parfunc")

    # ... Split all isotopes in separate keys
    levels_dict = dict_entries.pop("levels", {})
    for iso, levels_iso in levels_dict.items():
        config[dbname]["levels_iso{0}".format(iso)] = levels_iso

    if "levelsfmt" in dict_entries:
        config[dbname]["levelsfmt"] = dict_entries.pop("levelsfmt")

    # Allow some informative parameters
    # ... "download_url", "download_date", "wavenumber_min", "wavenumber_max" generated by fetch_hitemp
    for entry in ["download_url", "download_date", "wavenumber_min", "wavenumber_max"]:
        if entry in dict_entries:
            config[dbname][entry] = dict_entries.pop(entry)

    # Check nothing is left
    if dict_entries != {}:
        raise ValueError("Unexpected keys: {0}".format(list(dict_entries.keys())))

    # Write to ~/.radis
    # ... Note: what if there is a PermissionError here? Try/except pass?
    with open(CONFIG_PATH_OLD, "a") as configfile:
        configfile.write("\n")
        config.write(configfile)
    if verbose:
        print("Added {0} database in {1}".format(dbname, CONFIG_PATH_OLD))

    return


def diffDatabankEntries(dict_entries1, dict_entries2, verbose=True):
    r"""Compare two Databank entries under dict format (i.e: output of
    getDatabankEntries)

    Returns None if no differences are found, or the first different key
    """

    k = None

    try:
        verbose_compare = "if_different" if verbose else False
        assert len(dict_entries1) == len(dict_entries2)
        assert (
            compare_lists(
                list(dict_entries1.keys()),
                list(dict_entries2.keys()),
                verbose=verbose_compare,
            )
            == 1
        )
        for k in dict_entries1.keys():
            v1 = dict_entries1[k]
            v2 = dict_entries2[k]
            if k in ["info", "format", "parfunc", "parfuncfmt", "levelsfmt"]:
                assert v1.lower() == v2.lower()
            elif k in ["path"]:
                assert (
                    compare_lists(
                        [stdpath(path1).lower() for path1 in v1],
                        [stdpath(path2).lower() for path2 in v2],
                        verbose=verbose_compare,
                    )
                    == 1
                )
            elif k in ["levels"]:
                assert (
                    compare_dict(
                        v1,
                        v2,
                        compare_as_paths=list(v1.keys()),
                        verbose=verbose_compare,
                    )
                    == 1
                )
            else:
                raise ValueError("Unexpected key:", k)

        return None

    except AssertionError:
        if verbose:
            print("Key doesnt match:", k)
        return k


def printDatabankEntries_configformat(dbname, crop=200):
    r"""Print databank info.

    Parameters
    ----------
    dbname: str
        database name in :ref:`~/.radis <label_lbl_config_file>`
    crop: int
        if > 0, cutoff entries larger than that
    """
    warnings.warn(DeprecationWarning("config format changed to JSON in radis 0.9.29"))
    entries = getDatabankEntries_configformat(dbname)
    print(dbname, "\n-------")
    for k, v in entries.items():
        # Add extra arguments
        args = []
        if k == "levelszpe":
            args.append("cm-1")
        v = "{0}".format(v)
        if len(v) > crop and crop > 0:
            v = v[:crop] + "..."
        # Print item
        print(k, ":", v, *args)


def printDatabankList_configformat():
    """Print all databanks available in ~/.radis"""
    warnings.warn(DeprecationWarning("config format changed to JSON in radis 0.9.29"))
    try:
        print(
            "Databanks in {0}: ".format(CONFIG_PATH_OLD),
            ",".join(getDatabankList_configformat()),
        )
        for dbname in getDatabankList_configformat():
            print("\n")
            printDatabankEntries_configformat(dbname)
    except FileNotFoundError:
        print("No config file {0}".format(CONFIG_PATH_OLD))
        # it's okay


def getDatabankEntries(dbname, get_extra_keys=[], configpath=CONFIG_PATH_JSON):
    r"""Read :ref:`~/radis.json <label_lbl_config_file>` config file and returns a dictionary of entries.

    Parameters
    ----------
    dbname: str
        database name in :ref:`~/radis.json <label_lbl_config_file>` config file
    get_extra_keys: list
        read additional parameters on top of the usual Databank format keys :

    Notes
    -----
    Databank format::

        "MY-HITEMP-CO2": {              # your databank name: use this in calc_spectrum()
                                        # or SpectrumFactory.load_databank()
        "info": "HITEMP 2010 databank"  # whatever you want
        "path": [                       # no "", multipath allowed
            "D:\Databases\HITEMP-CO2\hitemp_07",
            "D:\Databases\HITEMP-CO2\hitemp_08",
            "D:\Databases\HITEMP-CO2\hitemp_09"
        ],
        "format": "hitemp"              # 'hitran' (HITRAN / HITEMP), 'cdsd-hitemp', 'cdsd-4000'
                                        # Databank text file format. More info in
                                        # SpectrumFactory.load_databank function.

        # Optional:

        "parfunc":                      # path or 'USE_HAPI'
                                        # path to tabulated partition functions. If
                                        # `USE_HAPI`, then HAPI (HITRAN Python
                                        interface) [1]_ is used to retrieve them (valid
                                        if your databank is HITRAN data). HAPI
                                        is embedded into RADIS. Check the version.

        "parfuncfmt":                   # 'cdsd'
                                        # format to read tabulated partition function
                                        # file. If `USE_HAPI` is given as `parfunc`
                                        # parameter then this line should not be used.

        "levels_iso1":                  # path to energy levels (needed for non-eq)
                                        # calculations.
        "levels_iso2":                  # etc
        "levels_iso4":                  # etc

        "levelsfmt":                    # 'cdsd'
        }                               # how to read the previous file.



    References
    ----------
    .. [1] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_


    """
    # Loading `~/radis.json` file
    _config = get_user_config(configpath)

    # Make sure it looks like a databank (path and format are given)
    try:
        # Getting the databases
        config = _config["database"]
        config[dbname]["path"]
        config[dbname]["format"]
    except KeyError:
        msg = (
            "{1}\nDBFORMAT\n{0}\n".format(DBFORMAT, dbname)
            + "No databank named {0} in `{1}`. ".format(dbname, configpath)
            + "Available databanks: {0}. ".format(
                getDatabankList(configpath=configpath)
            )
            + "See databank format above. More information in "
            + "https://radis.readthedocs.io/en/latest/lbl/lbl.html#configuration-file"
        )
        if dbname in ["hitran", "hitemp", "exomol"]:
            msg += f"\nIf querying {dbname.upper()} you may want to be using SpectrumFactory.fetch_databank({dbname}) instead of SpectrumFactory.load_databank"
        raise DatabankNotFound(msg)

    entries = config[dbname]

    # Path correction for all tests
    if type(entries["path"]) == str:
        entries["path"] = [entries["path"]]

    # Merge all isotope-dependant levels into one dict
    iso_list = [k for k in entries if k.startswith("levels_iso")]
    if len(iso_list) > 0:
        levels = {}
        for k in iso_list:
            iso = float(k[10:])
            levels[iso] = entries.pop(k)
        entries["levels"] = levels
    else:
        if "levels" in entries:
            raise SyntaxError(
                "in {0}: `levels` replaced with ".format(configpath)
                + "`levels_iso#` where # is the isotope number"
            )

    return entries


def getDatabankList(configpath=CONFIG_PATH_JSON):
    """Get all databanks available in :ref:`~/radis.json"""

    _config = get_user_config(configpath)

    # If `_config` is empty dictionary then return empty list
    if len(_config) == 0:
        return []

    # Getting the databases
    config = _config["database"]

    # Get databank path and format
    validdb = []
    for dbname in config:
        try:
            config[dbname]["path"]
            config[dbname]["format"]
        except configparser.NoSectionError:
            # not a db
            continue
        except configparser.NoOptionError:
            # not a db
            continue
        # looks like a db. add to db
        validdb.append(dbname)
    return validdb


def addDatabankEntries(dbname, dict_entries, verbose=True, configpath=CONFIG_PATH_JSON):
    r"""Add database dbname with entries from dict_entries.

    If database already exists in :ref:`~/radis.json <label_lbl_config_file>`, raises an error

    Examples
    --------
    ::

        addDatabankEntries("HITEMP2010-CO2",
            {
                "info": "HITEMP2020 CO2 lines with TIPS-2017 for partition functions (equilibrium) and RADIS for rovibrational energies (nonequilibrium) ",
                "path": ["PATH/TO/HITEMP/CO2/*.par"],
                "format": "hitran",
                "parfuncfmt": "hapi",
                "levelsfmt": "radis",
            })
    """

    # Get ~/radis.json if exists, else create it
    dbnames = getDatabankList(configpath)

    # Check database doesnt exist
    if dbname in dbnames:
        raise DatabaseAlreadyExists(
            f"Database {dbname} already exists in {configpath}. Cant add it"
        )

    # Loading `~/radis.json`
    with open(configpath, "r") as json_file:
        _config = json.load(json_file)

    try:
        # Accessing `database` key in file
        config = _config["database"]
    except KeyError:
        # Creating `database` key
        _config["database"] = {}
        config = _config["database"]

    # Adding entries in `config[dbname]`
    config[dbname] = {}

    if "info" in dict_entries:  # optional
        config[dbname]["info"] = dict_entries.pop("info")

    # ... Parse paths correctly
    if isinstance(dict_entries["path"], str):
        config[dbname]["path"] = dict_entries.pop("path")
    else:  # list
        # If list has a single path
        if len(dict_entries["path"]) == 1:
            config[dbname]["path"] = dict_entries.pop("path")[0]
        else:
            config[dbname]["path"] = dict_entries.pop("path")

    config[dbname]["format"] = dict_entries.pop("format")
    config[dbname]["parfuncfmt"] = dict_entries.pop("parfuncfmt")

    # Optional:
    # ... partition functions:
    if "parfunc" in dict_entries:
        config[dbname]["parfunc"] = dict_entries.pop("parfunc")

    # ... Split all isotopes in separate keys
    levels_dict = dict_entries.pop("levels", {})
    for iso, levels_iso in levels_dict.items():
        config[dbname]["levels_iso{0}".format(iso)] = levels_iso

    if "levelsfmt" in dict_entries:
        config[dbname]["levelsfmt"] = dict_entries.pop("levelsfmt")

    # Allow some informative parameters
    # ... "download_url", "download_date", "wavenumber_min", "wavenumber_max" generated by fetch_hitemp
    for entry in ["download_url", "download_date", "wavenumber_min", "wavenumber_max"]:
        if entry in dict_entries:
            config[dbname][entry] = dict_entries.pop(entry)

    # Check nothing is left
    if dict_entries != {}:
        raise ValueError("Unexpected keys: {0}".format(list(dict_entries.keys())))

    # Write to `~/radis.json`
    with open(configpath, "w") as configfile:
        json.dump(_config, configfile, indent=3)
        configfile.close()

    if verbose:
        print("Added {0} database in {1}".format(dbname, configpath))

    return


def printDatabankEntries(dbname, crop=200, configpath=CONFIG_PATH_JSON):
    r"""Print databank info.

    Parameters
    ----------
    dbname: str
        database name in :ref:`~/radis.json <label_lbl_config_file>`
    crop: int
        if > 0, cutoff entries larger than that
    """
    entries = getDatabankEntries(dbname, configpath=configpath)
    print("'{0}':".format(dbname))
    print(entries, "\n")


def printDatabankList(configpath=CONFIG_PATH_JSON):
    """Print all databanks available in ~/radis.json"""
    try:
        print(
            "Databanks in {0}: ".format(configpath),
            ",".join(getDatabankList(configpath=configpath)),
        )
        for dbname in getDatabankList(configpath=configpath):
            print("\n")
            printDatabankEntries(dbname, configpath=configpath)
    except FileNotFoundError:
        print("No config file {0}".format(configpath))
        # it's okay


# %% Test

if __name__ == "__main__":

    # from radis.test.misc.test_config import _run_testcases

    # _run_testcases(verbose=True)

    a = get_config()

    print(a.keys())
