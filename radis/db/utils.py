# -*- coding: utf-8 -*-
"""Created on Thu May 28 14:47:36 2015.

@author: Erwan
"""


import json
import os
from collections import OrderedDict
from os.path import abspath

from radis.misc.basics import is_number


def getFile(*relpath):
    r"""Converts the relative path of a database file in a the full path.
    Used by processing script not to worry about where the database is stored

    Examples of use::

        radis.db.getFile('CN','CN_Violet_vib.dat')
        radis.db.getFile('CN\CN_Violet_vib.dat')
        radis.db.getFile('CN\\CN_Violet_vib.dat')

    """

    #    return os.path.join(os.path.dirname(__file__), *relpath)
    from radis.misc.utils import getProjectRoot

    return os.path.join(getProjectRoot(), "db", *relpath)


def check_molecule_data_structure(fname, verbose=True):
    """Check that ``fname`` has a valid JSON structure for molecular data.

    Parameters
    ----------

    fname: str
        molecular data JSON file


    Notes
    -----

    Order in the json doesnt matter, however, the order of the
    electronic_level_names matters: the index of all levels should
    match the `index` key of these states.
    """

    with open(fname) as f:
        try:
            db = json.load(f)  # , object_pairs_hook=OrderedDict)
        except json.JSONDecodeError as err:
            raise json.JSONDecodeError(
                f"Error reading '{fname}' (line {err.lineno} col {err.colno}): \n{err.msg}",
                err.doc,
                err.pos,
            ) from err

    for molecule, db_molec in db.items():
        # ... Check number of isotopes is correct
        isotope_names = db_molec["isotopes_names"]
        isotopes = db_molec["isotopes"]
        if len(isotopes) != len(isotope_names):
            raise ValueError(
                f"In molecule {molecule}: isotope names "
                + f"({isotope_names}) dont match the number of isotopes ({list(isotopes.keys())})"
            )

        # ... Check number of electronic states is correct
        for isotope, db_iso in db_molec["isotopes"].items():
            elec_states_names = db_iso["electronic_levels_names"]
            elec_states = db_iso["electronic_level"]
            if len(elec_states_names) != len(elec_states):
                raise ValueError(
                    f"In molecule {molecule}, isotope {isotope}: electronic "
                    + f"levels names ({elec_states_names}) dont match the number of levels ({list(elec_states.keys())})"
                )

            # ... Check they are properly ordered
            for state, db_state in elec_states.items():
                if elec_states_names.index(db_state["name"]) + 1 != db_state["index"]:
                    # ... + 1 because Python index start at 0 and FORTRAN at 1 unless allocated with 0:N
                    raise ValueError(
                        f"In molecule {molecule}, isotope {isotope}: index of electronic "
                        + f"state {state} ({db_state['name']}): {db_state['index']} does not match the list of states: {elec_states_names}. "
                    )

    if verbose:
        print(f"Structure of {fname} looks correct")


def parse_doi(rovib_coeffs):
    if "#doi" in rovib_coeffs:
        doi = rovib_coeffs["#doi"]
    elif "#DOI" in rovib_coeffs:
        doi = rovib_coeffs["#DOI"]
    else:
        doi = None

    return doi


def _get_rovib_coefficients(
    molecule, isotope, electronic_state, jsonfile, remove_trailing_cm1=True
):
    """Returns all rovib coefficients for ``molecule``, ``isotope``,
    ``electronic_state`` by parsing a JSON file of molecule data.

    Parameters
    ----------

    molecule: str
        molecule name

    isotope: int
        isotope number

    electronic_state: str
        electronic state name

    jsonfile: str, or ``default``
        path to json file.

    Other Parameters
    ----------------

    remove_trailing_cm1: bool
        if ``True``, remove the trailing '_cm-1' in the spectroscopic coefficient name,
        if defined. Ex::

            ('we1_cm-1', 1333.93)
            >>> ('we1', 1333.93)

    See Also
    --------

    :py:func:`~radis.db.utils.get_dunham_coefficients`,
    :py:func:`~radis.db.utils.get_herzberg_coefficients`,
    """

    check_molecule_data_structure(jsonfile, verbose=False)

    with open(jsonfile) as f:
        try:
            db = json.load(f, object_pairs_hook=OrderedDict)
        except json.JSONDecodeError as err:
            raise json.JSONDecodeError(
                f"Error reading '{jsonfile}' (line {err.lineno} col {err.colno}): \n{err.msg}",
                err.doc,
                err.pos,
            ) from err

    # Get Dunham coefficients in 001 state (X)
    elec_state_names = db[molecule]["isotopes"][str(isotope)]["electronic_levels_names"]

    try:
        elec_state_index = elec_state_names.index(electronic_state)
    except ValueError:
        raise ValueError(
            f"{electronic_state} not in the electronic state list for {molecule}(iso={isotope}): {elec_state_names}"
        )
    else:
        elec_state_index = f"{elec_state_index + 1:03d}"  # 1 based index

    rovib_coeffs = db[molecule]["isotopes"][str(isotope)]["electronic_level"][
        elec_state_index
    ]

    if remove_trailing_cm1:

        def remove_cm1(coef):
            """Remove the trailing '_cm-1' in the spectroscopic coefficient
            name, if defined."""
            if coef.endswith("_cm-1"):
                coef = coef[:-5]
            return coef

        rovib_coeffs = {remove_cm1(k): v for (k, v) in rovib_coeffs.items()}

    return rovib_coeffs


def get_default_jsonfile(molecule):
    r"""Return full path of default jsonfile for spectroscopic constants for a
    molecule.

    These are stored in:

            radis\db\[molecule]\*.json

    and defined in radis/config.json
    """

    from radis import config

    name = config["spectroscopic_constants"][molecule]

    return abspath(getFile(f"{molecule}/{name}"))


def get_dunham_coefficients(
    molecule, isotope, electronic_state, jsonfile="default", return_doi=False
):
    r"""Returns Dunham coefficients ``Yij`` for ``molecule``, ``isotope``, ``electronic_state``
    by parsing a JSON file of molecule data.

    Dunham coefficients are identified as starting with ``Y`` (i.e. ``alpha_e``
    won't be recognized)

    Parameters
    ----------
    molecule: str
        molecule name
    isotope: int
        isotope number
    electronic_state: str
        electronic state name
    jsonfile: str, or ``default``
        path to json file. If ``default``, the default JSON file
        defined in radis/config.json is used. See :py:func:`~radis.db.utils.get_default_jsonfile`

    Other Parameters
    ----------------
    return_doi: bool
        if True, returns the value of key ``"#doi"`` or ``"#DOI"`` if it exists,
        and None if it doesn't.

    Returns
    -------
    dict: coeffs
    str: doi, if ``return_doi``

    See Also
    --------

    :py:func:`~radis.db.utils.get_herzberg_coefficients`

    """

    if jsonfile == "default":
        jsonfile = get_default_jsonfile(molecule)

    rovib_coeffs = _get_rovib_coefficients(
        molecule, isotope, electronic_state, jsonfile=jsonfile
    )

    # Only get Dunham coeffs, i.e, these that start with Y
    dunham_coeffs = {k: v for (k, v) in rovib_coeffs.items() if k.startswith("Y")}

    if len(dunham_coeffs) == 0:
        raise ValueError(
            f"No Dunham coefficients (Yij) found for {molecule} {electronic_state}(iso={isotope}) in {jsonfile}"
            + f"\n\nGot: {rovib_coeffs.keys()}"
        )

    if return_doi:
        return dunham_coeffs, parse_doi(rovib_coeffs)
    else:
        return dunham_coeffs


def get_herzberg_coefficients(
    molecule, isotope, electronic_state, jsonfile="default", return_doi=False
):
    r"""Returns spectroscopic coefficients with Herzberg conventions for
    ``molecule``, ``isotope``, ``electronic_state``  by parsing a JSON file of molecule data.

    Herzberg coefficients are the usual:

    .. math::

        \\omega_e, \\alpha_e, B_e, D_e, etc.

    The extensive list of parameters considered as Herzberg coefficients is found in
    :py:data:`~radis.db.conventions.herzberg_coefficients`

    Parameters
    ----------
    molecule: str
        molecule name
    isotope: int
        isotope number
    electronic_state: str
        electronic state name
    jsonfile: str, or ``default``
        path to json file. If ``default``, the ``molecules_data`` JSON file
        in the ``radis.db`` database is used. See :py:func:`~radis.db.utils.get_default_jsonfile`

    Other Parameters
    ----------------
    return_doi: bool
    if True, returns the value of key ``"#doi"`` or ``"#DOI"`` if it exists,
    and None if it doesn't.

    Returns
    -------
    dict: coeffs
    str: doi, if ``return_doi``

    See Also
    --------
    :py:data:`~radis.db.conventions.herzberg_coefficients`,
    :py:func:`~radis.db.utils.get_dunham_coefficients`

    """

    from radis.db.conventions import herzberg_coefficients

    if jsonfile == "default":
        jsonfile = get_default_jsonfile(molecule)

    rovib_coeffs = _get_rovib_coefficients(
        molecule, isotope, electronic_state, jsonfile=jsonfile
    )

    #    # Only get Herzberg coeffs, i.e, these that are defined in :py:data:`~radis.db.conventions.herzberg_coefficients`

    herzberg_coeffs = {
        k: v
        for (k, v) in rovib_coeffs.items()
        if ignore_trailing_number(k) in herzberg_coefficients
    }

    if len(herzberg_coeffs) == 0:
        raise ValueError(
            f"No Herzberg spectroscopic coefficients (we, Be, etc.) found for {molecule} {electronic_state}(iso={isotope}) in {jsonfile}"
            + f"\n\nGot: {rovib_coeffs.keys()}"
        )

    if return_doi:
        return herzberg_coeffs, parse_doi(rovib_coeffs)
    else:
        return herzberg_coeffs


def ignore_trailing_number(coef):
    """Used so that ``wexe1`` matches ``wexe`` as a well defined Herzberg
    coefficient."""

    if is_number(coef[-1]):
        coef = coef[:-1]
    return coef


def compare(databank, compare_with):
    """Simply to make databank case-insensitive!"""
    if isinstance(compare_with, str) and (
        (isinstance(databank, str) and databank.casefold() == compare_with.casefold())
        or (
            isinstance(databank, tuple)
            and databank[0].casefold() == compare_with.casefold()
        )
    ):
        return True
    if isinstance(databank, str) and databank.casefold() in compare_with:
        return True
    return False
