# -*- coding: utf-8 -*-
"""Registry of all HITRAN isotopes per molecule.

Derived from HITRAN 2020 [1]_ and ``isotope_name_dict`` in
:py:mod:`radis.db.molparam`.

Used to know which per-isotope database files to expect when
``isotope='all'`` is requested, without querying the HITRAN server
each time.

References
----------
.. [1] `HITRAN molecule metadata <https://hitran.org/docs/molec-meta/>`__

See Also
--------
:py:data:`radis.db.molparam.isotope_name_dict`
"""

from collections import defaultdict

from radis.db.molparam import isotope_name_dict


def _build_hitran_isotopes():
    """Build HITRAN_ISOTOPES from isotope_name_dict.

    Returns a dict  ``{molecule_name: [iso1, iso2, ...]}``
    sorted by isotope number.
    """
    from radis.db.classes import get_molecule

    mol_isos = defaultdict(list)
    for mol_id, iso_num in isotope_name_dict:
        mol_name = get_molecule(mol_id)
        mol_isos[mol_name].append(iso_num)
    # Sort isotope lists
    return {mol: sorted(isos) for mol, isos in mol_isos.items()}


HITRAN_ISOTOPES = _build_hitran_isotopes()
"""dict: Mapping of HITRAN molecule name -> list of valid isotope numbers.

Example::

    >>> from radis.db.hitran_isotopes import HITRAN_ISOTOPES
    >>> HITRAN_ISOTOPES["CO"]
    [1, 2, 3, 4, 5, 6]
    >>> HITRAN_ISOTOPES["CO2"]
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
"""


def get_hitran_isotopes(molecule):
    """Return the list of valid HITRAN isotope numbers for *molecule*.

    Parameters
    ----------
    molecule : str
        HITRAN molecule name, e.g. ``"CO"``, ``"CO2"``.

    Returns
    -------
    list of int

    Raises
    ------
    KeyError
        If *molecule* is not found in the HITRAN registry.
    """
    try:
        return HITRAN_ISOTOPES[molecule]
    except KeyError:
        raise KeyError(
            f"Molecule '{molecule}' not found in HITRAN isotope registry. "
            f"Available molecules: {sorted(HITRAN_ISOTOPES.keys())}"
        )
