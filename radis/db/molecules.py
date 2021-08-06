# -*- coding: utf-8 -*-
"""

Predefined molecules using RADIS built-in
:ref:`spectroscopic constants <label_db_spectroscopic_constants>`

See :py:data:`~radis.db.molecules.Molecules`

-------------------------------------------------------------------------------


"""

from radis.db.classes import ElectronicState
from radis.phys.convert import eV2cm

# %% Define some commonly used molecules


# CO
# --

# Define with default diatomic constants
CO_X_iso1 = ElectronicState(
    "CO",
    isotope=1,
    state="X",
    term_symbol="1Σ+",
    spectroscopic_constants="default",
    spectroscopic_constants_type="dunham",
    vmax=17,  # max level for Dunham's expansion
    vmax_morse=48,
    Ediss=eV2cm(11.16),
)
"""CO first isotope (:math:`^{16}O^{12}C`), electronic ground state"""
CO_X_iso2 = ElectronicState(
    "CO",
    isotope=2,
    state="X",
    term_symbol="1Σ+",
    spectroscopic_constants="default",
    spectroscopic_constants_type="dunham",
    vmax=17,  # max level for Dunham's expansion
    vmax_morse=48,
    Ediss=eV2cm(11.16),
)
"""CO 2nd isotope (:math:`^{16}O^{13}C`), electronic ground state"""
CO_X_iso3 = ElectronicState(
    "CO",
    isotope=3,
    state="X",
    term_symbol="1Σ+",
    spectroscopic_constants="default",
    spectroscopic_constants_type="dunham",
    vmax=17,  # max level for Dunham's expansion
    vmax_morse=48,
    Ediss=eV2cm(11.16),
)
"""CO 3rd isotope (:math:`^{18}O^{12}C`) electronic ground state"""

# CO2
# ---

CO2_X_626 = ElectronicState(
    "CO2",
    isotope=1,
    state="X",
    term_symbol="1Σu+",
    spectroscopic_constants="default",
    spectroscopic_constants_type="herzberg",
    Ediss=44600,
)
"""CO2 1st isotope (:math:`^{16}O^{12}C^{16}O`), electronic ground state"""
CO2_X_636 = ElectronicState(
    "CO2",
    isotope=2,
    state="X",
    term_symbol="1Σu+",
    spectroscopic_constants="default",
    spectroscopic_constants_type="herzberg",
    Ediss=44600,
)
"""CO2 2nd isotope (:math:`^{16}O^{13}C^{16}O`), electronic ground state"""
CO2_X_628 = ElectronicState(
    "CO2",
    isotope=3,
    state="X",
    term_symbol="1Σu+",
    spectroscopic_constants="default",
    spectroscopic_constants_type="herzberg",
    Ediss=44600,
)
"""CO2 3rd isotope (:math:`^{16}O^{12}C^{18}O`), electronic ground state"""
CO2_X_627 = ElectronicState(
    "CO2",
    isotope=4,
    state="X",
    term_symbol="1Σu+",
    spectroscopic_constants="default",
    spectroscopic_constants_type="herzberg",
    Ediss=44600,
)
"""CO2 4th isotope (:math:`^{16}O^{12}C^{17}O`), electronic ground state"""


# %% Dictionary of predefined molecules
# -----------

#  Molecule  Isotope  ElecState
#         :      :       :
Molecules = {
    "CO": {
        1: {"X": CO_X_iso1},
        2: {"X": CO_X_iso2},
        3: {"X": CO_X_iso3},
    },
    "CO2": {
        1: {"X": CO2_X_626},
        2: {"X": CO2_X_636},
        3: {"X": CO2_X_628},
        4: {"X": CO2_X_627},
    },
}
"""dict: list of Electronic states whose energy levels can be calculated with RADIS
built-in :ref:`spectroscopic constants <label_db_spectroscopic_constants>`.

For references refer to the definition of each molecule:
CO:

- :py:data:`~radis.db.molecules.CO_X_iso1`
- :py:data:`~radis.db.molecules.CO_X_iso2`
- :py:data:`~radis.db.molecules.CO_X_iso3`

CO2:

- :py:data:`~radis.db.molecules.CO2_X_626`
- :py:data:`~radis.db.molecules.CO2_X_636`
- :py:data:`~radis.db.molecules.CO2_X_628`
- :py:data:`~radis.db.molecules.CO2_X_627`

See Also
--------

:py:func:`~radis.db.molecules.getMolecule`

"""


def getMolecule(
    molecule, isotope=None, electronic_state=None, verbose=True
) -> ElectronicState:
    """Get an :py:class:`~radis.db.classes.ElectronicState` object in the
    RADIS :py:data:`~radis.db.molecules.Molecules` list, which use the defaults
    :ref:`spectroscopic constants <label_db_spectroscopic_constants>`.

    Parameters
    ----------
    molecule: str
        molecule name
    isotope: int, or ``None``
        isotope number. if None, only one isotope must exist in database. Else,
        an error is raised
    electronic_state: str
        if None, only one electronic state must exist in database. Else, an error
        is raised
    verbose: boolean
        if ``True``, print which electronic state we got

    Returns
    -------
    ElectronicState: an :py:class:`~radis.db.classes.ElectronicState` object.

    Examples
    --------

    .. minigallery:: radis.db.molecules.getMolecule
        :add-heading:

    See Also
    --------
    :py:data:`~radis.db.molecules.Molecules`

    """

    # molecule
    try:
        mol = Molecules[molecule]
    except KeyError as err:
        raise KeyError(
            "{0} is not defined in molecules with built-in ".format(molecule)
            + "spectroscopic constants. Choose one of: {0}".format(
                list(Molecules.keys())
            )
        ) from err

    # isotope
    if isotope is None:
        if len(list(mol.keys())) != 1:
            raise ValueError(
                "Please precise which isotope among: {0}".format(list(mol.keys()))
            )
        isotope = list(mol.keys())[0]
    try:
        iso = mol[isotope]
    except KeyError as err:
        raise KeyError(
            "Isotope {0} is not defined for molecule {1}. Choose one of: {2}".format(
                isotope, molecule, list(mol.keys())
            )
        ) from err

    # electronic state
    if electronic_state is None:
        if len(list(iso.keys())) != 1:
            raise ValueError(
                "Please choose which electronic state among: {0}".format(
                    list(iso.keys())
                )
            )
        electronic_state = list(iso.keys())[0]
    try:
        state = iso[electronic_state]
    except KeyError as err:
        raise KeyError(
            "{0} is not defined for molecule {1}(iso={2}). Choose one of: {3}".format(
                electronic_state, molecule, isotope, list(mol.keys())
            )
        ) from err

    # print name
    if verbose >= 2:
        print(r"Found {0} in RADIS database".format(state.get_fullname()))

    # Return
    return state


# %% Test


if __name__ == "__main__":

    from radis.test.db.test_molecules import _run_testcases

    print(("Testing molecules.py", _run_testcases()))
