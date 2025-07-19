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

# OH 
# --
OH_X_iso1 = ElectronicState(
    "OH",
    isotope=1, 
    state="X",
    term_symbol="2PI",  # Ground state
    g_e=4,    # Electronic degeneracy 2(2S+1) = 2(2*1/2+1) = 4
    spectroscopic_constants="default",
    spectroscopic_constants_type="dunham", 
    vmax=10,  # max vibrational level
    vmax_morse=20,  # max level for Morse potential
    Ediss=35593  # Dissociation energy in cm-1
)
"""OH first isotope (16O-1H), electronic ground state X²Π"""

OH_A_iso1 = ElectronicState(
    "OH", 
    isotope=1,
    state="A",
    term_symbol="2SIGMA+",  # First excited state
    g_e=2,    # Electronic degeneracy 2(2S+1) = 2(2*1/2+1) = 2
    spectroscopic_constants="default",
    spectroscopic_constants_type="dunham",
    Te=32402.38,  # Term energy in cm-1 (from ExoMol MoLLIST)
    Ediss=35593,  # Same dissociation energy as ground state
)
"""OH first isotope (16O-1H), first excited state A²Σ⁺"""

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
    "OH": {
        1: {"X": OH_X_iso1, "A": OH_A_iso1},
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
:py:class:`~radis.db.classes.ElectronicState`

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
    Get rovibrational energies using the default :ref:`spectroscopic constants <label_db_spectroscopic_constants>`:
    ::

        from radis import getMolecule

        # Here we get the energy of the v=6, J=3 level of the 2nd isotope of CO:

        CO = getMolecule("CO", 2, "X")
        print(CO.Erovib(6, 3))

    .. minigallery:: radis.getMolecule

    See Also
    --------
    :py:data:`~radis.db.molecules.Molecules`,
    :py:class:`~radis.db.classes.ElectronicState`

    """

    # molecule
    try:
        mol = Molecules[molecule]
    except KeyError as err:
        raise KeyError(
            f"{molecule} is not defined in molecules with built-in "
            + f"spectroscopic constants. Choose one of: {list(Molecules.keys())}"
        ) from err

    # isotope
    if isotope is None:
        if len(list(mol.keys())) != 1:
            raise ValueError(f"Please precise which isotope among: {list(mol.keys())}")
        isotope = list(mol.keys())[0]
    try:
        iso = mol[isotope]
    except KeyError as err:
        raise KeyError(
            f"Isotope {isotope} is not defined for molecule {molecule}. Choose one of: {list(mol.keys())}"
        ) from err

    # electronic state
    if electronic_state is None:
        if len(list(iso.keys())) != 1:
            raise ValueError(
                f"Please choose which electronic state among: {list(iso.keys())}"
            )
        electronic_state = list(iso.keys())[0]
    try:
        state = iso[electronic_state]
    except KeyError as err:
        raise KeyError(
            f"{electronic_state} is not defined for molecule {molecule}(iso={isotope}). Choose one of: {list(mol.keys())}"
        ) from err

    # print name
    if verbose >= 2:
        print(f"Found {state.get_fullname()} in RADIS database")

    # Return
    return state


# %% Test


if __name__ == "__main__":

    from radis.test.db.test_molecules import _run_testcases

    print(("Testing molecules.py", _run_testcases()))
