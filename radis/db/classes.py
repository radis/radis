# -*- coding: utf-8 -*-
"""Created on Tue Jul 18 16:15:03 2017.

@author: erwan

Summary
-------

Define :py:class:`~radis.db.classes.Molecule`, :py:class:`~radis.db.classes.Isotope`
and :py:class:`~radis.db.classes.ElectronicState` classes

ElectronicState has:

- Molecular parameters (read from HITRAN molparams.txt)
- Vibrational and rotational energies

Calculate Vibrational and Rotational energies for some database molecules

Molecules
-----------

Shortcut to some molecules ground states are defined at the end of the file
in the dictionary `Molecules`, with the following structure:

{molecule_name: {isotope_number: {electronic_state_name: ElectronicState object}}}

Import with::

    from radis.db.molecules import Molecules
    CO2_X = Molecules['CO2'][1]['X']            # 1 for first isotope

Or::

    from radis import getMolecule                # directly
    CO2_2 = getMolecule('CO2', 1, 'X')         # better KeyError messages


Routine Listing
---------------

- :func:`~radis.db.classes.get_molecule`
- :func:`~radis.db.classes.get_molecule_identifier`

-------------------------------------------------------------------------------
"""

import re
from os.path import dirname, exists, join

from radis.db.conventions import get_convention
from radis.db.utils import (
    get_default_jsonfile,
    get_dunham_coefficients,
    get_herzberg_coefficients,
    parse_doi,
)
from radis.levels.dunham import EvJ, Fv, Gv

# %% Hitran groups and classes
# As defined in Rothman et al, "The HITRAN 2004 molecular spectroscopic database"
# Tables 3 and 4
# Groups and classes are needed to compute energy in nonequilibrium mode.

# Groups (define local quanta)

HITRAN_GROUP1 = [
    "H2O",
    "O3",
    "SO2",
    "NO2",
    "HNO3",
    "H2CO",
    "HOCl",
    "H2O2",
    "COF2",
    "H2S",
    "HO2",
    "HCOOH",
    "ClONO2",
    "HOBr",
    "C2H4",
]
"""str: asymmetric rotors"""

HITRAN_GROUP2 = [
    "CO2",
    "N2O",
    "CO",
    "HF",
    "HCl",
    "HBr",
    "HI",
    "OCS",
    "N2",
    "HCN",
    "C2H2",
    "NO+",
]
"""str: diatomic and linear molecules"""

HITRAN_GROUP3 = ["SF6", "CH4"]
"""str: Spherical rotors"""

HITRAN_GROUP4 = ["CH3D", "CH3Cl", "C2H6", "NH3", "PH3", "CH3OH"]
"""str: symmetric rotors"""

HITRAN_GROUP5 = ["O2"]
"""str: Triplet-Sigma ground electronic states"""

HITRAN_GROUP6 = ["NO", "OH", "ClO"]
"""str: Doublet-Pi ground electronic states"""

# Classes (define global quanta)

HITRAN_CLASS1 = ["CO", "HF", "HCl", "HBr", "HI", "N2", "NO+"]
"""str: Diatomic molecules with ? """

HITRAN_CLASS2 = ["O2"]
"""str: Diatomic molecules with different electronic levels"""

HITRAN_CLASS3 = ["NO", "OH", "ClO"]
"""str: Diatomic molecules with doublet-Pi electronic state"""

HITRAN_CLASS4 = ["N2O", "OCS", "HCN"]
"""str: Linear triatomic"""

HITRAN_CLASS5 = ["CO2"]
"""str: Linear triatomic with large Fermi resonance"""

HITRAN_CLASS6 = ["H2O", "O3", "SO2", "NO2", "HOCl", "H2S", "HO2", "HOBr"]
"""str: Non-linear triatomic"""

HITRAN_CLASS7 = ["C2H2"]
"""str: Linear tetratomic"""

HITRAN_CLASS8 = ["NH3", "PH3"]
"""str: Pyramidal tetratomic"""

HITRAN_CLASS9 = ["H2CO", "H2O2", "COF2"]
"""str: Non-linear tetratomic"""

HITRAN_CLASS10 = [
    "CH4",
    "CH3Cl",
    "C2H6",
    "HNO3",
    "SF6",
    "HCOOH",
    "ClONO2",
    "C2H4",
    "CH3OH",
]
"""str: Pentatomic or greater polyatomic"""
# 0.9.22: Removed CH3D has no HITRAN identifier


# %% HITRAN ids

trans = {
    "1": "H2O",
    "2": "CO2",
    "3": "O3",
    "4": "N2O",
    "5": "CO",
    "6": "CH4",
    "7": "O2",
    "8": "NO",
    "9": "SO2",
    "10": "NO2",
    "11": "NH3",
    "12": "HNO3",
    "13": "OH",
    "14": "HF",
    "15": "HCl",
    "16": "HBr",
    "17": "HI",
    "18": "ClO",
    "19": "OCS",
    "20": "H2CO",
    "21": "HOCl",
    "22": "N2",
    "23": "HCN",
    "24": "CH3Cl",
    "25": "H2O2",
    "26": "C2H2",
    "27": "C2H6",
    "28": "PH3",
    "29": "COF2",
    "30": "SF6",
    "31": "H2S",
    "32": "HCOOH",
    "33": "HO2",
    "34": "O",
    "35": "ClONO2",
    "36": "NO+",
    "37": "HOBr",
    "38": "C2H4",
    "39": "CH3OH",
    "40": "CH3Br",
    "41": "CH3CN",
    "42": "CF4",
    "43": "C4H2",
    "44": "HC3N",
    "45": "H2",
    "46": "CS",
    "47": "SO3",
    "48": "C2N2",
    "49": "COCl2",
    "50": "SO",
    "51": "CH3F",
    "52": "GeH4",
    "53": "CS2",
    "54": "CH3I",
    "55": "NF3",
}
HITRAN_MOLECULES = list(trans.values())
""" str: list of [HITRAN-2020]_ molecules. """


def get_molecule_identifier(molecule_name):
    r"""
    For a given input molecular formula, return the corresponding
    :py:data:`~radis.db.classes.HITRAN_MOLECULES` identifier number [1]_.

    Parameters
    ----------
    molecular_formula : str
        The string describing the molecule.

    Returns
    -------
    M: int
        The HITRAN molecular identified number.

    References
    ----------
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__

    Function is from from https://github.com/nzhagen/hitran/blob/master/hitran.py

    """

    # Invert the dictionary.
    trans2 = {v: k for k, v in trans.items()}

    try:
        return int(trans2[molecule_name])
    except KeyError:
        raise NotImplementedError(
            "Molecule '{0}' not supported. Choose one of {1}".format(
                molecule_name, sorted(list(trans2.keys()))
            )
        )


def get_molecule(molecule_id):
    r"""
    For a given input molecular identifier, return the corresponding
    :py:data:`~radis.db.classes.HITRAN_MOLECULES` name [1]_.

    Parameters
    ----------
    molecular_id: str
        Hitran identifier of the molecule.


    References
    ----------
    .. [1] `HITRAN 1996, Rothman et al., 1998 <https://www.sciencedirect.com/science/article/pii/S0022407398000788>`__

    """

    # assert str
    id = "{:d}".format(int(molecule_id))

    try:
        return trans[id]
    except KeyError:
        raise NotImplementedError(
            "Molecule ID '{0}' unknown. Choose one of {1}".format(molecule_id, trans)
        )


# %% ExoMol molecules


from radis.db.molparam import get_extra_molparams

isotope_full_names = get_extra_molparams()["isotope_full_name"]

EXOMOL_ONLY_ISOTOPES_NAMES = {
    (molecule, int(iso_number)): iso_name
    for molecule, iso_dict in isotope_full_names.items()
    for iso_number, iso_name in iso_dict.items()
}
"""All :py:data:`~radis.db.classes.HITRAN_MOLECULES` are also converted to their ExoMol full-name format
in :py:func:`~radis.io.exomol.get_exomol_full_isotope_name`
"""
assert EXOMOL_ONLY_ISOTOPES_NAMES[("FeH", 1)] == "56Fe-1H"
assert EXOMOL_ONLY_ISOTOPES_NAMES[("SiO", 1)] == "28Si-16O"
assert EXOMOL_ONLY_ISOTOPES_NAMES[("CN", 1)] == "12C-14N"


EXOMOL_MOLECULES = [
    "AlCl",
    "AlF",
    "AlH",
    "AlO",
    "AsH3",
    "BeH",
    "C2",
    "C2H2",
    "C2H4",
    "CH",
    "CH3",
    "CH3Cl",
    "CH3F",
    "CH4",
    "CN",
    "CO",
    "CO2",
    "CP",
    "CS",
    "CaF",
    "CaH",
    "CaO",
    "CaOH",
    "CrH",
    "FeH",
    "H2",
    "H2CO",
    "H2CS",
    "H2O",
    "H2O2",
    "H2S",
    "H2_p",
    "H3O_p",
    "H3_p",
    "HCN",
    "HCl",
    "HF",
    "HNO3",
    "HeH_p",
    "KCl",
    "KF",
    "KOH",
    "LaO",
    "LiCl",
    "LiF",
    "LiH",
    "LiH_p",
    "MgF",
    "MgH",
    "MgO",
    "N2",
    "N2O",
    "NH",
    "NH3",
    "NO",
    "NS",
    "NaCl",
    "NaF",
    "NaH",
    "NaO",
    "NaOH",
    "NiH",
    "O2",
    "OCS",
    "OH",
    "OH_p",
    "PF3",
    "PH",
    "PH3",
    "PN",
    "PO",
    "PS",
    "SH",
    "SO2",
    "SO3",
    "ScH",
    "SiH",
    "SiH2",
    "SiH4",
    "SiN",
    "SiO",
    "SiO2",
    "SiS",
    "TiH",
    "TiO",
    "VO",
    "YO",
    "cis-P2H2",
    "trans-P2H2",
]
EXOMOL_ONLY_MOLECULES = sorted(
    list(set([M for M, iso in EXOMOL_ONLY_ISOTOPES_NAMES.keys()]))
)

for M in EXOMOL_ONLY_ISOTOPES_NAMES.keys():
    assert M[0] in EXOMOL_MOLECULES

# %% Molecule class


class Molecule(object):
    r"""Define a new molecule

    Parameters
    ----------
    name: str, int
        molecule name, or HITRAN identifier

    Other Parameters
    ----------------
    verbose: boolean
        more blabla

    See Also
    --------
    :py:class:`~radis.db.classes.Isotope`, :py:class:`~radis.db.classes.ElectronicState`,
    :py:func:`~radis.db.molecules.getMolecule`

    """

    def __init__(
        self, name, verbose=True
    ):  # @dev: no optional kwargs here (final stop)

        self.verbose = verbose

        # Get name and integer id
        if isinstance(name, str):
            self.name = name

            # Get name without parenthesis (without state) for HITRAN identification
            filtername = re.sub(r"[\(\[].*?[\)\]]", "", name)

            try:
                self.id = get_molecule_identifier(filtername)
            except NotImplementedError:  # Not an HITRAN molecule
                self.id = None
        elif type(name) == int:
            self.id = name
            self.name = get_molecule(name)
        else:
            raise ValueError("Wrong name type:", name)


# %% Isotope class


class Isotope(Molecule):
    r"""Create an isotope of a given :py:class:`~radis.db.classes.Molecule`

    Parameters
    ----------
    molecule_name: str, or int
        molecule name or HITRAN identifier
    isotope: int
        isotope identifier, sorted by decreasing abundance (cf HITRAN
        nomenclature)

    Other Parameters
    ----------------
    isotope_name: str
        (optional) isotope name. Default ''
    abundance: float
        isotopologue abundance. Default ``None``

    See Also
    --------
    :py:class:`~radis.db.classes.ElectronicState`, :py:class:`~radis.db.classes.Molecule`,
    :py:func:`~radis.db.molecules.getMolecule`
    """

    def __init__(
        self, molecule_name, isotope, isotope_name="", abundance=None, **kwargs
    ):

        super(Isotope, self).__init__(
            name=molecule_name, **kwargs
        )  # initialize Molecule

        # Check input
        if not isinstance(isotope, int):
            raise TypeError(
                "Wrong format for isotope {0}: expected int, got {1}".format(
                    isotope, type(isotope)
                )
            )

        # Store
        self.iso = isotope
        self.isotope_name = isotope_name
        self.abundance = abundance


# %% ElectronicState class

_term_symbols = {"Σ": "SIG", "Π": "PI", "Δ": "DEL"}


def _format_term_symbol(term_symbol):
    """standardized name of greek symbols."""
    for k, v in _term_symbols.items():
        term_symbol = term_symbol.replace(k, v)
    return term_symbol


def _print_term_symbol(term_symbol):
    """nice print of greek letters."""
    for k, v in _term_symbols.items():
        term_symbol = term_symbol.replace(v, k)
    return term_symbol


class ElectronicState(Isotope):
    r"""Define an electronic state of an :py:class:`~radis.db.classes.Isotope`
    of a :py:class:`~radis.db.classes.Molecule`

    Parameters
    ----------
    molecule_name: str, int
        molecule name, or HITRAN identifier
    isotope: int
        isotope number of the molecule (sorted by abundance on Earth, HITRAN
        convention).
    state: str
        Electronic state name
    term_symbol: str
        Term symbol. Default ''
    g_e: int
        Degeneracy. Default ``None``
    spectroscopic_constants: str, ``'default'``, or dict
        filename of spectroscopic constants under Herzberg or Dunham format.

        Expected in the file:

            Yij: cm-1
                rovibrational coefficients in Dunham convention

        or

            wexe, Be, etc. : cm-1
                rovibrational coefficients in Herzberg convention

        or

            Te: cm-1
                electronic energy. Default ``None`` if not given

        If ``default``, the constants defined in
        :ref:`spectroscopic constants <label_db_spectroscopic_constants>` are used.

        One can also directly give a dictionary directly::

            spectroscopic_constants = {'we':2300, 'wexe':20}

    spectroscopic_constants_type: ``'herzberg'``, ``'dunham'``
        convention for spectroscopic constants. Default ``'herzberg'``
    Erovib: function
        Rovibrational energy of the molecule. Accessible with, typically, Mol.Erovib(v,J)
    Ehaj: function
        To return vibrational energy in harmonic and anharmonic component,
        plus rotational energy. Used to build Treanor distributions
    vmax, vmax_morse: int, or ``None``
        maximum vibrational number (required for partition function calculation)
        to be calculated with Dunham expansion, and Morse potential
        If None, number will be infered from dissociation energy
    Jmax: int, or ``None``
        maximum rotational number (required for partition function calculation)
        If None, number will be infered from dissociation energy
    Ediss: cm-1
        dissociation energy. Required for partition function calculation
        if neither vmax nor Jmax are given
    kwargs: **dict
        forwarded to parent class

    Examples
    --------

    .. minigallery:: radis.db.classes.ElectronicState
        :add-heading:


    See Also
    --------
    :py:class:`~radis.db.classes.Isotope`, :py:class:`~radis.db.classes.Molecule`,
    :py:func:`~radis.db.molecules.getMolecule`

    """

    #'''
    # TODO
    # To allow for more flexibility, ElectronicState should also inherit from some
    # built-in Symmetry / species class, Sigma / Pi / Delta etc... where the energy
    # is calculated appropriately.  (ex: Rigid Rotor, Symmetric Top etc.)

    # (major code rewritting ahead...)
    # That should be part of RADIS extension to non Sigma electronic states.

    def __init__(
        self,
        molecule_name,
        isotope,
        state,
        term_symbol="",
        g_e=None,
        spectroscopic_constants="default",
        spectroscopic_constants_type="herzberg",
        Erovib=None,
        Ehaj=None,
        vmax=None,
        vmax_morse=None,
        Jmax=None,
        Ediss=None,
        **kwargs
    ):

        super(ElectronicState, self).__init__(
            molecule_name=molecule_name, isotope=isotope, **kwargs
        )  # initialize Isotope

        # Initialize electronic state specific parameters
        self.state = state
        self.term_symbol = term_symbol
        self.g_e = g_e

        self.doi = (
            None  #: str: stores DOI of spectroscopic constants if given in database
        )
        self._parse_rovib_constants(
            spectroscopic_constants, spectroscopic_constants_type
        )

        #        self.Erovib = None   #: func: overwritten by Erovib if given, or by default Dunham developments if spectroscopic_constants are given
        self.Ehaj = None  #: func: overwritten by Erovib if given, or by default Dunham developments if spectroscopic_constants are given
        self._assign_E(Erovib, Ehaj)
        self.state = state
        self.vmax = vmax
        self.vmax_morse = vmax_morse
        self.Jmax = Jmax
        self.Ediss = Ediss

    def _parse_rovib_constants(
        self, spectroscopic_constants, spectroscopic_constants_type
    ):
        r"""Parse spectroscopic constants

        Stores :py:attr:`~radis.db.classes.ElectronicState.Te` and
        :py:attr:`~radis.db.classes.ElectronicState.re` as electronic state
        attributes, and the rest under :py:attr:`~radis.db.classes.ElectronicState.rovib_constants`

        Parameters
        ----------
        spectroscopic_constants: str, or ``'default'``, or dict
            filename of spectroscopic constants under Herzberg or Dunham format.

            Expected in the file:

                Yij: cm-1
                    rovibrational coefficients in Dunham convention

            or

                wexe, Be, etc. : cm-1
                    rovibrational coefficients in Herzberg convention

            or

                Te: cm-1
                    electronic energy. Default ``None`` if not given

            If ``default``, the constants defined in
            :ref:`spectroscopic constants <label_db_spectroscopic_constants>` are used.

            One can also directly give a dictionary directly::

                spectroscopic_constants = {'we':2300, 'wexe':20}

        spectroscopic_constants_type: ``'herzberg'``, ``'dunham'``
            convention for spectroscopic constants. Default ``'herzberg'``

        Returns
        -------
        None:
            but constants are stored under :py:attr:`~radis.db.classes.ElectronicState.rovib_constants`,
            and store json file in :py:attr:`~radis.db.classes.ElectronicState.jsonfile`

        """

        if isinstance(spectroscopic_constants, dict):
            rovib_constants = spectroscopic_constants
            jsonfile = None

            doi = parse_doi(rovib_constants)

        else:  # file

            # Get file name
            if spectroscopic_constants == "default":
                jsonfile = get_default_jsonfile(self.name)
            elif exists(spectroscopic_constants):  # absolute path
                jsonfile = spectroscopic_constants
            else:  # assume a json file stored in the default folder
                jsonfile = join(
                    dirname(get_default_jsonfile(self.name)), spectroscopic_constants
                )

            # Parse file
            if spectroscopic_constants_type == "dunham":
                rovib_constants, doi = get_dunham_coefficients(
                    self.name,
                    self.iso,
                    self.get_statename_utf(),
                    jsonfile=jsonfile,
                    return_doi=True,
                )
            elif spectroscopic_constants_type == "herzberg":
                rovib_constants, doi = get_herzberg_coefficients(
                    self.name,
                    self.iso,
                    self.get_statename_utf(),
                    jsonfile=jsonfile,
                    return_doi=True,
                )
            else:
                raise ValueError(
                    "Unexpected spectroscopic constant type: {0}".format(
                        spectroscopic_constants_type
                    )
                )

        # Clean keys
        # In particular, remove trailing '_cm-1' if given in dict or database
        import re

        rovib_constants = {
            re.sub("_cm-1$", "", k): v for (k, v) in rovib_constants.items()
        }

        # Get specific keys
        self.Te = rovib_constants.pop("Te", None)  # default None
        self.re = rovib_constants.pop("re", None)

        # Store
        self.rovib_constants = rovib_constants
        self.jsonfile = jsonfile
        self.doi = doi

    # Default method to calculate energy

    def _assign_E(self, Erovib, Ehaj):
        r"""Finds appropriate Electrorovibrational energy function for this molecule ,
        based on which convention is used for spectroscopic coefficients
        (Dunham, or Herzberg)

        Replaces the :py:meth:`~radis.db.classes.ElectronicState.Erovib` template
        with the appropriate function for this molecule / isotope / electronic state.

        Default methods only work for diatomic molecules. For polyatomic molecules,
        you should give the energy function directly (see how CO2 is defined in
        radis.db.molecules.py)

        Returns
        -------

        None:
            but assigns the ``self.Erovib`` function to the molecule

        """
        # Check input
        if Erovib is None and Ehaj is not None:
            raise ValueError(
                "If giving Ehaj (harmonic and anharmonic component calculation) "
                + "you must also supply Erovib (total rovibrational energy)"
            )

        if Erovib is not None:
            # overwrite energy calculation
            self.Erovib = Erovib
            self.Ehaj = Ehaj
            if self.verbose >= 2:
                print(
                    "{0}: overwritting Energy calculation with {1}".format(
                        self.get_fullname(), Erovib
                    )
                )
        else:
            # Autofind which energy model to use for the given molecule
            c = self.rovib_constants
            if len(c) == 0:
                pass  # keep the default function, that will raise an error on first call.
                if self.verbose >= 2:
                    print("{0}: No rovibconstants found".format(self.get_fullname()))
            else:
                convention = get_convention(c)  # Herzberg or Dunham

                if self.name in HITRAN_CLASS1:
                    if convention == "herzberg":
                        self.Erovib = self._E_Herzberg
                        self.Ehaj = None  # NotImplemented
                        if self.verbose >= 2:
                            print(
                                "{0}: using Herzberg coefficients for Energy calculation".format(
                                    self.get_fullname()
                                )
                            )
                    else:
                        self.Erovib = self._E_Dunham
                        self.Ehaj = None  # NotImplemented
                        if self.verbose >= 2:
                            print(
                                "{0}: using Dunham coefficients for Energy calculation".format(
                                    self.get_fullname()
                                )
                            )
                        # TODO: reformat these methods as external functions, but keep docstrings
                elif self.name in HITRAN_CLASS5:
                    if convention == "herzberg":
                        from radis.levels.energies_co2 import EvJ_co2, EvJah_co2

                        self._Erovib = EvJ_co2
                        self._Ehaj = EvJah_co2
                        self.Erovib = self._Erovib_default_coefs
                        #                        self.Erovib.__doc__ == EvJ_co2.__doc__
                        self.Ehaj = self._Ehaj_default_coefs
                    #                        self.Ehaj.__doc__ == EvJah_co2.__doc__
                    else:
                        raise NotImplementedError(
                            "Only Herzberg convention spectroscopic "
                            + "constants defined for {0}. ".format(self.name)
                            + "You still define your own Erovib function "
                            + "and overwrite the ElectronicState class one."
                        )

    def Erovib(self, v, J, remove_ZPE=True):
        r"""Calculate rovibrational energy of molecule

        .. math::

            E(v,J) = G(v) + F(v,J)

        Works with Herzberg convention for the spectroscopic coefficients used
        in :func:`~radis.levels.dunham.Gv`, :func:`~radis.levels.dunham.Fv` ::

        Parameters
        ----------

        v: int
            vibrational state
            For polyatomic molecules, use ``v1, v2, v3, etc.`` instead of ``v``.

        J: int
            rotational state

        remove_ZPE: boolean
            if ``True``, removes energy of v=0,J=0 vibrational level (zero-point-energy)
            Default ``True``

        Returns
        -------

        energy of state in cm-1

        Notes
        -----

        Method is overwritten on molecule creation for polyatomic molecules.
        Refer to the correct method below:

        - CO2: :py:func:`~radis.levels.energies_co2.EvJ_co2`


        See Also
        --------

        :func:`~radis.levels.dunham.Gv`, :func:`~radis.levels.dunham.Fv`


        """

        raise NotImplementedError(
            "Rovibrational energy not implemented for {0}".format(self.get_fullname)
        )  # @dev: see radis.db.classes.ElectronicState._assign_E

    def _Erovib_default_coefs(self, *args, **kwargs):
        r"""Calculate rovibrational energy of molecule

        .. math::

            E(v,J) = G(v) + F(v,J)

        Works with Herzberg convention for the spectroscopic coefficients used
        in :func:`~radis.levels.dunham.Gv`, :func:`~radis.levels.dunham.Fv` ::

        Parameters
        ----------

        v: int
            vibrational state
            For polyatomic molecules, use ``v1, v2, v3, etc.`` instead of ``v``.

        J: int
            rotational state

        remove_ZPE: boolean
            if ``True``, removes energy of v=0,J=0 vibrational level (zero-point-energy)
            Default ``True``

        Returns
        -------

        energy of state in cm-1

        Notes
        -----

        Method is overwritten on molecule creation for polyatomic molecules.
        Refer to the correct method below:

        - CO2: :py:func:`~radis.levels.energies_co2.EvJ_co2`


        See Also
        --------

        :func:`~radis.levels.dunham.Gv`, :func:`~radis.levels.dunham.Fv`

        """
        assert "coeff_dict" not in kwargs
        kwargs.update({"coeff_dict": self.rovib_constants})

        return self._Erovib(*args, **kwargs)

    def _Ehaj_default_coefs(self, *args, **kwargs):
        """Call default _Ehaj() function with default rovib_constants
        coefficients."""

        assert "coeff_dict" not in kwargs
        kwargs.update({"coeff_dict": self.rovib_constants})

        return self._Ehaj(*args, **kwargs)

    def _E_Dunham(self, v, J, remove_ZPE=True):
        # TODO: move in levels.energies as was done for CO2

        c = self.rovib_constants

        if remove_ZPE:
            ZPE = EvJ(0, 0, **c)
        else:
            ZPE = 0

        return EvJ(v, J, **c) - ZPE

    def _E_Herzberg(self, v, J, remove_ZPE=True):
        # TODO: move in levels.energies as was done for CO2
        r"""Calculate rovibrational energy of molecule

        .. math::

            E(v,J) = G(v) + F(v,J)

        Works with Herzberg convention for the spectroscopic coefficients used
        in :func:`~radis.levels.dunham.Gv`, :func:`~radis.levels.dunham.Fv` ::

        Only works for diatomic molecules. Method should be overwritten on molecule
        creation for other molecules

        Parameters
        ----------

        v: int
            vibrational state

        J: int
            rotational state

        remove_ZPE: boolean
            if ``True``, removes energy of v=0,J=0 vibrational level (zero-point-energy)
            Default ``True``

        Returns
        -------

        energy of state in cm-1

        See Also
        --------

        :func:`~radis.levels.dunham.Gv`, :func:`~radis.levels.dunham.Fv`

        """

        try:

            # Get from constants (hardcoded, see reference above)
            c = self.rovib_constants

            # Vibrational constants
            we = c["we"]  # mandatory
            wexe = c.get("wexe", 0)  # optional: 0 is default value
            weye = c.get("weye", 0)
            weze = c.get("weze", 0)
            weae = c.get("weae", 0)
            webe = c.get("webe", 0)

            # Rotational constants
            Be = c["Be"]
            De = c.get("De", 0)

            alpha_e = c.get("alpha_e", 0)
            beta_e = c.get("beta_e", 0)
            gamma_e = c.get("gamma_e", 0)
            delta_e = c.get("delta_e", 0)
            pi_e = c.get("pi_e", 0)

            He = c.get("He", 0)
            eta_e = c.get("eta_e", 0)

            # Energies
            #            Te = c['Te']      # electronic energy
            G = Gv(v, we, wexe, weye, weze, weae, webe)  # vibrational energy
            F = Fv(
                v,
                J,
                Be,
                De,
                alpha_e,
                beta_e,
                gamma_e,
                delta_e,
                pi_e=pi_e,
                He=He,
                eta_e=eta_e,
            )  # rotational energy

            if remove_ZPE:
                ZPE = Gv(0, we, wexe, weye, weze)
            else:
                ZPE = 0

            E = G + F - ZPE  # cm-1

        except KeyError as err:
            raise KeyError(
                "Mandatory spectroscopic constant `{0}` ".format(err.args[0])
                + "not defined for electronic state {0}".format(self.get_fullname())
                + ". Check your ElectronicState definition"
            )

        return E

    # %% Informative methods

    def get_Morse_inc(self):
        r"""Get Morse potential energy increment correction for given molecule

        Examples
        --------

        ::

            inc = ElecState.get_Morse_inc()
            Delta_E(vi, vi+1) = Delta_E0 - (vi+1-vi0)*inc

        for the energy gap between vibrational level vi and vi+1,
        where vi0 is the last vibrational level calculated with Dunham
        expansion, and  and Delta_E0 = DeltaE(vi0-1, vi0)

        See Also
        --------

        :py:func:`~radis.phys.morse.morse_increment`

        """

        from radis.phys.morse import morse_increment

        c = self.rovib_constants

        if get_convention(c) == "dunham":

            # convert to get wexe, weae etc in Herzberg convention
            from radis.db.conventions import herzberg2dunham

            def convert(coef, default=0):
                if coef in c:
                    sign, Yij = herzberg2dunham[coef]
                    return sign * Yij
                else:
                    return default  # default

            we = convert("we")  # required
            wexe = convert("wexe", 0)  # optional: 0 is default value
            weye = convert("weye", 0)
            weze = convert("weze", 0)
            weae = convert("weae", 0)
            webe = convert("webe", 0)
            wece = convert("wece", 0)

            if we is None:
                raise KeyError(
                    "Mandatory spectroscopic constant `we` "
                    + "not defined for electronic state {0}".format(self.get_fullname())
                    + ". Check your ElectronicState definition. Given "
                    + "constants: {0}".format(list(c.keys()))
                )

        else:
            try:
                we = c["we"]  # required
                wexe = c.get("wexe", 0)  # optional: 0 is default value
                weye = c.get("weye", 0)
                weze = c.get("weze", 0)
                weae = c.get("weae", 0)
                webe = c.get("webe", 0)
                wece = c.get("wece", 0)

            except KeyError as err:
                raise KeyError(
                    "Mandatory spectroscopic constant `{0}` ".format(err.args[0])
                    + "not defined for electronic state {0}".format(self.get_fullname())
                    + ". Check your ElectronicState definition. Given "
                    + "constants: {0}".format(list(c.keys()))
                )

        return morse_increment(
            self.Ediss,
            we,
            wexe=wexe,
            weye=weye,
            weze=weze,
            weae=weae,
            webe=webe,
            wece=wece,
        )

    def get_fullname(self):
        return r"{0}({1}{2})(iso{3})".format(
            self.name, self.state, _print_term_symbol(self.term_symbol), self.iso
        )

    def get_statename_utf(self):
        r"""Returns Electronic state name in UTF format, used in the standard
        RADIS JSON databases. Examples::

           'X', '1Σ+'
           >>> X1SIG+

           'C', '3Πu'
           >>> C3PIu

        """

        return self.state + _format_term_symbol(self.term_symbol)


# %% Test


if __name__ == "__main__":

    from radis.test.db.test_molecules import _run_testcases

    print(("Testing molecules.py", _run_testcases()))
