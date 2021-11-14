# -*- coding: utf-8 -*-
"""
Test that RADIS molecule lists are always up-to-date

Note: Radis HITRAN molecules remain hardcoded so we do not require an internet
connection; fetching is only done while running tests
"""

from urllib.request import HTTPError, urlopen

import pytest
from bs4 import BeautifulSoup
from tqdm import tqdm


def fetch_ExoMol_molecules():
    """ Get list of molecules in ExoMol by scrapping the ExoMol page """

    url = "https://exomol.com/data/molecules/"
    try:
        response = urlopen(url).read()
    except HTTPError as err:
        raise ValueError(f"HTTPError opening url={url}") from err

    soup = BeautifulSoup(
        response, features="lxml"
    )  # make soup that is parse-able by bs

    # All others
    rows = soup.find_all(
        "a", {"class": "list-group-item link-list-group-item molecule_link"}
    )
    molecules = [r.get_attribute_list("href")[0] for r in rows]

    return molecules


@pytest.mark.needs_connection
def test_ExoMol_molecules_list(verbose=True, *args, **kwargs):
    """ Test that ExoMol molecule list in RADIS remains up to date"""

    from radis.db.classes import EXOMOL_MOLECULES
    from radis.misc.basics import compare_lists

    molecules = fetch_ExoMol_molecules()

    if verbose:
        print("ExoMol molecules, fetched online ")
        print(molecules)

        print("Comparing Radis hardcoded ExoMol molecules list to the ExoMol website: ")

    assert (
        compare_lists(
            EXOMOL_MOLECULES,
            molecules,
            l1_str="Radis molecules",
            l2_str="Fetched from ExoMol website",
            print_index=True,
        )
        == 1
    )


def fetch_ExoMol_isotopes(verbose=True, *args, **kwargs):
    """ Get dictionary of isotope names in ExoMol by scrapping the ExoMol pages """

    isotopes_full_names = {}
    molecules = fetch_ExoMol_molecules()

    if verbose:
        print("Parsing ExoMol molecules", flush=True)
    for molecule in tqdm(molecules, disable=not verbose):

        url = f"https://exomol.com/data/molecules/{molecule}"
        try:
            response = urlopen(url).read()
        except HTTPError as err:
            raise ValueError(f"HTTPError opening url={url}") from err

        soup = BeautifulSoup(
            response, features="lxml"
        )  # make soup that is parse-able by bs

        # Recommended database
        rows = soup.find_all("a", {"class": "list-group-item link-list-group-item"})
        isotopes = [r.get_attribute_list("href")[0] for r in rows]

        isotopes_full_names[molecule] = isotopes
    return isotopes_full_names


#%%


def fetch_HITRAN_molecules():

    import re

    url = "https://hitran.org/lbl/"
    try:
        response = urlopen(url).read()
    except HTTPError as err:
        raise ValueError(f"HTTPError opening url={url}") from err

    soup = BeautifulSoup(
        response, features="lxml"
    )  # make soup that is parse-able by bs

    # All others
    rows = soup.find_all("td", {"class": "molecule-formula"})
    molecules = [re.sub(r"<[^>]*>", "", str(r)) for r in rows]

    return molecules


@pytest.mark.needs_connection
def test_HITRAN_molecules_list(verbose=True, *args, **kwargs):
    """ Test that HITRAN molecule list in RADIS remains up to date"""

    from radis.db.classes import HITRAN_MOLECULES
    from radis.misc.basics import compare_lists

    molecules = fetch_HITRAN_molecules()

    if verbose:
        print("HITRAN molecules, fetched online ")
        print(molecules)

        print("Comparing Radis hardcoded HITRAN molecules list to the HITRAN website: ")

    assert (
        compare_lists(
            HITRAN_MOLECULES,
            molecules,
            l1_str="Radis molecules",
            l2_str="Fetched from HITRAN website",
            print_index=True,
        )
        == 1
    )


#%%


def generate_molparam_for_non_HITRAN_species():
    """Simple estimator of molecular parameters when not given in HITRAN.

    isotopic abundance  : in a simple approach, we multiply isotopic abundances
    of all elements of the molecule .  We typically get 0.04% error for main isotope,
    3-10% error in general , up to 80% in the worst cases, so always the same order of magnitude,
    (when compared with species that exist in HITRAN)

    molar mass : we sum the mass of all elements


    Output ``EXTRA_ABUNDANCES_DICT, EXTRA_MOLAR_MASS_DICT, EXTRA_ISOTOPE_FULLNAME_DICT``
    should be pasted in ``default_radis.json``
    """

    #%%
    import numpy as np

    exomol_molecules = fetch_ExoMol_molecules()
    # hitran_molecules = fetch_HITRAN_molecules()
    isotopes_full_names = fetch_ExoMol_isotopes()

    #%%
    from mendeleev import element  # may create ImportError, not added in requirement.

    Ia_dict = {}
    Mm_dict = {}
    for atom in [
        "H",
        "C",
        "O",
        "N",
        "S",
        "F",
        "Cl",
        "Br",
        "I",
        "P",
        "Mg",
        "Na",
        "Ni",
        "Al",
        "Cr",
        "Ca",
        "Be",
        "Ti",
        "Fe",
        "Li",
        "Sc",
        "Si",
        "V",
        "Y",
        "K",
        "As",
        "He",
    ]:
        el = element(atom)
        Ia_dict[atom] = {}
        Mm_dict[atom] = {}
        for iso in el.isotopes:
            if iso.abundance is None:
                Ia_dict[atom][iso.mass_number] = 1e-6  # arbitrary
            else:
                Ia_dict[atom][iso.mass_number] = iso.abundance
            Mm_dict[atom][iso.mass_number] = iso.mass
    # add missing:
    Ia_dict["Al"][26] = 1e-6  # arbitrary
    Mm_dict["Al"][26] = 26.9815385 - 1.006  # doesn't exist by default
    Mm_dict["H"][3] = 2.0141017781 + 1.006  # 3H  'None' by default
    Mm_dict["C"][14] = 13.003354835 + 1.006  # 14C

    import re

    _parse_molecule = re.compile("(\d+)+([A-Z][a-z]?)(\d+)?")

    def get_isotopic_atomic_composition(molecule):

        atoms = {}
        for mass_number, atom, number in _parse_molecule.findall(molecule):
            if number == "":
                number = 1
            if atom in atoms:
                atoms[atom][int(mass_number)] = int(number)
            else:
                atoms[atom] = {int(mass_number): int(number)}
        return atoms

        # out = dict(_parse_molecule.findall(molecule))
        # for k, v in out.items():
        #     if v == '':
        #         out[k] = 1
        #     else:
        #         out[k] = int(v)
        # return out

    # print(get_atoms_in("CO2"))
    print(get_isotopic_atomic_composition("12C-16O2"))

    def estimate_terrestrial_abundance(molecule):
        """``molecule`` given in ExoMol full isotope format"""
        abundance = 1
        atoms_isotopes = get_isotopic_atomic_composition(molecule)
        for atom in atoms_isotopes:
            for Z, N_Z in atoms_isotopes[atom].items():
                abundance *= Ia_dict[atom][Z] ** N_Z
                # if Z in Ia_dict[atom]:
                # else:
                #     abundance *= (1e-6)**N_Z  # arbitrary (low)
        return abundance

    def get_molar_mass(molecule):
        """``molecule`` given in ExoMol full isotope format"""
        atoms_isotopes = get_isotopic_atomic_composition(molecule)
        molar_mass = sum(
            [
                Mm_dict[atom][Z] * N_Z
                for atom in atoms_isotopes
                for Z, N_Z in atoms_isotopes[atom].items()
            ]
        )
        return molar_mass

    assert get_molar_mass("12C-16O2") == 43.98982924

    # Test : compare with terrestrial abundance when known:

    from radis.db.classes import get_molecule
    from radis.db.molparam import MolParams

    molpar = MolParams()

    for _, r in molpar.df.reset_index().iterrows():
        print(get_molecule(r.id), "-", r.isotope_name_exomol)
        if not r.isotope_name_exomol:
            raise
        estimate = estimate_terrestrial_abundance(r.isotope_name_exomol)
        print(
            "... Estim: {0:.3e}".format(
                estimate_terrestrial_abundance(r.isotope_name_exomol)
            )
        )
        print(
            "... Real:  {0:.3e}".format(r.abundance),
            "         ({0:.2%}% error)".format(
                abs(estimate - r.abundance) / r.abundance
            ),
        )

    #%% Generate dictionary of terrestrial abuneances for all molecules not in HITRAN
    EXTRA_ABUNDANCES_DICT = {}
    EXTRA_MOLAR_MASS_DICT = {}
    EXTRA_ISOTOPE_FULLNAME_DICT = {}

    print("Check if isotopologues in ExoMol are sorted as in HITRAN")
    for M in exomol_molecules:
        print(M)
        abundances_naive = np.array(
            [
                estimate_terrestrial_abundance(isotope_fullname)
                for isotope_fullname in isotopes_full_names[M]
            ]
        )
        b_sort_naive = np.argsort(abundances_naive)[::-1]  # most abundant isotope first
        abundances_naive = abundances_naive[b_sort_naive]
        molar_mass = [
            get_molar_mass(isotope_fullname)
            for isotope_fullname in np.array(isotopes_full_names[M])[b_sort_naive]
        ]
        # sorted_naive = np.array(isotopes_full_names[M])[b]
        for i, isotope_fullname in enumerate(
            np.array(isotopes_full_names[M])[b_sort_naive]
        ):
            # check if isotopes are sorted by abundances
            try:
                sorted_equally = (
                    molpar.get(M, i + 1, key="isotope_name_exomol") == isotope_fullname
                )
                # sorted_naively = molpar.get(M, i+1, key="isotope_name_exomol") == isotope_fullname
                print(
                    "... ",
                    isotope_fullname,
                    "\t",
                    "ok"
                    if sorted_equally
                    else "⚠️ DIFFERENT ({0})".format(
                        molpar.get(M, i + 1, key="isotope_name_exomol")
                    ),
                )
            except (KeyError, NotImplementedError):
                print("... ", isotope_fullname, "not in HITRAN")
                try:  # add isotope i+1
                    EXTRA_ABUNDANCES_DICT[M][str(i + 1)] = abundances_naive[i]
                    EXTRA_MOLAR_MASS_DICT[M][str(i + 1)] = molar_mass[i]
                    EXTRA_ISOTOPE_FULLNAME_DICT[M][str(i + 1)] = isotope_fullname
                except KeyError:
                    EXTRA_ABUNDANCES_DICT[M] = {str(i + 1): abundances_naive[i]}
                    EXTRA_MOLAR_MASS_DICT[M] = {str(i + 1): molar_mass[i]}
                    EXTRA_ISOTOPE_FULLNAME_DICT[M] = {str(i + 1): isotope_fullname}

    EXTRA_ABUNDANCES_DICT = {
        k: EXTRA_ABUNDANCES_DICT[k] for k in sorted(EXTRA_ABUNDANCES_DICT)
    }
    EXTRA_MOLAR_MASS_DICT = {
        k: EXTRA_MOLAR_MASS_DICT[k] for k in sorted(EXTRA_MOLAR_MASS_DICT)
    }
    EXTRA_ISOTOPE_FULLNAME_DICT = {
        k: EXTRA_ISOTOPE_FULLNAME_DICT[k] for k in sorted(EXTRA_ISOTOPE_FULLNAME_DICT)
    }

    return EXTRA_ABUNDANCES_DICT, EXTRA_MOLAR_MASS_DICT, EXTRA_ISOTOPE_FULLNAME_DICT


#%%

if __name__ == "__main__":

    test_ExoMol_molecules_list()
    test_HITRAN_molecules_list()
