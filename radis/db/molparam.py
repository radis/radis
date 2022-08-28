# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
"""


import re

import pandas as pd

from radis.db.utils import getFile

# Isotope definitions in HITRAN
isotope_name_dict = {
    (1, 1): "H2(16O)",
    (1, 2): "H2(18O)",
    (1, 3): "H2(17O)",
    (1, 4): "HD(16O)",
    (1, 5): "HD(18O)",
    (1, 6): "HD(17O)",
    (1, 7): "D2(16O)",
    (2, 1): "(12C)(16O)2",
    (2, 2): "(13C)(16O)2",
    (2, 3): "(16O)(12C)(18O)",
    (2, 4): "(16O)(12C)(17O)",
    (2, 5): "(16O)(13C)(18O)",
    (2, 6): "(16O)(13C)(17O)",
    (2, 7): "(12C)(18O)2",
    (2, 8): "(17O)(12C)(18O)",
    (2, 9): "(12C)(17O)2",
    (2, 10): "(13C)(18O)2",
    (2, 11): "(18O)(13C)(17O)",
    (2, 12): "(13C)(17O)2",
    (3, 1): "(16O)3",
    (3, 2): "(16O)(16O)(18O)",
    (3, 3): "(16O)(18O)(16O)",
    (3, 4): "(16O)(16O)(17O)",
    (3, 5): "(16O)(17O)(16O)",
    (4, 1): "(14N)2(16O)",
    (4, 2): "(14N)(15N)(16O)",
    (4, 3): "(15N)(14N)(16O)",
    (4, 4): "(14N)2(18O)",
    (4, 5): "(14N)2(17O)",
    (5, 1): "(12C)(16O)",
    (5, 2): "(13C)(16O)",
    (5, 3): "(12C)(18O)",
    (5, 4): "(12C)(17O)",
    (5, 5): "(13C)(18O)",
    (5, 6): "(13C)(17O)",
    (6, 1): "(12C)H4",
    (6, 2): "(13C)H4",
    (6, 3): "(12C)H3D",
    (6, 4): "(13C)H3D",
    (7, 1): "(16O)2",
    (7, 2): "(16O)(18O)",
    (7, 3): "(16O)(17O)",
    (8, 1): "(14N)(16O)",
    (8, 2): "(15N)(16O)",
    (8, 3): "(14N)(18O)",
    (9, 1): "(32S)(16O)2",
    (9, 2): "(34S)(16O)2",
    (10, 1): "(14N)(16O)2",
    (11, 1): "(14N)H3",
    (11, 2): "(15N)H3",
    (12, 1): "H(14N)(16O)3",
    (12, 2): "H(15N)(16O)3",
    (13, 1): "(16O)H",
    (13, 2): "(18O)H",
    (13, 3): "(16O)D",
    (14, 1): "H(19F)",
    (14, 2): "D(19F)",
    (15, 1): "H(35Cl)",
    (15, 2): "H(37Cl)",
    (15, 3): "D(35Cl)",
    (15, 4): "D(37Cl)",
    (16, 1): "H(79Br)",
    (16, 2): "H(81Br)",
    (16, 3): "D(79Br)",
    (16, 4): "D(81Br)",
    (17, 1): "H(127I)",
    (17, 2): "D(127I)",
    (18, 1): "(35Cl)(16O)",
    (18, 2): "(37Cl)(16O)",
    (19, 1): "(16O)(12C)(32S)",
    (19, 2): "(16O)(12C)(34S)",
    (19, 3): "(16O)(13C)(32S)",
    (19, 4): "(16O)(12C)(33S)",
    (19, 5): "(18O)(12C)(32S)",
    (20, 1): "H2(12C)(16O)",
    (20, 2): "H2(13C)(16O)",
    (20, 3): "H2(12C)(18O)",
    (21, 1): "H(16O)(35Cl)",
    (21, 2): "H(16O)(37Cl)",
    (22, 1): "(14N)2",
    (22, 2): "(14N)(15N)",
    (23, 1): "H(12C)(14N)",
    (23, 2): "H(13C)(14N)",
    (23, 3): "H(12C)(15N)",
    (24, 1): "(12C)H3(35Cl)",
    (24, 2): "(12C)H3(37Cl)",
    (25, 1): "H2(16O)2",
    (26, 1): "(12C)2H2",
    (26, 2): "(12C)(13C)H2",
    (26, 3): "(12C)2HD",
    (27, 1): "(12C)2H6",
    (27, 2): "(12C)H3(13C)H3",
    (28, 1): "(31P)H3",
    (29, 1): "(12C)(16O)(19F)2",
    (29, 2): "(13C)(16O)(19F)2",
    (30, 1): "(32S)(19F)6",
    (31, 1): "H2(32S)",
    (31, 2): "H2(34S)",
    (31, 3): "H2(33S)",
    (32, 1): "H(12C)(16O)(16O)H",
    (33, 1): "H(16O)2",
    (34, 1): "(16O)",
    (35, 1): "(35Cl)(16O)(14N)(16O)2",
    (35, 2): "(37Cl)(16O)(14N)(16O)2",
    (36, 1): "(14N)(16O)+",
    (37, 1): "H(16O)(79Br)",
    (37, 2): "H(16O)(81Br)",
    (38, 1): "(12C)2H4",
    (38, 2): "(12C)H2(13C)H2",
    (39, 1): "(12C)H3(16O)H",
    (40, 1): "(12C)H3(79Br)",
    (40, 2): "(12C)H3(81Br)",
    (41, 1): "(12C)H3(12C)(14N)",
    (42, 1): "(12C)(19F)4",
    (43, 1): "(12C)4H2",
    (44, 1): "H(12C)3(14N)",
    (45, 1): "H2",
    (45, 2): "HD",
    (46, 1): "(12C)(32S)",
    (46, 2): "(12C)(34S)",
    (46, 3): "(13C)(32S)",
    (46, 4): "(12C)(33S)",
    (47, 1): "(32S)(16O)3",
    (48, 1): "(12C)2(14N)2",
    (49, 1): "(12C)(16O)(35Cl)2",
    (49, 2): "(12C)(16O)(35Cl)(37Cl)",
    (49, 2): "(12C)(16O)(35Cl)(37Cl)",
    (50, 1): "(32S)(16O)",
    (50, 2): "(34S)(16O)",
    (50, 3): "(32S)(18O)",
    (51, 1): "(12C)H3(19F)",
    (52, 1): "(74Ge)H4",
    (52, 2): "(72Ge)H4",
    (52, 3): "(70Ge)H4",
    (52, 4): "(73Ge)H4",
    (52, 5): "(76Ge)H4",
    (53, 1): "(12C)(32S)2",
    (53, 2): "(32S)(12C)(34S)",
    (53, 3): "(32S)(12C)(33S)",
    (53, 4): "(13C)(32S)2",
    (54, 1): "(12C)H3(127I)",
    (55, 1): "(14N)(19F)3",
}
# TODO : add test in test_molparams comparing dictionary above to
# latest values from HAPI; to make sure we never have any problem
# even if HITRAN eventually changes the conventions and labels.


#%%
# Read Extra parameters taken from molparams_extra.json
# TODO : Refactor : have either .json either .txt, but not both

from json import JSONDecodeError
from os.path import exists, join

import hjson

from radis.misc.utils import getProjectRoot

# %% Functions to parse radis/config.json

MOLPARAMS_EXTRA_PATH = join(getProjectRoot(), "db", "molparams_extra.json")
assert exists(MOLPARAMS_EXTRA_PATH)


def get_extra_molparams(path=MOLPARAMS_EXTRA_PATH):
    """Read the extra molecular parameters in
    :py:attr:`~radis.db.molparam.MOLPARAMS_EXTRA_PATH` to complement the
    HITRAN molecular parameters of molparam.txt"""

    with open(path) as f:
        try:
            molparams_dict = hjson.load(f)
            # read with Hjson to allow comments in MOLPARAMS_EXTRA_PATH
        except JSONDecodeError as err:
            raise JSONDecodeError(
                "Error reading '{0}' (line {2} col {3}): \n{1}".format(
                    path, err.msg, err.lineno, err.colno
                ),
                err.doc,
                err.pos,
            ) from err

    return molparams_dict


#%%


class MolParams(object):
    """Easy access to molecular parameters taken from HITRAN `molparam.txt`.

    .. note::
        Starting from Radis 0.13.1 MolParams can also read extra parameters
        from a .json file

    Parameters
    ----------
    file: str
        if None the one in RADIS is taken. See https://github.com/radis/radis/blob/master/radis/db/molparam.txt
    extra_file_json: str
        if not ``None``, read extra molecular parameters from a Json file, for
        molecules not implemented in HITRAN  `molparam.txt`.
        See :py:attr:`~radis.db.molparam.MOLPARAMS_EXTRA_PATH`

    Examples
    --------
    Get earth abundance of CO2, isotope 1::

        from radis.db.molparam import MolParams
        molpar = MolParams()
        molpar.get(2, 1, 'abundance')         # 2 for CO2, 1 for isotope 1

    .. minigallery:: radis.db.molparam.MolParams
    :add-heading:

    Note
    ----
    Isotope number was derived manually assuming the isonames were ordered in the database
    The isotope name (ex: CO2 626) is kept for comparison if ever needed

    References
    ----------
    http://hitran.org/media/molparam.txt
    """

    __slots__ = ["df", "terrestrial_abundances", "extra_molparams"]

    def __init__(
        self,
        file=None,
        terrestrial_abundances=True,
        extra_file_json=MOLPARAMS_EXTRA_PATH,
    ):
        if file is None:
            file = getFile("molparam.txt")

        df = pd.read_csv(file, comment="#", delim_whitespace=True)
        df = df.set_index(["id", "iso"])

        df["isotope_name_exomol"] = _add_exomol_name(df)

        self.df = df
        self.terrestrial_abundances = terrestrial_abundances
        self.extra_molparams = {}

        if extra_file_json is not None:
            self.extra_molparams.update(get_extra_molparams(extra_file_json))

        # ------
        try:  # Add hints (Python >3.6 only)
            self.get.__annotations__["key"] = list(df.keys())
        except AttributeError:
            pass  # old Python version

    def _get_molecule_list(self):
        """Returns list of [molecule_id, molecule_name] in loaded molparam dataframe

        Examples
        --------
        ::

            [1,  'H2O',  2,  'CO2',  3,  'O3', 4, 'N2O',  etc.  ]
        """
        from radis.db.classes import get_molecule

        return [
            item
            for m in self.df.reset_index(inplace=False)["id"].unique()
            for item in [m, get_molecule(m)]
        ]

    def get(self, molecule, isotope, key):
        r"""Get attribute of molecule, isotope.

        Parameters
        ----------
        molecule: str or int
            molecule name or molecule id
        iso: int
            isotope number # (sorted by terrestrial abundance)
        key: ``'abundance'``, ``'mol_mass'``, ``'isotope_name'``, ``'isotope_name_exomol'``
            parameter

        Examples
        --------
        Get explicit name of an isotope::

            from radis.db.molparam import MolParams
            mp = MolParams()
            print(mp.get("CO2", 1, "isotope_name"))   # >> (12C)(16O)2
            print(mp.get("CO2", 2, "isotope_name"))   # >> (13C)(16O)2


        .. minigallery:: radis.db.molparam.MolParams.get
            :add-heading:

        """
        if key == "mol_mass":  # backward compatibility >=0.13.1
            key = "molar_mass"

        if key not in self.df.columns:
            raise ValueError(
                f"Unknown molecular parameter : {key}. Use one of {self.df.columns}"
            )

        M = molecule
        try:
            try:
                float(M)
            except ValueError:
                # M is not a number (i.e. not the HITRAN molecule id)
                from radis.db.classes import get_molecule_identifier

                M = get_molecule_identifier(M)
        except NotImplementedError:
            # Molecule not in HITRAN . Try to read from extra parameters
            if key in self.extra_molparams and M in self.extra_molparams[key]:
                try:
                    return self.extra_molparams[key][M][str(isotope)]
                except KeyError:

                    raise NotImplementedError(
                        f"{key} for isotope `{str(isotope)}` for molecule {M} not supported"
                        + "\nChoose one of{0}".format(
                            list(self.extra_molparams[key][M].keys())
                        )
                        + f"\nYou can manually add it in the {MOLPARAMS_EXTRA_PATH} file"
                        + "{'molparams':{'molar_mass':{'molecule':{'ISOTOPE':...}}}}; or in the `SpectrumFactory.molpar.extra_molparams[key][molecule][isotope]` dictionary. Please also report on GitHub so we can update !"
                    )
            else:
                raise NotImplementedError(
                    f"{key} for molecule '{M}' not supported"
                    + "\nChoose one of{0}".format(self._get_molecule_list())
                    + "\nor `\n{0}".format(
                        list(self.extra_molparams.get(key, {}).keys())
                    )
                    + f"\nYou can manually add it in the {MOLPARAMS_EXTRA_PATH} file"
                    + "{'molparams':{'molar_mass':{'molecule':{'ISOTOPE':...}}}}; or in the `SpectrumFactory.molpar.extra_molparams[key][molecule][isotope]` dictionary. Please also report on GitHub so we can update !"
                )

        else:
            return self.df.loc[(M, isotope), key]


def _add_exomol_name(df):
    r"""Convert HITRAN isotope name from :py:data:`~radis.db.molparam.isotope_name_dict`
    in ExoMol format
    """

    df["isotope_name"] = df.index.map(isotope_name_dict)

    # add exomol full name:
    def replace_parenthesis(s):
        s = str(s)
        # add parenthesis for hydrogen : Eg. H2 -> (1H)2,   H->(1H)
        s = re.sub(r"H([1-9])*", r"(1H)\1", s)
        # add parenthesis for deuterium : Eg. D2 -> (2H)2,   D->(2H)
        s = re.sub(r"D([1-9])*", r"(2H)\1", s)

        # Cut elements and replace parenthesis with - :
        elements = re.findall(r"\(.*?\)[1-9]*", s)
        elements = [e.replace(r"(", "").replace(")", "") for e in elements]
        s = "-".join(elements)

        return s

    return df.isotope_name.apply(replace_parenthesis)


if __name__ == "__main__":
    from radis.test.io.test_exomol import test_exomol_parsing_functions

    test_exomol_parsing_functions()

    from radis.test.db.test_molparams import test_molparams

    test_molparams()
