# -*- coding: utf-8 -*-
"""Created on Fri Apr 28 10:00:25 2017.

@author: erwan

-------------------------------------------------------------------------------
"""


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
}


class MolParams:
    def __init__(self, file=None):
        """Easy access to molecular parameters taken from HITRAN molparam.txt.

        Parameters
        ----------
        file: str
            if None the one in RADIS is taken

        Examples
        --------
        Get earth abundance of CO2, isotope 1::

            molpar = Molparams()
            molpar.get(2, 1, 'abundance')         # 2 for CO2, 1 for isotope 1

        Note
        ----
        Isotope number was derived manually assuming the isonames were ordered in the database
        The isotope name (ex: CO2 626) is kept for comparison if ever needed

        References
        ----------

        http://hitran.org/media/molparam.txt
        """

        if file is None:
            file = getFile("molparam.txt")

        df = pd.read_csv(file, comment="#", delim_whitespace=True)
        df = df.set_index(["id", "iso"])

        df["isotope_name_exomol"] = _add_exomol_name(df)

        self.df = df

        # ------
        try:  # Add hints (Python >3.6 only)
            self.get.__annotations__["key"] = list(df.keys())
        except AttributeError:
            pass  # old Python version

    def get(self, M, I, key):
        """Get attribute of molecule, isotope.

        Parameters
        ----------
        M: int
            molecule id
        I: int
            molecule isotope #
        key: ``'abundance'``, ``'mol_mass'``
            parameter

        """
        try:
            float(M)
        except:
            # not a number
            from radis.db.classes import get_molecule_identifier

            M = get_molecule_identifier(M)
        return self.df.loc[(M, I), key]


def _add_exomol_name(df):
    """Convert HITRAN isotope name from :py:data:`~radis.db.molparam.isotope_name_dict`
    in ExoMol format
    """

    # Add in molparams

    df["isotope_name"] = df.index.map(isotope_name_dict)

    # add exomol full name:
    # deal with Hydrogen : H2, HD
    isotope_fullname = (
        df.isotope_name.str.replace("H2", "(1H2)", regex=False)
        .str.replace("HD", "(1H)(2H)", regex=False)
        .str.replace("H3", "(1H3)", regex=False)
    )
    # Convert:
    return (
        isotope_fullname.str.replace(")(", "-", regex=False)
        .str.replace("(", "", regex=False)
        .str.replace(")", "", regex=False)
    )


def _test(verbose=True, *args, **kwargs):
    molpar = MolParams()

    b = True

    b *= molpar.get(2, 2, "abundance") == 0.0110574
    if verbose:
        print("CO2-636 abundance:", molpar.get(2, 2, "abundance"))

    if verbose:
        print("Testing molparams.py:", b)

    return True


if __name__ == "__main__":

    _test()
