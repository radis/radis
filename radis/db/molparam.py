# -*- coding: utf-8 -*-
"""Created on Fri Apr 28 10:00:25 2017.

@author: erwan

-------------------------------------------------------------------------------
"""


import pandas as pd

from radis.db.utils import getFile


class MolParams:
    def __init__(self, file=None):
        """Easy access to molecular parameters taken from HITRAN molparam.txt.

        Parameters
        ----------

        file: str
            if None the one in RADIS is taken

        Examples
        --------

        Get abundance of CO2, isotope 1::

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

        self.df = df

        # ------
        try:  # Add hints (Python >3.6 only)
            self.get.__annotations__["key"] = list(df.keys())
        except AttributeError:
            pass  # old Python version

    def get(self, M, I, key):
        """
        Parameters
        ----------

        M: int
            molecule id
            # TODO: allow name here with get_molecule()

        I: int
            molecule isotope #

        key: ``'abundance'``, ``'mol_mass'``
            parameter

        """
        return self.df.loc[(M, I), key]


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
