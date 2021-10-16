# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 00:47:30 2021

@author: erwan
"""

from radis.db.molparam import MolParams


def test_molparams(verbose=True, *args, **kwargs):

    molpar = MolParams()

    assert molpar.get(2, 2, "abundance") == 0.0110574


if __name__ == "__main__":
    test_molparams()
