# -*- coding: utf-8 -*-
"""
Test names and labels
"""

import pytest

from radis.lbl.labels import (
    vib_lvl_name_cdsd_pc,
    vib_lvl_name_cdsd_pcN,
    vib_lvl_name_hitran_class1,
    vib_lvl_name_hitran_class5,
    vib_lvl_name_hitran_class5_short,
)
from radis.misc.printer import printm


@pytest.mark.fast
def test_vibrational_levels_labelling(verbose=True, *args, **kwargs):

    # If verbose: print some
    if verbose:
        printm("Some vibrational level formats:")

        printm("... CO, HITRAN format (v):\t\t", vib_lvl_name_hitran_class1(10))

        printm(
            "... CO2, HITRAN format (v1,v2,l2,v3):\t\t",
            vib_lvl_name_hitran_class5(2, 1, 1, 3),
        )
        printm(
            "... CO2, HITRAN format, short (v1v2l2v3):\t",
            vib_lvl_name_hitran_class5_short(2, 1, 1, 3),
        )

        printm("... CO2, CDSD format (p,c):\t\t", vib_lvl_name_cdsd_pc(14, 1))
        printm("... CO2, CDSD format (p,c,N):\t\t", vib_lvl_name_cdsd_pcN(14, 1, 1))

    # Do the tests

    assert vib_lvl_name_hitran_class1(10) == "(10)"

    assert vib_lvl_name_hitran_class5(2, 1, 1, 3) == "(2,1,1,3)"
    assert vib_lvl_name_hitran_class5_short(2, 1, 1, 3) == "21`1`3"

    assert vib_lvl_name_cdsd_pc(14, 1) == "(14,1)"
    assert vib_lvl_name_cdsd_pcN(14, 1, 1) == "(14,1,1)"

    return


if __name__ == "__main__":

    test_vibrational_levels_labelling()
