import numpy as np
import pandas as pd
import pytest

from radis.api.kuruczapi import load_pf_Barklem2016, read_kurucz
from radis.test.utils import getTestFile


@pytest.mark.fast
def test_read_kurucz():
    df = read_kurucz(
        getTestFile("gf4000.all")
    )  # from: http://kurucz.harvard.edu/atoms/4000/gf4000.all
    a = df.loc[
        5,
        [
            "orig_wavelen",
            "loggf",
            "species",
            "El",
            "jl",
            "labellower",
            "Eu",
            "ju",
            "labelupper",
            "gamRad",
            "gamSta",
            "gamvdW",
            "iso",
        ],
    ]  # .to_numpy()
    b = np.array(
        [
            504403.3610,
            -7.259,
            "40.00",
            35456.250,
            2.0,
            "d3p w3P",
            35476.070,
            3.0,
            "5s4F6s e5F",
            8.41,
            -5.27,
            -7.53,
            0,
        ],
        dtype="object",
    )
    assert (a == b).all()
    a = df.loc[
        7,
        [
            "orig_wavelen",
            "loggf",
            "species",
            "El",
            "jl",
            "labellower",
            "Eu",
            "ju",
            "labelupper",
            "gamRad",
            "gamSta",
            "gamvdW",
            "iso",
        ],
    ]  # .to_numpy()
    b = np.array(
        [
            328425.5787,
            -6.296,
            "40.00",
            21943.740,
            2.0,
            "d4 a5D",
            21974.180,
            1.0,
            "d2sp z3S",
            7.22,
            -5.99,
            -7.79,
            0,
        ],
        dtype="object",
    )
    assert (a == b).all()


@pytest.mark.fast
def test_barklem_pf():
    # Testing load_pf_Barklem2016 function
    pfTdat, pfdat = load_pf_Barklem2016()
    assert isinstance(pfTdat, pd.Series)
    assert isinstance(pfdat, pd.DataFrame)
    assert len(pfTdat) == 42
    assert len(pfdat) == 284


if __name__ == "__main__":
    test_read_kurucz()
    test_barklem_pf()
