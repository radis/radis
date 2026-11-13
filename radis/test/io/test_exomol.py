# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 22:44:39 2021

@author: erwan
"""

from decimal import Decimal
from fractions import Fraction

import astropy.units as u
import numpy as np
import pytest

from radis.io.exomol import get_exomol_database_list, get_exomol_full_isotope_name

conditions = {
    "wmin": 2002 / u.cm,
    "wmax": 2300 / u.cm,
    "molecule": "CO",
    "isotope": "1",
    "pressure": 1.01325,  # bar
    "mole_fraction": 0.1,
    "path_length": 1,  # cm
    "verbose": True,
}


@pytest.mark.fast
@pytest.mark.parametrize(
    "fetch_kwargs",
    [{}, {"isotope": 1}, {"isotope": "1"}, {"isotope": np.int64(1)}],
    ids=["default", "int", "str", "np.int64"],
)
def test_fetch_exomol_valid_isotope(fetch_kwargs, monkeypatch):
    """Check that valid ``isotope`` values pass validation in fetch_exomol
    and are converted to int (offline).

    See https://github.com/radis/radis/issues/1045"""

    import radis.io.exomol
    from radis.io.exomol import fetch_exomol

    class _Stop(Exception):
        pass

    received = []

    def _fake_get_exomol_full_isotope_name(molecule, isotope):
        received.append(isotope)
        raise _Stop

    monkeypatch.setattr(
        radis.io.exomol,
        "get_exomol_full_isotope_name",
        _fake_get_exomol_full_isotope_name,
    )

    with pytest.raises(_Stop):
        fetch_exomol("SiO", "EBJT", **fetch_kwargs)
    assert received == [1]
    assert type(received[0]) is int


@pytest.mark.fast
@pytest.mark.parametrize(
    "isotope",
    ["1,2", "abc", None, 1.5, 1.0, True, Decimal("1.7"), Fraction(3, 2)],
    ids=repr,
)
def test_fetch_exomol_invalid_isotope(isotope, monkeypatch):
    """Check that invalid ``isotope`` values raise before anything is fetched
    in fetch_exomol (offline).

    See https://github.com/radis/radis/issues/1045"""

    import radis.io.exomol
    from radis.io.exomol import fetch_exomol

    def _fake_get_exomol_full_isotope_name(molecule, isotope):
        raise AssertionError(f"isotope={isotope!r} passed validation")

    monkeypatch.setattr(
        radis.io.exomol,
        "get_exomol_full_isotope_name",
        _fake_get_exomol_full_isotope_name,
    )

    with pytest.raises(ValueError, match="must be an integer"):
        fetch_exomol("SiO", "EBJT", isotope=isotope)


@pytest.mark.needs_connection
def test_exomol_parsing_functions(verbose=True, *args, **kwargs):
    """Test functions used to parse ExoMol website"""

    assert get_exomol_full_isotope_name("H2O", 1) == "1H2-16O"
    assert get_exomol_full_isotope_name("H2O", 4) == "1H-2H-16O"
    assert get_exomol_full_isotope_name("CO2", 3) == "16O-12C-18O"
    assert get_exomol_full_isotope_name("H2O", 1) == "1H2-16O"
    assert get_exomol_full_isotope_name("NH3", 1) == "14N-1H3"
    assert get_exomol_full_isotope_name("H2S", 1) == "1H2-32S"
    assert get_exomol_full_isotope_name("FeH", 1) == "56Fe-1H"
    assert get_exomol_full_isotope_name("CH4", 1) == "12C-1H4"
    assert get_exomol_full_isotope_name("C2H4", 1) == "12C2-1H4"
    assert get_exomol_full_isotope_name("C2H2", 1) == "12C2-1H2"

    #  Note : this may change if new databases are added ; test would have to
    # be updated in that case.
    databases, recommended = get_exomol_database_list("CH4", "12C-1H4")
    assert sorted(databases) == sorted(["xsec-MM", "YT10to10", "YT34to10", "MM"])
    assert recommended == "MM"

    # Test that these databases are found
    KNOWN_EXOMOL_DATABASE_NAMES = {
        "H2O": {"1H2-16O": ["POKAZATEL"]},
        "NH3": {"14N-1H3": ["CoYuTe"]},
        "H2S": {"1H2-32S": ["AYT2"]},
        "FeH": {"56Fe-1H": ["MoLLIST"]},
    }
    for molecule, val in KNOWN_EXOMOL_DATABASE_NAMES.items():
        for isotope_name, lookup_databases in val.items():
            known_databases, _ = get_exomol_database_list(molecule, isotope_name)
            for dat in lookup_databases:
                assert dat in known_databases


@pytest.mark.needs_connection
@pytest.mark.needs_HITRAN_credentials
def test_calc_exomol_vs_hitemp(verbose=True, plot=False, *args, **kwargs):
    """Auto-fetch and calculate a CO spectrum from the ExoMol database
    (HITEMP lineset) and compare to HITEMP (on the HITRAN server)

    https://github.com/radis/radis/pull/320#issuecomment-884508206
    """

    from radis import SpectrumFactory

    sf = SpectrumFactory(**conditions)

    # ExoMol
    sf.fetch_databank(
        source="exomol",
        database="Li2015",
        broadf=False,
        broadf_download=False,  # accelerates the test!
    )
    s_exomol = sf.eq_spectrum(Tgas=1000, path_length=1)

    sf.fetch_databank(
        source="hitemp",
    )
    s_hitemp = sf.eq_spectrum(Tgas=1000, path_length=1, name="HITEMP (Air broadened)")

    if plot:
        s_exomol.plot(
            lw=3,
        )
        s_hitemp.plot(lw=1, nfig="same")
        import matplotlib.pyplot as plt

        plt.legend()

    # Broadening coefficients are different but areas under the lines should be the same:
    assert np.isclose(
        s_exomol.get_integral("abscoeff"), s_hitemp.get_integral("abscoeff"), rtol=0.001
    )


@pytest.mark.needs_connection
def test_exomol_parallel_vs_sequential(tmp_path, verbose=True, *args, **kwargs):
    """Test that parallel fetch_exomol produces same DataFrame as sequential"""
    import pandas as pd

    from radis.io.exomol import fetch_exomol

    molecule = "CO"
    isotope = 1
    wmin = 2000
    wmax = 2010

    def get_df(parallel, path):
        return fetch_exomol(
            molecule=molecule,
            isotope=isotope,
            load_wavenum_min=wmin,
            load_wavenum_max=wmax,
            parallel=parallel,
            clean_cache_files=False,
            verbose=False,
            local_databases=path,
        )

    df_seq = get_df(parallel=False, path=tmp_path / "seq")
    df_par = get_df(parallel=True, path=tmp_path / "par")
    pd.testing.assert_frame_equal(df_seq, df_par)



@pytest.mark.fast
def test_exomol_broadening_partners(tmp_path):
    """Diluents are matched case-insensitively to ExoMol .broad files; air and self
    are always loaded separately, and diluents without a file are skipped"""
    from types import SimpleNamespace

    from radis.io.exomol import _exomol_broadening_partners

    partners = ["H2", "He", "air", "self", "CO2"]
    files = {p: tmp_path / f"12C-16O__{p}.broad" for p in partners}
    for p in ("H2", "He", "air"):
        files[p].write_text("a0 0.07 0.5 0\n")
    mdb = SimpleNamespace(broad_partners=partners, broad_files=files)

    assert _exomol_broadening_partners({"H2": 0.85, "He": 0.15}, mdb) == ["H2", "He"]
    assert _exomol_broadening_partners({"h2": 0.5, "air": 0.5}, mdb) == ["H2"]
    assert _exomol_broadening_partners("CO2", mdb) == []  # no file
    assert _exomol_broadening_partners("air", mdb) == []
    assert _exomol_broadening_partners(None, mdb) == []


@pytest.mark.needs_connection
def test_fetch_exomol_diluent_broadening(tmp_path):
    """H2 and He .broad files give gamma_h2 / gamma_he columns, and all ExoMol
    broadening coefficients are converted from cm-1/bar to cm-1/atm"""
    from radis.api.exomolapi import read_broad
    from radis.io.exomol import BAR_PER_ATM, fetch_exomol

    df = fetch_exomol(
        molecule="CO",
        database="Li2015",
        isotope=1,
        load_wavenum_min=2000,
        load_wavenum_max=2010,
        diluent={"H2": 0.85, "He": 0.15},
        clean_cache_files=False,
        verbose=False,
        local_databases=tmp_path,
    )
    path = tmp_path / "CO" / "12C-16O" / "Li2015"
    for partner, col in (("H2", "h2"), ("He", "he"), ("air", None)):
        bdat = read_broad(path / f"12C-16O__{partner}.broad", output="pandas")
        gamma_col, n_col = ("airbrd", "Tdpair") if col is None else (f"gamma_{col}", f"n_{col}")
        if partner == "air" and (bdat["code"] != "a0").any():
            continue  # only check a0 (J'') recipes here
        table = bdat.set_index("jlower")
        line = df[df["jl"].isin(table.index)].iloc[0]
        assert np.isclose(line[gamma_col], table.loc[line["jl"], "alpha_ref"] * BAR_PER_ATM)
        assert np.isclose(line[n_col], table.loc[line["jl"], "n_Texp"])


if __name__ == "__main__":
    # test_exomol_parsing_functions()
    test_calc_exomol_vs_hitemp()
