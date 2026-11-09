import numpy as np
import pytest

from radis.api.exomolapi import _map_m0_parameter, check_code_level


@pytest.mark.fast
@pytest.mark.parametrize("bdat_list", [["a1"], ["a0", "a1"], ["a1", "a0"]])
def test_check_bdat_a1(bdat_list):
    bdat = {}
    bdat["code"] = bdat_list
    assert check_code_level(bdat) == "a1"


@pytest.mark.fast
@pytest.mark.parametrize("bdat_list", [["a0"]])
def test_check_bdat_a0(bdat_list):
    bdat = {}
    bdat["code"] = bdat_list
    assert check_code_level(bdat) == "a0"


@pytest.mark.fast
@pytest.mark.parametrize("bdat_list", [["a0", "a1", "a2"]])
def test_check_bdat_no_code_level(bdat_list):
    """a2 is not a valid code level"""
    bdat = {}
    bdat["code"] = bdat_list
    assert check_code_level(bdat) == None


@pytest.mark.fast
def test_map_m0_parameter_replaces_missing_values():
    import pandas as pd

    values = pd.Series([1, 2, 3])
    mapped, missing_count = _map_m0_parameter(values, {1: 0.1, 2: np.nan}, 0.07)

    assert missing_count == 2
    assert mapped.tolist() == [0.1, 0.07, 0.07]


@pytest.mark.fast
def test_make_j2b_beyond_table_uses_last_jlower():
    """J'' above the .broad table get the values of its largest J'' (as PyExoCross),
    J'' missing inside the table keep the defaults"""
    import pandas as pd

    from radis.api.exomolapi import make_j2b
    from radis.misc.warning import AccuracyWarning

    bdat = pd.DataFrame(
        {
            "code": ["a0", "a0", "a0"],
            "alpha_ref": [0.09, 0.08, 0.06],
            "n_Texp": [0.6, 0.55, 0.45],
            "jlower": [0, 1, 3],
        }
    )
    with pytest.warns(AccuracyWarning):
        j2alpha_ref, j2n_Texp = make_j2b(
            bdat, alpha_ref_default=0.07, n_Texp_default=0.5, jlower_max=6
        )
    assert j2alpha_ref.tolist() == [0.09, 0.08, 0.07, 0.06, 0.06, 0.06, 0.06]
    assert j2n_Texp.tolist() == [0.6, 0.55, 0.5, 0.45, 0.45, 0.45, 0.45]


@pytest.mark.fast
def test_broadening_from_recipes_layers_a0_m0():
    """Files mixing recipes (e.g. NH3: a0, m0, m1) use m0 where it has an entry
    (|m| = J'' for P and Q, J'' + 1 for R branch lines) and a0 elsewhere"""
    import pandas as pd

    from radis.api.exomolapi import broadening_from_recipes
    from radis.misc.warning import AccuracyWarning

    bdat = pd.DataFrame(
        {
            "code": ["a0", "a0", "a0", "m0", "m0", "m1"],
            "alpha_ref": [0.09, 0.08, 0.06, 0.11, 0.10, 0.5],
            "n_Texp": [0.6, 0.55, 0.45, 0.7, 0.65, 0.9],
            "jlower": [0, 1, 2, 1, 2, 1],
            "jupper": [np.nan] * 5 + [0],
        }
    )
    #               P(1)  Q(1)  R(1)  R(0)  P(2)  R(2)  Q(5)
    jlower = np.array([1, 1, 1, 0, 2, 2, 5])
    jupper = np.array([0, 1, 2, 1, 1, 3, 5])
    with pytest.warns(AccuracyWarning, match="m1"):
        alpha_ref, n_Texp = broadening_from_recipes(bdat, jlower, jupper)
    # |m|:            1     1     2     1     2     3(a0)  5 (a0, beyond: J''=2)
    assert alpha_ref.tolist() == [0.11, 0.11, 0.10, 0.11, 0.10, 0.06, 0.06]
    assert n_Texp.tolist() == [0.7, 0.7, 0.65, 0.7, 0.65, 0.45, 0.45]

    # m0 alone: |m| beyond the table carry its largest |m|
    m0 = bdat[bdat["code"] == "m0"]
    alpha_ref, _ = broadening_from_recipes(m0, jlower, jupper)
    assert alpha_ref.tolist() == [0.11, 0.11, 0.10, 0.11, 0.10, 0.10, 0.10]


@pytest.mark.fast
@pytest.mark.parametrize("engine", ["pytables", "vaex"])
@pytest.mark.parametrize("with_trans_wavenumber", [True, False])
def test_pickup_gE_line_positions_from_states(engine, with_trans_wavenumber):
    """nu_lines = E_upper - E_lower from the states file, even when the .trans
    wavenumber column disagrees (e.g. computed before MARVELisation)"""
    import warnings

    import pandas as pd

    from radis.api.exomolapi import pickup_gE

    states = pd.DataFrame(
        {
            "i": [1, 2, 3],
            "E": [0.0, 1431.977929, 2800.678343],
            "g": [4, 4, 4],
            "J": [0.5, 0.5, 0.5],
        }
    )
    trans = pd.DataFrame({"i_upper": [2, 3], "i_lower": [1, 2], "A": [1.0, 2.0]})
    if with_trans_wavenumber:
        trans["nu_lines"] = [1432.066979, 1368.738823]  # calculated energies
    if engine == "vaex":
        vaex = pytest.importorskip("vaex")
        states, trans = vaex.from_pandas(states), vaex.from_pandas(trans)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = pickup_gE(states, trans, {}, engine=engine)
    if engine == "vaex":
        out = out.to_pandas_df()
    out = out.sort_values("i_upper")

    np.testing.assert_allclose(out["nu_lines"], [1431.977929, 1368.700414])
    assert any("differ from the states-file" in str(x.message) for x in w) == (
        with_trans_wavenumber
    )
