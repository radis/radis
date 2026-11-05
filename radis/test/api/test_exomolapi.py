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
