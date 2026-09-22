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
