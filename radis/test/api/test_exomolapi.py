import pandas as pd
import pytest

from radis.api.exomolapi import check_code_level


@pytest.mark.fast
def test_check_bdat_a1():
    for code in [["a1"], ["a0", "a1"], ["a1", "a0"]]:
        bdat = pd.DataFrame()
        bdat["code"] = code
        assert check_code_level(bdat) == "a1"


@pytest.mark.fast
def test_check_bdat_a0():
    bdat = pd.DataFrame()
    bdat["code"] = ["a0"]
    assert check_code_level(bdat) == "a0"


@pytest.mark.fast
def test_check_bdat_no_code_level():
    """a2 is not a valid code level"""
    bdat = pd.DataFrame()
    bdat["code"] = ["a0", "a1", "a2"]
    assert check_code_level(bdat) == None
