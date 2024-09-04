import pytest

from radis.api.exomolapi import check_code_level


@pytest.mark.parametrize("bdat_list", [["a1"], ["a0", "a1"], ["a1", "a0"]])
def test_check_bdat_a1(bdat_list):
    bdat = {}
    bdat["code"] = bdat_list
    assert check_code_level(bdat) == "a1"


@pytest.mark.parametrize("bdat_list", [["a0"]])
def test_check_bdat_a0(bdat_list):
    bdat = {}
    bdat["code"] = bdat_list
    assert check_code_level(bdat) == "a0"


@pytest.mark.parametrize("bdat_list", [["a0", "a1", "a2"]])
def test_check_bdat_no_code_level(bdat_list):
    """a2 is not a valid code level"""
    bdat = {}
    bdat["code"] = bdat_list
    assert check_code_level(bdat) == None
