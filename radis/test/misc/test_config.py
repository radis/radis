# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 15:38:18 2018

@author: erwan
"""


import pytest

from radis.misc.config import get_config, printDatabankList


@pytest.mark.fast
def test_json_config_file(*args, **kwargs):
    """Test that it's readable"""
    get_config()


@pytest.mark.fast
def test_databanks(*args, **kwargs):
    """Test radis.json is readable"""
    printDatabankList()


def _run_testcases(verbose=True, *args, **kwargs):
    """Test array functions"""

    test_json_config_file()
    test_databanks()

    return True


if __name__ == "__main__":
    print("test config files:", _run_testcases())
