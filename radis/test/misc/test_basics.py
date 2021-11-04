# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 15:38:18 2018

@author: erwan
"""


import pytest


@pytest.mark.fast
def test_flatten(*args, **kwargs):
    from radis.misc.basics import flatten

    assert flatten(*[600, 300, (600, 600, 700)]) == [600, 300, 600, 600, 700]


@pytest.mark.fast
def test_round_off(*args, **kwargs):
    from radis.misc.basics import round_off

    assert round_off(0.024231242) == 0.024
    assert round_off(0.0019932133) == 0.002
    assert round_off(0.0000000000002193) == 0


@pytest.mark.fast
def test_str2bool(*args, **kwargs):
    import pandas as pd

    from radis.misc.basics import str2bool

    df = pd.Series(["1.0", "true", "True", 1])
    assert df.dtypes == object

    # Test
    assert df.map(str2bool).all()

    df = pd.Series(["0.0", "false", "False", 0])
    assert df.dtypes == object

    # Test
    assert not df.map(str2bool).any()
