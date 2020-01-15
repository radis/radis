# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 15:38:18 2018

@author: erwan
"""


from __future__ import print_function, absolute_import, division, unicode_literals
import pytest


@pytest.mark.fast
def test_flatten(*args, **kwargs):
    from radis.misc.basics import flatten

    assert flatten(*[600, 300, (600, 600, 700)]) == [600, 300, 600, 600, 700]
