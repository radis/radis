# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 10:35:24 2018

@author: erwan
"""

from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from radis.misc.arrays import is_sorted, is_sorted_backward


def test_is_sorted(*args, **kwargs):

    a = np.arange(10)

    assert is_sorted(a)
    assert is_sorted_backward(a[::-1])

    assert not is_sorted_backward(a)
    assert not is_sorted(a[::-1])


def _run_testcases(verbose=True, *args, **kwargs):
    """ Test array functions

    """

    test_is_sorted()

    return True


if __name__ == "__main__":
    print("test_arrays functions:", _run_testcases())
