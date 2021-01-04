# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 10:35:24 2018

@author: erwan
"""

import numpy as np

from radis.misc.arrays import (
    autoturn,
    bining,
    calc_diff,
    centered_diff,
    find_first,
    find_nearest,
    is_sorted,
    is_sorted_backward,
    logspace,
)


def test_is_sorted(*args, **kwargs):

    a = np.arange(10)

    assert is_sorted(a)
    assert is_sorted_backward(a[::-1])

    assert not is_sorted_backward(a)
    assert not is_sorted(a[::-1])


def test_find_first(*args, **kwargs):

    a = np.arange(10)
    assert find_first(a, -1) == 0
    assert find_first(a, 0) == 1
    assert find_first(a, 5) == 6
    assert find_first(a, 8) == 9
    assert find_first(a, 9) == 0
    assert find_first(a, 20) == 0

    assert not find_first(a, -1) == -1
    assert not find_first(a, 0) == 0
    assert not find_first(a, 5) == 5
    assert not find_first(a, 10) == 10


def test_bining(*args, **kwargs):

    a = np.arange(20).reshape(4, 5)
    assert (bining(a) == np.array([2, 7, 12, 17])).all()
    assert (bining(a, ymin=1) == np.array([2.5, 7.5, 12.5, 17.5])).all()
    assert (bining(a, ymax=3) == np.array([1, 6, 11, 16])).all()
    assert (bining(a, ymin=1, ymax=3) == np.array([1.5, 6.5, 11.5, 16.5])).all()


def test_calc_diff(*args, **kwargs):
    t1 = np.arange(5)
    t2 = np.arange(5)

    v1 = np.array([10, 12, 14, 16, 18])
    v2 = np.array([10, 12, 14, 16, 18])

    t_res1, v_res1 = calc_diff(t1, v1, t2, v2)
    t_res2, v_res2 = calc_diff(t1, v1, t1[::-1], v2)
    t_res3, v_res3 = calc_diff(t2[::-1], v1, t2, v2)

    assert (t_res1 == np.array([1, 2, 3])).all()
    assert (v_res1 == np.array([0, 0, 0])).all()

    assert (t_res2 == np.array([1, 2, 3])).all()
    assert (v_res2 == np.array([-4, 0, 4])).all()

    assert (t_res3 == np.array([1, 2, 3])).all()
    assert (v_res3 == np.array([4, 0, -4])).all()


def test_autoturn(*args, **kwargs):
    dat = np.arange(20)
    dat.resize(2, 10)
    dat_rot = np.transpose(dat)

    assert (autoturn(dat, key=0) == dat).all()
    assert (autoturn(dat, key=1) == dat_rot).all()
    assert (autoturn(dat) == dat).all()


def test_centered_diff(*args, **kwargs):
    a = np.arange(10)
    ones = np.ones_like(a)
    zeros = np.zeros_like(a)

    assert (centered_diff(a) == ones).all()
    assert not (centered_diff(a) == zeros).all()
    assert (centered_diff(ones) == zeros).all()
    assert len(centered_diff(a)) == len(a)


def test_logspace(*args, **kwargs):
    dat1 = logspace(1, 100, 10)
    dat2 = logspace(17, 250, 37)
    dat3 = logspace(5, 19, 3)

    dats = [dat1, dat2, dat3]

    assert (dats[0][0] - 1) <= 1e-6 and (dats[0][9] - 100) <= 1e-6
    assert (dats[1][0] - 17) <= 1e-6 and (dats[1][36] - 250) <= 1e-6
    assert (dats[2][0] - 5) <= 1e-6 and (dats[2][2] - 19) <= 1e-6

    for dat in dats:
        for i in range(2, len(dat)):
            assert (dat[i] / dat[i - 1] - dat[i - 1] / dat[i - 2]) <= 1e-6


def test_find_nearest(*args, **kwargs):
    a = np.arange(10)
    b = np.ones(5)

    x1, y1 = find_nearest(a, a, True)
    assert (x1 == a).all()
    assert len(y1) == len(a)
    assert (y1 == np.array([True for i in a])).all()

    x2, y2 = find_nearest(a, b, True)
    assert (x2 == b).all()
    assert len(x2) == len(b)
    assert len(y2) == len(a)
    assert (y2 == np.array([i == 1 for i in a])).all()

    x3, y3 = find_nearest(np.array([10]), np.array([-10, 0, 10]), True)
    assert (x3 == np.array([10, 10, 10])).all()
    assert len(y3) == 1
    assert y3 == [True]

    x4, y4 = find_nearest(np.array([1, 4, 3]), np.array([]), True)
    assert (x4 == np.array([])).all()
    assert len(y4) == 3
    assert (y4 == [False for _ in range(3)]).all()

    assert (find_nearest(a, a[::-1]) == a[::-1]).all()
    assert (
        find_nearest(np.array([1.5, 2.0, 2.5, 3.0]), np.array([-10, 10, 2.25]))
        == np.array([1.5, 3.0, 2.0])
    ).all()

    assert (find_nearest(np.array([1, 3]), np.array([2])) == np.array([1])).all()
    assert (find_nearest(np.array([3, 1]), np.array([2])) == np.array([3])).all()


if __name__ == "__main__":
    import pytest

    pytest.main(["test_arrays.py", "-s"])  # -s for showing console output
