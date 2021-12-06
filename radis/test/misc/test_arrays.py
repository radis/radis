# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 10:35:24 2018

@author: erwan
"""

from time import perf_counter

import numpy as np
import pytest

from radis import get_residual
from radis.lbl.factory import SpectrumFactory
from radis.misc.arrays import (
    add_at,
    arange_len,
    autoturn,
    bining,
    calc_diff,
    centered_diff,
    find_first,
    find_nearest,
    first_nonnan_index,
    is_sorted,
    is_sorted_backward,
    last_nonnan_index,
    logspace,
    non_zero_values_around,
)
from radis.test.utils import setup_test_line_databases


@pytest.mark.fast
def test_arange_len(*args, **kwargs):

    # Positive arrays
    wmin, wmax, wstep = (380, 700, 0.1)
    assert arange_len(wmin, wmax, wstep) == len(np.arange(wmin, wmax, wstep))
    wmin, wmax, wstep = (380, 700, 0.31)
    assert arange_len(wmin, wmax, wstep) == len(np.arange(wmin, wmax, wstep))

    # Negative arrays
    wmin, wmax, wstep = (-100, 0, 0.1)
    assert arange_len(wmin, wmax, wstep) == len(np.arange(wmin, wmax, wstep))

    wmin, wmax, wstep = (-100, 0, 0.31)
    assert arange_len(wmin, wmax, wstep) == len(np.arange(wmin, wmax, wstep))

    # Centered arrays
    wmin, wmax, wstep = (-100, 100, 0.1)
    assert arange_len(wmin, wmax, wstep) == len(np.arange(wmin, wmax, wstep))
    wmin, wmax, wstep = (-100, 100, 0.31)
    assert arange_len(wmin, wmax, wstep) == len(np.arange(wmin, wmax, wstep))

    wmin, wmax = (-100, 100)
    for wstep in np.random.rand(100):
        assert arange_len(wmin, wmax, wstep) == len(np.arange(wmin, wmax, wstep))


@pytest.mark.fast
def test_is_sorted(*args, **kwargs):

    a = np.arange(10)

    assert is_sorted(a)
    assert is_sorted_backward(a[::-1])

    assert not is_sorted_backward(a)
    assert not is_sorted(a[::-1])


@pytest.mark.fast
def test_nonnan_index(*args, **kwargs):

    a = np.arange(1000) * 0.2

    assert first_nonnan_index(a) == 0  # is None
    assert last_nonnan_index(a) == 999  # is None

    a[:10] = np.nan
    a[-60:] = np.nan

    assert first_nonnan_index(a) == 10
    assert last_nonnan_index(a) == 939

    a[:] = np.nan

    assert first_nonnan_index(a) == None
    assert last_nonnan_index(a) == None  # len(a)-1


@pytest.mark.fast
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


@pytest.mark.fast
def test_bining(*args, **kwargs):

    a = np.arange(20).reshape(4, 5)
    assert (bining(a) == np.array([2, 7, 12, 17])).all()
    assert (bining(a, ymin=1) == np.array([2.5, 7.5, 12.5, 17.5])).all()
    assert (bining(a, ymax=3) == np.array([1, 6, 11, 16])).all()
    assert (bining(a, ymin=1, ymax=3) == np.array([1.5, 6.5, 11.5, 16.5])).all()


@pytest.mark.fast
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


@pytest.mark.fast
def test_autoturn(*args, **kwargs):
    dat = np.arange(20)
    dat.resize(2, 10)
    dat_rot = np.transpose(dat)

    assert (autoturn(dat, key=0) == dat).all()
    assert (autoturn(dat, key=1) == dat_rot).all()
    assert (autoturn(dat) == dat).all()


@pytest.mark.fast
def test_centered_diff(*args, **kwargs):
    a = np.arange(10)
    ones = np.ones_like(a)
    zeros = np.zeros_like(a)

    assert (centered_diff(a) == ones).all()
    assert not (centered_diff(a) == zeros).all()
    assert (centered_diff(ones) == zeros).all()
    assert len(centered_diff(a)) == len(a)


@pytest.mark.fast
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


@pytest.mark.fast
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


@pytest.mark.fast
def test_cython_add_at(*args, **kwargs):
    """
    Compare the workings of the Cython compiled add_at() function
    versus the numpy add.at() function with bogus data.
    """

    # First check if Python was able to import the Cython version of add at:
    from radis_cython_extensions import add_at as cython_add_at

    assert add_at == cython_add_at

    # Compare output LDM's between the two additions:
    Nv = 300000
    NG = 4
    NL = 16
    Ni = 1000000

    I = np.random.rand(Ni).astype(np.float32)
    k = np.random.randint(Nv, size=Ni, dtype=np.int32)
    l = np.random.randint(NG, size=Ni, dtype=np.int32)
    m = np.random.randint(NL, size=Ni, dtype=np.int32)

    LDM1 = np.zeros((Nv, NG, NL), dtype=np.float32)
    t0 = perf_counter()
    np.add.at(LDM1, (k, l, m), I)
    print("Numpy add.at(): ", perf_counter() - t0)

    LDM2 = np.zeros((Nv, NG, NL), dtype=np.float32)
    t0 = perf_counter()
    cython_add_at(LDM2, k, l, m, I)
    print("Cython add_at(): ", perf_counter() - t0)

    print("Residual: ", np.sum(np.abs(LDM1 - LDM2)))
    assert np.allclose(LDM1, LDM2)


def test_cython_add_at_spectra(*args, **kwargs):
    """
    Test if the Cython add_at() produces the same spectra as
    with numpy add.at().
    """

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    # Conditions
    wstep = 0.005
    wmin = 2100  # cm-1
    wmax = 2200  # cm-1

    T = 1200  # K
    p = 0.1  # bar

    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=1,
        path_length=1,  # doesnt change anything
        wstep=wstep,
        pressure=p,
        isotope="1",
        verbose=False,
        warnings={
            "MissingSelfBroadeningWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "HighTemperatureWarning": "ignore",
            "OutOfRangeLinesWarning": "ignore",
            "GaussianBroadeningWarning": "ignore",
            "CollisionalBroadeningWarning": "ignore",
            "AccuracyWarning": "ignore",
        },
    )
    sf.load_databank("HITRAN-CO-TEST")

    sf.use_cython = False
    s_numpy = sf.eq_spectrum(Tgas=T, name="numpy")
    s_numpy.apply_slit(0.5, "nm")
    assert sf.misc.add_at_used == "numpy"

    sf.use_cython = True
    s_cython = sf.eq_spectrum(Tgas=T, name="cython")
    s_cython.apply_slit(0.5, "nm")
    assert sf.misc.add_at_used == "cython"

    res = get_residual(s_numpy, s_cython, "transmittance")
    assert res < 2e-4


def test_non_zero_values_around(*args, **kwargs):
    """return a boolean array of same size as ``a`` where each position ``i``
    is ``True`` if there are non-zero points less than ``n`` index position
    away from ``a[i]``, and ``False`` if all points in ``a`` are 0 ``n``  index
    position away from from ``a[i]``
    """

    a = np.array(
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0], dtype=np.float64
    )

    # n = 1
    b = np.array([0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1])
    out = non_zero_values_around(a, 1)
    print(a, "width 1")
    print(b)
    print(np.array(out, dtype=int))
    assert (out == b).all()

    # n = 2
    b = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    out = non_zero_values_around(a, 2)
    print(a, "width 2")
    print(b)
    print(np.array(out, dtype=int))
    assert (out == b).all()

    # n = 2
    a = np.array(
        [0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype=np.float64
    )
    b = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1])
    out = non_zero_values_around(a, 2)
    print(a, "width 2")
    print(b)
    print(np.array(out, dtype=int))
    assert (out == b).all()

    # n = 1
    a = np.array(
        [0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype=np.float64
    )
    b = np.array([1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1])
    out = non_zero_values_around(a, 1)
    print(a, "width 1")
    print(b)
    print(np.array(out, dtype=int))
    assert (out == b).all()

    # Also test non_zero_ranges_in_array, boolean_array_from_ranges
    from radis.misc.arrays import boolean_array_from_ranges, non_zero_ranges_in_array

    b = np.array([0, 0, 1, 1, 0, 1, 0, 1], dtype=bool)
    assert (non_zero_ranges_in_array(b) == np.array([(2, 4), (5, 6), (7, 8)])).all()

    L = np.array([[2, 4], [5, 6], [7, 8]], dtype=np.int64)
    assert (
        boolean_array_from_ranges(L, 8)
        == np.array([0, 0, 1, 1, 0, 1, 0, 1], dtype=bool)
    ).all()


if __name__ == "__main__":

    pytest.main(["test_arrays.py", "-s"])  # -s for showing console output
