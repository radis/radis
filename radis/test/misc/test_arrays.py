# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 10:35:24 2018

@author: erwan
"""

# from time import perf_counter

import numpy as np
import pytest

# from radis import get_residual
# from radis.lbl.factory import SpectrumFactory
from radis.misc.arrays import (  # add_at,
    anynan,
    arange_len,
    array_allclose,
    autoturn,
    bining,
    calc_diff,
    centered_diff,
    count_nans,
    evenly_distributed,
    evenly_distributed_fast,
    find_first,
    find_nearest,
    first_nonnan_index,
    is_sorted,
    is_sorted_backward,
    last_nonnan_index,
    logspace,
    nantrapz,
    non_zero_values_around,
    norm,
    norm_on,
    numpy_add_at,
    scale_to,
)

# from radis.test.utils import setup_test_line_databases


@pytest.mark.fast
def test_norm(*args, **kwargs):
    """Test the norm() function for array normalization."""
    # Test basic max normalization
    a = np.array([1, 2, 3, 4, 5])
    result = norm(a)
    assert result[-1] == 1.0  # max value normalized to 1
    assert np.allclose(result, a / 5)

    # Test with negative values
    a_neg = np.array([-5, -2, 0, 2, 5])
    result_neg = norm(a_neg)
    assert result_neg[-1] == 1.0

    # Test normalization with another array (normby)
    a = np.array([1, 2, 3, 4, 5])
    normby = np.array([2, 4, 6, 8, 10])
    result = norm(a, normby=normby)
    assert np.allclose(result, a / 10)

    # Test mean normalization
    a = np.array([1, 2, 3, 4, 5])
    result_mean = norm(a, how="mean")
    assert np.allclose(result_mean, a / np.mean(np.abs(a)))

    # Test with NaN values
    a_nan = np.array([1, np.nan, 3, 4, 5])
    result_nan = norm(a_nan)
    assert result_nan[-1] == 1.0  # nanmax should ignore NaN

    # Test invalid method raises error
    with pytest.raises(ValueError):
        norm(a, how="invalid")


@pytest.mark.fast
def test_norm_on(*args, **kwargs):
    """Test the norm_on() function for normalization on a specific range."""
    a = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=float)
    w = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=float)

    # Test normalization on a range
    result = norm_on(a, w, wmin=2, wmax=5)
    assert np.max(np.abs(result)) > 0  # normalized array should have values

    # Test with no range specified (normalizes by full array max)
    result_full = norm_on(a, w)
    assert np.allclose(result_full, a / 10)

    # Test mean normalization on range
    result_mean = norm_on(a, w, wmin=2, wmax=5, how="mean")
    assert result_mean is not None


@pytest.mark.fast
def test_scale_to(*args, **kwargs):
    """Test the scale_to() function for scaling arrays."""
    a = np.array([1, 2, 3, 4, 5])
    b = np.array([10, 20, 30, 40, 50])

    # Test basic scaling (k=1)
    result = scale_to(a, b)
    assert np.max(np.abs(result)) == np.max(np.abs(b))

    # Test with scaling factor k=2
    result_k2 = scale_to(a, b, k=2)
    assert np.max(np.abs(result_k2)) == 2 * np.max(np.abs(b))

    # Test with negative values
    a_neg = np.array([-5, -2, 0, 2, 5])
    b_neg = np.array([-10, -5, 0, 5, 10])
    result_neg = scale_to(a_neg, b_neg)
    assert np.max(np.abs(result_neg)) == np.max(np.abs(b_neg))


@pytest.mark.fast
def test_array_allclose(*args, **kwargs):
    """Test the array_allclose() function."""
    a = np.array([1, 2, 3])
    b = np.array([1, 2, 3])
    c = np.array([1, 2])  # different size
    d = np.array([1, 2, 4])  # different values

    # Test equal arrays
    assert array_allclose(a, b)

    # Test different size arrays
    assert not array_allclose(a, c)

    # Test different values
    assert not array_allclose(a, d)

    # Test with tolerance
    e = np.array([1.0, 2.0, 3.0001])
    assert array_allclose(a, e, atol=1e-3)
    assert not array_allclose(a, e, atol=1e-5)

    # Test with NaN values (equal_nan=True by default)
    f = np.array([1, np.nan, 3])
    g = np.array([1, np.nan, 3])
    assert array_allclose(f, g)

    # Test with equal_nan=False
    assert not array_allclose(f, g, equal_nan=False)


@pytest.mark.fast
def test_nantrapz(*args, **kwargs):
    """Test the nantrapz() function for integration ignoring NaN."""
    # Simple case without NaN
    I = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    w = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    result = nantrapz(I, w)
    expected = np.trapezoid(I, w)
    assert np.isclose(result, expected)

    # Case with NaN values
    I_nan = np.array([1.0, np.nan, 3.0, 4.0, 5.0])
    w_nan = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    result_nan = nantrapz(I_nan, w_nan)
    # Should integrate only non-NaN values
    assert not np.isnan(result_nan)


@pytest.mark.fast
def test_evenly_distributed(*args, **kwargs):
    """Test the evenly_distributed() function."""
    # Evenly distributed array
    w_even = np.linspace(0, 10, 100)
    assert evenly_distributed(w_even)

    # Not evenly distributed array
    w_uneven = np.array([0, 1, 2, 4, 8, 16])
    assert not evenly_distributed(w_uneven)

    # Test with custom tolerance
    w_almost_even = np.linspace(0, 10, 100)
    w_almost_even[50] += 1e-4
    assert evenly_distributed(w_almost_even, atolerance=1e-3)
    assert not evenly_distributed(w_almost_even, atolerance=1e-5)


@pytest.mark.fast
def test_evenly_distributed_fast(*args, **kwargs):
    """Test the evenly_distributed_fast() function."""
    # Evenly distributed array
    w_even = np.linspace(0, 10, 100)
    assert evenly_distributed_fast(w_even)

    # Not evenly distributed (different first and last steps)
    w_uneven = np.array([0, 1, 2, 3, 5])  # last step is 2, first step is 1
    assert not evenly_distributed_fast(w_uneven)


@pytest.mark.fast
def test_anynan(*args, **kwargs):
    """Test the anynan() function."""
    # Array with NaN
    a_nan = np.array([1, 2, np.nan, 4, 5])
    assert anynan(a_nan)

    # Array without NaN
    a_no_nan = np.array([1, 2, 3, 4, 5])
    assert not anynan(a_no_nan)


@pytest.mark.fast
def test_count_nans(*args, **kwargs):
    """Test the count_nans() function."""
    # No NaN
    a_no_nan = np.array([1, 2, 3, 4, 5])
    assert count_nans(a_no_nan) == 0

    # Some NaN
    a_some_nan = np.array([1, np.nan, 3, np.nan, 5])
    assert count_nans(a_some_nan) == 2

    # All NaN
    a_all_nan = np.array([np.nan, np.nan, np.nan])
    assert count_nans(a_all_nan) == 3

    # Single NaN
    a_single_nan = np.array([np.nan])
    assert count_nans(a_single_nan) == 1


@pytest.mark.fast
def test_numpy_add_at(*args, **kwargs):
    """Test the numpy_add_at() function."""
    # Create a 3D array
    LDM = np.zeros((5, 3, 3), dtype=np.float64)
    k = np.array([0, 1, 2, 0])
    l = np.array([0, 1, 2, 0])
    m = np.array([0, 1, 2, 0])
    I = np.array([1.0, 2.0, 3.0, 4.0])

    numpy_add_at(LDM, k, l, m, I)

    # Check that values were added correctly
    assert LDM[0, 0, 0] == 5.0  # 1.0 + 4.0
    assert LDM[1, 1, 1] == 2.0
    assert LDM[2, 2, 2] == 3.0
    assert LDM[3, 0, 0] == 0.0  # untouched


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


# @pytest.mark.fast
# def test_cython_add_at(*args, **kwargs):
##    """
##    Compare the workings of the Cython compiled add_at() function
##    versus the numpy add.at() function with bogus data.
##    """
##
##    # First check if Python was able to import the Cython version of add at:
##    from radis_cython_extensions import add_at as cython_add_at
##
##    assert add_at == cython_add_at
##
##    # Compare output LDM's between the two additions:
##    Nv = 300000
##    NG = 4
##    NL = 16
##    Ni = 1000000
##
##    I = np.random.rand(Ni).astype(np.float32)
##    k = np.random.randint(Nv, size=Ni, dtype=np.int32)
##    l = np.random.randint(NG, size=Ni, dtype=np.int32)
##    m = np.random.randint(NL, size=Ni, dtype=np.int32)
##
##    LDM1 = np.zeros((Nv, NG, NL), dtype=np.float32)
##    t0 = perf_counter()
##    np.add.at(LDM1, (k, l, m), I)
##    print("Numpy add.at(): ", perf_counter() - t0)
##
##    LDM2 = np.zeros((Nv, NG, NL), dtype=np.float32)
##    t0 = perf_counter()
##    cython_add_at(LDM2, k, l, m, I)
##    print("Cython add_at(): ", perf_counter() - t0)
##
##    print("Residual: ", np.sum(np.abs(LDM1 - LDM2)))
##    assert np.allclose(LDM1, LDM2)
# assert True


# def test_cython_add_at_spectra(*args, **kwargs):
##    """
##    Test if the Cython add_at() produces the same spectra as
##    with numpy add.at().
##    """
##
##    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there
##
##    # Conditions
##    wstep = 0.005
##    wmin = 2100  # cm-1
##    wmax = 2200  # cm-1
##
##    T = 1200  # K
##    p = 0.1  # bar
##
##    sf = SpectrumFactory(
##        wavenum_min=wmin,
##        wavenum_max=wmax,
##        mole_fraction=1,
##        path_length=1,  # doesnt change anything
##        wstep=wstep,
##        pressure=p,
##        isotope="1",
##        verbose=False,
##        warnings={
##            "MissingSelfBroadeningWarning": "ignore",
##            "NegativeEnergiesWarning": "ignore",
##            "HighTemperatureWarning": "ignore",
##            "OutOfRangeLinesWarning": "ignore",
##            "GaussianBroadeningWarning": "ignore",
##            "CollisionalBroadeningWarning": "ignore",
##            "AccuracyWarning": "ignore",
##        },
##    )
##    sf.load_databank("HITRAN-CO-TEST")
##
##    sf.use_cython = False
##    s_numpy = sf.eq_spectrum(Tgas=T, name="numpy")
##    s_numpy.apply_slit(0.5, "nm")
##    assert sf.misc.add_at_used == "numpy"
##
##    sf.use_cython = True
##    s_cython = sf.eq_spectrum(Tgas=T, name="cython")
##    s_cython.apply_slit(0.5, "nm")
##    assert sf.misc.add_at_used == "cython"
##
##    res = get_residual(s_numpy, s_cython, "transmittance")
##    assert res < 2e-4
# assert True


def test_non_zero_values_around(*args, **kwargs):
    """return a boolean array of same size as ``a`` where each position ``i``
    is ``True`` if there are non-zero points less than ``n`` index position
    away from ``a[i]``, and ``False`` if all points in ``a`` are 0 ``n``  index
    position away from ``a[i]``
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
