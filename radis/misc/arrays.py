# -*- coding: utf-8 -*-
"""
Description
------------

Functions to deal with numpy arrays:

- :py:func:`~radis.misc.arrays.norm`
- :py:func:`~radis.misc.arrays.norm_on`
- :py:func:`~radis.misc.arrays.scale_to`
- :py:func:`~radis.misc.arrays.array_allclose`
- :py:func:`~radis.misc.arrays.nantrapz`
- :py:func:`~radis.misc.arrays.arange_len`
- :py:func:`~radis.misc.arrays.calc_diff`
- :py:func:`~radis.misc.arrays.find_nearest`
- :py:func:`~radis.misc.arrays.find_first`
- :py:func:`~radis.misc.arrays.autoturn`
- :py:func:`~radis.misc.arrays.centered_diff`
- :py:func:`~radis.misc.arrays.evenly_distributed`
- :py:func:`~radis.misc.arrays.anynan`
- :py:func:`~radis.misc.arrays.first_nonnan_index`
- :py:func:`~radis.misc.arrays.last_nonan_index`
- :py:func:`~radis.misc.arrays.is_sorted`
- :py:func:`~radis.misc.arrays.is_sorted_backward`
- :py:func:`~radis.misc.arrays.bining`
- :py:func:`~radis.misc.arrays.count_nans`
- :py:func:`~radis.misc.arrays.logspace`
- :py:func:`~radis.misc.arrays.numpy_add_at`





-------------------------------------------------------------------------------

"""


from math import ceil

import numba
import numpy as np
from numba import bool_, float64, int32, int64
from numpy import hstack
from scipy.interpolate import interp1d

# Normalize


def norm(a, normby=None, how="max"):
    """Normalize a numpy array with its maximum. Or normalize it with another
    vector. Works if array contains nans.

    Parameters
    ----------

    normby: array, or None
        if array, norm with this other array's maximum. If None, normalize with
        its maximum.
    """
    if normby is not None:
        normby = np.abs(normby)
    else:
        normby = np.abs(a)

    if how == "max":
        return a / np.nanmax(normby)
    elif how == "mean":
        return a / np.nanmean(normby)
    else:
        raise ValueError("Unknown normalisation method")


def norm_on(a, w, wmin=None, wmax=None, how="max"):
    """Normalize `a` on a specific range of `w`

    Parameters
    ----------

    a: array
        array
    w: array
        x-axis array

    Other Parameters
    ----------------

    wmin, wmax: float
        crop range
    how: 'mean', 'max'
        how to normalize

    Returns
    -------

    a_norm: array
        normalized array
    """
    imin = np.argmin((w - wmin) ** 2) if wmin else None
    imax = np.argmin((w - wmax) ** 2) if wmax else None
    if imin is not None and imax is not None:
        if imin > imax:
            imin, imax = imax, imin
        elif imin == imax:
            imax += 1
    return norm(a, a[imin:imax], how=how)


def scale_to(a, b, k=1):
    """Scale function a to k*b."""
    return a * k * max(np.abs(b)) / max(np.abs(a))


def array_allclose(a, b, rtol=1e-5, atol=1e-8, equal_nan=True):
    """Returns wheter a and b are all close (element wise). If not the same
    size, returns False (instead of crashing like the numpy version). Cf
    numpy.allclose docs for more information.

    Parameters
    ----------

    a, b: arrays
    rtol: float
    atol: float
    equal_nan: bool
        whether to consider Nan's as equal. Contrary to the numpy version this
        one is set to True by default
    """

    if len(a) != len(b):
        return False

    return np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)


def nantrapz(I, w, dx=1.0, axis=-1):
    """Returns :py:func:`~numpy.trapz` (I, w) discarding nan."""
    b = ~np.isnan(I)
    return np.trapz(I[b], w[b], dx=dx, axis=axis)


# %%
# ==============================================================================
# Numpy Function
# ==============================================================================


def arange_len(wmin, wmax, wstep) -> int:
    """Returns len of a :py:func:`numpy.arange` ``(wmin, max, wstep)`` array,
    accounting for floating point errors

    Note: :py:func:`numpy.arange` is useful to maintain the input ``wstep``.
    If you don't have this requirement, you better use :py:func:`numpy.linspace`
    directly.
    """
    return ceil((wmax - wmin) / wstep)


def calc_diff(t1, v1, t2, v2):
    """Substract two vectors that may have slightly offset abscisses
    interpolating the correct values.

    Parameters
    ----------

    t1, v1: array_like
        first vector and its abscisses
    t2, v2: array_like
        second vector and its abscisses


    Returns
    -------

    tdiff, vdiff: array_like
        substracted vector and its abscisses
    """

    t1, v1, t2, v2 = list(map(np.array, (t1, v1, t2, v2)))

    # Deal with inversely sorted case
    if t1[-1] < t1[0]:
        t1, v1 = t1[::-1], v1[::-1]
    if t2[-1] < t2[0]:
        t2, v2 = t2[::-1], v2[::-1]

    # Get the overlapping range
    b = np.logical_and(t2 > t1[0], t2 < t1[-1])

    tdiff = t2[b]
    v2 = v2[b]

    # Interpolate the correct values
    f = interp1d(t1, v1)
    v1 = f(tdiff)

    # Finally, substract:
    vdiff = v1 - v2

    return tdiff, vdiff


def find_nearest(array, searched, return_bool=False):
    """Return the closest elements in array for each element in 'searched'
    array. In case of multiple elements in `array` having equal difference with
    `searched` element, one with least index is returned. Also returns a
    boolean array with indices of elements occuring in output list set to true.

    Examples
    --------

    ::

        from numpy import array
        find_nearest(array([1,2,3,4]), array([2.1,2]), True)

        >>> (array([2, 2]), array([False, True, False, False]))

        find_nearest(np.array([1,2,3,4]), np.array([2.6,2]), True)

        >>> (array([3, 2]), array([False,  True,  True, False]))

        find_nearest(np.array([1, 3]), np.array([2]))

        >>> array([1])

        find_nearest(np.array([3, 1]), np.array([2]))

        >>> array([3])
    """
    if len(array) == 0:
        raise ValueError("Array to be searched cannot be empty")

    b = np.zeros_like(array, dtype=bool)

    def find_nearest(array, value):
        idx = (np.abs(array - value)).argmin()
        return idx, array[idx]

    nearest_els = []
    try:
        for s in searched:
            idx, el = find_nearest(array, s)
            b[idx] = True
            nearest_els.append(el)
    except:
        idx, el = find_nearest(array, searched)
        b[idx] = True
        nearest_els.append(el)

    if return_bool:
        out = nearest_els, b
    else:
        out = nearest_els

    return out


def find_first(arr, treshold):
    """Return the index of the first element of the array arr whose value is
    more than the treshold."""

    return np.argmax(arr > treshold)


def autoturn(data, key=-1):
    """Turns array data. key value:

    - ``0`` don't transpose
    - ``1`` : transpose
    - ``-1`` : auto : make sure the vectors are along the longest dimension
    """

    if key == 0:
        return data
    elif key == 1:
        return data.transpose()
    elif key == -1:
        if data.shape[0] == max(data.shape):
            # Probably columns :
            return data.transpose()
        else:
            return data


def centered_diff(w):
    """Return w[i+1]-w[i-1]/2, same size as w.

    Similar to :py:func:`numpy.diff`, but does not change the array
    size.
    """
    dw = np.diff(w)
    return (hstack((dw, dw[-1])) + hstack((dw[0], dw))) / 2


def evenly_distributed(w, tolerance=1e-5):
    """Make sure array `w` is evenly distributed.

    Parameters
    ----------

    w : numpy array
        array to test
    tolerance: float
        absolute tolerance

    Returns
    -------

    out: bool
        ``True`` or ``False`` if ``w`` is evenly distributed.
    """
    mean_step = np.diff(w).mean()
    return (np.abs((np.diff(w) - mean_step)) > tolerance).sum() == 0


def anynan(a):
    """Returns whether ``a`` has at least one :py:attr:`~numpy.nan`

    Fastest implementation for arrays with >10^4 elements
    https://stackoverflow.com/a/45011547/5622825
    """
    return np.isnan(np.dot(a, a))


@numba.njit
def first_nonnan_index(a):
    """Returns index of first non-nan value in ``a``

    Returns None is all values are :py:attr:`~numpy.nan`

    See Also
    --------

    :func:`~radis.misc.arrays.last_nonnan_index`
    """
    for i in range(a.size):
        if not np.isnan(a[i]):
            return i
    return None


@numba.njit
def last_nonnan_index(a):
    """Returns index of first non-nan value in ``a``

    Returns None is all values are :py:attr:`~numpy.nan`

    See Also
    --------

    :func:`~radis.misc.arrays.first_nonnan_index`
    """
    for i in range(a.size - 1, 0, -1):
        if not np.isnan(a[i]):
            return i
    return None


@numba.njit
def is_sorted(a):
    """Returns whether ``a`` is sorted in ascending order.

    From B.M. answer on StackOverflow: https://stackoverflow.com/a/47004533/5622825

    See Also
    --------

    :func:`~radis.misc.arrays.is_sorted_backward`
    """
    for i in range(a.size - 1):
        if a[i + 1] < a[i]:
            return False
    return True


@numba.njit
def is_sorted_backward(a):
    """Returns whether ``a`` is sorted in descending order.

    See Also
    --------

    :func:`~radis.misc.arrays.is_sorted`
    """
    for i in range(a.size - 1):
        if a[i + 1] > a[i]:
            return False
    return True


def bining(I, ymin=None, ymax=None, axis=1):
    """Averages a I multi-dimensional array (typically an image) along the y
    axis bining(I) corresponds to I.mean(axis=1) Nan are not taken into
    account.

    Parameters
    ----------

    I: numpy array
        intensity
    ymin: int [0-I.shape[1]]
        If None, 0 is used. Default ``None``.
    ymax: int [0-I.shape[1]]
        If None, I.shape[1] is used. Default ``None``.
    axis: int
        Default 1
    """
    I = np.array(I)  # convert to array in case it's a Pandas dataframe for instance
    if ymin is None:
        ymin = 0
    if ymax is None:
        ymax = I.shape[axis]
    if ymin < 0:
        print("Warning in bining. ymin ({0}) < 0".format(ymin))
    if ymax > I.shape[axis]:
        print(
            "Warning in bining. ymax ({0}) > yaxis length ({1})".format(
                ymax, I.shape[axis]
            )
        )
    return np.nanmean(I[:, ymin:ymax], axis=axis)


def count_nans(a):
    """Nan are good but only in India."""

    return np.isnan(a).sum()


def logspace(xmin, xmax, npoints):
    """Returns points from xmin to xmax regularly distributed on a logarithm
    space.

    Numpy's :py:func:`numpy.logspace` does the same from 10**xmin to 10**xmax
    """

    return np.logspace(np.log10(xmin), np.log10(xmax), npoints)


def numpy_add_at(LDM, k, l, m, I):
    """Add the linestrengths on the LDM grid.

    Uses the numpy implementation of :py:func:`~numpy.add.at`, which
    add arguments element-wise.

    Parameters
    ----------
    LDM : ndarray
        LDM grid
    k, l, m : array
        index
    I : array
        intensity to add

    Returns
    -------
    add: ndarray
        linestrengths distributed over the LDM grid

    Notes
    -----
    Cython version implemented in https://github.com/radis/radis/pull/234

    See Also
    --------
    :py:func:`numpy.add.at`

    """
    # print('Numpy add.at()')
    return np.add.at(LDM, (k, l, m), I)


## Cython functions:
try:
    import radis_cython_extensions as rcx
except (ModuleNotFoundError):
    add_at = numpy_add_at
    #  or use radis.misc.utils.NotInstalled() ?
# EP: Also got (after switching back to former Cython version) :
#     File "radis\cython\radis_cython_extensions.pyx", line 1, in init radis_cython_extensions
#     ValueError: numpy.ndarray size changed, may indicate binary incompatibility. Expected 88 from C header, got 80 from PyObject
# we may consider escaping this.
else:
    add_at = rcx.add_at


# @numba.njit
# def non_zero_values_around2(a, n):
#     """return a boolean array of same size as ``a`` where each position ``i``
#     is ``True`` if there are non-zero points less than ``n`` index position
#     away from ``a[i]``, and ``False`` if all points in ``a`` are 0 ``n``  index
#     position away from from ``a[i]``
#     """
#     # # Cleaner version, but about 2x slower on real-life test cases:
#     # #
#     # %timeit non_zero_values_around(LDM[:, l, m], truncation)
#     # 2.92 ms ± 99.3 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)

#     # %timeit non_zero_values_around2(LDM[:, l, m], truncation)
#     # 5.11 ms ± 110 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
#     # --------------

#     # non_zero_values_around2()
#     b = np.zeros(len(a)+2*n+1,  dtype=np.bool_)
#     for pos in np.nonzero(a)[0]:
#         b[pos:pos+2*n+1] = True
#     return b[n:len(a)+n]


@numba.njit(
    bool_[:](float64[:], int64),
    cache=True,
)
def non_zero_values_around(a, n):
    """return a boolean array of same size as ``a`` where each position ``i``
    is ``True`` if there are non-zero points less than ``n`` index position
    away from ``a[i]``, and ``False`` if all points in ``a`` are 0 ``n``  index
    position away from from ``a[i]``
    """

    # build the list
    L = []
    pos = -1  # found position of non-zero point
    start = -1  # start position of subrange

    for i, ai in enumerate(a):
        if ai != 0:
            pos = i
            if start == -1:
                if i == 0:
                    start = 0
                else:
                    start = i - n
        else:
            if start == -1:
                continue
            elif pos == i - n - 1:
                L.append((start, i - 1))
                start = -1
                pos = -1

    if start != -1:
        L.append((start, len(a) - 1))

    # turn it into a boolean array
    b = np.zeros(len(a), dtype=np.bool_)
    for start, end in L:
        b[max(0, start) : end + 1] = 1

    return b


# @numba.njit(
#     numba.types.List(numba.types.containers.UniTuple(int64, 2))(int32[:], int64),
#     cache=True,
# )
# def non_zero_ranges_around(a, n):
#     """return a list of coordinates corresponding to non-zero values
#     where each position ``i``
#     is ``True`` if there are non-zero points less than ``n`` index position
#     away from ``a[i]``, and ``False`` if all points in ``a`` are 0 ``n``  index
#     position away from from ``a[i]``

#     Parameters
#     ----------
#     a : array of int, line centers
#     n : int or array of int, width
#     """

#     # build the list
#     L = []
#     pos = -1  # found position of non-zero point
#     start = -1  # start position of subrange

#     for i, ai in enumerate(a):
#         if ai != 0:
#             pos = i
#             if start == -1:
#                 if i == 0:
#                     start = 0
#                 else:
#                     start = i - n
#         else:
#             if start == -1:
#                 continue
#             elif pos == i - n - 1:
#                 L.append((start, i - 1))
#                 start = -1
#                 pos = -1

#     if start != -1:
#         L.append((start, len(a) - 1))

#     return L


# @numba.njit(
#     bool_[:](int32[:], int64),
#     cache=True,
# )
# def non_zero_values_around(a, n):
#     """return a boolean array of same size as ``a`` where each position ``i``
#     is ``True`` if there are non-zero points less than ``n`` index position
#     away from ``a[i]``, and ``False`` if all points in ``a`` are 0 ``n``  index
#     position away from from ``a[i]``
#     """
#     L = non_zero_ranges_around(a, n)

#     # turn it into a boolean array
#     b = np.zeros(len(a), dtype=np.bool_)
#     for start, end in L:
#         b[max(0, start) : end + 1] = 1

#     return b


@numba.njit(
    (int64[:, :])(bool_[:]),
    cache=True,
)
def non_zero_ranges_in_array(b):
    """return a list of coordinates corresponding to non-zero values
    in boolean array ``b``

    Parameters
    ----------
    b : boolean array

    Examples
    --------
    ::

        b = np.array([0,0,1,1,0,1,0,1], dtype=bool)
        non_zero_ranges_in_array(b)
        >> ([2, 5, 7],
            [4, 6, 8])
    """

    # build the list
    L = []
    start = -1  # start position of subrange

    for i, bi in enumerate(b):
        if bi:
            if start == -1:
                start = i
        else:
            if start == -1:
                continue
            else:
                L.append((start, i))
                start = -1

    if start != -1:
        L.append((start, len(b)))

    return np.array(L)


@numba.njit(
    (bool_[:])(int64[:, :], int64),
    cache=True,
)
def boolean_array_from_ranges(ranges, n):
    """return a boolean array of length ``n`` where (``L[i][0]``, ``L[i][1]``)
    give the ranges set to ``True``

    Examples
    --------
    ::

        L = np.array([[2,4], [5,6], [7,8]])
        boolean_array_from_coordinates(*L.T, 8)
        >>> np.array([0, 0, 1, 1, 0, 1, 0, 1], dtype=bool)

    """
    # build the list
    b = np.zeros(n, dtype=bool_)
    for i in range(len(ranges)):
        b[ranges[i][0] : ranges[i][1]] = True

    return b


@numba.njit(
    numba.types.Tuple((int64[:, :], float64[:]))(
        int32[:], float64[:], float64[:], float64[:], int64, int64
    ),
    cache=True,
)
def sparse_add_at(ki0, Iv0, Iv1, weight, max_range, truncation_pts):
    """Returns (start, stop) of non-zero ranges, and intensity ``I`` (length
    ``truncation_pts``) for all line intensities ``Iv0`` and ``Iv1`` at position ``ki0``
    and ``ki0+1`` with weights ``weight``

    .. warning::
        ``ki0`` must be sorted.
    """
    assert is_sorted(ki0)
    I = np.zeros(max_range)
    L = [
        (0, 0) for k in range(0)
    ]  # empty list; helps Numba infer type. See https://numba.pydata.org/numba-doc/0.41.0/user/troubleshoot.html#my-code-has-an-untyped-list-problem
    n = truncation_pts

    # initiate start
    pos = ki0[0]
    start = max(pos - n, 0)
    end = max_range  # will be rewritten in the first call of the loop

    # Loop over all lines :
    for i in range(len(ki0)):
        pos = ki0[i]

        # Adjust range:
        if pos > end:
            # close the previous range:
            L.append((start, end))
            # reopen new range as `[pos-n, pos+n]`:
            start = max(pos - n, 0)
            end = min(
                pos + 2 + n, max_range
            )  # TODO:  +2 because line is added on ki0 and ki0+1
        else:
            # we're in the middle of a range:
            # keep the start, move the end
            end = min(pos + 2 + n, max_range)

        # Sum intensity  (equivalent of "add-at")
        I[pos] += Iv0[i] * weight[i]
        I[pos + 1] += Iv1[i] * weight[i]
    # final close:
    if len(L) == 0 or L[-1] != (start, end):
        L.append((start, end))

    return np.array(L), I


if __name__ == "__main__":
    import pytest

    pytest.main(["../test/misc/test_arrays.py", "-s"])  # -s for showing console output
