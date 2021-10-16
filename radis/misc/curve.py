# -*- coding: utf-8 -*-
"""Created on Mon Jan  8 16:29:18 2018.

@author: erwan

Operations on Curves, where a curve is a (w, I) array tuple

Similar to OriginPro's Simple Curve Math operators


Routine Listing
---------------

- :func:`~radis.misc.curve.curve_distance`
- :func:`~radis.misc.curve.curve_add`
- :func:`~radis.misc.curve.curve_substract`
- :func:`~radis.misc.curve.curve_multiply`
- :func:`~radis.misc.curve.curve_divide`


-------------------------------------------------------------------------------
"""

import warnings

import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist


def curve_distance(w1, I1, w2, I2, discard_out_of_bounds=True):
    r"""Get a regularized euclidian distance from curve (w1, I1) to curve (w2, I2)

    .. math::

        D(w_1)[i] = \sqrt{ \sum_j (\hat{I_1}[i]  - \hat{I_2}[j] )^2 + (\hat{w_1}[i] - \hat{w_2}[j])^2}

    Where values are normalized as:

    .. math::

        \hat{A} = \\frac{A}{max(A) - min(A)}

    This regularized Euclidian distance minimizes the effect of a small shift in between the two
    curves in case of stiff curves (like a spectrum bandhead can be)

    No interpolation needed neither.

    Distances for out of bounds values is set to nan

    .. warning ::

        This is a distance on both the waverange and the intensity axis.
        It may be used to compensate for a small offset in your experimental
        spectrum (due to wavelength calibration, for instance) but can lead
        to wrong fits easily. Plus, it is very cost-intensive!


    Parameters
    ----------
    w1, I1: array
        range and values for first curve
    w2, I2: array
        range and values for 2nd curve
    discard_out_of_bounds: boolean
        if True, distance for out of bound values is set to nan. Else, it
        will be the distance from the last point.


    Returns
    -------

    w1, Idist: array
        minimal distance from I1 to I2, for each point in (w1, I1)
    """

    norm_w1 = np.max(w1) - np.min(w1)
    norm_w2 = np.max(w2) - np.min(w2)
    norm_I1 = np.max(I1) - np.min(I1)
    norm_I2 = np.max(I2) - np.min(I2)

    dist = cdist(
        np.array((w1 / norm_w1, I1 / norm_I1)).T,
        np.array((w2 / norm_w2, I2 / norm_I2)).T,
    ).min(axis=1)

    # discard out of bound values
    if discard_out_of_bounds:
        b = np.logical_or(w1 < w2.min(), w1 > w2.max())
        dist[b] = np.nan

    return w1, dist


def curve_add(w1, I1, w2, I2, is_sorted=False, kind="linear"):
    """Add curve (w2, I2) from (w1, I1) Linearly interpolates if the two ranges
    dont match. Fills out of bound parameters with nan.

    Similar to OriginPro's "Simple Curve Math Substract"

    .. note::
        for higher accuracy, choose (w2, I2) has the curve with the highest
        resolution


    Parameters
    ----------
    w1, I1: array
        range and values for first curve
    w2, I2: array
        range and values for 2nd curve
    is_sorted: boolean
        (optional) if True, doesnt sort input arrays
    kind: str
        interpolation kind. Default 'linear'. See scipy.interpolate.interp1d


    Returns
    -------

    w1, Iadd: array
        sum ``I1 + I2`` interpolated on the first range ``w1``
    """

    I2_interp = _curve_interpolate(w1, I1, w2, I2, is_sorted=is_sorted, kind=kind)

    return w1, I1 + I2_interp


def curve_substract(w1, I1, w2, I2, is_sorted=False, kind="linear"):
    """Substracts curve (w2, I2) from (w1, I1) Linearly interpolates if the two
    ranges dont match. Fills out of bound parameters with nan.

    Similar to OriginPro's "Simple Curve Math Substract"

    .. note::
        for higher accuracy, choose (w2, I2) has the curve with the highest
        resolution


    Parameters
    ----------
    w1, I1: array
        range and values for first curve
    w2, I2: array
        range and values for 2nd curve
    is_sorted: boolean
        (optional) if True, doesnt sort input arrays
    kind: str
        interpolation kind. Default 'linear'. See scipy.interpolate.interp1d


    Returns
    -------

    w1, Idiff: array
        difference ``I1 - I2`` interpolated on the first range ``w1``
    """

    I2_interp = _curve_interpolate(w1, I1, w2, I2, is_sorted=is_sorted, kind=kind)

    return w1, I1 - I2_interp


def curve_multiply(w1, I1, w2, I2, is_sorted=False, kind="linear"):
    """Multiply curve (w2, I2) with (w1, I1) Linearly interpolates if the two
    ranges dont match. Fills out of bound parameters with nan.

    Similar to OriginPro's "Simple Curve Math Substract"

    .. note::
        for higher accuracy, choose (w2, I2) has the curve with the highest
        resolution


    Parameters
    ----------
    w1, I1: array
        range and values for first curve
    w2, I2: array
        range and values for 2nd curve
    is_sorted: boolean
        (optional) if True, doesnt sort input arrays
    kind: str
        interpolation kind. Default 'linear'. See scipy.interpolate.interp1d


    Returns
    -------
    w1, Iproduct: array
        product ``I1 * I2`` interpolated on the first range ``w1``
    """

    I2_interp = _curve_interpolate(w1, I1, w2, I2, is_sorted=is_sorted, kind=kind)

    return w1, I1 * I2_interp


def curve_divide(w1, I1, w2, I2, is_sorted=False, kind="linear", interpolation=1):
    """Divides curve (w1, I1) by (w2, I2) Linearly interpolates if the two
    ranges dont match. Fills out of bound parameters with nan.

    Similar to OriginPro's "Simple Curve Math Substract"

    .. note::
        for higher accuracy, choose (w2, I2) has the curve with the highest
        resolution


    Parameters
    ----------
    w1, I1: array
        range and values for first curve
    w2, I2: array
        range and values for 2nd curve

    Other Parameters
    ----------------
    is_sorted: boolean
        (optional) if ``True``, assumes that both input arrays are sorted already.
        Default ``False``.
    kind: str
        interpolation kind. Default 'linear'. See scipy.interpolate.interp1d
    interpolation: int, optional
        If 1, interpolate on w1, I1. Else, on w2, I2.
        Default 1

    Returns
    -------
    w1, Idiv: array
        Division ``I1 / I2`` interpolated on the first or second range according to reverseInterpolation
    """
    if interpolation == 1:
        I2_interp = _curve_interpolate(w1, I1, w2, I2, is_sorted=is_sorted, kind=kind)
        output = w1, I1 / I2_interp
    else:
        I1_interp = _curve_interpolate(w2, I2, w1, I1, is_sorted=is_sorted, kind=kind)
        output = w2, I1_interp / I2

    if np.isnan(output[1]).any():
        warnings.warn(
            "Presence of NaN in curve_divide!\nThink about interpolation=2", UserWarning
        )

    return output


def _curve_interpolate(w1, I1, w2, I2, is_sorted=False, kind="linear"):
    """Resample ``I2`` on ``w1``

    For the interpolation to work, this requires to first sort ``w1`` and ``w2``,
    interpolate, then revert to the initial ``w1``

    .. note::
        for higher accuracy, choose (w2, I2) has the curve with the highest
        resolution

    Returns
    -------

    I2_inter: np.array
        ``I2`` interpolated on ``w1``.
    """
    if (abs((w1[1] - w1[0])) < 0.5 * abs((w2[1] - w2[0]))) and (
        abs((w1[-1] - w1[-2])) < 0.5 * abs((w2[-1] - w2[-2]))
    ):
        # w1 looks like the spectra with highest accuracy
        warnings.warn(
            "First spectrum has more resolution than 2nd. Reverse your spectra in interpolation/comparison for a better accuracy"
        )

    # First sort both
    if not is_sorted:  # TODO : add first check with `radis.misc.is_sorted(a) `
        b = np.argsort(w1)
        w1, I1 = w1[b], I1[b]
        # Get reverse sorting array:
        bm1 = np.zeros_like(w1, dtype=np.int64)
        bm1[b] = np.arange(len(w1))

        b = np.argsort(w2)
        w2, I2 = w2[b], I2[b]

    # resample if needed
    if not np.array_equal(w1, w2):
        f = interp1d(w2, I2, kind=kind, bounds_error=False)
        I2 = f(w1)

    # revert to initial array
    if not is_sorted:
        out = I2[bm1]
    else:
        out = I2

    return out
