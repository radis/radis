# -*- coding: utf-8 -*-
"""

Summary
-------

Signal processing functions

Resampling  & smoothing.


-------------------------------------------------------------------------------


"""

import numpy as np
from numpy import abs, isnan, linspace, nan, trapz, zeros_like
from scipy.interpolate import splev, splrep
from scipy.linalg import solveh_banded

from radis.misc.arrays import (
    anynan,
    first_nonnan_index,
    is_sorted,
    is_sorted_backward,
    last_nonnan_index,
)
from radis.misc.debug import printdbg


def resample(
    xspace,
    vector,
    xspace_new,
    k=1,
    ext="error",
    energy_threshold=5e-3,
    print_conservation=True,
):
    """Resample (xspace, vector) on a new space (xspace_new) of evenly
    distributed data and whose bounds are taken as the same as `xspace`.

    Uses spline interpolation to create the intermediary points. Number of points
    is the same as the initial xspace, times a resolution factor. Verifies energy
    conservation on the intersecting range at the end.


    Parameters
    ----------
    xspace: array
        space on which vector was generated
    vector: array
        quantity to resample
    xspace_new: array
        space on which to resample

    Other Parameters
    ----------------
    resfactor: array
        xspace vector to resample on
    k: int
        order of spline interpolation. 3: cubic, 1: linear. Default 1.
    ext: 'error', 'extrapolate', 0, 1
        Controls the value returned for elements of xspace_new not in the interval
        defined by xspace. If 'error', raise a ValueError. If 'extrapolate', well,
        extrapolate. If '0' or 0, then fill with 0. If 1, fills with 1.
        Default 'error'.
    energy_threshold: float or ``None``
        if energy conservation (integrals on the intersecting range) is above
        this threshold, raise an error. If ``None``, dont check for energy conservation
        Default 5e-3 (0.5%)
    print_conservation: boolean
        if True, prints energy conservation


    Returns
    -------
    array: resampled vector on evenly spaced array. Number of element is conserved.

    Notes
    -----
    Note that depending upon the from_space > to_space operation, sorting may
    be reversed.

    Examples
    --------
    Resample a :class:`~radis.spectrum.spectrum.Spectrum` radiance
    on an evenly spaced wavenumber space::

        w_nm, I_nm = s.get('radiance')
        w_cm, I_cm = resample_even(nm2cm(w_nm), I_nm)

    See Also
    --------
    :py:meth:`~radis.spectrum.spectrum.Spectrum.resample`
    """

    if len(xspace) != len(vector):
        raise ValueError(
            "vector and xspace should have the same length. "
            + "Got {0}, {1}".format(len(vector), len(xspace))
        )

    # Check reversed (interpolation requires objects are sorted)
    if is_sorted(xspace):
        reverse = False
    elif is_sorted_backward(xspace):
        reverse = True
    else:
        raise ValueError(
            "Resampling requires wavespace to be sorted. It is not! Use .sort()? "
        )

    if reverse:
        xspace = xspace[::-1]
        xspace_new = xspace_new[::-1]
        vector = vector[::-1]

    # translate ext in FITPACK syntax for splev
    if ext == "extrapolate":
        ext_fitpack = 0  # splev returns extrapolated value
    elif ext in [0, "0", 1, "1", nan, "nan"]:
        ext_fitpack = 1  # splev returns 0  (fixed in post-processing)
    elif ext == "error":
        ext_fitpack = 2  # splev raises ValueError
    else:
        raise ValueError("Unexpected value for `ext`: {0}".format(ext))

    # Handle nans in vector :
    # FITPACK will fail if there are nans in the vector being interpolated.
    # here we trim the sides; and raise an error if there are nans in the middle
    m = first_nonnan_index(vector)
    M = last_nonnan_index(vector)
    if m is None:
        raise ValueError("Vector being resample has only nans. Check your data")
    nanmask = False
    if m > 0 or M < len(vector) - 1:  # else, no change has to be made
        nanmask = zeros_like(vector, dtype=bool)
        nanmask[m : M + 1] = True
        xspace = xspace[nanmask]
        vector = vector[nanmask]
    if anynan(vector):
        raise ValueError(
            "Vector being resample has {0} nans. Interpolation will fail".format(
                isnan(vector).sum()
            )
        )

    # Resample the slit function on the spectrum grid
    try:
        tck = splrep(xspace, vector, k=k)
    except ValueError:
        # Probably error on input data. Print it before crashing.
        print("\nValueError - Input data below:")
        print("-" * 5)
        print(xspace)
        print(vector)
        print("Check plot 101 too")
        import matplotlib.pyplot as plt

        plt.figure(101).clear()
        plt.plot(xspace, vector)
        plt.xlabel("xspace")
        plt.ylabel("vector")
        plt.title("ValueError")
        raise
    vector_new = splev(xspace_new, tck, ext=ext_fitpack)

    # ... get masks
    b = (xspace >= xspace_new.min()) * (xspace <= xspace_new.max())
    b_new = (xspace_new >= xspace.min()) * (xspace_new <= xspace.max())

    # fix filling for out of boundary values
    if ext in [1, "1"]:
        vector_new[~b_new] = 1
        if __debug__:
            printdbg(
                "Filling with 1 on w<{0}, w>{1} ({2} points)".format(
                    xspace.min(), xspace.max(), (1 - b_new).sum()
                )
            )
    elif ext in [nan, "nan"]:
        vector_new[~b_new] = nan
        if __debug__:
            printdbg(
                "Filling with nans on w<{0}, w>{1} ({2} points)".format(
                    xspace.min(), xspace.max(), (1 - b_new).sum()
                )
            )

    # Check energy conservation:

    # ... calculate energy
    energy0 = abs(trapz(vector[b], x=xspace[b]))
    energy_new = abs(trapz(vector_new[b_new], x=xspace_new[b_new]))
    if energy_new == 0:  # deal with particular case of energy = 0
        if energy0 == 0:
            energy_ratio = 1
        else:
            energy_ratio = 0
    else:  # general case
        energy_ratio = energy0 / energy_new
    if energy_threshold:
        if abs(energy_ratio - 1) > energy_threshold:
            import matplotlib.pyplot as plt

            plt.figure(101).clear()
            plt.plot(xspace, vector, "-ok", label="original")
            plt.plot(xspace_new, vector_new, "-or", label="resampled")
            plt.xlabel("xspace")
            plt.ylabel("vector")
            plt.legend()
            raise ValueError(
                "Error in resampling: "
                + "energy conservation ({0:.5g}%) below tolerance level ({1:.5g}%)".format(
                    (1 - energy_ratio) * 100, energy_threshold * 100
                )
                + ". Check graph 101. "
                + "Increasing energy_threshold is possible but not recommended"
            )
    if print_conservation:
        print("Resampling - Energy conservation: {0:.5g}%".format(energy_ratio * 100))

    # Reverse again
    if reverse:
        # xspace_new = xspace_new[::-1]
        vector_new = vector_new[::-1]

    return vector_new


def resample_even(
    xspace,
    vector,
    resfactor=2,
    k=1,
    ext="error",
    energy_threshold=1e-3,
    print_conservation=True,
):
    """Resample (xspace, vector) on a new space (xspace_new) of evenly
    distributed data and whose bounds are taken as the same as `xspace`.

    Uses spline interpolation to create the intermediary points. Number of points
    is the same as the initial xspace, times a resolution factor. Verifies energy
    conservation at the end.


    Parameters
    ----------

    xspace: array
        space on which vector was generated

    vector: array
        quantity to resample

    resfactor: float
        increase of resolution. If 1, output vector has the same number of
        points as the input vector. Default 2.

    k: int
        order of spline interpolation. 3: cubic, 1: linear. Default 1.

    ext: 'error', 'extrapolate', 0
        Controls the value returned for elements of xspace_new not in the interval
        defined by xspace. If 'error', raise a ValueError. If 'extrapolate', well,
        extrapolate. If '0' or 0, then fill with 0. Default 'error'.

    energy_threshold: float
        if energy conservation (integrals) is above this threshold, raise an
        error

    print_conservation: boolean
        if True, prints energy conservation


    Returns
    -------

    xspace_new: array
        evenly spaced mapping of xspace (same min, same max)

    vector_new: array
        resampled vector on evenly spaced array. Number of element is conserved.

    Note that depending upon the from_space > to_space operation, sorting may
    be reversed.


    Examples
    --------

    Resample a :class:`~radis.spectrum.spectrum.Spectrum` radiance
    on an evenly spaced wavenumber space::

        w_nm, I_nm = s.get('radiance')
        w_cm, I_cm = resample_even(nm2cm(w_nm), I_nm)
    """

    # Get new evenly space array
    xspace_new = linspace(xspace[0], xspace[-1], len(vector) * resfactor)

    # Resample
    vector_new = resample(
        xspace,
        vector,
        xspace_new,
        k=k,
        ext=ext,
        energy_threshold=energy_threshold,
        print_conservation=print_conservation,
    )

    return xspace_new, vector_new


def als_baseline(
    intensities,
    asymmetry_param=0.05,
    smoothness_param=1e6,
    max_iters=10,
    conv_thresh=1e-5,
    verbose=False,
):
    """Computes the asymmetric least squares baseline

    Parameters
    ----------
    intensities: array_like
        vector to smooth
    asymmetry_param: float
        value will shift the baseline fit below or above your average line. To cancel gaussian noise,
        you would want the value to be ``0.5``. To remove a positive peak would want it to be close
        to ``0``. To remove a negative peaks (from a transmittance signal for instance) you would
        want it to be close to ``1``. Default ``0.05``
    smoothness_param: float
        Relative importance of smoothness of the predicted response:
        the higher, the smoother the baseline. Suggested range: ``1e2`` to ``1e8``. Default ``1e6``

    Other Parameters
    ----------------
    max_iters: int
        number of iterations
    conv_thresh: float
        convergence
    verbose: boolean


    Examples
    --------
    ::

        from neq.math.smooth import als_baseline
        I = als_baseline(I, smoothness_param=1e5, aymmetry_param=0.1)

    References
    ----------

    * https://zanran_storage.s3.amazonaws.com/www.science.uva.nl/ContentPages/443199618.pdf
    * http://www.science.uva.nl/~hboelens/publications/draftpub/Eilers_2005.pdf

    Notes
    -----

    Implementation:

    - uses :func:`~neq.math.smooth.WhittakerSmoother` internally
    - kudos to https://gist.github.com/perimosocordiae

    See Also
    --------

    :func:`~radis.misc.signal.WhittakerSmoother`

    """
    smoother = WhittakerSmoother(intensities, smoothness_param, deriv_order=2)
    # Rename p for concision.
    p = asymmetry_param
    # Initialize weights.
    w = np.ones(intensities.shape[0])
    for i in range(max_iters):
        z = smoother.smooth(w)
        mask = intensities > z
        new_w = p * mask + (1 - p) * (~mask)
        conv = np.linalg.norm(new_w - w)
        if verbose:
            print(i + 1, conv)
        if conv < conv_thresh:
            break
        w = new_w
    else:
        print("ALS did not converge in %d iterations" % max_iters)
    return z


class WhittakerSmoother(object):
    """

    References
    ----------
    - kudos to https://gist.github.com/perimosocordiae

    See Also
    --------

    :func:`~radis.misc.signal.als_baseline`
    """

    def __init__(self, signal, smoothness_param, deriv_order=1):
        self.y = signal
        assert deriv_order > 0, "deriv_order must be an int > 0"
        # Compute the fixed derivative of identity (D).
        d = np.zeros(deriv_order * 2 + 1, dtype=int)
        d[deriv_order] = 1
        d = np.diff(d, n=deriv_order)
        n = self.y.shape[0]
        k = len(d)
        s = float(smoothness_param)

        # Here be dragons: essentially we're faking a big banded matrix D,
        # doing s * D.T.dot(D) with it, then taking the upper triangular bands.
        diag_sums = np.vstack(
            [
                np.pad(s * np.cumsum(d[-i:] * d[:i]), ((k - i, 0),), "constant")
                for i in range(1, k + 1)
            ]
        )
        upper_bands = np.tile(diag_sums[:, -1:], n)
        upper_bands[:, :k] = diag_sums
        for i, ds in enumerate(diag_sums):
            upper_bands[i, -i - 1 :] = ds[::-1][: i + 1]
        self.upper_bands = upper_bands

    def smooth(self, w):
        foo = self.upper_bands.copy()
        foo[-1] += w  # last row is the diagonal

        return solveh_banded(foo, w * self.y, overwrite_ab=True, overwrite_b=True)


if __name__ == "__main__":
    from radis.test.spectrum.test_spectrum import (
        test_resampling_function,
        test_resampling_nan_function,
    )

    test_resampling_function()
    test_resampling_nan_function()
