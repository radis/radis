"""
Gaussian and Lorentzian Width Envelope Pre-computation
======================================================

Compute tight min/max bounds of Gaussian and Lorentzian half-widths
across all spectral lines, so that LDM grid dimensions can be determined
*before* per-line broadening widths are evaluated.

Algorithms moved from :py:mod:`radis.gpu.params` so both GPU and CPU
code paths can reuse them.
"""

import json

import numpy as np

#  Serialization helpers for HDF5 metadata


def envelope_to_metadata(g_envelope, l_envelope):
    """Serialize envelopes to a flat dict storable in HDF5 metadata.

    Gaussian bounds are stored as plain floats.  Lorentzian piecewise-linear
    segments are stored as JSON strings (arrays are not HDF5-attribute-safe).

    Parameters
    ----------
    g_envelope : tuple (float, float)
    l_envelope : tuple of two tuples ``((A, B, X), (A, B, X))``

    Returns
    -------
    dict
    """
    meta = {
        "envelope_g_min": g_envelope[0],
        "envelope_g_max": g_envelope[1],
    }
    l_min, l_max = l_envelope
    meta["envelope_l_min"] = json.dumps(
        {
            "A": l_min[0].tolist(),
            "B": l_min[1].tolist(),
            "X": l_min[2].tolist(),
        }
    )
    meta["envelope_l_max"] = json.dumps(
        {
            "A": l_max[0].tolist(),
            "B": l_max[1].tolist(),
            "X": l_max[2].tolist(),
        }
    )
    return meta


#  Gaussian envelope


def compute_gaussian_envelope(log_2vMm):
    """Compute the Gaussian width envelope from per-line ``log_2vMm`` values.

    Parameters
    ----------
    log_2vMm : array_like, shape (N_lines,)
        ``log(v0) + 0.5 * log(8 k ln2 / (c² Mm))`` for each line.

    Returns
    -------
    tuple (float, float)
        ``(log_2vMm_min, log_2vMm_max)``
    """
    return (float(np.min(log_2vMm)), float(np.max(log_2vMm)))


#  Lorentzian envelope


def compute_lorentzian_envelope(na, gamma_arr):
    """Compute piecewise-linear min/max envelopes for the Lorentzian width.

    The Lorentzian half-width in log-space is a piecewise-linear function
    of ``log(T_ref / T)`` (denoted ``log_rT`` elsewhere).  This function
    determines the piecewise segments via a convex-hull sweep over the
    ``(na, log(gamma))`` parameter space.

    Parameters
    ----------
    na : array_like, shape (N_lines,)
        Temperature-dependence exponent for each line.
    gamma_arr : array_like, shape (N_broadeners, N_lines)
        Broadening half-widths at reference conditions for each broadener
        (e.g. row 0 = self-broadening, row 1 = air-broadening).

    Returns
    -------
    tuple of two tuples
        ``((A_min, B_min, X_min), (A_max, B_max, X_max))`` where each
        inner tuple represents the piecewise-linear envelope segments.

        * ``A`` – slopes  (the ``na`` values at each breakpoint)
        * ``B`` – offsets  (``log(gamma)`` at each breakpoint)
        * ``X`` – breakpoints in ``log_rT`` space (first is ``-inf``,
          last is ``+inf``)

        To evaluate at a given ``log_rT``::

            # find segment i such that X[i-1] <= log_rT < X[i]
            log_wL = log_rT * A[i] + B[i] + log(2 * p)
    """
    result = []
    for minmax in (np.min, np.max):

        # Combine (na, gamma) per line, taking min or max across broadeners
        na_gamma_arr = np.zeros((na.size, 2), dtype=np.float32)
        na_gamma_arr[:, 0] = na
        na_gamma_arr[:, 1] = minmax(gamma_arr, axis=0)

        # Deduplicate lines with identical (na, gamma) pairs
        unique_lines = np.unique(na_gamma_arr.reshape(2 * na.size).view(np.int64))
        unique_lines = unique_lines.view(np.float32).reshape(unique_lines.size, 2)

        # For lines sharing the same na, keep only the extreme gamma
        test_dict = {}
        for na_i, gamma_i in unique_lines:
            try:
                test_dict[na_i] = minmax((gamma_i, test_dict[na_i]))
            except KeyError:
                test_dict[na_i] = gamma_i

        # Build the piecewise-linear envelope via a convex-hull sweep
        keys = sorted(test_dict.keys(), reverse=(minmax == np.min))
        A = [keys[0]]
        B = [np.log(test_dict[keys[0]])]
        X = [-np.inf]

        for key in keys[1:]:
            if key == A[-1]:
                # Same slope as the last segment — dominated, skip
                continue
            for i in range(len(X)):
                xi = (np.log(test_dict[key]) - B[i]) / (A[i] - key)
                if xi >= X[i]:
                    if i < len(X) - 1:
                        if xi < X[i + 1]:
                            break
                    else:
                        break

            while X[i] == xi:
                i -= 1

            A = A[: i + 1] + [key]
            B = B[: i + 1] + [np.log(test_dict[key])]
            X = X[: i + 1] + [xi]

        X = X[1:] + [np.inf]
        result.append((np.array(A), np.array(B), np.array(X)))

    return tuple(result)
