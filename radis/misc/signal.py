# -*- coding: utf-8 -*-
"""

Summary
-------

Signal processing functions 


-------------------------------------------------------------------------------


"""

from __future__ import absolute_import, division, print_function, unicode_literals

from numpy import trapz, abs, linspace, isnan, nan
from scipy.interpolate import splev, splrep
from warnings import warn
from radis.misc.debug import printdbg
from radis.misc.arrays import is_sorted, is_sorted_backward


def resample(
    xspace,
    vector,
    xspace_new,
    k=1,
    ext="error",
    energy_threshold=1e-3,
    print_conservation=True,
):
    """ Resample (xspace, vector) on a new space (xspace_new) of evenly distributed
    data and whose bounds are taken as the same as `xspace`. 

    Uses spline interpolation to create the intermediary points. Number of points
    is the same as the initial xspace, times a resolution factor. Verifies energy 
    conservation on the intersecting range at the end.


    Parameters    
    ----------

    xspace: array
        space on which vector was generated 

    vector: array
        quantity to resample

    resfactor: array
        xspace vector to resample on

    k: int
        order of spline interpolation. 3: cubic, 1: linear. Default 1.

    ext: 'error', 'extrapolate', 0, 1 
        Controls the value returned for elements of xspace_new not in the interval 
        defined by xspace. If 'error', raise a ValueError. If 'extrapolate', well,
        extrapolate. If '0' or 0, then fill with 0. If 1, fills with 1. 
        Default 'error'.

    energy_threshold: float
        if energy conservation (integrals on the intersecting range) is above 
        this threshold, raise an error. If None, dont check for energy conservation
        Default 1e-3 (0.1%)

    print_conservation: boolean
        if True, prints energy conservation


    Returns
    -------

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
        raise ValueError("Resampling requires wavespace to be sorted. It is not!")

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

    if isnan(vector).sum() > 0:
        raise ValueError(
            "Resampled vector has {0} nans. Interpolation will fail".format(
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
    """ Resample (xspace, vector) on a new space (xspace_new) of evenly distributed
    data and whose bounds are taken as the same as `xspace`. 

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


# %% Test routines


def _test(verbose=True, debug=False, plot=True, warnings=True, *args, **kwargs):
    """ Test procedures


    Parameters    
    ----------

    debug: boolean
        swamps the console namespace with local variables. Default ``False``

    """

    from radis.test.utils import getTestFile
    from radis.phys.convert import nm2cm, cm2nm
    import matplotlib.pyplot as plt
    from numpy import loadtxt, linspace

    # Test even resampling

    w_nm, I_nm = loadtxt(getTestFile("spectrum.txt")).T
    w_cm, I_cm = resample_even(
        nm2cm(w_nm),
        I_nm,
        resfactor=2,
        energy_threshold=1e-3,
        print_conservation=verbose,
    )

    if plot:
        plt.figure()
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity")
        plt.plot(w_nm, I_nm, "-ok", label="original")
        plt.plot(cm2nm(w_cm), I_cm, "-or", label="resampled")
        plt.legend()

    # Test resampling

    w_crop = linspace(376, 381, 100)
    I_crop = resample(w_nm, I_nm, w_crop, energy_threshold=0.01)

    if plot:
        plt.figure()
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity")
        plt.plot(w_nm, I_nm, "-ok", label="original")
        plt.plot(w_crop, I_crop, "-or", label="resampled")
        plt.legend()

    if debug:
        globals().update(locals())

    if warnings:
        warn("Testing resampling: no quantitative tests defined yet")

    return True  # no standard tests yet


if __name__ == "__main__":
    _test(debug=False)
