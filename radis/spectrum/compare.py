# -*- coding: utf-8 -*-
"""
Summary
-------

Functions to compare and plot comparison of :class:`~radis.spectrum.spectrum.Spectrum`
objects


Routine Listings
----------------

- :func:`~radis.spectrum.compare.get_diff`
- :func:`~radis.spectrum.compare.get_distance`
- :func:`~radis.spectrum.compare.get_ratio`
- :func:`~radis.spectrum.compare.get_residual`
- :func:`~radis.spectrum.compare.get_residual_integral`
- :func:`~radis.spectrum.compare.plot_diff`


-------------------------------------------------------------------------------



"""

from warnings import warn

import numpy as np

from radis.misc.arrays import array_allclose
from radis.misc.basics import compare_dict, compare_lists
from radis.misc.curve import curve_distance, curve_divide, curve_substract
from radis.spectrum.spectrum import Spectrum
from radis.spectrum.utils import cast_waveunit, format_xlabel, make_up, make_up_unit

# %% ======================================================================
# External functions
# ----------------
# XXX =====================================================================


def get_diff(
    s1: Spectrum,
    s2: Spectrum,
    var,
    wunit="default",
    Iunit="default",
    resample=True,
    diff_window=0,
):
    # type: (Spectrum, Spectrum, str, str, str, str, bool, int) -> np.array, np.array
    """Get the difference between 2 spectra.

    Basically returns w1, I1 - I2 where (w1, I1) and (w2, I2) are the values of
    s1 and s2 for variable var. (w2, I2) is linearly interpolated if needed.

        .. math::

            dI = I_1 - I_2


    Parameters
    ----------
    s1, s2: Spectrum objects
        2 spectra to compare.
    var: str
        spectral quantity (ex: ``'radiance'``, ``'transmittance'``...)
    wunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
        waveunit to compare in: wavelength air, wavenumber, wavelength vacuum
    Iunit: str
        if ``'default'`` use s1 unit for variable var
    medium: 'air', 'vacuum', default'
        propagating medium to compare in (if in wavelength)

    Other Parameters
    ----------------
    resample: bool
        if not ``True``, wavelength must be equals. Else, resample ``s2`` on
        ``s1`` if needed.
    diff_window: int
        If non 0, calculates diff by offsetting s1 by ``diff_window`` number of
        units on either side, and returns the minimum. Compensates for experimental
        errors on the w axis. Default 0. (look up code for more details...)

    Returns
    -------
    w1, Idiff: array
        difference interpolated on the wavespace range of the first Spectrum

    Notes
    -----

    Uses :func:`~radis.misc.curve.curve_substract` internally

    See Also
    --------

    :func:`~radis.spectrum.compare.get_ratio`,
    :func:`~radis.spectrum.compare.get_distance`,
    :func:`~radis.spectrum.compare.get_residual`,
    :func:`~radis.spectrum.compare.get_residual_integral`,
    :func:`~radis.spectrum.compare.plot_diff`,
    :meth:`~radis.spectrum.spectrum.compare_with`
    """
    # TODO : ensure objects are sorted and if so; forward to curve_substract()

    var, wunit, Iunit = get_default_units(s1, s2, var=var, wunit=wunit, Iunit=Iunit)

    # Get data
    # ----
    w1, I1 = s1.get(var, wunit=wunit, Iunit=Iunit, trim_nan=True)
    w2, I2 = s2.get(var, wunit=wunit, Iunit=Iunit, trim_nan=True)

    if not resample:
        if not array_allclose(w1, w2):
            raise AssertionError("Wavespace are not the same: use Spectrum.resample()")

    # Compute : I1 - I2 (on same range)
    w1, Idiff = curve_substract(w1, I1, w2, I2)
    if diff_window:
        # allow fluctuation from diff_window units. Kinda compensates
        # for experimental errors on x-axis
        diff_list = []
        I2_interp = I1 - Idiff
        for i in range(-diff_window, diff_window + 1):
            diff_list.append(I2_interp - np.roll(I1, i))
        # get minimum in abs value
        diff_list = np.array(diff_list)
        b = np.abs(diff_list).argmin(axis=0)
        Idiff = np.choose(b, diff_list)  # keeps sign
    return w1, Idiff


def get_ratio(
    s1: Spectrum, s2: Spectrum, var, wunit="default", Iunit="default", resample=True
):
    # type: (Spectrum, Spectrum, str, str, str, str, bool) -> np.array, np.array
    """Get the ratio between two spectra Basically returns w1, I1 / I2 where
    (w1, I1) and (w2, I2) are the values of s1 and s2 for variable var. (w2,
    I2) is linearly interpolated if needed.

        .. math::

            R = I_1 / I_2


    Parameters
    ----------
    s1, s2: Spectrum objects
        :py:class:`~radis.spectrum.spectrum.Spectrum`
    var: str
        spectral quantity
    wunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
        waveunit to compare in: wavelength air, wavenumber, wavelength vacuum
    Iunit: str
        if ``'default'`` use s1 unit for variable var

    Notes
    -----
    Uses :func:`~radis.misc.curve.curve_divide` internally

    See Also
    --------

    :func:`~radis.spectrum.compare.get_diff`,
    :func:`~radis.spectrum.compare.get_distance`,
    :func:`~radis.spectrum.compare.get_residual`,
    :func:`~radis.spectrum.compare.get_residual_integral`,
    :func:`~radis.spectrum.compare.plot_diff`,
    :meth:`~radis.spectrum.spectrum.compare_with`
    """

    var, wunit, Iunit = get_default_units(s1, s2, var=var, wunit=wunit, Iunit=Iunit)

    # Get data
    # ----
    w1, I1 = s1.get(var, wunit=wunit, Iunit=Iunit)
    w2, I2 = s2.get(var, wunit=wunit, Iunit=Iunit)

    if not resample:
        if not array_allclose(w1, w2):
            raise AssertionError("Wavespace are not the same: use Spectrum.resample()")

    # Compute: I1 / I2 (on same range)
    return curve_divide(w1, I1, w2, I2)


def get_distance(
    s1: Spectrum, s2: Spectrum, var, wunit="default", Iunit="default", resample=True
):
    # type: (Spectrum, Spectrum, str, str, str, str, bool) -> np.array, np.array
    """Get a regularized Euclidian distance between two spectra ``s1`` and
    ``s2``

    This regularized Euclidian distance minimizes the effect of a small shift in
    wavelength between the two spectra

    .. math::

        D(w_1)[i] = \sqrt{ \sum_j (\hat{I_1}[i]  - \hat{I_2}[j] )^2 + (\hat{w_1}[i] - \hat{w_2}[j])^2}

    Where values are normalized as:

    .. math::

        \hat{A} = \\frac{A}{max(A) - min(A)}

    If waveranges dont match, ``s2`` is interpolated over ``s1``.

    .. warning ::

        This is a distance on both the waverange and the intensity axis.
        It may be used to compensate for a small offset in your experimental
        spectrum (due to wavelength calibration, for instance) but can lead
        to wrong fits easily. Plus, it is very cost-intensive! Better use
        :func:`~radis.spectrum.compare.get_residual` for an automatized procedure.


    Parameters
    ----------
    s1, s2: Spectrum objects
        :py:class:`~radis.spectrum.spectrum.Spectrum`
    var: str
        spectral quantity
    wunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
        waveunit to compare in: wavelength air, wavenumber, wavelength vacuum
    Iunit: str
        if ``'default'`` use s1 unit for variable var
    medium: 'air', 'vacuum', default'
        propagating medium to compare in (if in wavelength)

    Notes
    -----
    Uses :func:`~radis.misc.curve.curve_distance` internally

    See Also
    --------
    :func:`~radis.spectrum.compare.get_diff`,
    :func:`~radis.spectrum.compare.get_ratio`,
    :func:`~radis.spectrum.compare.get_residual`,
    :func:`~radis.spectrum.compare.get_residual_integral`,
    :func:`~radis.spectrum.compare.plot_diff`,
    :meth:`~radis.spectrum.spectrum.compare_with`
    """
    # TODO: normalize with Imax, wmax

    var, wunit, Iunit = get_default_units(s1, s2, var=var, wunit=wunit, Iunit=Iunit)

    # Get data
    # ----
    w1, I1 = s1.get(var, wunit=wunit, Iunit=Iunit)
    w2, I2 = s2.get(var, wunit=wunit, Iunit=Iunit)

    if not resample:
        if not array_allclose(w1, w2):
            raise AssertionError("Wavespace are not the same: use Spectrum.resample()")

    # Compute: euclidian distance from w1, I1
    return curve_distance(w1, I1, w2, I2, discard_out_of_bounds=True)


def get_residual(
    s1: Spectrum,
    s2: Spectrum,
    var,
    norm="L2",
    ignore_nan=False,
    diff_window=0,
    normalize=False,
    normalize_how="max",
    wunit="default",
) -> float:
    # type: (Spectrum, Spectrum, str, bool, int) -> np.array, np.array
    """Returns L2 norm of ``s1`` and ``s2``

    For ``I1``, ``I2``, the values of variable ``var`` in ``s1`` and ``s2``,
    respectively, residual is calculated as:

    For ``L2`` norm:

        .. math::

            res = \\frac{\\sqrt{\\sum_i {(s_1[i]-s_2[i])^2}}}{N}.

    For ``L1`` norm:

        .. math::

            res = \\frac{\\sqrt{\\sum_i {|s_1[i]-s_2[i]|}}}{N}.


    Parameters
    ----------
    s1, s2: :class:`~radis.spectrum.spectrum.Spectrum` objects
        if not on the same range, ``s2`` is resampled on ``s1``.
    var: str
        spectral quantity
    norm: 'L2', 'L1'
        which norm to use

    Other Parameters
    ----------------
    ignore_nan: boolean
        if ``True``, ignore ``nan`` in the difference between s1 and s2 (ex: out of bound)
        when calculating residual. Default ``False``. Note: ``get_residual`` will still
        fail if there are ``nan`` in initial Spectrum.
    normalize: bool, or tuple
        if ``True``, normalize the two spectra before calculating the residual.
        If a tuple (ex: ``(4168, 4180)``), normalize on this range only. The unit
        is that of the first Spectrum by default (use ``wunit`` to change). Ex::

            s_exp   # in 'nm'
            s_calc  # in 'cm-1'
            get_residual(s_exp, s_calc, normalize=(4178, 4180))  # interpreted as 'nm'

    normalize_how: ``'max'``, ``'area'``, ``'mean'``
        how to normalize. ``'max'`` is the default but may not be suited for very
        noisy experimental spectra. ``'area'`` will normalize the integral to 1.
        ``'mean'`` will normalize by the mean amplitude value
    wunit: str
        used if normalized is a range

    Notes
    -----

    0 values for I1 yield nans except if I2 = I1 = 0

    when s1 and s2 dont have the size wavespace range, they are automatically
    resampled through get_diff on 's1' range

    Implementation of ``L2`` norm::

        np.sqrt((dI**2).sum())/len(dI)

    Implementation of ``L1`` norm::

        np.abs(dI).sum()/len(dI)

    See Also
    --------
    :func:`~radis.spectrum.compare.get_diff`,
    :func:`~radis.spectrum.compare.get_ratio`,
    :func:`~radis.spectrum.compare.get_distance`,
    :func:`~radis.spectrum.compare.plot_diff`,
    :func:`~radis.spectrum.compare.get_residual_integral`,
    :meth:`~radis.spectrum.spectrum.compare_with`
    """

    if normalize:
        if isinstance(normalize, tuple):
            wrange = normalize
        else:
            wrange = ()
        var, wunit, Iunit = get_default_units(s1, s2, var=var, wunit=wunit)
        s1 = s1.take(var).normalize(
            wrange=wrange, normalize_how=normalize_how, wunit=wunit
        )
        s2 = s2.take(var).normalize(
            wrange=wrange, normalize_how=normalize_how, wunit=wunit
        )

    var, wunit, Iunit = get_default_units(s1, s2, var=var, wunit=wunit)
    wdiff, dI = get_diff(
        s1,
        s2,
        var,
        resample=True,
        diff_window=diff_window,
        wunit=wunit,
        Iunit=Iunit,
    )

    if ignore_nan:
        b = np.isnan(dI)
        wdiff, dI = wdiff[~b], dI[~b]

    if norm == "L2":
        output = np.sqrt((dI ** 2).sum()) / len(dI)
    elif norm == "L1":
        output = (np.abs(dI)).sum() / len(dI)
    else:
        raise ValueError("unexpected value for norm")

    if np.isnan(output):
        warn(
            'NaN output in residual. You should use "ignore_nan=True". Read the help of get_residual.',
            UserWarning,
        )

    return output


def get_residual_integral(s1: Spectrum, s2: Spectrum, var, ignore_nan=False) -> float:
    # type: (Spectrum, Spectrum, str, bool) -> float
    """Returns integral of the difference between two spectra s1 and s2,
    relatively to the integral of spectrum s1.

    Compared to :func:`~radis.spectrum.compare.get_residual`, this tends to
    cancel the effect of the gaussian noise of an experimental spectrum.

    ::

        res = trapz(I2-I1, w1) / trapz(I1, w1)

    Note: when the considered variable is ``transmittance`` or ``transmittance_noslit``,
    the *upper* integral is used (up to 1) to normalize the integral difference

    ::

        res = trapz(I2-I1, w1) / trapz(1-I1, w1)

    Parameters
    ----------
    s1, s2: Spectrum objects
        :py:class:`~radis.spectrum.spectrum.Spectrum`
    var: str
        spectral quantity

    Other Parameters
    ----------------
    ignore_nan: boolean
        if ``True``, ignore nan in the difference between s1 and s2 (ex: out of bound)
        when calculating residual. Default ``False``. Note: ``get_residual_integral``
        will still fail if there are nan in initial Spectrum.

    Notes
    -----
    For I1, I2, the values of 'var' in s1 and s2, respectively, residual
    is calculated as::

        res = trapz(I2-I1, w1) / trapz(I1, w1)

    0 values for I1 yield nans except if I2 = I1 = 0

    when s1 and s2 dont have the size wavespace range, they are automatically
    resampled through get_diff on 's1' range


    See Also
    --------
    :func:`~radis.spectrum.compare.get_diff`,
    :func:`~radis.spectrum.compare.get_ratio`,
    :func:`~radis.spectrum.compare.get_distance`,
    :func:`~radis.spectrum.compare.get_residual`,
    :func:`~radis.spectrum.compare.plot_diff`,
    :meth:`~radis.spectrum.spectrum.compare_with`
    """
    var, wunit, Iunit = get_default_units(s1, s2, var=var)

    # Get data
    # ----
    w1, I1 = s1.get(var, wunit=wunit, Iunit=Iunit)
    wdiff, dI = get_diff(s1, s2, var, wunit=wunit, Iunit=Iunit, resample=True)

    # mask for 0
    if ignore_nan:
        b = np.isnan(dI)
        wdiff, dI = wdiff[~b], dI[~b]
        b = np.isnan(I1)
        w1, I1 = w1[~b], I1[~b]

    if var in ["transmittance", "transmittance_noslit"]:
        norm = 1 - np.trapz(I1, w1)
    else:
        norm = np.trapz(I1, w1)

    return np.abs(np.trapz(dI, wdiff)) / norm


def get_default_units(
    s1: Spectrum, s2: Spectrum, var=None, wunit="default", Iunit="default"
):
    """Get ``wunit``, Iunit for var; compatible with both spectra s1 and s2

    Parameters
    ----------
    s1, s2: Spectrum objects
    var: str, or None
        spectral quantity to plot (ex: ``'abscoeff'``). If None,
        plot the first one in the Spectrum from ``'radiance'``,
        ``'radiance_noslit'``, ``'transmittance'``, etc.
    wunit: ``'default'``, ``'nm'``, ``'cm-1'``, ``'nm_vac'``
        wavespace unit:  wavelength air, wavenumber, wavelength vacuum.
        If ``'default'``, use first spectrum wunit
    Iunit: str
        if ``'default'``, use first spectrum unit

    Returns
    -------
    var: str
        variable (evaluated if input was ``None``)
    wunit, Iunit: str
        units (evaluated if input was ``'default'``)
    """

    # Get defaults
    # ---
    if var is None:  # if nothing is defined, try these first:
        params = s1.get_vars()
        if "radiance" in params:
            var = "radiance"
        elif "radiance_noslit" in params:
            var = "radiance_noslit"
        elif "transmittance" in params:
            var = "transmittance"
        elif "transmittance_noslit" in params:
            var = "transmittance_noslit"
        else:
            # or plot the first variable we find
            var = list(params)[0]
            if var.replace("_noslit", "") in params:
                var = var.replace("_noslit", "")
    # ... create variable if it doesnt exist
    if var not in s1.get_vars():
        s1.get(var)
    if var not in s2.get_vars():
        s2.get(var)
    if Iunit == "default":
        try:
            Iunit = s1.units[var]
        except KeyError:  # unit not defined in dictionary
            raise KeyError(
                "Iunit not defined in spectrum for variable {0}. ".format(var)
                + "Cant use default unit. Specify unit in s.units['{0}'].".format(var)
            )
    if wunit == "default":
        wunit = s1.get_waveunit()
    wunit = cast_waveunit(wunit)

    return var, wunit, Iunit


def _get_wdiff_Idiff(
    s1: Spectrum,
    s2: Spectrum,
    var=None,
    wunit="default",
    Iunit="default",
    resample=True,
    method="diff",
    diff_window=0,
    normalize=False,
    verbose=False,
):
    """Get difference/ratio between two spectra.

    Properly make sure all units are the same.

    Calls :py:func:`radis.spectrum.compare.get_diff`,  :py:func:`radis.spectrum.compare.get_ratio`
    or  :py:func:`radis.spectrum.compare.get_distance`

    Called by :py:func:`radis.spectrum.compare.plot_diff`

    Parameters
    ----------
    all:
        see :py:func:`radis.spectrum.compare.plot_diff`

    Returns
    -------
    wdiff, Iunit: list of arrays
        difference or ratio (depending on ``method``)
    methods: list of str
        methods
    """

    # Get all comparison methods:
    # ------------------------
    if isinstance(method, list):
        methods = method
    else:
        methods = [method]

    for method in methods:
        if diff_window != 0 and method != "diff":
            raise NotImplementedError("diff_window with method {0}".format(method))

    # Get data
    # ----
    if normalize:
        # copy before modifying directly in spectrum
        s1 = s1.copy(copy_lines=False)
        s2 = s2.copy(copy_lines=False)
        w1, I1 = s1.get(var, wunit=wunit, copy=False)
        w2, I2 = s2.get(var, wunit=wunit, copy=False)
        ratio = np.nanmax(I1) / np.nanmax(I2)
        I1 /= np.nanmax(I1)
        I2 /= np.nanmax(I2)

        if verbose:
            print(("Rescale factor: " + str(ratio)))

    wdiffs, Idiffs = [], []
    for method in methods:
        if not normalize:
            if method == "distance":
                wdiff, Idiff = get_distance(
                    s1, s2, var=var, wunit=wunit, Iunit=Iunit, resample=resample
                )
            elif method == "diff":
                wdiff, Idiff = get_diff(
                    s1,
                    s2,
                    var=var,
                    wunit=wunit,
                    Iunit=Iunit,
                    diff_window=diff_window,
                    resample=resample,
                )
            elif method == "ratio":
                wdiff, Idiff = get_ratio(
                    s1, s2, var=var, wunit=wunit, Iunit=Iunit, resample=resample
                )
            else:
                raise ValueError("Unknown comparison method: {0}".format(method))
            wdiffs.append(wdiff)
            Idiffs.append(Idiff)
        else:
            if method == "distance":
                raise NotImplementedError(
                    "{0} was not implemented yet for normalized spectra".format(method)
                )
            elif method == "diff":
                if not resample:
                    raise NotImplementedError()
                wdiff, Idiff = curve_substract(w1, I1, w2, I2)
            elif method == "ratio":
                wdiff, Idiff = get_ratio(
                    s1, s2, var=var, wunit=wunit, Iunit=Iunit, resample=resample
                )
            else:
                raise ValueError("Unknown comparison method: {0}".format(method))
            wdiffs.append(wdiff)
            Idiffs.append(Idiff)
    return s1, s2, wdiffs, Idiffs, methods


def plot_diff(
    s1: Spectrum,
    s2: Spectrum,
    var=None,
    wunit="default",
    Iunit="default",
    resample=True,
    method="diff",
    diff_window=0,
    show_points=False,
    label1=None,
    label2=None,
    figsize=None,
    title=None,
    nfig=None,
    normalize=False,
    verbose=True,
    save=False,
    show=True,
    show_residual=False,
    lw_multiplier=1,
    diff_scale_multiplier=1,
    discard_centile=0,
    plot_medium="vacuum_only",
    legendargs={"loc": "best"},
    show_ruler=False,
):
    """Plot two spectra, and the difference between them. ``method=`` allows
    you to plot the absolute difference, ratio, or both.

    If waveranges dont match, ``s2`` is interpolated over ``s1``.


    Parameters
    ----------
    s1, s2: Spectrum objects
    var: str, or None
        spectral quantity to plot (ex: ``'abscoeff'``). If None,
        plot the first one in the Spectrum from ``'radiance'``,
        ``'radiance_noslit'``, ``'transmittance'``, etc.
    wunit: ``'default'``, ``'nm'``, ``'cm-1'``, ``'nm_vac'``
        wavespace unit:  wavelength air, wavenumber, wavelength vacuum.
        If ``'default'``, use first spectrum wunit
    Iunit: str
        if ``'default'``, use first spectrum unit
    method: ``'distance'``, ``'diff'``, ``'ratio'``, or list of them.
        If ``'diff'``, plot difference of the two spectra.
        If ``'distance'``, plot Euclidian distance (note that units are meaningless then)
        If ``'ratio'``, plot ratio of two spectra
        Default ``'diff'``.

        .. warning::
            with ``'distance'``, calculation scales as ~N^2 with N the number
            of points in a spectrum (against ~N with ``'diff'``). This can quickly
            override all memory.

        Can also be a list::

            method=['diff', 'ratio']

    normalize: bool
        Normalize the spectra to be ploted

    Other Parameters
    ----------------
    diff_window: int
        If non 0, calculates diff by offsetting s1 by ``diff_window`` number of
        units on either side, and returns the minimum. Kinda compensates for experimental
        errors on the w axis. Default 0. (look up code to understand...)
    show_points: boolean
        if ``True``, make all points appear with 'o'
    label1, label2: str
        curve names
    figsize
        figure size
    nfig: int, str
        figure number of name
    title: str
        title
    verbose: boolean
        if ``True``, plot stuff such as rescale ratio in normalize mode. Default ``True``
    save: str
        Default is ``False``. By default won't save anything, type the path of the
        destination if you want to save it (format in the name).
    show: Bool
        Default is ``True``. Will show the plots : bad if there are more than 20.
    show_residual: bool
        if ``True``, calculates and shows on the graph the residual in L2 norm.
        See :func:`~radis.spectrum.compare.get_residual`. ``diff_window`` is
        used in the residual calculation too. ``normalize`` has no effect.
    diff_scale_multiplier: float
        dilate the diff plot scale. Default ``1``
    discard_centile: int
        if not ``0``, discard the firsts and lasts centile when setting the limits
        of the diff window. Example::

            discard_centile=1     #  --> discards the smallest 1% and largest 1%
            discard_centile=10    #  --> discards the smallest 10% and largest 10%

        Useful to remove spikes in a ratio, for instance.
        Note that this does not change the values of the residual. It's just
        a plot feature.
        Default ``0``

    plot_medium: bool, ``'vacuum_only'``
        if ``True`` and ``wunit`` are wavelengths, plot the propagation medium
        in the xaxis label (``[air]`` or ``[vacuum]``). If ``'vacuum_only'``,
        plot only if ``wunit=='nm_vac'``. Default ``'vacuum_only'``
        (prevents from inadvertently plotting spectra with different propagation
        medium on the same graph).
    legendargs: dict
        format arguments forwarded to the legend
    show_ruler: bool
        if `True`, add a ruler tool to the Matplotlib toolbar.

        .. warning::
            still experimental in 0.9.30 ! Try it, feedback welcome !

    Returns
    -------
    fig: figure
        fig
    [ax0, ax1]: axes
        spectra and difference axis


    Examples
    --------
    Simple use::

        from radis import plot_diff
        plot_diff(s10, s50)                # s10, s50 are two spectra

    Advanced use, plotting the total power in the label, and getting the figure
    and axes handle to edit them afterwards::

        Punit = 'mW/cm2/sr'
        fig, axes = plot_diff(s10, s50, 'radiance_noslit', figsize=(18,6),
              label1='brd 10 cm-1, P={0:.2f} {1}'.format(s10.get_power(unit=Punit),Punit),
              label2='brd 50 cm-1, P={0:.2f} {1}'.format(s50.get_power(unit=Punit),Punit)
              )
        # modify fig, axes..

    See an example in :ref:`label_spectrum_howto_compare`, which produces the output below:

    .. image:: https://radis.readthedocs.io/en/latest/_images/cdsd4000_vs_hitemp_3409K.svg

    If you wish to plot in a logscale, it can be done in the following way::

        fig, [ax0, ax1] = plot_diff(s0, s1, normalize=False, verbose=False)
        ylim0 = ax0.get_ybound()
        ax0.set_yscale("log")
        ax0.set_ybound(ylim0)

    See Also
    --------
    :func:`~radis.spectrum.compare.get_diff`,
    :func:`~radis.spectrum.compare.get_ratio`,
    :func:`~radis.spectrum.compare.get_distance`,
    :func:`~radis.spectrum.compare.get_residual`,
    :meth:`~radis.spectrum.spectrum.compare_with`
    """
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.widgets import MultiCursor

    from radis.misc.plot import fix_style, set_style

    if (not show) and (
        not save
    ):  # I added this line to avoid calculus in the case there is nothing to do (Minou)
        if verbose:
            print("plot_diff : Nothing to do")
        return None, None

    var, wunit, Iunit = get_default_units(s1, s2, var=var, wunit=wunit, Iunit=Iunit)

    s1, s2, wdiffs, Idiffs, methods = _get_wdiff_Idiff(
        s1,
        s2,
        var=var,
        method=method,
        wunit=wunit,
        Iunit=Iunit,
        resample=resample,
        diff_window=diff_window,
        normalize=normalize,
        verbose=verbose,
    )

    # Plot
    # ----

    # Format units
    xlabel = format_xlabel(wunit, plot_medium)

    # Init figure
    set_style()
    fig = plt.figure(num=nfig, figsize=figsize)
    gs = gridspec.GridSpec(1 + len(methods), 1, height_ratios=[3] + [1] * len(methods))
    ax0 = fig.add_subplot(gs[0])
    ax0.ticklabel_format(useOffset=False)
    ax1 = []
    for i in range(len(methods)):
        ax1i = fig.add_subplot(gs[i + 1])
        ax1i.get_shared_x_axes().join(ax0, ax1i)
        ax1i.ticklabel_format(useOffset=False)
        ax1.append(ax1i)

    # Plotting style
    if show_points:
        style = "-o"
    else:
        style = "-"

    # Get labels and names
    if label1 is None:
        label1 = s1.get_name()
    if label2 is None:
        label2 = s2.get_name()
    # Max label length:
    if len(label1) > 60:
        label1 = label1[:58] + "..."
    if len(label2) > 60:
        label2 = label2[:58] + "..."

    # Plot compared spectra
    # ... note: if 'normalize', s1 & s2 have been edited by _get_wdiff_Idiff
    ax0.plot(
        *s1.get(var, wunit=wunit, Iunit=Iunit),
        style,
        color="k",
        lw=3 * lw_multiplier,
        label=label1
    )
    ax0.plot(
        *s2.get(var, wunit=wunit, Iunit=Iunit),
        style,
        color="r",
        lw=1 * lw_multiplier,
        label=label2
    )

    # cosmetic changes
    Iunit = make_up_unit(Iunit, var)

    ax0.tick_params(labelbottom=False)
    if label1 is not None or label2 is not None:
        ax0.legend(**legendargs)

    # Start to 0
    if var in ["radiance_noslit", "radiance", "abscoeff", "absorbance"]:
        ax0.set_ylim(bottom=0)

    # plot difference (sorted)
    for ax1i, wdiff, Idiff in zip(ax1, wdiffs, Idiffs):
        b = np.argsort(wdiff)
        ax1i.plot(wdiff[b], Idiff[b], style, color="k", lw=1 * lw_multiplier)

    # Write labels
    ax1[-1].set_xlabel(make_up(xlabel))
    if normalize:
        fig.text(
            0.02,
            0.5,
            ("{0} (norm.)".format(make_up(var))),
            va="center",
            rotation="vertical",
        )
    else:
        fig.text(
            0.02,
            0.5,
            ("{0} ({1})".format(make_up(var), Iunit)),
            va="center",
            rotation="vertical",
        )

    # Set limits of 'diff' window
    for i, method in enumerate(methods):
        if method == "diff":
            # symmetrize error scale:
            # auto-zoom on min, max, but discard first and last centile (case of spikes / divergences)
            Idiff = Idiffs[i]
            if discard_centile:
                Idiff_sorted = np.sort(Idiff[~np.isnan(Idiff)])
                ymax = max(
                    abs(Idiff_sorted[-int(discard_centile * len(Idiff_sorted) // 100)]),
                    abs(Idiff_sorted[int(discard_centile * len(Idiff_sorted) // 100)]),
                )
            else:
                ymax = np.nanmax(abs(Idiff))
            ax1[i].set_ylim(-ymax * diff_scale_multiplier, ymax * diff_scale_multiplier)
        elif method == "distance":
            if discard_centile:
                raise NotImplementedError(
                    "discard_centile not implemented for method=distance"
                )
            _, ymax = ax1[i].get_ylim()
            ax1[i].set_ylim(0, ymax * diff_scale_multiplier)
        elif method == "ratio":
            # auto-zoom on min, max, but discard first and last centile (case of spikes / divergences)
            Idiff = Idiffs[i]
            if discard_centile:
                Idiff_sorted = np.sort(Idiff[~np.isnan(Idiff)])
                ymin = Idiff_sorted[int(discard_centile * len(Idiff_sorted) // 100)]
                ymax = Idiff_sorted[-int(discard_centile * len(Idiff_sorted) // 100)]
            else:
                ymin = np.nanmin(Idiff)
                ymax = np.nanmax(Idiff)
            ax1[i].set_ylim(
                bottom=((ymin - 1) * diff_scale_multiplier + 1),
                top=((ymax - 1) * diff_scale_multiplier + 1),
            )
    #            ymax = max(abs(Idiff_sorted[len(Idiff_sorted)//100]-1),
    #                       abs(Idiff_sorted[len(-Idiff_sorted)//100]-1))
    #            ax1[i].set_ylim(ymax*diff_scale_multiplier+1, -ymax*diff_scale_multiplier+1)

    if title:
        fig.suptitle(title)
    # Fix format
    fix_style(ax=ax0)
    for ax1i in ax1:
        fix_style(ax=ax1i)
    plt.tight_layout()
    if title:
        plt.subplots_adjust(left=0.15, top=0.92)
    else:
        plt.subplots_adjust(left=0.15)

    # Plot difference text
    for i, method in enumerate(methods):
        if method == "diff":
            difftext = "diff"
            ax1[i].axhline(y=0, color="grey", zorder=-1)
        elif method == "distance":
            difftext = "distance"
            ax1[i].axhline(y=0, color="grey", zorder=-1)
        elif method == "ratio":
            difftext = "ratio"
            ax1[i].axhline(y=1, color="grey", zorder=-1)

        # Show residualget_residual
        if show_residual:
            difftext += " (residual={0:.2g})".format(
                get_residual(
                    s1, s2, var=var, norm="L2", ignore_nan=True, diff_window=diff_window
                )
            )
        pos = ax1[i].get_position()
        fig.text(0.09, pos.ymax + 0.02, difftext)

    # Add tools
    # ... Add cursors
    axes = [ax0] + ax1
    fig.cursors = MultiCursor(
        fig.canvas,
        axes,
        color="r",
        lw=1 * lw_multiplier,
        alpha=0.2,
        horizOn=False,
        vertOn=True,
    )

    # ... Add
    if show_ruler:
        from radis.tools.plot_tools import add_ruler

        add_ruler(fig, wunit=wunit, Iunit=Iunit, ax=axes[0])

    if show:
        plt.show()
    if save:
        fig.savefig(save)
        if not show:
            plt.close(fig)  # to avoid memory load

    # Return graphs
    return fig, axes


def averageDistance(s1: Spectrum, s2: Spectrum, var="radiance"):
    """Return the average distance between two spectra.
    It's important to note for fitting that if averageDistance(s1, s2)==0
    then s1 = s2

    Parameters
    ----------

    s1, s2: Spectrum objects
        spectra to be compared

    var: str, optional
        spectral quantity (ex: 'radiance', 'transmittance'...)

    Returns
    -------

    distance: float
        Average distance as in the following equation:

        .. math::

            dist = \\frac{\\sqrt{\\sum_i {(s1_i-s2_i)^2}}}{N}.

    """

    warn(
        FutureWarning(
            "will be deprecated in future versions. Use get_residual(s1, s2, norm=L2) instead"
        )
    )

    distArray = get_diff(s1, s2, var)
    distArray_Y2 = distArray[1] * distArray[1]
    N = np.size(distArray[1])
    distance = np.sqrt(np.sum(distArray_Y2)) / N
    return distance


# %% Function to compare Spectra including conditions and lines


def compare_spectra(
    first: Spectrum,
    other: Spectrum,
    spectra_only=False,
    plot=True,
    wunit="default",
    verbose=True,
    rtol=1e-5,
    ignore_nan=False,
    ignore_outliers=False,
    ignore_conditions=["calculation_time"],
    normalize=False,
    **kwargs
) -> bool:
    """Compare Spectrum with another Spectrum object.

    Parameters
    ----------
    first: type Spectrum
        a Spectrum to be compared
    other: type Spectrum
        another Spectrum to compare with
    spectra_only: boolean, or str
        if ``True``, only compares spectral quantities (in the same waveunit)
        and not lines or conditions. If str, compare a particular quantity
        name. If False, compare everything (including lines and conditions
        and populations). Default ``False``
    plot: boolean
        if ``True``, use plot_diff to plot all quantities for the 2 spectra
        and the difference between them. Default ``True``.
    wunit: 'nm', 'cm-1', 'default'
        in which wavespace to compare (and plot). If default, natural wavespace
        of first Spectrum is taken
    rtol: float
        relative difference to use for spectral quantities comparison.
        absolute difference ``atol`` is set to 0.
    ignore_nan: boolean
        if ``True``, nans are ignored when comparing spectral quantities
    ignore_outliers: boolean, or float
        if not False, outliers are discarded. i.e, output is determined by::

            out = (~np.isclose(I, Ie, rtol=rtol, atol=0)).sum()/len(I) < ignore_outliers
    ignore_conditions: list
        do not compare the metadata of these keys
    normalize: bool
        Normalize the spectra to be plotted

    Other Parameters
    ----------------
    kwargs: dict
        arguments are forwarded to :func:`~radis.spectrum.compare.plot_diff`

    Returns
    -------
    equals: boolean
        return True if spectra are equal (respective to tolerance defined by
        rtol and other input conditions)


    Examples
    --------
    Compare two Spectrum objects, or specifically the transmittance::

        s1.compare_with(s2)
        s1.compare_with(s2, 'transmittance')


    Note that you can also simply use `s1 == s2`, that uses
    :meth:`~radis.spectrum.spectrum.Spectrum.compare_with` internally::

        s1 == s2       # will return True or False
    """

    # Check inputs
    if not 0 <= ignore_outliers < 1:
        raise ValueError("ignore_outliers should be < 1, or False")
    if not isinstance(other, Spectrum):
        raise TypeError(
            "2nd object is not a Spectrum: got class {0}".format(other.__class__)
        )
    if isinstance(spectra_only, str):  # case where we compare all quantities
        if not spectra_only in first.get_vars():
            raise ValueError(
                "{0} is not a spectral quantity in our Spectrum ({1})".format(
                    spectra_only, first.get_vars()
                )
            )
        if not spectra_only in other.get_vars():
            raise ValueError(
                "{0} is not a spectral quantity in the other Spectrum ({1})".format(
                    spectra_only, other.get_vars()
                )
            )
    if verbose:  # print conditions
        what = spectra_only if isinstance(spectra_only, str) else "all quantities"
        msg = "compare {0} with rtol={1}".format(what, rtol)
        if ignore_nan:
            msg += ", ignore_nan"
        if ignore_outliers:
            msg += ", ignore_outliers={0}".format(ignore_outliers)
        print(msg)
    if not plot and len(kwargs) > 0:
        raise ValueError("Unexpected argument: {0}".format(kwargs))

    if wunit == "default":
        wunit = first.get_waveunit()

    def _compare_dataframes(df1, df2, name):
        """

        Parameters
        ----------

        df1, df2: pandas Dataframe
            lines, or vib/rovib levels dataframes

        name: str
            for error message
        """

        from pandas.util.testing import assert_frame_equal

        try:
            assert_frame_equal(
                df1.sort_index(axis=0).sort_index(axis=1),
                df2.sort_index(axis=0).sort_index(axis=1),
                check_names=True,
                check_column_type=False,  # solves problem in Python 2/3 dataframes (unicode/str)
            )
            out = True

        except AssertionError as err:
            if verbose:
                print("Comparing ", name)
                print(err.args[0])
            out = False

        return out

    def _compare_arrays(I, Ie):
        """Compare spectral quantities I and Ie."""

        if ignore_nan:
            b = ~(np.isnan(I) + np.isnan(Ie))
            I = I[b]
            Ie = Ie[b]

        # in all cases, we ignore if there are nans at the same
        # on both arrays (happens with convolved values in 'valid' mode)
        equal_nan = True

        if ignore_outliers:
            out = (
                ~np.isclose(I, Ie, rtol=rtol, atol=0, equal_nan=equal_nan)
            ).sum() / len(I) < ignore_outliers
        else:
            out = np.allclose(I, Ie, rtol=rtol, atol=0, equal_nan=equal_nan)

        return bool(out)

    def _display_difference(q, q0):
        error = np.nanmax(abs(q / q0 - 1))
        avgerr = np.nanmean(abs(q / q0 - 1))
        print(
            "...",
            k,
            "don't match (up to {0:.3}% diff.,".format(error * 100)
            + " average {0:.3f}%)".format(avgerr * 100),
        )

    b = True
    if isinstance(spectra_only, str):  # compare this quantity
        vars = [spectra_only]
    else:  # compare all quantities
        b = set(first.get_vars()) == set(other.get_vars())
        if not b and verbose:
            print(
                "... list of quantities dont match: {0} vs {1}".format(
                    first.get_vars(), other.get_vars()
                )
            )
        vars = [k for k in first.get_vars() if k in other.get_vars()]

    if spectra_only:
        # Compare spectral arrays
        # -----------
        for k in vars:
            w, q = first.get(k, wunit=wunit)
            w0, q0 = other.get(k, wunit=wunit)
            if len(w) != len(w0):
                print(
                    "Wavespaces have different length (for {0}: {1} vs {2})".format(
                        k, len(w), len(w0)
                    )
                )
                print("We interpolate one spectrum on the other one")

                from scipy.interpolate import interp1d

                if len(q) > len(q0):
                    f = interp1d(w, q, kind="cubic")
                    new_q = f(w0)
                    b1 = _compare_arrays(new_q, q0)
                    if not b1 and verbose:
                        _display_difference(new_q, q0)
                else:
                    f = interp1d(w0, q0, kind="cubic")
                    new_q0 = f(w)
                    b1 = _compare_arrays(q, new_q0)
                    if not b1 and verbose:
                        _display_difference(q, new_q0)

            else:  # no need to interpolate
                b1 = np.allclose(w, w0, rtol=rtol, atol=0)
                b1 *= _compare_arrays(q, q0)

                if not b1 and verbose:
                    _display_difference(q, q0)
            b *= b1

            if plot:
                try:
                    plot_diff(
                        first,
                        other,
                        var=k,
                        wunit=wunit,
                        normalize=normalize,
                        verbose=verbose,
                        **kwargs
                    )
                except:
                    print("... couldn't plot {0}".format(k))

    else:
        # Compare spectral variables
        # -----------
        for k in vars:
            w, q = first.get(k, wunit=wunit)
            w0, q0 = other.get(k, wunit=wunit)
            if len(w) != len(w0):
                print(
                    "Wavespaces have different length (for {0}: {1} vs {2})".format(
                        k, len(w), len(w0)
                    )
                )
                b1 = False
            else:
                b1 = np.allclose(w, w0, rtol=rtol, atol=0)
                b1 *= _compare_arrays(q, q0)
                if not b1 and verbose:
                    error = np.nanmax(abs(q / q0 - 1))
                    avgerr = np.nanmean(abs(q / q0 - 1))
                    print(
                        "...",
                        k,
                        "don't match (up to {0:.3}% diff.,".format(error * 100)
                        + " average {0:.3f}%)".format(avgerr * 100),
                    )
            b *= b1

            if plot:
                try:
                    plot_diff(
                        first,
                        other,
                        var=k,
                        wunit=wunit,
                        normalize=normalize,
                        verbose=verbose,
                        **kwargs
                    )
                except:
                    print("... there was an error while plotting {0}".format(k))

        # Compare conditions and units
        # -----------
        verbose_dict = "if_different" if verbose else False
        b1 = (
            compare_dict(
                first.conditions,
                other.conditions,
                verbose=verbose_dict,
                ignore_keys=ignore_conditions,
            )
            == 1
        )
        b2 = compare_dict(first.cond_units, other.cond_units, verbose=verbose_dict) == 1
        b3 = compare_dict(first.units, other.units, verbose=verbose_dict) == 1
        if not b1 and verbose:
            print("... conditions don't match")
        if not b2 and verbose:
            print("... conditions units don't match")
        if not b3 and verbose:
            print("... units don't match")
        b *= b1 * b2 * b3

        # Compare lines
        # -----------
        if first.lines is None and other.lines is None:
            b4 = True
        elif first.lines is None:
            b4 = False
        elif other.lines is None:
            b4 = False
        else:
            b4 = _compare_dataframes(first.lines, other.lines, "lines")
        if not b4 and verbose:
            print("... lines dont match")
        b *= b4

        # Compare populations
        # -----------
        if first.populations is None and other.populations is None:
            b5 = True
        elif first.populations is None:
            b5 = False
        elif other.populations is None:
            b5 = False
        else:
            # Compare keys in populations
            b5 = True
            if (
                compare_lists(
                    first.populations, other.populations, verbose="if_different"
                )
                == 1
            ):
                # same molecules, compare isotopes
                for molecule, isotopes in first.populations.items():
                    if (
                        compare_lists(
                            isotopes,
                            other.populations[molecule],
                            verbose="if_different",
                        )
                        == 1
                    ):
                        # same isotopes, compare electronic states
                        for isotope, states in isotopes.items():
                            if (
                                compare_lists(
                                    states,
                                    other.populations[molecule][isotope],
                                    verbose="if_different",
                                )
                                == 1
                            ):
                                # same electronic states, compare levels + other information
                                for state, content in states.items():
                                    for k, v in content.items():
                                        if k in ["vib", "rovib"]:
                                            b5 *= _compare_dataframes(
                                                v,
                                                other.populations[molecule][isotope][
                                                    state
                                                ][k],
                                                "populations of {0}({1})(iso{2})".format(
                                                    molecule, state, isotope
                                                ),
                                            )
                                        else:
                                            b5 *= (
                                                v
                                                == other.populations[molecule][isotope][
                                                    state
                                                ][k]
                                            )
                            else:
                                b5 = False
                                if verbose:
                                    print(
                                        "{0}(iso{1}) states are different (see above)".format(
                                            molecule, isotope
                                        )
                                    )
                    else:
                        b5 = False
                        if verbose:
                            print(
                                "{0} isotopes are different (see above)".format(
                                    molecule
                                )
                            )
            else:
                b5 = False
                if verbose:
                    print("Molecules are different (see above)")
        if not b5 and verbose:
            print("... populations dont match (see detail above)")
        b *= b5

        # Compare slit
        # -----------

        if len(first._slit) == len(other._slit) == 0:
            # no slit anywhere
            b6 = True
        elif len(first._slit) != len(other._slit):
            b6 = False
            if verbose:
                print("A spectrum has slit function array but the other doesnt")
        else:
            ws, Is = first.get_slit()
            ws0, Is0 = other.get_slit()
            if len(ws) != len(ws0):
                if verbose:
                    print("Slits have different length")
                b6 = False
            else:
                b6 = np.allclose(ws, ws0, rtol=rtol, atol=0)
                b6 *= _compare_arrays(Is, Is0)
                if not b6 and verbose:
                    print("Slit functions dont match")
        b *= b6

    return bool(b)


if __name__ == "__main__":

    from radis.test.spectrum.test_compare import _run_testcases

    print(("Test compare.py: ", _run_testcases(verbose=True)))
