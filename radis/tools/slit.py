#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Routine Listing
---------------

This module contains...

a function to convolve with a slit function, that can resample and correct
for slit dilation:

- :func:`~radis.tools.slit.convolve_with_slit`

This function is called directly by the Spectrum method
:py:meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`

predefined slit functions generators (triangular, gaussian, etc... see :data:`~radis.tools.slit.SLIT_SHAPES`)
and experimental slit function importer:

- :func:`~radis.tools.slit.get_slit_function`
- :func:`~radis.tools.slit.import_experimental_slit`
- :func:`~radis.tools.slit.triangular_slit`
- :func:`~radis.tools.slit.trapezoidal_slit`
- :func:`~radis.tools.slit.gaussian_slit`

functions to manipulate slit functions (plot, get_FWHM, get_effective_FWHM,
recenter_slit, crop_slit):

- :func:`~radis.tools.slit.get_FWHM`
- :func:`~radis.tools.slit.get_effective_FWHM`
- :func:`~radis.tools.slit.plot_slit`
- :func:`~radis.tools.slit.recenter_slit`
- :func:`~radis.tools.slit.crop_slit`


-------------------------------------------------------------------------------



"""


from warnings import warn

import numpy as np
from numpy import exp
from numpy import log as ln
from numpy import sqrt, trapz
from scipy.interpolate import splev, splrep
from scipy.signal import oaconvolve

from radis.misc.arrays import anynan, evenly_distributed, evenly_distributed_fast
from radis.misc.basics import is_float
from radis.misc.debug import printdbg
from radis.misc.signal import resample_even
from radis.misc.warning import SlitDispersionWarning
from radis.phys.convert import cm2nm, dcm2dnm, dnm2dcm, nm2cm
from radis.spectrum.spectrum import cast_waveunit

SLIT_SHAPES = ["triangular", "trapezoidal", "gaussian"]
"""list : list of predefined slit shapes
"""

# %% Get slit function


def get_slit_function(
    slit_function,
    unit="nm",
    norm_by="area",
    shape="triangular",
    center_wavespace=None,
    return_unit="same",
    wstep=None,
    plot=False,
    resfactor=2,
    verbose=True,
    auto_recenter_crop=True,
    *args,
    **kwargs,
):
    """Import or generate slit function in correct wavespace
    Give a file path to import, or a float / tuple to generate arbitrary shapes

    Warning with units: read about unit and return_unit parameters.

    See :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit` and
    :func:`~radis.tools.slit.convolve_with_slit` for more info


    Parameters
    ----------

    slit_function: tuple, or str
        If float:
            generate slit function with FWHM of `slit_function` (in `unit`)
        If .txt file:
            import experimental slit function (in `unit`): format must be 2-columns
            with wavelengths and intensity (doesn't have to be normalized)

    unit: 'nm' or 'cm-1'
        unit of slit_function FWHM, or unit of imported file

    norm_by: ``'area'``, ``'max'``, or ``None``
        how to normalize. ``'area'`` conserves energy. With ``'max'`` the slit is normalized
        at peak so that the maximum is one.

        .. note::

            ``'max'`` changes the unit of the spectral array, e.g. from
            ``'mW/cm2/sr/µm'`` to ``'mW/cm2/sr')``

        ``None`` doesnt normalize. Default ``'area'``

    shape: ``'triangular'``, ``'trapezoidal'``, ``'gaussian'``
        which shape to use when generating a slit. Default ``'triangular'``.
        Should be one of :py:data:`~radis.tools.slit.SLIT_SHAPES`.

    center_wavespace: float, or None
        center of slit when generated (in unit). Not used if slit is imported.

    return_unit: ``'nm'``, ``'cm-1'``, or ``'same'``
        if not ``'same'``, slit is converted to the given wavespace.
        Default ``'same'``.

    wstep: float
        which discretization step to use (in return_unit) when generating a slit
        function. Not used if importing


    Other Parameters
    ----------------

    resfactor: int
        resolution increase when resampling from nm to cm-1, or the other way
        round. Default 2.

    energy_threshold: float
         tolerance fraction. Only used when importing experimental slit as the
         theoretical slit functions are directly generated in spectrum wavespace
         Default 1e-3 (0.1%)

    auto_recenter_crop: bool
        if ``True``, recenter slit and crop zeros on the side when importing
        an experimental slit. Default ``True``.
        See :func:`~radis.tools.slit.recenter_slit`, :func:`~radis.tools.slit.crop_slit`


    Returns
    -------

    wslit, Islit: array
        wslit is in `return_unit` . Islit is normalized according to  `norm_by`


    Examples
    --------

    >>> wslit, Islit = get_slit_function(1, 'nm', shape='triangular',
    center_wavespace=600, wstep=0.01)

    Returns a triangular slit function of FWHM = 1 nm, centered on 600 nm, with
    a step of 0.01 nm

    >>> wslit, Islit = get_slit_function(1, 'nm', shape='triangular',
    center_wavespace=600, return_unit='cm-1', wstep=0.01)

    Returns a triangular slit function expressed in cm-1, with a FWHM = 1 nm
    (converted in equivalent width in cm-1 at 600 nm), centered on 600 nm, with
    a step of 0.01 cm-1 (!)


    Notes
    -----

    In ``norm_by='max'`` mode, slit is normalized by slit max. In RADIS, this is done
    in the spectrum wavespace (to avoid errors that would be caused by interpolating
    the spectrum).

    A problem arise if the spectrum wavespace is different from the slit wavespace:
    typically, slit is in ``'nm'`` but a spectrum calculated by RADIS is stored in
    ``'cm-1'``: in that case, the convoluted spectrum is divided by ``int(Islit*dν)``
    instead of ``int(Islit*dλ)``. The output unit is then ``[radiance]*[spec_unit]``
    instead of ``[radiance]*[slit_unit]``, i.e. typically, ``[mW/cm2/sr/nm]*[cm-1]``
    instead of ``[mW/cm2/sr/nm]*[nm]=[mW/cm2/sr]``

    While this remains true if the units are taken into account, this is not the
    expected behaviour. For instance, Specair users are used to having a FWHM(nm)
    factor between spectra convolved with slit normalized by max and slit normalized
    by area.

    The ``norm_by='max'`` behaviour adds a correction factor
    ```/int(Islit*dλ)/int(Islit*dν)``` to maintain an output spectrum in ``[radiance]*[slit_unit]``

    See Also
    --------

    :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`,
    :func:`~radis.tools.slit.convolve_with_slit`

    """

    if "waveunit" in kwargs:
        assert return_unit == "same"  # default
        return_unit = kwargs.pop("waveunit")
        warn(DeprecationWarning("waveunit renamed return_unit"))
    if "slit_unit" in kwargs:
        assert unit == "nm"  # default
        unit = kwargs.pop("slit_unit")
        warn(DeprecationWarning("slit_unit renamed unit"))

    energy_threshold = kwargs.pop("energy_threshold", 1e-3)  # type: float
    # tolerance fraction
    # when resampling (only used in experimental slit as the)
    # theoretical slit functions are directly generated in
    # spectrum wavespace

    def check_input_gen():
        if center_wavespace is None:
            raise ValueError(
                "center_wavespace has to be given when generating " + "slit function"
            )
        if wstep is None:
            raise ValueError("wstep has to be given when generating " + "slit function")

    # Cast units
    if return_unit == "same":
        return_unit = unit
    unit = cast_waveunit(unit)
    return_unit = cast_waveunit(return_unit)
    # in norm_by=max mode, used to keep units in [Iunit]*return_unit in [Iunit]*unit
    scale_slit = 1
    # not used in norm_by=area mode

    # Generate Slit
    # -------------
    # First get the slit in return_unit space
    if is_float(
        slit_function
    ):  # Generate slit function (directly in return_unit space)

        check_input_gen()

        # ... first get FWHM in return_unit  (it is in `unit` right now)
        FWHM = slit_function
        if return_unit == "cm-1" and unit == "nm":
            # center_wavespace ~ nm, FWHM ~ nm
            FWHM = dnm2dcm(FWHM, center_wavespace)  # wavelength > wavenumber
            center_wavespace = nm2cm(center_wavespace)
            if norm_by == "max":
                scale_slit = slit_function / FWHM  # [unit/return_unit]
        elif return_unit == "nm" and unit == "cm-1":
            # center_wavespace ~ cm-1, FWHM ~ cm-1
            FWHM = dcm2dnm(FWHM, center_wavespace)  # wavenumber > wavelength
            center_wavespace = cm2nm(center_wavespace)
            if norm_by == "max":
                scale_slit = slit_function / FWHM  # [unit/return_unit]
        else:
            pass  # correct unit already
        # Now FWHM is in 'return_unit'

        # ... now, build it (in our wavespace)
        if __debug__:
            printdbg(
                "get_slit_function: {0} FWHM {1:.2f}{2}, center {3:.2f}{2}, norm_by {4}".format(
                    shape, FWHM, return_unit, center_wavespace, norm_by
                )
            )

        if shape == "triangular":
            wslit, Islit = triangular_slit(
                FWHM,
                wstep,
                center=center_wavespace,
                bplot=plot,
                norm_by=norm_by,
                wunit=return_unit,
                scale=scale_slit,
                *args,
                **kwargs,
            )

        # Insert other slit shapes here
        # ...

        elif shape == "gaussian":
            wslit, Islit = gaussian_slit(
                FWHM,
                wstep,
                center=center_wavespace,
                bplot=plot,
                norm_by=norm_by,
                wunit=return_unit,
                scale=scale_slit,
                *args,
                **kwargs,
            )

        elif shape == "trapezoidal":
            raise TypeError("A (top, base) tuple must be given with a trapezoidal slit")

        else:
            raise TypeError(
                "Slit function ({0}) not in known slit shapes: {1}".format(
                    shape, SLIT_SHAPES
                )
            )

    elif isinstance(slit_function, tuple):

        check_input_gen()

        try:
            top, base = slit_function
        except:
            raise TypeError("Wrong format for slit function: {0}".format(slit_function))
        if shape == "trapezoidal":
            pass
        elif shape == "triangular":  # it's the default
            warn("Triangular slit given with a tuple: we used trapezoidal slit instead")
            shape = "trapezoidal"
        else:
            raise TypeError("A (top, base) tuple must be used with a trapezoidal slit")

        # ... first get FWHM in our wavespace unit
        if return_unit == "cm-1" and unit == "nm":
            # center_wavespace ~ nm, FWHM ~ nm
            top = dnm2dcm(top, center_wavespace)  # wavelength > wavenumber
            base = dnm2dcm(base, center_wavespace)  # wavelength > wavenumber
            center_wavespace = nm2cm(center_wavespace)
            if norm_by == "max":
                scale_slit = sum(slit_function) / (top + base)  # [unit/return_unit]
        elif return_unit == "nm" and unit == "cm-1":
            # center_wavespace ~ cm-1, FWHM ~ cm-1
            top = dcm2dnm(top, center_wavespace)  # wavenumber > wavelength
            base = dcm2dnm(base, center_wavespace)  # wavenumber > wavelength
            center_wavespace = cm2nm(center_wavespace)
            if norm_by == "max":
                scale_slit = sum(slit_function) / (top + base)  # [unit/return_unit]
        else:
            pass  # correct unit already

        FWHM = (top + base) / 2

        # ... now, build it (in our wavespace)
        if __debug__:
            printdbg(
                "get_slit_function: {0}, FWHM {1:.2f}{2}, center {3:.2f}{2}, norm_by {4}".format(
                    shape, FWHM, return_unit, center_wavespace, norm_by
                )
            )

        wslit, Islit = trapezoidal_slit(
            top,
            base,
            wstep,
            center=center_wavespace,
            bplot=plot,
            norm_by=norm_by,
            wunit=return_unit,
            scale=scale_slit,
            *args,
            **kwargs,
        )

    # Or import it from file or numpy input
    # ------------
    elif isinstance(slit_function, str) or isinstance(
        slit_function, np.ndarray
    ):  # import it
        if __debug__:
            printdbg(
                "get_slit_function: {0} in {1}, norm_by {2}, return in {3}".format(
                    slit_function, unit, norm_by, return_unit
                )
            )

        wslit, Islit = import_experimental_slit(
            slit_function,
            norm_by=norm_by,  # norm is done later anyway
            wunit=unit,
            verbose=verbose,
            auto_crop=auto_recenter_crop,
            auto_recenter=auto_recenter_crop,
            bplot=False,  # we will plot after resampling
            *args,
            **kwargs,
        )

        # ... get unit
        # Normalize
        if norm_by == "area":  # normalize by the area
            #        I_slit /= np.trapz(I_slit, x=w_slit)
            Iunit = "1/{0}".format(unit)
        elif norm_by == "max":  # set maximum to 1
            Iunit = ""
        elif norm_by is None:
            Iunit = None
        else:
            raise ValueError(
                "Unknown normalization type: `norm_by` = {0}".format(norm_by)
            )

        # ... check it looks correct
        unq, counts = np.unique(wslit, return_counts=True)
        dup = counts > 1
        if dup.sum() > 0:
            raise ValueError(
                "Not all wavespace points are unique: slit function "
                + "format may be wrong. Duplicates for w={0}".format(unq[dup])
            )

        # ... resample if needed
        if return_unit == "cm-1" and unit == "nm":  # wavelength > wavenumber
            wold, Iold = wslit, Islit
            wslit, Islit = resample_even(
                nm2cm(wslit),
                Islit,
                resfactor=resfactor,
                energy_threshold=energy_threshold,
                print_conservation=True,
            )
            scale_slit = trapz(Iold, wold) / trapz(Islit, wslit)  # [unit/return_unit]
            renormalize = True
        elif return_unit == "nm" and unit == "cm-1":  # wavenumber > wavelength
            wold, Iold = wslit, Islit
            wslit, Islit = resample_even(
                cm2nm(wslit),
                Islit,
                resfactor=resfactor,
                energy_threshold=energy_threshold,
                print_conservation=True,
            )
            scale_slit = trapz(Iold, wold) / trapz(Islit, wslit)  # [unit/return_unit]
            renormalize = True
        else:  # return_unit == unit
            renormalize = False
        # Note: if wstep dont match with quantity it's alright as it gets
        # interpolated in the `convolve_with_slit` function

        # re-Normalize if needed (after changing units)
        if renormalize:
            if __debug__:
                printdbg("get_slit_function: renormalize")
            if norm_by == "area":  # normalize by the area
                Islit /= abs(np.trapz(Islit, x=wslit))
                Iunit = "1/{0}".format(return_unit)
            elif norm_by == "max":  # set maximum to 1
                Islit /= abs(np.max(Islit))
                Islit *= scale_slit
                Iunit = ""
                if scale_slit != 1:
                    Iunit += "x{0}".format(scale_slit)
            elif norm_by is None:
                Iunit = None
            else:
                raise ValueError(
                    "Unknown normalization type: `norm_by` = {0}".format(norm_by)
                )

        if plot:  # (plot after resampling / renormalizing)
            # Plot slit
            plot_slit(wslit, Islit, wunit=return_unit, Iunit=Iunit)

    else:
        raise TypeError(
            "Unexpected type for slit function: {0}".format(type(slit_function))
        )

    # Final shape
    if len(wslit) < 5:
        raise ValueError(
            "Slit should have at least 5 points. Got {0} only. ".format(len(wslit))
            + "Reduce `wstep=`?"
        )

    return wslit, Islit


try:  # Python >3.6 only
    get_slit_function.__annotations__["shape"] = SLIT_SHAPES
except AttributeError:
    pass  # old Python version

# %% Convolve with slit function


def convolve_with_slit(
    w,
    I,
    w_slit,
    I_slit,
    mode="valid",
    k=1,
    bplot=False,
    verbose=True,
    assert_evenly_spaced=True,
    wunit="",
    waveunit=None,
):
    """Convolves spectrum (w,I) with instrumental slit function (w_slit, I_slit)

        Returns a convolved spectrum on a valid range.

        Check Notes section below for more details about its implementation.

        This function is called directly from a :py:class:`~radis.spectrum.spectrum.Spectrum`
        object with the :py:meth:`~radis.spectrum.spectrum.Spectrum.apply_slit` method.


    .. warning::

        By default, convoluted spectra are thinner than non convoluted spectra,
        as spectra are cropped to remove boundary effects. See ``mode=`` to change
        this behaviour.

    Parameters
    ----------

    w, I: array
        theoretical spectrum (wavespace unit (nm / cm-1))

    w_slit, I_slit: array
        instrumental slit function  (wavespace unit (nm / cm-1))
        It is recommended to truncate the input slit function to its minimum spectral
        extension (see Notes).

        .. warning::

            Both wavespaces have to be the same!

    mode: ``'valid'``, ``'same'``

        ``'same'``:
        returns output of same length as initial spectra,
        but boundary effects are still visible.

        ``'valid'``:
        returns output of length len(spectra) - len(slit) + 1, assuming len(I_slit) < len(I),
        for which lines outside of the calculated range have no impact.

        Default ``'valid'``.

    Other Parameters
    ----------------
    k: int
        order of spline interpolation. 3: cubic, 1: linear. Default 1.
    bplot: boolean
        if ``True``, plots the convolve slit function (for debugging)
    verbose: boolean
        more blabla
    assert_evenly_spaced: boolean, or ``'resample'``
        for the convolution to be accurate, ``w`` should be evenly spaced. If
        ``assert_evenly_spaced=True``, then we check this is the case, and raise
        an error if arrays is not evenly spaced. If ``'resample'``, then we resample
        ``w`` and ``I`` if needed. Recommended, but it takes some time.
    wunit: ``'nm'``, ``'cm-1'``
        just for printing messages. However, ``w`` and ``w_slit`` should be in the
        same wavespace.

    Returns
    -------

    w_conv, I_conv: arrays
        convoluted quantity. w_conv is defined on a valid range
        (sides are cut) if ``mode='valid'``

    Notes
    -----
    The discrete convolution operation is defined as

    .. math:: (I * I_{slit})[n] = \\sum_{m = -\\infty}^{\\infty} I[m] I_{slit}[n - m]

    This numerical discrete convolution doesn't take care of respecting the physical distances.
    Therefore we have to make sure that the ``I`` and ``I_slit`` are expressed using the same
    spectral spacing (e.g. that ``w`` and ``w`_slit`` have the same step in nm between
    two adjacent points). The minimal step between values of ``w`` is used to determine
    this reference spacing.

    If the slit function covers a spectral interval much larger than its
    Full Width a Half Maximum (FWHM), it will potentially add non-zero elements in the sum
    and mix spectral features that are much further away from each other than the FWHM, as
    it would have been expected.
    For an experimental slit function, this also means introducing more noise and artificial
    broadening. It is thus recommended to truncate the input slit function to its minimum
    useful spectral range.

    For instance, if you have an instrumental slit function that looks like so, we recommend
    to cut it where the vertical lines are drawn.

    ::

        I_slit  ^
                |     |   _   |
                |     |  / \\  |
                |     | /   \\ |
               0|_____|/     \\|________
              --|-----         --------->
                    cut       cut        w_slit

    When we apply a slit function, we should consider the contributions of spectral features outside
    of the initial spectral range (i.e. ``w[0]`` - FWHM, ``w[-1]`` + FWHM) that would matter to reproduce
    the entire spectrum. Those points are probably not taken into account during the simulation and
    it could be interesting to extend the spectral range of ``w`` when possible to gauge this contribution.

    Since, we don't have access to this extra points, we truncated the output spectrum to the spectral
    range we can trust (no missing information). The ``valid`` mode operates this truncation by cutting
    the edge of the spectrum ever a spectral width corresponding to the spectral width of the provided
    slit function. This is meant to remove incorrect points on the edge. We recommend using this mode.
    You can however access those removed points using the ``same`` mode (at your own risk).

    Implementation of :func:`~radis.tools.slit.convolve_with_slit` is done in 5 steps:

    - Check input slit function: compare input slit spectral range to spectrum spectral range
      if ranges are comparable, check FHWM and send a warning.
    - Interpolate the slit function on the spectrum grid, resample it if not
      evenly spaced (in order to match physical distances)
    - Check slit aspect, plot slit if asked for
    - Convolve!
    - Remove boundary effects

    See Also
    --------

    :func:`~radis.tools.slit.get_slit_function`,
    :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`

    """

    # 1. Check input
    # --------------

    # Deprecated input:
    if waveunit is not None:
        warn(
            "`waveunit=` parameter in convolve_with_slit is now named `wunit=`",
            DeprecationWarning,
        )
        wunit = waveunit

    # Assert slit function is thin enough

    w_range = abs(w[-1] - w[0])
    w_slit_range = abs(w_slit[-1] - w_slit[0])
    slit_FWHM = get_FWHM(w_slit, I_slit)

    if w_slit_range > 10 * slit_FWHM:
        warn(
            f"FWHM >> spectral range!\n"
            + f"Slit Function is provided with a spectral range {w_slit_range:.1f} nm that is much "
            + f"wider than its estimated FWHM {slit_FWHM:.1f}. We recommend truncating the input "
            + f"slit function. See radis documentation about convolve_with_slit for more information."
        )

    if w_range < w_slit_range:
        if mode == "valid":
            raise AssertionError(
                f"Slit function is provided with a spectral range ({w_slit_range:.1f} {waveunit}) larger than "
                + f"the spectral range of the spectrum ({w_range:.1f} {waveunit}) : the output spectrum will therefore be empty "
                + f"as boundary effects of the convolution are automatically discarded. If you still want to apply the slit "
                + f"despite potential boundary effects, set `mode='same'`. However, we recommend truncating the input slit function if possible, or compute the spectrum on a larger spectral range. "
                + f"See radis documentation about convolve_with_slit for more information."
            )
        else:
            warn(
                f"Slit spectral range > Spectrum spectral range!\n"
                + f"Slit function is provided with a spectral range ({w_slit_range:.1f} nm) larger than "
                + f"the spectral range of the spectrum ({w_range:.1f} nm). We recommend truncating the input "
                + f"slit function. See radis documentation about convolve_with_slit for more information."
            )

    # 2. Interpolate the slit function on the spectrum grid, resample it if not
    #    evenly spaced
    # --------------

    # check first & last spacing, always :
    is_evenly_spaced = evenly_distributed_fast(w, rtolerance=1e-3)
    wstep = w[1] - w[0]

    # ... check with details
    if is_evenly_spaced and assert_evenly_spaced:
        wstep = abs(np.diff(w)).min()  # spectrum wavelength spacing
        is_evenly_spaced = evenly_distributed(w, atolerance=wstep * 1e-3)

    if not is_evenly_spaced:
        if assert_evenly_spaced == "resample":
            w, I = resample_even(w, I, resfactor=2, print_conservation=True)
            wstep = abs(np.diff(w)).min()  # new spectrum wavelength spacing

        raise ValueError(
            "Spectrum is not evenly spaced. Cannot apply slit function. Please resample with `convolve_with_slit(..., convolve_with_slit='resample')` ; or if using a Spectrum object `s`, with `s.resample_even()`"
        )

    # ... Check that the slit is not reversed (interpolation requires objects are sorted)
    reverse = w_slit[-1] < w_slit[0]
    if reverse:
        w_slit = w_slit[::-1]
        I_slit = I_slit[::-1]

    if not np.allclose(np.diff(w_slit), wstep):
        if verbose >= 2:
            # numerical errors can be produced
            print("interpolating slit function over spectrum grid")
        try:
            tck = splrep(w_slit, I_slit, k=k)
        except ValueError:
            # Probably error on input data. Print it before crashing.
            print("\nValueError - Input data below:")
            print("-" * 5)
            print(w_slit)
            print(I_slit)
            print("Check figure")
            plot_slit(w_slit, I_slit, wunit=wunit)
            raise

        # can be a bug here if wstep has the wrong sign.
        w_slit_interp = np.arange(w_slit[0], w_slit[-1] + wstep, wstep)
        # more test needed.
        I_slit_interp = splev(w_slit_interp, tck)
    else:
        w_slit_interp = w_slit
        I_slit_interp = I_slit

    # 3. Check aspect
    # -------------

    # Check no nan
    if anynan(I_slit):
        raise ValueError("Slit has nan value")

    # check slit is positive
    if not (I_slit >= 0).all():
        plot_slit(w_slit, I_slit, wunit="")
        raise ValueError("Slit is partially negative. Check Figure")

    # Plot slit if asked for
    if bplot:
        plot_slit(w_slit, I_slit, wunit="")

    # 4. Convolve!
    # --------------

    # We actually do not use mode valid in np.convolve,
    # instead we use mode=same and remove the same boundaries from I and W in remove_boundary()
    I_conv = oaconvolve(I, I_slit_interp, mode="same") * wstep

    # 5. Remove boundary effects
    # --------------

    w_conv, I_conv = remove_boundary(
        w, I_conv, mode, len_I=len(I), len_I_slit_interp=len(I_slit_interp)
    )

    # reverse back if needed
    # Todo: add test case for that
    if reverse:
        w_slit = w_slit[::-1]
        I_slit = I_slit[::-1]

    return w_conv, I_conv


# %% Slit function methods


def get_FWHM(w, I, return_index=False):
    """Calculate full width half maximum (FWHM) by comparing amplitudes


    Parameters
    ----------

    w, I: arrays

    return_index: boolean
        if True, returns indexes for half width boundaries. Default ``False``

    Returns
    -------

    FWHM: float
        full width at half maximum

    [xmin, xmax]: int


    See Also
    --------

    :py:func:`~radis.tools.slit.get_effective_FWHM`,
    :py:func:`~radis.tools.slit.offset_dilate_slit_function`,
    :py:func:`~radis.tools.slit.normalize_slit`,
    :py:func:`~radis.tools.slit.remove_boundary`,
    :py:func:`~radis.tools.slit.plot_slit`,
    :py:func:`~radis.tools.slit.recenter_slit`,
    :py:func:`~radis.tools.slit.crop_slit`

    """
    # TODO: Linearly interpolate at the boundary? insignificant for large number of points

    half = I.max() / 2
    I_offset = I - half
    upper = np.argwhere(I_offset >= 0)

    xmin = upper.min()
    xmax = upper.max()

    def get_zero(x1, y1, x2, y2):
        """
        Linear interpolation on data to get zero
        Return location of zero by solving for f(x) = (y2-y1)/(x2-x1) * (x-x1) + y1 = 0
        """
        return x1 - y1 * (x2 - x1) / (y2 - y1)

    wmin = get_zero(w[xmin - 1], I_offset[xmin - 1], w[xmin], I_offset[xmin])
    wmax = get_zero(w[xmax], I_offset[xmax], w[xmax + 1], I_offset[xmax + 1])

    if return_index:
        return abs(wmax - wmin), xmin, xmax
    else:
        return abs(wmax - wmin)


def get_effective_FWHM(w, I):
    """Calculate full-width-at-half-maximum (FWHM) of a triangular
    slit of same area and height 1


    Parameters
    ----------

    w, I: arrays

    Returns
    -------

    fwhm: float
        effective FWHM

    See Also
    --------

    :py:func:`~radis.tools.slit.get_FWHM`,
    :py:func:`~radis.tools.slit.offset_dilate_slit_function`,
    :py:func:`~radis.tools.slit.normalize_slit`,
    :py:func:`~radis.tools.slit.remove_boundary`,
    :py:func:`~radis.tools.slit.plot_slit`,
    :py:func:`~radis.tools.slit.recenter_slit`,
    :py:func:`~radis.tools.slit.crop_slit`

    """

    Imax = I.max()

    area = abs(np.trapz(I, w))

    return area / Imax


def offset_dilate_slit_function(
    w_slit_nm, I_slit, w_nm, slit_dispersion, threshold=0.01, verbose=True
):
    """Offset the slit wavelengths ``w_slit_nm`` to the center of range ``w_nm``, and
    dilate them with ``slit_dispersion(w0)``

    Parameters
    ----------

    w_slit_nm: np.array
        wavelength

        .. warning::
            slit_dispersion is assumed to be in wavelength. Warning if you're
            applying this to an array stored in wavenumbers: don't forget to
            convert.

    I_slit: np.array
        slit intensity

    w_nm: np.array
        wavelength whose center is used to center the slit function

    slit_dispersion: func of (lambda), or ``None``
        spectrometer reciprocal function : dλ/dx(λ)
        If not None, then the slit_dispersion function is used to correct the
        slit function for the whole range. Can be important if slit function
        was measured far from the measured spectrum  (e.g: a slit function
        measured at 632.8 nm will look broader at 350 nm because the spectrometer
        dispersion is higher at 350 nm. Therefore it should be corrected)
        Default ``None``. For more details see :func:`~radis.tools.slit.convolve_with_slit`

    Other Parameters
    ----------------

    threshold: float
        if not ``None``, check that slit dispersion is about constant (< ``threshold`` change)
        on the calculated range. Default 0.01 (1%)

    Returns
    -------

    w_slit, I_slit: np.array
        dilated wavelength (or wavenumber), slit intensity (renormalized)

    See Also
    --------

    :py:func:`~radis.tools.slit.get_FWHM`,
    :py:func:`~radis.tools.slit.get_effective_FWHM`,
    :py:func:`~radis.tools.slit.normalize_slit`,
    :py:func:`~radis.tools.slit.remove_boundary`,
    :py:func:`~radis.tools.slit.plot_slit`,
    :py:func:`~radis.tools.slit.recenter_slit`,
    :py:func:`~radis.tools.slit.crop_slit`

    """
    w0 = w_nm[len(w_nm) // 2]
    wslit0 = w_slit_nm[len(w_slit_nm) // 2]
    slit_disp_0 = slit_dispersion(wslit0)
    range_nm = abs(w_slit_nm[-1] - w_slit_nm[0]) / 2
    w_slit_init = w_slit_nm

    # Check that slit dispersion is about constant (<1% change) on the calculated range
    if threshold:
        if (
            not 1 - threshold
            <= slit_dispersion(
                w_nm.max() - range_nm / slit_disp_0 * slit_dispersion(w_nm.max())
            )
            / slit_dispersion(
                w_nm.min() + range_nm / slit_disp_0 * slit_dispersion(w_nm.min())
            )
            <= 1 + threshold
        ):
            warn(
                "Slit dispersion changes slightly ({2:.2f}%) between {0:.3f} and {1:.3f}nm".format(
                    w_nm.min(),
                    w_nm.max(),
                    abs(slit_dispersion(w_nm.max()) / slit_dispersion(w_nm.min()) - 1)
                    * 100,
                )
                + ". Consider splitting your spectrum",
                SlitDispersionWarning,
            )

    # Offset slit and correct for dispersion
    w_slit_nm = w0 + slit_dispersion(w0) / slit_dispersion(wslit0) * (
        w_slit_nm - wslit0
    )

    if verbose > 2:
        print(
            "{0:.2f} to {1:.2f}nm: slit function FWHM changed from {2:.2f} to {3:.2f}".format(
                wslit0,
                w0,
                get_effective_FWHM(w_slit_init, I_slit),
                get_effective_FWHM(w_slit_nm, I_slit),
            )
        )

    return w_slit_nm, I_slit


def normalize_slit(w_slit, I_slit, norm_by="area"):
    """Normalize slit function with different normalization modes. Warning,
    some change units after convolution!

    Parameters
    ----------

    w_slit, I_slit: np.array
        wavelength and slit intensity

    norm_by: ``'area'``, ``'max'``, or ``None``
        how to normalize. ``'area'`` conserves energy. With ``'max'`` the slit is normalized
        at peak so that the maximum is one.

        .. note::

            ``'max'`` changes the unit of the spectral array, e.g. from
            ``'mW/cm2/sr/µm'`` to ``'mW/cm2/sr')``

        ``None`` doesnt normalize. Default ``'area'``

    Returns
    -------

    w_slit, I_slit: np.array
        wavelength, and normalized slit intensity

    See Also
    --------

    :py:func:`~radis.tools.slit.get_FWHM`,
    :py:func:`~radis.tools.slit.get_effective_FWHM`,
    :py:func:`~radis.tools.slit.offset_dilate_slit_function`,
    :py:func:`~radis.tools.slit.remove_boundary`,
    :py:func:`~radis.tools.slit.plot_slit`,
    :py:func:`~radis.tools.slit.recenter_slit`,
    :py:func:`~radis.tools.slit.crop_slit`

    """

    # Renormalize
    # ---------

    if norm_by == "area":  # normalize by the area
        I_slit = I_slit / abs(np.trapz(I_slit, x=w_slit))
    elif norm_by == "max":  # set maximum to 1
        I_slit = I_slit / abs(np.max(I_slit))
    elif norm_by is None:
        pass
    else:
        raise ValueError("Unknown normalization type: `norm_by` = {0}".format(norm_by))

    return w_slit, I_slit


def remove_boundary(
    w, I_conv, mode, len_I=None, len_I_slit_interp=None, crop_left=None, crop_right=None
):
    """Crop convoluted array to remove boundary effects

    Parameters
    ----------

    w, I_conv: numpy array
        wavelength and intensity of already convoluted quantity
        (ex: ``radiance``)

    mode: ``'valid'``, ``'same'``, ``''crop'``
        ``'same'`` returns output of same length as initial spectra,
        but boundary effects are still visible. ``'valid'`` returns
        output of length len(spectra) - len(slit) + 1, for
        which lines outside of the calculated range have
        no impact. ``'crop'`` just removes ``crop_wings`` points on
        the side.

    Other Parameters
    ----------------

    len_I, len_I_slit_interp: int
        length of initial quantity (ex: ``radiance_noslit``) and length
        of slit function intensity, before convolution.
        Needed to determine the valid range if ``mode='valid'``

    crop_left, crop_right: int
        number of points to discard on each side if ``mode='crop'``
        Values are replaced with ``nan``

    Returns
    -------

    w_conv, I_conv: numpy arrays
        cropped waverange and quantity

    Notes
    -----

    .. note::
        new in 0.9.30 : remove_boundary replaces off-range values with ``nan``
        but keeps the same array size.

    See Also
    --------

    :py:func:`~radis.tools.slit.get_FWHM`,
    :py:func:`~radis.tools.slit.get_effective_FWHM`,
    :py:func:`~radis.tools.slit.offset_dilate_slit_function`,
    :py:func:`~radis.tools.slit.normalize_slit`,
    :py:func:`~radis.tools.slit.plot_slit`,
    :py:func:`~radis.tools.slit.recenter_slit`,
    :py:func:`~radis.tools.slit.crop_slit`

    """

    # Remove boundary effects by adding nans; keep the same x-axis
    if mode == "valid":
        la = min(len_I, len_I_slit_interp)
        a = int((la - 1) / 2)
        b = int((la) / 2)
        # I_conv = I[a:-b]    # former version : we would change the array size
        # w_conv = w[a:-b]
        I_conv[:a] = np.nan
        I_conv[-b:] = np.nan
        w_conv = w
    elif mode == "same":
        I_conv = I_conv
        w_conv = w
    elif mode == "crop":
        if crop_right == 0:
            _crop_right = None
        else:
            _crop_right = -crop_right
        # l = len(I_conv)
        # I_conv = I_conv[crop_left:_crop_right]
        # w_conv = w[crop_left:_crop_right]
        # assert len(I_conv) == l - crop_left - crop_right
        I_conv[:crop_left] = np.nan
        I_conv[_crop_right:] = np.nan
    else:
        raise ValueError("Unexpected mode: {0}".format(mode))

    return w_conv, I_conv


def plot_slit(
    w,
    I=None,
    wunit="",
    plot_unit="same",
    Iunit=None,
    warnings=True,
    ls="-",
    title=None,
    waveunit=None,
):
    """Plot slit, calculate and display FWHM, and calculate effective FWHM.
    FWHM is calculated from the limits of the range above the half width,
    while FWHM is the equivalent width of a triangular slit with the same area


    Parameters
    ----------
    w, I: arrays    or   (str, None)
        if str, open file directly
    waveunit: ``'nm'``, ``'cm-1'`` or ``''``
        unit of input w
    plot_unit: ``'nm'``, ``'cm-1'`` or ``'same'``
        change plot unit (and FWHM units)
    Iunit: str, or None
        give Iunit
    warnings: boolean
        if True, test if slit is correctly centered and output a warning if it
        is not. Default ``True``

    Returns
    -------
    fix, ax: matplotlib objects
        figure and ax

    Examples
    --------

    .. minigallery:: radis.plot_slit

    See Also
    --------

    :py:func:`~radis.tools.slit.get_FWHM`,
    :py:func:`~radis.tools.slit.get_effective_FWHM`,
    :py:func:`~radis.tools.slit.offset_dilate_slit_function`,
    :py:func:`~radis.tools.slit.normalize_slit`,
    :py:func:`~radis.tools.slit.remove_boundary`,
    :py:func:`~radis.tools.slit.recenter_slit`,
    :py:func:`~radis.tools.slit.crop_slit`

    """

    # Deprecated input:
    if waveunit is not None:
        warn(
            "`waveunit=` parameter in convolve_with_slit is now named `wunit=`",
            DeprecationWarning,
        )
        wunit = waveunit

    import matplotlib.pyplot as plt

    from radis.misc.plot import set_style

    set_style()

    try:
        from radis.plot.toolbar import add_tools

        add_tools()  # includes a Ruler to measure slit
    except:
        pass

    # Check input
    if isinstance(w, str) and I is None:
        w, I = np.loadtxt(w).T
    assert len(w) == len(I)
    if anynan(I):
        warn("Slit function has nans")
        w = w[~np.isnan(I)]
        I = I[~np.isnan(I)]
    assert len(I) > 0

    # cast units
    wunit = cast_waveunit(wunit, force_match=False)
    plot_unit = cast_waveunit(plot_unit, force_match=False)
    if plot_unit == "same":
        plot_unit = wunit

    # Convert wavespace unit if needed
    elif wunit == "cm-1" and plot_unit == "nm":  # wavelength > wavenumber
        w = cm2nm(w)
        wunit = "nm"
    elif wunit == "nm" and plot_unit == "cm-1":  # wavenumber > wavelength
        w = nm2cm(w)
        wunit = "cm-1"
    elif wunit == plot_unit:  # same units
        pass
    elif plot_unit == "":  # ???
        pass
    else:
        raise ValueError("Unknown plot unit: {0}".format(plot_unit))

    # Recalculate FWHM
    FWHM, xmin, xmax = get_FWHM(w, I, return_index=True)
    FWHM_eff = get_effective_FWHM(w, I)

    # Get labels
    if plot_unit == "nm":
        xlabel = "Wavelength (nm)"
    elif plot_unit == "cm-1":
        xlabel = "Wavenumber (cm-1)"
    elif plot_unit == "":
        xlabel = "Wavespace"
    else:
        raise ValueError("Unknown unit for plot_unit: {0}".format(plot_unit))
    ylabel = "Slit function"
    if Iunit is not None:
        ylabel += " ({0})".format(Iunit)

    fig, ax = plt.subplots()
    ax.plot(w, I, "o", color="lightgrey")
    ax.plot(
        w,
        I,
        "k",
        ls=ls,
        label="FWHM: {0:.3f} {1}".format(FWHM, plot_unit)
        + "\nEff. FWHM: {0:.3f} {1}".format(FWHM_eff, plot_unit)
        + "\nArea: {0:.3f}".format(abs(np.trapz(I, x=w))),
    )

    # Vertical lines on center, and FWHM
    plt.axvline(w[len(w) // 2], ls="-", lw=2, color="lightgrey")  # center
    plt.axvline(
        w[(xmin + xmax) // 2], ls="--", color="k", lw=0.5
    )  # maximum (should be center)
    plt.axvline(w[xmin], ls="--", color="k", lw=0.5)  # FWHM min
    plt.axvline(w[xmax], ls="--", color="k", lw=0.5)  # FWHM max
    plt.axhline(I.max() / 2, ls="--", color="k", lw=0.5)  # half maximum

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.legend(loc="best", prop={"size": 16})

    if title:
        plt.title(title)

    # extend axis:
    fig.tight_layout()
    xlmin, xlmax = ax.get_xlim()
    ax.set_xlim((xlmin - 0.5, xlmax + 0.5))

    if warnings:
        if w[(xmin + xmax) // 2] != w[len(w) // 2]:
            warn(
                "Slit function doesnt seem centered: center measured with FWHM"
                + " is not the array center (shift: {0:.3f}{1}): This can induce offsets!".format(
                    abs(w[(xmin + xmax) // 2] - w[len(w) // 2]), wunit
                )
            )

        if I[0] != 0 or I[-1] != 0:
            warn("Slit function should have zeros on both sides")

    return fig, ax


def recenter_slit(w_slit, I_slit, verbose=True):
    """Recenters the slit on the maximum calculated from the two
    FWHM limits. To recenter, zeros are added on the shorter side (zero padding)

    See Also
    --------

    :py:func:`~radis.tools.slit.get_FWHM`,
    :py:func:`~radis.tools.slit.get_effective_FWHM`,
    :py:func:`~radis.tools.slit.offset_dilate_slit_function`,
    :py:func:`~radis.tools.slit.normalize_slit`,
    :py:func:`~radis.tools.slit.remove_boundary`,
    :py:func:`~radis.tools.slit.plot_slit`,
    :py:func:`~radis.tools.slit.crop_slit`

    """

    _, xmin, xmax = get_FWHM(w_slit, I_slit, return_index=True)
    xcenter = (xmax + xmin) / 2
    offset = len(w_slit) / 2 - xcenter
    add_zeros = int(2 * offset)
    if add_zeros > 0:
        # add zeros on the left (assume regular step)
        wstep = w_slit[1] - w_slit[0]
        w_left = np.linspace(
            w_slit[0] - wstep * (add_zeros + 1), w_slit[0] - wstep, add_zeros
        )
        w_slit = np.hstack((w_left, w_slit))
        I_slit = np.hstack((np.zeros(add_zeros), I_slit))
        if verbose:
            print("... added zeros to the slit to center it")
    elif add_zeros < 0:
        add_zeros = -add_zeros
        # add zeros on the right (assume regular step)
        wstep = w_slit[-1] - w_slit[-2]
        w_right = np.linspace(
            w_slit[-1] + wstep, w_slit[-1] + wstep * (add_zeros + 1), add_zeros
        )
        w_slit = np.hstack((w_slit, w_right))
        I_slit = np.hstack((I_slit, np.zeros(add_zeros)))
        if verbose:
            print("... added zeros to the slit to center it")
    else:
        pass

    return w_slit, I_slit


def crop_slit(w_slit, I_slit, verbose=True):
    """Removes unnecessary zeros on the side for a faster convolution.
    (remove as many zeros on the left as on the right).

    See Also
    --------

    :py:func:`~radis.tools.slit.get_FWHM`,
    :py:func:`~radis.tools.slit.get_effective_FWHM`,
    :py:func:`~radis.tools.slit.offset_dilate_slit_function`,
    :py:func:`~radis.tools.slit.normalize_slit`,
    :py:func:`~radis.tools.slit.remove_boundary`,
    :py:func:`~radis.tools.slit.plot_slit`,
    :py:func:`~radis.tools.slit.recenter_slit`

    """

    nzeros_index = np.argwhere(I_slit != 0)
    zeros_left = nzeros_index.min()
    zeros_right = len(I_slit) - nzeros_index.max() - 1
    remove = 0
    if zeros_left > 1 and zeros_right > zeros_left:
        remove = zeros_left - 1
    elif zeros_right > 1 and zeros_left > zeros_right:
        remove = zeros_right - 1
    if remove > 0:
        w_slit = w_slit[remove:-remove]
        I_slit = I_slit[remove:-remove]
        if verbose:
            print("... removed {0} zeros to the slit on each side".format(remove))

    return w_slit, I_slit


# %% Slit function models


def import_experimental_slit(
    slit,
    norm_by="area",
    bplot=False,
    wunit="nm",
    auto_recenter=True,
    auto_crop=True,
    center=None,
    scale=1,
    verbose=True,
    waveunit=None,
):
    """Import instrumental slit function and normalize it


    Parameters
    ----------

    slit: str or numpy object (wslit, Islit)
        slit function spectrum.
        If ``str``, the slit is loaded from the corresponding file with::

            w_slit, I_slit = np.loadtxt(slit).T

    norm_by: ``None``, ``'area'``, ``'max'``
        normalisation type. ``'area'`` divides by the area (conserves energy / units
        are unchanged). ``'max'`` divides by the peak value (similar to what is
        done in Specair / changes units). ``None`` does not normalize the
        experimental slit function. Default ``'area'``

    bplot: boolean
        plot normalized slit function (for debugging). Default ``False``

    wunit: ``'nm'``, ``'cm-1'``
        used for plot only. Slit function is generated assuming you use the
        correct wavespace. No conversions are made here. Default ``'nm'``

    Other Parameters
    ----------------

    auto_recenter: bool
        if True, recenters the slit on the maximum calculated from the two
        FWHM limits. To recenter, zeros are added on the shorter side (zero padding)
        Default ``True``

    auto_crop: bool
        If True, remove unnecessary zeros on the side for a faster convolution.
        (remove as many zeros on the left as on the right). Default ``True``

    center: float, None
        normally the slit instrumental slit function is centered on where it
        was measured. If not None, center overwrites the center and
        offsets the slit to the given value. Default ``None``

    scale: float
        multiply slit by an arbitrary factor. Default 1.

    verbose: boolean
        Display messages

    Returns
    -------

    w_slit, I_slit: numpy arrays
        wavelength and intensity of normalized slit function

    See Also
    --------

    :py:func:`~radis.tools.slit.triangular_slit`,
    :py:func:`~radis.tools.slit.trapezoidal_slit`,
    :py:func:`~radis.tools.slit.gaussian_slit`


    """

    # Deprecated input:
    if waveunit is not None:
        warn(
            "`waveunit=` parameter in convolve_with_slit is now named `wunit=`",
            DeprecationWarning,
        )
        wunit = waveunit

    # import
    if isinstance(slit, str):
        try:
            w_slit, I_slit = np.loadtxt(slit).T
        except ValueError:
            w_slit, I_slit = np.loadtxt(slit)
    # numpy input
    elif isinstance(slit, np.ndarray):
        a, b = np.shape(slit)
        if a == 2:
            w_slit, I_slit = slit
        elif b == 2:
            w_slit, I_slit = slit.T
        else:
            raise ValueError("Wrong format of slit_function. Should be 2 x n")
    else:
        raise TypeError("Unexpected type for slit function: {0}".format(type(slit)))

    assert np.shape(w_slit) == np.shape(I_slit)

    # check sides are zeros
    if not I_slit[0] == 0 and I_slit[-1] == 0:
        raise ValueError("Slit function must be null on each side. Fix it")

    # recenter if asked for
    if auto_recenter:
        w_slit, I_slit = recenter_slit(w_slit, I_slit, verbose=verbose)
        # note: if auto_crop is true we may be adding zeros just to remove them
        # right way. Fix that. Someday.

    # remove unecessary zeros
    if auto_crop:  # (note that we know I_slit has zeros on each side already)
        w_slit, I_slit = crop_slit(w_slit, I_slit, verbose=verbose)

    # Offset slit if needed
    if center is not None:
        w_slit += center - w_slit[len(w_slit) // 2]

    # Normalize
    if norm_by == "area":  # normalize by the area
        #        I_slit /= np.trapz(I_slit, x=w_slit)
        I_slit /= abs(np.trapz(I_slit, x=w_slit))
        Iunit = "1/{0}".format(waveunit)
    elif norm_by == "max":  # set maximum to 1
        I_slit /= abs(np.max(I_slit))
        Iunit = ""
    elif norm_by is None:
        Iunit = None
    else:
        raise ValueError("Unknown normalization type: `norm_by` = {0}".format(norm_by))

    # scale
    I_slit *= scale
    #    if Iunit is not None and scale != 1:
    #        Iunit += 'x{0:.2f}'.format(scale)

    # Plot slit
    if bplot:
        plot_slit(w_slit, I_slit, wunit=wunit, Iunit=Iunit)

    return w_slit, I_slit


def triangular_slit(
    FWHM,
    wstep,
    center=0,
    norm_by="area",
    bplot=False,
    wunit="",
    scale=1,
    footerspacing=0,
    waveunit=None,  # Deprecated
):
    r""" Generate normalized slit function


    Parameters
    ----------

    FWHM: (nm)
        full-width at half maximum

    wstep: (nm)
        wavelength step

    center: (nm)
        center wavelength for the wavelength axs of the slit function

    norm_by: ``'area'``, ``'max'``
        normalisation type. ``'area'`` conserves energy. ``'max'`` is what is
        done in Specair and changes units. Default ``'area'``

    bplot: boolean
        plot normalized slit function (for debugging). Default ``False``

    wunit: '', 'nm', 'cm-1'
        used for plot only. Slit function is generated assuming you use the
        correct wavespace. No conversions are made here.

    scale: float
        multiply slit by an arbitrary factor. Default 1.

    footerspacing: int
        spacing (footer) on left and right. Default 10.


    Returns
    -------

    w, I: numpy arrays
        Slit function of shape::

                       ^
                      / \
                 _ _ /   \ _ _
                :   :  :  :   :
               -a-f -a 0  a  a+f


                with FWHM = (a+1/2)*wstep, `f` spacing on left & right


    Notes
    -----

    slit is generated with an odd number of elements and centered on its
    maximum. This may result in sightly different FWHM if wstep is too large.
    However, slit is corrected so that effective FWHM is preserved.

    See Also
    --------

    :py:func:`~radis.tools.slit.import_experimental_slit`,
    :py:func:`~radis.tools.slit.trapezoidal_slit`,
    :py:func:`~radis.tools.slit.gaussian_slit`

    """

    # Deprecated input:
    if waveunit is not None:
        warn(
            "`waveunit=` parameter in convolve_with_slit is now named `wunit=`",
            DeprecationWarning,
        )
        wunit = waveunit

    # Build first half
    slope = 1 / (FWHM - wstep)
    a = int(FWHM / wstep)
    I = 1 - np.arange(0, a) * slope * wstep  # slope

    # Zeros
    if FWHM % wstep:  # add one extra zero when not a divider
        footerspacing += 1
    f = int(footerspacing)  # add zeros (footer)
    I = np.hstack((I, np.zeros(f)))
    w = np.linspace(0, len(I) * wstep, len(I))

    # Mirror to get second half
    I = np.hstack((I[1:][::-1], I))
    w = np.hstack((-w[1:][::-1], w)) + center

    # Normalize
    if norm_by == "area":  # normalize by the area
        I /= np.trapz(I, x=w)
        Iunit = "1/{0}".format(wunit)
    elif norm_by == "max":  # set maximum to 1
        I /= np.max(I)
        Iunit = ""
    elif norm_by is None:
        Iunit = None
    else:
        raise ValueError("Unknown normalization type: `norm_by` = {0}".format(norm_by))

    # Scale
    I *= scale
    #    if Iunit is not None and scale != 1:
    #        Iunit += 'x{0:.2f}'.format(scale)

    # Plot slit
    if bplot:
        plot_slit(w, I, wunit=wunit, Iunit=Iunit)

    return w, I


# trapezoidal instrumental broadening function of base base nm and top top nm


def trapezoidal_slit(
    top,
    base,
    wstep,
    center=0,
    norm_by="area",
    bplot=False,
    wunit="",
    scale=1,
    footerspacing=0,
    waveunit=None,
):
    r""" Build a trapezoidal slit. Remember that FWHM = (top + base) / 2


    Parameters
    ----------
    top: (nm)
        top of the trapeze
    base: (nm)
        base of the trapeze
    wstep: (nm)
        wavelength step

    Other Parameters
    ----------------
    center: (nm)
        center wavelength for the wavelength axs of the slit function
    norm_by: ``'area'``, ``'max'``
        normalisation type. ``'area'`` conserves energy. ``'max'`` is what is
        done in Specair and changes units. Default ``'area'``
    bplot: boolean
        plot normalized slit function (for debugging). Default ``False``
    waveunit: '', 'nm', 'cm-1'
        used for plot only. Slit function is generated assuming you use the
        correct wavespace. No conversions are made here.
    scale: float
        multiply slit by an arbitrary factor. Default 1.
    footerspacing: int
        spacing (footer) on left and right. Default 10.


    Returns
    -------

    w, I: numpy arrays
        Slit function of shape::

                     _______
                    /       \
                   / :     : \
             _ _ _/  :     :  \____
             0    :  :     :  :   0
             :   -b  :     :  b   :
           -b-f     -t     t     b+f


            with `t`and `b` the half-top and half-base in number of elements, `f` spacing
            on left & right. FWHM = (t+b+1)*wstep


    Notes
    -----

    slit is generated with an odd number of elements and centered on its
    maximum. This may result in sightly different FWHM if wstep is too large.
    However, slit is corrected so that effective FWHM is preserved.

    See Also
    --------

    :py:func:`~radis.tools.slit.import_experimental_slit`,
    :py:func:`~radis.tools.slit.triangular_slit`,
    :py:func:`~radis.tools.slit.gaussian_slit`

    """

    if top > base:
        top, base = base, top

    FWHM = (base + top) / 2
    b = 2 * int(top / wstep // 2) + 1  # number of points on top (even)

    # Build first half
    slope = 1 / (FWHM - b * wstep)
    a = int(FWHM / wstep) - b  # number of points in slope
    I = 1 - np.arange(0, a + 1) * slope * wstep  # slope

    if len(I) == 0:
        I = np.ones(1)

    if I[-1] < 0:
        I[-1] = 0
    elif I[-1] > 0:
        # add one extra zero
        footerspacing += 1

    # Zeros
    #    if abs((base-top) % wstep) < wstep:  # add one extra zero when a divider
    #        footerspacing += 1

    f = int(footerspacing)  # add zeros (footer)
    I = np.hstack((I, np.zeros(f)))
    #    w = np.linspace(0, len(I)*wstep, len(I))

    # Mirror to get second half, add top
    I = np.hstack((I[1:][::-1], np.ones(b), I[1:]))
    w = wstep * np.linspace(-len(I) / 2, len(I) / 2, len(I)) + center

    # Normalize
    if norm_by == "area":  # normalize by the area
        I /= np.trapz(I, x=w)
        Iunit = "1/{0}".format(wunit)
    elif norm_by == "max":  # set maximum to 1
        I /= np.max(I)
        Iunit = ""
    elif norm_by is None:
        Iunit = None
    else:
        raise ValueError("Unknown normalization type: `norm_by` = {0}".format(norm_by))

    # scale
    I *= scale
    #    if Iunit is not None and scale != 1:
    #        Iunit += 'x{0:.2f}'.format(scale)

    # Plot slit
    if bplot:
        plot_slit(w, I, wunit=wunit, Iunit=Iunit)

    return w, I


def gaussian_slit(
    FWHM,
    wstep,
    center=0,
    norm_by="area",
    bplot=False,
    wunit="",
    calc_range=4,
    scale=1,
    footerspacing=0,
    waveunit=None,  # Deprecated
):
    r""" Generate normalized slit function


    Parameters
    ----------
    FWHM: (nm)
        full-width at half maximum
    wstep: (nm)
        wavelength step

    Other Parameters
    ----------------
    center: (nm)
        center wavelength for the wavelength axs of the slit function
    norm_by: ``'area'``, ``'max'``
        normalisation type. ``'area'`` conserves energy. ``'max'`` is what is
        done in Specair and changes units. Default ``'area'``
    bplot: boolean
        plot normalized slit function (for debugging). Default ``False``
    wunit: '', 'nm', 'cm-1'
        used for plot only. Slit function is generated assuming you use the
        correct wavespace. No conversions are made here.
    calc_range: float   (number of sigma)
        function broadening range expressed in number of standard deviation
    scale: float
        multiply slit by an arbitrary factor. Default 1.
    footerspacing: int
        spacing (footer) on left and right. Default 10


    Returns
    -------
    w, I: numpy arrays
        Slit function of shape::

                     .-.
                    /   \
                   /     \
             _ _.-'       `-._ _
             0    [gaussian]   0
             :  :     :     :  :
           -c-f -c    0     c  c+f          (c : calc_range)


            with `c` cutoff in number of elements, `f` spacing on left & right


    Notes
    -----

    slit is generated with an odd number of elements and centered on its
    maximum. This may result in sightly different FWHM if wstep is too large.

    See Also
    --------

    :py:func:`~radis.tools.slit.import_experimental_slit`,
    :py:func:`~radis.tools.slit.triangular_slit`,
    :py:func:`~radis.tools.slit.trapezoidal_slit`


    """

    f = int(footerspacing)  # spacing (footer) on left and right
    sigma = FWHM / 2 / sqrt(2 * ln(2))

    # half-base in number of elements (even number)
    a = 2 * int(calc_range * sigma // wstep / 2)

    # 2 sigma: gaussian non calculated residual: 5%
    # 3 sigma: gaussian non calculated residual: 1%

    w0 = wstep * np.linspace(-a, a, 2 * a + 1)  # centered
    Igauss = exp(-(w0**2) / (2 * sigma**2))

    I = np.hstack((np.zeros(f), Igauss, np.zeros(f)))
    w = wstep * np.linspace(-(a + f), (a + f), 2 * a + 2 * f + 1) + center

    # Normalize
    if norm_by == "area":  # normalize by the area
        I /= np.trapz(I, x=w)
        Iunit = "1/{0}".format(wunit)
    elif norm_by == "max":  # set maximum to 1
        I /= np.max(I)
        Iunit = ""
    elif norm_by is None:
        Iunit = None
    else:
        raise ValueError("Unknown normalization type: `norm_by` = {0}".format(norm_by))

    # scale
    I *= scale
    #    if Iunit is not None and scale != 1:
    #        Iunit += 'x{0:.2f}'.format(scale)

    # Plot slit
    if bplot:
        plot_slit(w, I, wunit=wunit, Iunit=Iunit)

    return w, I


# %% Test
if __name__ == "__main__":
    from radis.test.tools.test_slit import _run_testcases

    print("Testing slit.py: ", _run_testcases())
