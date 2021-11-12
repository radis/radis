# -*- coding: utf-8 -*-
"""
Summary
-------
Classes to deal with multi-slabs configurations:

- :func:`~radis.los.slabs.MergeSlabs` for several species at the same spatial position
- :func:`~radis.los.slabs.SerialSlabs` to add several spectra along the line-of-path
One Slab is just a :class:`~radis.spectrum.spectrum.Spectrum` object
Examples
--------
See more examples in the
`RADIS line-of-sight module <https://radis.readthedocs.io/en/latest/los/index.html>`__
-------------------------------------------------------------------------------
"""
# Todo:
#
# - transport emisscoeff too
# - emisscoeff default unit


from warnings import warn

import numpy as np
from numpy import abs, allclose, arange, diff

from radis.misc.arrays import anynan
from radis.misc.basics import in_all, intersect, merge_lists
from radis.misc.debug import printdbg
from radis.spectrum.spectrum import Spectrum

# %% Slabs / Multi-layers / Radiative Transfer Equation (RTE)
# ----------------------------------------------------------------------


def SerialSlabs(*slabs, **kwargs) -> Spectrum:
    # type: (*Spectrum, **dict) -> Spectrum
    r"""Adds several slabs along the line-of-sight.
    If adding two slabs only, you can also use::

        s1>s2
    Serial spectrum ``1>2`` is calculated with Eqn (4.2) of the [RADIS-2018]_ article,
    generalized to N slabs :

    .. math::
        I_{\lambda, 1>2} = I_{\lambda, 1} \tau_{\lambda, 2} + I_{\lambda, 2}

        \tau_{\lambda, 1+2} = \tau_{\lambda, 1} \cdot \tau_{\lambda, 2}

    where

        .. math:: I_{\lambda}, \tau_{\lambda}

    are the radiance and transmittance of the two slabs ``1`` and ``2``.
    Radiance and transmittance are calculated if not given in the
    initial slabs (if possible).

    Parameters
    ----------
    slabs: list of Spectra, each representing a slab
        line-of-sight::

            slabs       [0]     [1]  ............... [n]
                         :       :                    :         \====
            light        *   ->  *        ->          *    ->    )===  observer
                                                                /====
    resample_wavespace: ``'never'``, ``'intersect'``, ``'full'``
        what to do when spectra have different wavespaces:

        - If ``'never'``, raises an error
        - If ``'intersect'``, uses the intersection of all ranges, and resample
          spectra on the most resolved wavespace.
        - If ``'full``', uses the overlap of all ranges, resample spectra on the
          most resolved wavespace, and fill missing data with 0 emission and 0
          absorption

        Default ``'never'``

    out: ``'transparent'``, ``'nan'``, ``'error'``
        what to do if resampling is out of bounds:

        - ``'transparent'``: fills with transparent medium.
        - ``'nan'``: fills with nan.
        - ``'error'``: raises an error.

        Default ``'nan'``

    Other Parameters
    ----------------
    verbose: bool
        if ``True``, more blabla. Default ``False``
    modify_inputs: False
        if ``True``, slabs wavelengths/wavenumbers are modified directly when
        they are resampled. This avoids making a copy so it is slightly faster.
        Default ``False``.

        .. note::
            for large number of slabs (in radiative transfer calculations) you
            surely want to use this option !

    Returns
    -------
    Spectrum: object representing total emission and total transmittance as
        observed at the output (slab[n+1]). Conditions and units are transported too,
        unless there is a mismatch then conditions are dropped (and units mismatch
        raises an error because it doesnt make sense)

    Examples
    --------
    Add s1 and s2 along the line of sight: s1 --> s2::

        s1 = calc_spectrum(...)
        s2 = calc_spectrum(...)
        s3 = SerialSlabs(s1, s2)

    The last line is equivalent to::

        s3 = s1>s2

    .. minigallery:: radis.SerialSlabs

    See Also
    --------
    :func:`~radis.los.slabs.MergeSlabs`,
    :ref:`See more examples in Line-of-Sight module <label_los_index>`

    """
    # TODO: rewrite with 'recompute' list like in MergeSlabs ?

    if "resample_wavespace" in kwargs:
        warn(DeprecationWarning("'resample_wavespace' replaced with 'resample'"))
        kwargs["resample"] = kwargs.pop("resample_wavespace")
    if "out_of_bounds" in kwargs:
        warn(DeprecationWarning("'out_of_bounds' replaced with 'out'"))
        kwargs["out"] = kwargs.pop("out_of_bounds")

    # Check inputs, get defaults
    resample_wavespace = kwargs.pop("resample", "never")  # default 'never'
    out_of_bounds = kwargs.pop("out", "nan")  # default 'nan'
    verbose = kwargs.pop("verbose", False)  # type: bool
    modify_inputs = kwargs.pop("modify_inputs", False)  # type: bool
    if len(kwargs) > 0:
        raise ValueError("Unexpected input: {0}".format(list(kwargs.keys())))
    if resample_wavespace not in ["never", "intersect", "full"]:
        raise ValueError(
            "resample should be one of: {0}".format(
                ", ".join(["never", "intersect", "full"])
            )
        )

    if len(slabs) == 0:
        raise ValueError("Empty list of slabs")

    elif len(slabs) == 1:
        if not isinstance(slabs[0], Spectrum):
            raise TypeError(
                "SerialSlabs takes an unfolded list of Spectrum as "
                + "argument: *list (got {0})".format(type(slabs[0]))
            )
        return slabs[0]

    else:
        # recursively calculate serial slabs
        slabs = list(slabs)

        # Recursively deal with the rest of Spectra --> call it s
        sn = slabs.pop(-1)  # type: Spectrum
        _check_valid(sn)  # check it is a spectrum
        s = SerialSlabs(
            *slabs,
            resample=resample_wavespace,
            out=out_of_bounds,
            modify_inputs=modify_inputs
        )

        # Now calculate sn and s in Serial
        quantities = {}
        unitsn = sn.units

        # make sure we use the same wavespace type (even if sn is in 'nm' and s in 'cm-1')
        # also make sure we use the same units
        waveunit = s.get_waveunit()

        # Make all our slabs copies with the same wavespace range
        # (note: wavespace range may be different for different quantities, but
        # equal for all slabs)
        s, sn = resample_slabs(
            waveunit, resample_wavespace, out_of_bounds, modify_inputs, s, sn
        )
        try:
            w = s._q["wavespace"]
        except KeyError:
            raise KeyError(
                "Cannot calculate the RTE if non convoluted quantities "
                + "are not defined. Got: {0}".format(s.get_vars())
            )

        # Get all data
        # -------------

        I, In, T, Tn = None, None, None, None

        # To make it easier, the radiative transfer equation is solved with 'radiance_noslit' and
        # 'transmittance_noslit' only. Here we first try to get these quantities:

        # ... get sn quantities
        try:
            sn.update("transmittance_noslit", verbose=verbose)
        except ValueError:
            pass
        else:
            Tn = sn.get(
                "transmittance_noslit",
                wunit=waveunit,
                Iunit=unitsn["transmittance_noslit"],
                copy=False,
            )[1]
        try:
            sn.update("radiance_noslit", verbose=verbose)
        except ValueError:
            pass
        else:
            In = sn.get(
                "radiance_noslit",
                wunit=waveunit,
                Iunit=unitsn["radiance_noslit"],
                copy=False,
            )[1]
        # ... get s quantities
        try:
            s.update("transmittance_noslit", verbose=verbose)
        except ValueError:
            pass
        else:
            T = s.get(
                "transmittance_noslit",
                wunit=waveunit,
                Iunit=unitsn["transmittance_noslit"],
                copy=False,
            )[1]
        try:
            s.update("radiance_noslit", verbose=verbose)
        except ValueError:
            pass
        else:
            I = s.get(
                "radiance_noslit",
                wunit=waveunit,
                Iunit=unitsn["radiance_noslit"],
                copy=False,
            )[1]

        # Solve radiative transfer equation
        # ---------------------------------

        if I is not None and In is not None:
            # case where we may use SerialSlabs just to compute the products of all transmittances
            quantities["radiance_noslit"] = (w, I * Tn + In)

        if T is not None:  # note that we dont need the transmittance in the inner
            # slabs to calculate the total radiance
            quantities["transmittance_noslit"] = (w, Tn * T)

        # Update conditions
        # -----------------

        # Get conditions (if they're different, fill with 'N/A')
        conditions = intersect(s.conditions, sn.conditions)
        conditions["waveunit"] = waveunit
        if "thermal_equilibrium" in conditions:
            conditions["thermal_equilibrium"] = False
        # Add extensive parameters :
        for cond in [
            "path_length",
            "calculation_time",
            "total_lines",
            "lines_calculated",
            "lines_cutoff",
        ]:  # sum of all
            if cond in s.conditions and cond in sn.conditions:
                conditions[cond] = s.conditions[cond] + sn.conditions[cond]
                if cond in s.cond_units and cond in sn.cond_units:
                    assert s.cond_units == sn.cond_units
        cond_units = intersect(s.cond_units, sn.cond_units)

        # Update references
        # -----------------
        # (just add everything)
        references = {**s.references, **sn.references}

        # name
        name = _serial_slab_names(s, sn)

        return Spectrum(
            quantities=quantities,
            conditions=conditions,
            cond_units=cond_units,
            units=unitsn,
            name=name,
            references=references,
            warnings=False,  # we already know waveranges are properly spaced, etc.
        )


def _serial_slab_names(s, sn):
    # type: (Spectrum, Spectrum) -> Spectrum
    name_s = s.get_name()
    name_sn = sn.get_name()
    if "//" in name_s and not ">>" in name_s:
        name_s = "({0})".format(name_s)
    if "//" in name_sn:
        name_sn = "({0})".format(name_sn)
    return "{0}>>{1}".format(name_s, name_sn)


def _check_valid(s):
    # type: (Spectrum) -> bool
    """Check s is a valid Spectrum object. Raises an error if not.

    Valid if:

    - is a Spectrum

    Also print a warning if:

    - quantities used for solving the LOS have nan
    """

    if not isinstance(s, Spectrum):
        raise TypeError(
            "All inputs must be Spectrum objects (got: {0})".format(type(s))
        )
    for k, v in s._q.items():
        if (
            k in ["transmittance_noslit", "radiance_noslit", "abscoeff", "emisscoeff"]
        ) and anynan(v):
            warn(
                "Nans detected in Spectrum object for multi-slab operation. "
                + "Results may be wrong!"
            )

    return True


def _has_quantity(quantity, *slabs):
    # type: (str, *Spectrum) -> bool
    slabs = list(slabs)
    b = True
    for s in slabs:
        b *= quantity in slabs
    return b


def resample_slabs(
    waveunit, resample_wavespace, out_of_bounds="nan", modify_inputs=False, *slabs
):
    # type: (str, str, str, *Spectrum) -> *Spectrum
    """Resample slabs on the same wavespace: if the range are differents,
    depending on the mode we may fill with optically thin media, or raise an
    error

    Parameters
    ----------
    waveunit: ``'nm'``, ``'cm-1'``
        which wavespace we're working in
    resample_wavespace: 'never', 'intersect', 'full'
        what to do when spectra have different wavespaces:

        - If 'never', raises an error
        - If 'intersect', uses the intersection of all ranges, and resample
          spectra on the most resolved wavespace.
        - If 'full', uses the overlap of all ranges, resample spectra on the
          most resolved wavespace, and fill missing data with 0 emission and 0
          absorption

        Default 'never'
    out_of_bounds: 'transparent', 'nan', 'error'
        what to do if resampling is out of bounds:

        - 'transparent': fills with transparent medium.
        - 'nan': fills with nan.
        - 'error': raises an error.

        Default ``'nan'``
    *slabs: list of Spectrum objects

    Other Parameters
    ----------------
    modify_inputs: False
        if ``True``, slabs are modified directly when they are resampled. This
        avoids making a copy so is slightly faster. Default ``False``.

    Returns
    -------
    slabs: list of Spectrum objects
        resampled copies of inputs Spectra. All now have the same wavespace
    """

    def same_wavespace(wl):
        try:
            assert all([(w == wl[0]).all() for w in wl[1:]])
        except AssertionError:
            # ok if wavespace is the same within some tolerance
            #   @dev: testing with == first is 10x faster and works in most cases.
            if all([allclose(w, wl[0]) for w in wl[1:]]):
                return True
            return False
        except:  # ex:  different lengths
            return False
        else:
            return True

    # Work on copies
    if not modify_inputs:
        slabs = [s.copy(copy_lines=False) for s in slabs]

    # Get all keys
    keys = merge_lists([s.get_vars() for s in slabs])

    for k in keys:
        # We ensure all quantities are resampled. Note that .resample() deals
        # with all quantities so usually spectra are fully corrected after the
        # first iteration. However, it also deals with cases where quantities
        # have different wavespaces (ex: radiance and radiance_noslit are not
        # defined on the same range)
        slabsk = [s for s in slabs if k in s.get_vars()]  # note that these are
        # references to the actual Spectrum copy
        wl = [s.get(k, wunit=waveunit, copy=False)[0] for s in slabsk]
        if not same_wavespace(wl):
            # resample slabs if allowed
            if resample_wavespace == "never":
                raise ValueError(
                    "All wavelengths/wavenumbers must be the same for "
                    + "multi slabs configurations. Consider using "
                    + "`resample='intersect'` or "
                    + "`resample='full'`"
                )
            elif resample_wavespace == "full":
                # ... get bounds
                wmin = min([w.min() for w in wl])  # minimum of all
                wmax = max([w.max() for w in wl])  # maximum of all
                dw = min([abs(diff(w)).min() for w in wl])  # highest density
                wnew = arange(wmin, wmax + dw, dw)
                if wnew[-1] > wmax:  # sometimes arange doesnt work as expected
                    wnew = wnew[:-1]
                for s in slabsk:
                    s.resample(
                        wnew,
                        unit=waveunit,
                        out_of_bounds=out_of_bounds,
                        inplace=True,  # we already copied 'if not modify_inputs'
                    )
                # note: s.resample() fills with 0 when out of bounds
            elif resample_wavespace == "intersect":
                # ... get bounds
                wmin = max([w.min() for w in wl])  # maximum of all
                wmax = min([w.max() for w in wl])  # minimum of all
                dw = min([abs(diff(w)).min() for w in wl])  # highest density
                wnew = arange(wmin, wmax + dw, dw)
                if wnew[-1] > wmax:  # sometimes arange doesnt work as expected
                    wnew = wnew[:-1]
                if len(wnew) == 0:
                    raise ValueError("Intersect range is empty")
                for s in slabsk:
                    s.resample(
                        wnew,
                        unit=waveunit,
                        out_of_bounds=out_of_bounds,
                        inplace=True,  # we already copied 'if not modify_inputs'
                    )
                # note: s.resample() fills with 0 when out of bounds
        # Now all our slabs have the same wavespace

    return slabs


def MergeSlabs(*slabs, **kwargs) -> Spectrum:
    # type: (*Spectrum, **dict) -> Spectrum
    r"""Combines several slabs into one. Useful to calculate multi-gas slabs.
    Linear absorption coefficient is calculated as the sum of all linear absorption
    coefficients, and the RTE is recalculated to get the total radiance.
    You can also simply use::

        s1//s2
    Merged spectrum ``1+2`` is calculated with Eqn (4.3) of the [RADIS-2018]_ article,
    generalized to N slabs :

    .. math::

        j_{\lambda, 1+2} = j_{\lambda, 1} + j_{\lambda, 2}

        k_{\lambda, 1+2} = k_{\lambda, 1} + k_{\lambda, 2}

    where

    .. math:: j_{\lambda}, k_{\lambda}

    are the emission coefficient and absorption coefficient of the two slabs ``1`` and ``2``.
    Emission and absorption coefficients are calculated if not given in the
    initial slabs (if possible).

    Parameters
    ----------
    slabs: list of Spectra, each representing a slab
        ``path_length`` must be given in Spectrum conditions, and equal for all
        spectra.

        line-of-sight::

            slabs
                        [0]        \====
            light       [1]  ->     )===  observer
                        [n]        /====

    Other Parameters
    ----------------
    kwargs input:
    resample: ``'never'``, ``'intersect'``, ``'full'``
        what to do when spectra have different wavespaces:

        - If ``'never'``, raises an error
        - If ``'intersect'``, uses the intersection of all ranges, and resample
          spectra on the most resolved wavespace.
        - If ``'full'``, uses the overlap of all ranges, resample spectra on the
          most resolved wavespace, and fill missing data with 0 emission and 0
          absorption

        Default ``'never'``
    out: ``'transparent'``, ``'nan'``, ``'error'``
        what to do if resampling is out of bounds:

        - ``'transparent'``: fills with transparent medium.
        - ``'nan'``: fills with nan.
        - ``'error'``: raises an error.

        Default ``'nan'``
    optically_thin: boolean
        if ``True``, merge slabs in optically thin mode. Default ``False``
    verbose: boolean
        if ``True``, print messages and warnings. Default ``False``
    modify_inputs: False
        if ``True``, slabs are modified directly when they are resampled. This
        avoids making a copy so is slightly faster. Default ``False``.

    Returns
    -------
    Spectrum: object representing total emission and total transmittance as
        observed at the output. Conditions and units are transported too,
        unless there is a mismatch then conditions are dropped (and units mismatch
        raises an error because it doesnt make sense)

    Examples
    --------
    Merge two spectra calculated with different species (physically correct
    only if broadening coefficients dont change much)::

        from radis import calc_spectrum, MergeSlabs
        s1 = calc_spectrum(...)
        s2 = calc_spectrum(...)
        s3 = MergeSlabs(s1, s2)

    The last line is equivalent to::

        s3 = s1//s2

    Load a spectrum precalculated on several partial spectral ranges, for a same
    molecule (i.e, partial spectra are optically thin on the rest of the spectral
    range)::

        from radis import load_spec, MergeSlabs
        spectra = []
        for f in ['spec1.spec', 'spec2.spec', ...]:
            spectra.append(load_spec(f))
        s = MergeSlabs(*spectra, resample='full', out='transparent')
        s.update()   # Generate missing spectral quantities
        s.plot()

    .. minigallery:: radis.MergeSlabs

    See Also
    --------
    :func:`~radis.los.slabs.SerialSlabs`
    :ref:`See more examples in Line-of-Sight module <label_los_index>`

    """

    # Deprecation warnings
    if "resample_wavespace" in kwargs:
        warn(DeprecationWarning("'resample_wavespace' replaced with 'resample'"))
        kwargs["resample"] = kwargs.pop("resample_wavespace")
    if "out_of_bounds" in kwargs:
        warn(DeprecationWarning("'out_of_bounds' replaced with 'out'"))
        kwargs["out"] = kwargs.pop("out_of_bounds")

    # Check inputs, get defaults
    # inputs (Python 2 compatible)
    resample_wavespace = kwargs.pop("resample", "never")  # default 'never'
    out_of_bounds = kwargs.pop("out", "nan")  # default 'nan'
    optically_thin = kwargs.pop("optically_thin", False)  # default False
    verbose = kwargs.pop("verbose", False)  # type: bool
    kwargs.pop("debug", False)  # type: bool
    modify_inputs = kwargs.pop("modify_inputs", False)  # type: bool
    if len(kwargs) > 0:
        raise ValueError("Unexpected input: {0}".format(list(kwargs.keys())))

    # Check inputs
    if resample_wavespace not in ["never", "intersect", "full"]:
        raise ValueError(
            "'resample' should be one of: {0}".format(
                ", ".join(["never", "intersect", "full"])
            )
        )

    if len(slabs) == 0:
        raise ValueError("Empty list of slabs")

    elif len(slabs) == 1:
        if not isinstance(slabs[0], Spectrum):
            raise TypeError(
                "MergeSlabs takes an unfolded list of Spectrum as "
                + "argument: (got {0})".format(type(slabs[0]))
            )
        return slabs[0]

    else:  # calculate serial slabs

        slabs = list(slabs)

        #        # Check all items are valid Spectrum objects
        for s in slabs:
            _check_valid(s)

        # Check all path_lengths are defined and they exist
        try:
            path_lengths = [s.conditions["path_length"] for s in slabs]
        except KeyError:
            raise ValueError(
                "path_length must be defined for all slabs in MergeSlabs. "
                + "Set it with `s.conditions['path_length']=`. "
            )
        if not all([L == path_lengths[0] for L in path_lengths[1:]]):
            raise ValueError(
                "path_length must be equal for all MergeSlabs inputs"
                + "  (got {0})".format(path_lengths)
            )

        # make sure we use the same wavespace type (even if sn is in 'nm' and s in 'cm-1')
        waveunit = slabs[0].get_waveunit()
        # Make all our slabs copies with the same wavespace range
        # (note: wavespace range may be different for different quantities, but
        # equal for all slabs)
        slabs = resample_slabs(
            waveunit, resample_wavespace, out_of_bounds, modify_inputs, *slabs
        )
        w_noconv = slabs[0]._get_wavespace()

        # %% Update conditions of the Merged spectrum
        # -------------------------------------------
        conditions = slabs[0].conditions
        conditions["waveunit"] = waveunit
        cond_units = slabs[0].cond_units
        units0 = slabs[0].units
        # Define conditions as intersection of everything (N/A if unknown)
        # ... this will only keep intensive parameters (same for all)
        for s in slabs[1:]:
            conditions = intersect(conditions, s.conditions)
            cond_units = intersect(cond_units, s.cond_units)
            # units = intersect(units0, s.units)  # we're actually using [slabs0].units insteads
        # ... Add extensive parameters
        for cond in ["molecule"]:  # list of all
            if in_all(cond, [s.conditions for s in slabs]):
                conditions[cond] = set([s.conditions[cond] for s in slabs])
        for cond in [
            "thermal_equilibrium",
        ]:
            if in_all(cond, [s.conditions for s in slabs]):
                conditions[cond] = all([s.conditions[cond] for s in slabs])
        for cond in [
            "calculation_time",
            "total_lines",
            "lines_calculated",
            "lines_cutoff",
        ]:  # sum of all
            if in_all(cond, [s.conditions for s in slabs]):
                conditions[cond] = sum([s.conditions[cond] for s in slabs])
        # ... TODO @dev: create a large list/dictionary outside of SerialSlabS/MergeSlabs
        # ... with how to deal with all every condition (sum, list, intersect, etc.)
        # ... Example :
        # {"calculation_time":sum,
        #  "lines_calculated":sum,
        #  "lines_cutoff":sum,
        #  "lines_in_continuum":sum,
        #  "molecule":{'MergeSlabs':set},
        #  "path_length":{'SerialSlabs':sum}
        #  "isotope":{'MergeSlabs':dict},  # make a dict, same for mole fractions?
        #  "mole_fractions":{'MergeSlabs':dict},  # make a dict, same for mole fractions?
        #  }

        # Update references
        # -----------------
        # (just add everything)
        references = slabs[0].references
        for s in slabs[1:]:
            references.update(s.references)

        # %% Get quantities that should be calculated
        # Try to keep all the quantities of the initial slabs:
        requested = merge_lists([s.get_vars() for s in slabs])
        recompute = requested[:]  # copy
        if "radiance_noslit" in requested and not optically_thin:
            recompute.append("emisscoeff")
            recompute.append("abscoeff")
        if "abscoeff" in recompute and "path_length" in conditions:
            recompute.append("absorbance")
            recompute.append("transmittance_noslit")

        # To make it easier, we start from abscoeff and emisscoeff of all slabs
        # Let's recompute them all
        # TODO: if that changes the initial Spectra, maybe we should just work on copies
        for s in slabs:
            if "abscoeff" in recompute and not "abscoeff" in list(s._q.keys()):
                s.update("abscoeff", verbose=False)
                # that may crash if Spectrum doesnt have the correct inputs.
                # let update() handle that
            if "emisscoeff" in recompute and not "emisscoeff" in list(s._q.keys()):
                s.update("emisscoeff", verbose=False)
                # same

        # %% Calculate total emisscoeff and abscoeff
        added = {}

        # ... absorption coefficient (cm-1)
        if "abscoeff" in recompute:
            # TODO: deal with all cases
            if __debug__:
                printdbg("... merge: calculating abscoeff k=sum(k_i)")
            abscoeff_eq = np.sum(
                [
                    s.get("abscoeff", wunit=waveunit, Iunit=units0["abscoeff"])[1]
                    for s in slabs
                ],
                axis=0,
            )
            assert len(w_noconv) == len(abscoeff_eq)
            added["abscoeff"] = (w_noconv, abscoeff_eq)

        # ... emission coefficient
        if "emisscoeff" in recompute:
            if __debug__:
                printdbg("... merge: calculating emisscoeff j=sum(j_i)")
            emisscoeff_eq = np.sum(
                [
                    s.get("emisscoeff", wunit=waveunit, Iunit=units0["emisscoeff"])[1]
                    for s in slabs
                ],
                axis=0,
            )
            assert len(w_noconv) == len(emisscoeff_eq)
            added["emisscoeff"] = (w_noconv, emisscoeff_eq)

        # name
        name = "//".join([s.get_name() for s in slabs])

        # TODO: check units are consistent in all slabs inputs
        s = Spectrum(
            quantities=added,
            conditions=conditions,
            cond_units=cond_units,
            units=units0,
            name=name,
            references=references,
        )

        # %% Calculate all quantities from emisscoeff and abscoeff

        if "emissivity_noslit" in requested and (
            "thermal_equilibrium" not in s.conditions or s.is_at_equilibrium() != True
        ):
            requested.remove("emissivity_noslit")
            if __debug__:
                printdbg(
                    "... merge: all slabs are not proven to be at equilibrium. "
                    + "Emissivity was not calculated"
                )

        # Add the rest of the spectral quantities afterwards:
        s.update(
            [k for k in requested if k not in ["emisscoeff", "abscoeff"]],
            optically_thin=optically_thin,
            verbose=verbose,
        )

        return s


# %% Tests


if __name__ == "__main__":
    from radis.test.los.test_slabs import _run_testcases

    print("Testing merge slabs: ", _run_testcases(verbose=True))
