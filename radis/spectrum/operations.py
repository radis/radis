#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Handy functions for manipulation of :class:`~radis.spectrum.spectrum.Spectrum` objects

Routine Listing
---------------

Operators:

- :py:func:`~radis.spectrum.operations.multiply` : simply use ``s*2``
- :py:func:`~radis.spectrum.operations.add_constant` : simply use ``s+1``
- :py:func:`~radis.spectrum.operations.add_array` : simply use ``s+a``
- :py:func:`~radis.spectrum.operations.add_spectra` : simply use ``s1+s2``
- :py:func:`~radis.spectrum.operations.substract_spectra` : simply use ``s1-s2``

Note that these operators are purely algebraic and should not be used in place
of the line-of-sight functions, i.e, :py:func:`~radis.los.slabs.SerialSlabs` (``>``)
and :py:func:`~radis.los.slabs.MergeSlabs` (``//``)

Functions to manipulate one spectrum:

- :py:func:`~radis.spectrum.operations.crop`
- :py:func:`~radis.spectrum.operations.offset`

Play with baselines:

- :py:func:`~radis.spectrum.operations.get_baseline`
- :py:func:`~radis.spectrum.operations.sub_baseline`

Functions to discard all but one spectral quantity:

- :py:func:`~radis.spectrum.operations.Transmittance`
- :py:func:`~radis.spectrum.operations.Transmittance_noslit`
- :py:func:`~radis.spectrum.operations.Radiance`
- :py:func:`~radis.spectrum.operations.Radiance_noslit`

Keeps all spectral quantities, but make emission equal to 0 (useful when
calculating contributions of line of sight slabs):

- :py:func:`~radis.spectrum.operations.PerfectAbsorber`


Examples
--------

Most of these functions are implemented with the standard operators. Ex::

    ((s_exp - 0.1)*10).plot()   # works for a Spectrum s_exp

-------------------------------------------------------------------------------


"""


from warnings import warn

import astropy.units as u
from numpy import hstack, ones_like

from radis.phys.convert import (
    air2vacuum,
    cm2nm,
    cm2nm_air,
    dcm2dnm,
    dcm2dnm_air,
    dnm2dcm,
    dnm_air2dcm,
    nm2cm,
    nm_air2cm,
    vacuum2air,
)
from radis.phys.units import Unit
from radis.spectrum.spectrum import Spectrum

# %% Filter Spectra


def Transmittance(s: Spectrum) -> Spectrum:
    """Returns a new Spectrum with only the ``transmittance`` component of ``s``

    Parameters
    ----------
    s: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object

    Returns
    -------
    s_tr: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object, with only the ``transmittance``,
        ``absorbance`` and/or ``abscoeff`` part of ``s``, where ``radiance_noslit`` ,
        ``emisscoeff`` and ``emissivity_noslit`` (if they exist) have been set to 0

    Examples
    --------
    This function is useful to use :ref:`Spectrum algebra <label_spectrum_algebra>`
    operations::

        s = calc_spectrum(...)   # contains emission & absorption arrays
        tr = Transmittance(s)    # contains 'radiance_noslit' array only
        tr -= 0.1    # arithmetic operation is applied to Transmittance only

    Equivalent to::

        rad = s.take('transmittance')

    See Also
    --------
    :py:func:`~radis.spectrum.operations.Transmittance_noslit`,
    :py:func:`~radis.spectrum.operations.Radiance_noslit`,
    :py:func:`~radis.spectrum.operations.Radiance`,
    :py:meth:`~radis.spectrum.Spectrum.take`

    """

    return s.copy(copy_lines=True, quantity="transmittance")


def Transmittance_noslit(s: Spectrum) -> Spectrum:
    """Returns a new Spectrum with only the ``transmittance_noslit`` component of ``s``

    Parameters
    ----------
    s: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object

    Returns
    -------
    s_tr: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object, with only
        ``transmittance_noslit`` defined

    Examples
    --------
    This function is useful to use :ref:`Spectrum algebra <label_spectrum_algebra>`
    operations::

        s = calc_spectrum(...)   # contains emission & absorption arrays
        tr = Transmittance_noslit(s) # contains 'radiance_noslit' array only
        tr -= 0.1    # arithmetic operation is applied to Transmittance_noslit only

    Equivalent to::

        rad = s.take('transmittance_noslit')

    See Also
    --------
    :py:func:`~radis.spectrum.operations.Transmittance`
    :py:func:`~radis.spectrum.operations.Radiance_noslit`,
    :py:func:`~radis.spectrum.operations.Radiance`,
    :py:meth:`~radis.spectrum.Spectrum.take`

    """

    return s.copy(copy_lines=True, quantity="transmittance_noslit")


def Radiance(s: Spectrum) -> Spectrum:
    """Returns a new Spectrum with only the ``radiance`` component of ``s``

    Parameters
    ----------
    s: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object

    Returns
    -------
    s_tr: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object, with only ``radiance``
        defined

    Examples
    --------
    This function is useful to use :ref:`Spectrum algebra <label_spectrum_algebra>`
    operations::

        s = calc_spectrum(...)   # contains emission & absorption arrays
        rad = Radiance(s)        # contains radiance array only
        rad -= 0.1    # arithmetic operation is applied to Radiance only

    Equivalent to::

        rad = s.take('radiance')

    See Also
    --------
    :py:func:`~radis.spectrum.operations.Radiance_noslit`
    :py:func:`~radis.spectrum.operations.Transmittance_noslit`,
    :py:func:`~radis.spectrum.operations.Transmittance`,
    :py:meth:`~radis.spectrum.Spectrum.take`

    """

    return s.copy(copy_lines=True, quantity="radiance")


def Radiance_noslit(s: Spectrum) -> Spectrum:
    """Returns a new Spectrum with only the ``radiance_noslit`` component of ``s``

    Parameters
    ----------
    s: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object

    Returns
    -------
    s_tr: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object, with only ``radiance_noslit``
        defined

    Examples
    --------
    This function is useful to use :ref:`Spectrum algebra <label_spectrum_algebra>`
    operations::

        s = calc_spectrum(...)   # contains emission & absorption arrays
        rad = Radiance_noslit(s) # contains 'radiance_noslit' array only
        rad -= 0.1    # arithmetic operation is applied to Radiance_noslit only

    Equivalent to::

        rad = s.take('radiance_noslit')

    See Also
    --------
    :py:func:`~radis.spectrum.operations.Radiance`,
    :py:func:`~radis.spectrum.operations.Transmittance_noslit`,
    :py:func:`~radis.spectrum.operations.Transmittance`,
    :py:meth:`~radis.spectrum.Spectrum.take`

    """

    return s.copy(copy_lines=True, quantity="radiance_noslit")


# Useful for determining line-of-sight contributions:


def PerfectAbsorber(s: Spectrum) -> Spectrum:
    """Makes a new Spectrum with the same transmittance/absorbance as Spectrum
    ``s``, but with radiance set to 0.
    Useful to get contribution of different slabs in line-of-sight
    calculations (see example).

    .. note:

        formerly named "Transmittance", but "Transmittance(s)" wouldnt
        return the Transmittance exactly

    Parameters
    ----------

    s: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object

    Returns
    -------

    s_tr: Spectrum
        :class:`~radis.spectrum.spectrum.Spectrum` object, with only the ``transmittance``,
        ``absorbance`` and/or ``abscoeff`` part of ``s``, where ``radiance_noslit`` ,
        ``emisscoeff`` and ``emissivity_noslit`` (if they exist) have been set to 0

    Examples
    --------

    Let's say you have a total line of sight::

        s_los = s1 > s2 > s3

    If you now want to get the contribution of ``s2`` to the line-of-sight emission,
    you can do::

        (s2 > PerfectAbsorber(s3)).plot('radiance_noslit')

    And the contribution of ``s1`` would be::

        (s1 > PerfectAbsorber(s2>s3)).plot('radiance_noslit')

    See more examples in :ref:`Line-of-Sight module <label_los_index>`

    """

    s_tr = s.copy()
    for k in ["radiance_noslit", "emisscoeff", "emissivity_noslit"]:
        if k in s_tr._q:
            s_tr._q[k] *= 0

    for k in ["radiance", "emissivity"]:
        if k in s_tr._q:
            s_tr._q[k] *= 0

    # Deactivate equilibrium conditions (so Kirchoff's law cannot be used anymore)
    s_tr.conditions["thermal_equilibrium"] = False

    s_tr.name = "PureAbsorber({0})".format(s.get_name())

    return s_tr


# %% Change wavelength


def crop(s: Spectrum, wmin=None, wmax=None, wunit=None, inplace=False) -> Spectrum:
    # type: (Spectrum, float, float, str, str, bool) -> Spectrum
    """Crop spectrum to ``wmin-wmax`` range in ``wunit``

    Parameters
    ----------

    s: Spectrum object
        object to crop

    wmin, wmax: float, or None
        boundaries of spectral range (in ``wunit``)

        wunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
            which waveunit to use for ``wmin, wmax``. If ``default``:
            use the default Spectrum wavespace defined with
            :meth:`~radis.spectrum.spectrum.Spectrum.get_waveunit`.

    Other Parameters
    ----------------

    inplace: bool
        if ``True``, modifiy ``s`` directly. Else, returns a copy.

    Returns
    -------

    s_crop: Spectrum
        a cropped Spectrum.
        if using ``inplace``, then ``s_crop`` and ``s`` are still the same object

    Examples
    --------

    ::

        crop(s, 420, 480, 'nm', 'air')

    Or in ``cm-1``::

        crop(s, 2000, 2300, 'cm-1')

    """

    # Check inputs
    if wmin is None and wmax is None:
        raise ValueError("Choose at least `wmin=` or `wmax=`")
    if wunit is None:
        raise ValueError("Please precise unit for wmin and wmax with `unit=`")
    assert wunit in ["nm", "cm-1"]
    if (wmin is not None and wmax is not None) and wmin >= wmax:
        raise ValueError(
            "wmin should be < wmax (Got: {0:.2f}, {1:.2f})".format(wmin, wmax)
        )

    if not inplace:
        s = s.copy()

    # Convert wmin, wmax to Spectrum wavespace  (stored_waveunit)
    # (deal with cases where wavelength are given in 'air' or 'vacuum')
    # TODO @dev: rewrite with wunit='cm-1', 'nm_air', 'nm_vac'
    stored_waveunit = s.get_waveunit()
    wmin0, wmax0 = wmin, wmax
    if stored_waveunit == "cm-1":
        # convert wmin, wmax to wavenumber
        if wunit == "nm":
            if wmax0:
                wmin = nm_air2cm(wmax0)  # note: min/max inverted
            if wmin0:
                wmax = nm_air2cm(wmin0)  # note: min/max inverted
        elif wunit == "nm_vac":
            if wmax0:
                wmin = nm2cm(wmax0)  # note: min/max inverted
            if wmin0:
                wmax = nm2cm(wmin0)  # note: min/max inverted
        elif wunit == "cm-1":
            pass
        else:
            raise ValueError(wunit)
    elif stored_waveunit == "nm":
        # convert wmin, wmax to wavelength air
        if wunit == "nm":
            pass
        elif wunit == "nm_vac":
            if wmin0:
                wmin = vacuum2air(wmin0)
            if wmax0:
                wmax = vacuum2air(wmax0)
        elif wunit == "cm-1":
            if wmax0:
                wmin = cm2nm_air(wmax0)  # note: min/max inverted
            if wmin0:
                wmax = cm2nm_air(wmin0)  # note: min/max inverted
        else:
            raise ValueError(wunit)
    elif stored_waveunit == "nm_vac":
        # convert wmin, wmax to wavelength vacuum
        if wunit == "nm":
            if wmin0:
                wmin = air2vacuum(wmin0)
            if wmax0:
                wmax = air2vacuum(wmax0)
        elif wunit == "nm_vac":
            pass
        elif wunit == "cm-1":
            if wmax0:
                wmin = cm2nm(wmax0)  # note: min/max inverted
            if wmin0:
                wmax = cm2nm(wmin0)  # note: min/max inverted
        else:
            raise ValueError(wunit)
    else:
        raise ValueError(stored_waveunit)

    # Crop
    if len(s._q) > 0:
        b = ones_like(s._q["wavespace"], dtype=bool)
        if wmin:
            b *= wmin <= s._q["wavespace"]
        if wmax:
            b *= s._q["wavespace"] <= wmax
        for k, v in s._q.items():
            s._q[k] = v[b]

    return s


# %% Algebric operations on Spectra


def _get_unique_var(s, var, inplace):
    """Returns the unique spectral quantity in the Spectrum ``s``. If there are more,
    raises an error.
    """
    # If var is undefined, get it if there is no ambiguity
    if var is None:
        var = s._get_unique_var()
    # If inplace, assert there is only one var
    if inplace and len(s.get_vars()) > 1:
        raise ValueError(
            "Cant modify inplace one spectral quantity of a Spectrum "
            + "that contains several ({0}). Use 'inplace=False'".format(s.get_vars())
        )
    return var


def multiply(s, coef, unit=None, var=None, inplace=False):
    """Multiply s[var] by the float 'coef'

    Parameters
    ----------
    s: Spectrum object
        The spectrum to multiply.
    coef: float
        Coefficient of the multiplication.
    unit: str, or `~astropy.units.core.Unit`
        unit for ``coef``. If ``None``, ``coef`` is considered to be
        adimensioned. Else, the spectrum `~radis.spectrum.spectrum.Spectrum.units` is multiplied.
    var: str, or ``None``
        'radiance', 'transmittance', ... If ``None``, get the unique spectral
        quantity of ``s`` or raises an error if there is any ambiguity
    inplace: bool
        if ``True``, modifies ``s`` directly. Else, returns a copy.
        Default ``False``

    Returns
    -------
    s : Spectrum
        Spectrum object where intensity of s['var'] is multiplied by coef
        If ``inplace=True``, ``s`` has been modified directly.

    """
    # Check input
    var = _get_unique_var(s, var, inplace)

    # Case where a is dimensioned
    if isinstance(coef, u.quantity.Quantity):
        if unit is not None:
            raise ValueError(
                "Cannot use unit= when giving a dimensioned array ({0})".format(
                    coef.unit
                )
            )
        unit = coef.unit
        coef = coef.value

    if not inplace:
        s = s.copy(quantity=var)

    # Multiply inplace       ( @dev: we have copied already if needed )
    w, I = s.get(var, wunit=s.get_waveunit(), copy=False)
    I *= coef  # @dev: updates the Spectrum directly because of copy=False

    # Convert Spectrum unit
    if unit is not None:
        Iunit = s.units[var]
        s.units[var] = (Unit(Iunit) * Unit(unit)).to_string()

    return s


def add_constant(s, cst, unit=None, var=None, inplace=False):
    """Return a new spectrum with a constant added to s[var].
    Equivalent to::

        s + constant

    Parameters
    ----------
    s: Spectrum objects
        Spectrum you want to modify
    cst: float
        Constant to add.
    unit: str
        unit for ``cst``. If ``None``, uses the default unit in ``s`` for
        variable ``var``.
    var: str, or ``None``
        'radiance', 'transmittance', ... If ``None``, get the unique spectral
        quantity of ``s`` or raises an error if there is any ambiguity
    inplace: bool
        if ``True``, modifies ``s`` directly. Else, returns a copy.
        Default ``False``

    Returns
    -------
    s : Spectrum
        Spectrum object where cst is added to intensity of s['var']
        If ``inplace=True``, ``s`` has been modified directly.

    Notes
    -----

    Use only for rough work. If you want to work properly with spectrum
    objects, see :py:meth:`~radis.los.slabs.MergeSlabs`.
    """
    # Check input
    var = _get_unique_var(s, var, inplace)

    # Note: using dimensioned value for `cst` will make it an array,
    # processed by add_array

    # Convert to Spectrum unit
    if unit is not None:
        Iunit = s.units[var]
        if unit != Iunit:
            from radis.phys.convert import conv2

            cst = conv2(cst, unit, Iunit)

    if not inplace:
        s = s.copy(quantity=var)

    # Add inplace       ( @dev: we have copied already if needed )
    w, I = s.get(var, wunit=s.get_waveunit(), copy=False)
    I += cst
    # @dev: updates the Spectrum directly because of copy=False

    #    s.name = '{0}+{1}'.format(s.get_name(), cst)

    return s


def add_array(s, a, unit=None, var=None, inplace=False):
    """Return a new spectrum with a constant added to s[var].
    Equivalent to::

        s + array

    Parameters
    ----------
    s: Spectrum objects
        Spectrum you want to modify
    a: numpy array, or `~astropy.units.quantity.Quantity`
        array to add. Must have the same length as variable ``var`` in Spectrum
        ``s``. Can be dimensioned with :py:mod:`~astropy.units`.
    unit: str
        unit for ``a``. If ``None``, uses the default unit in ``s`` for
        variable ``var``.
    var: str, or ``None``
        'radiance', 'transmittance', ... If ``None``, get the unique spectral
        quantity of ``s`` or raises an error if there is any ambiguity
    inplace: bool
        if ``True``, modifies ``s`` directly. Else, returns a copy.
        Default ``False``

    Returns
    -------
    s : Spectrum
        Spectrum object where array ``a`` is added to intensity of s['var']
        If ``inplace=True``, ``s`` has been modified directly.

    Notes
    -----
    Use only for rough work. If you want to work properly with spectrum
    objects, see MergeSlabs.

    Examples
    --------

    Add Gaussian noise to your Spectrum (assuming there is only one spectral
    quantity defined)::

        s += np.random.normal(0,1,len(s.get_wavelength()))

    """
    # Check input
    var = _get_unique_var(s, var, inplace)

    # Case where a is dimensioned
    if isinstance(a, u.quantity.Quantity):
        if unit is not None:
            raise ValueError(
                "Cannot use unit= when giving a dimensioned array ({0})".format(a.unit)
            )
        unit = a.unit
        a = a.value

    # Convert to Spectrum unit
    if unit is not None:
        Iunit = s.units[var]
        if unit != Iunit:
            from radis.phys.convert import conv2

            a = conv2(a, unit, Iunit)

    if not inplace:
        s = s.copy(quantity=var)

    # Add inplace       ( @dev: we have copied already if needed )
    w, I = s.get(var, wunit=s.get_waveunit(), copy=False)
    I += a
    # @dev: updates the Spectrum directly because of copy=False

    return s


def sub_baseline(s, left, right, unit=None, var=None, inplace=False):
    """Return a new spectrum with a baseline substracted to s[var]

    Parameters
    ----------

    s: Spectrum objects
        Spectrum you want to modify
    left: Float
        Constant to substract on the left of the spectrum.
    right: Float
        Constant to substract on the right of the spectrum.
    unit: str
        unit for ``cst``. If ``None``, uses the default unit in ``s`` for
        variable ``var``.
    var: str
        'radiance', 'transmittance', ...  If ``None``, get the unique spectral
        quantity of ``s`` or raises an error if there is any ambiguity
    inplace: bool
        if ``True``, modifies ``s`` directly. Else, returns a copy.
        Default ``False``

    Returns
    -------

    s: Spectrum
        Spectrum object where the baseline was substracted to intensity of s['var']
        If ``inplace=True``, ``s`` has been modified directly.

    Notes
    -----
    Use only for rough work.

    See Also
    --------

    :py:func:`~radis.spectrum.operations.get_baseline`
    """
    import numpy as np

    # Check input
    var = _get_unique_var(s, var, inplace)

    # Case where left, right are dimensioned
    if isinstance(left, u.quantity.Quantity) or isinstance(right, u.quantity.Quantity):
        if unit is not None:
            raise ValueError(
                "Cannot use unit= when giving a dimensioned array ({0})".format(
                    left.unit
                )
            )
        if not isinstance(left, u.quantity.Quantity) or not isinstance(
            right, u.quantity.Quantity
        ):
            raise ValueError(
                "Both left and right arguments must be dimensioned with the same unit"
            )
        if left.unit != right.unit:
            raise ValueError(
                "Both left and right arguments must be dimensioned with the same unit"
            )
        unit = left.unit
        left = left.value
        right = right.value

    # Convert to Spectrum unit
    if unit is not None:
        Iunit = s.units[var]
        if unit != Iunit:
            from radis.phys.convert import conv2

            left = conv2(left, unit, Iunit)
            right = conv2(right, unit, Iunit)

    if not inplace:
        s = s.copy(quantity=var)

    # @EP:

    # Substract inplace       ( @dev: we have copied already if needed )
    w, I = s.get(var, wunit=s.get_waveunit(), copy=False)
    I -= np.linspace(left, right, num=np.size(I))
    # @dev: updates the Spectrum directly because of copy=False

    return s


def add_spectra(s1, s2, var=None, force=False):
    """Return a new spectrum with ``s2`` added to ``s1``.
    Equivalent to::

        s1 + s2

    .. warning::
        we are just algebrically adding the quantities. If you want to merge
        spectra while preserving the radiative transfer equation, see
        :func:`~radis.los.slabs.MergeSlabs` and :func:`~radis.los.slabs.SerialSlabs`

    Parameters
    ----------

    s1, s2: Spectrum objects
        Spectrum you want to substract
    var: str
        quantity to manipulate: 'radiance', 'transmittance', ... If ``None``,
        get the unique spectral quantity of ``s1``, or the unique spectral
        quantity of ``s2``, or raises an error if there is any ambiguity

    Returns
    -------

    s: Spectrum
        Spectrum object with the same units and waveunits as ``s1``

    See Also
    --------

    :func:`~radis.los.slabs.MergeSlabs`,
    :func:`~radis.spectrum.operations.substract_spectra`

    """

    # Get variable
    if var is None:
        try:
            var = _get_unique_var(
                s2, var, inplace=False
            )  # unique variable of 2nd spectrum
        except KeyError:
            var = _get_unique_var(
                s1, var, inplace=False
            )  # if doesnt exist, unique variable of 1st spectrum
            # if it fails, let it fail
    # Make sure it is in both Spectra
    if var not in s1.get_vars():
        raise KeyError("Variable {0} not in Spectrum {1}".format(var, s1.get_name()))
    if var not in s2.get_vars():
        raise KeyError("Variable {0} not in Spectrum {1}".format(var, s1.get_name()))

    if var in ["transmittance_noslit", "transmittance"] and not force:
        raise ValueError(
            "It does not make much physical sense to sum transmittances. Are "
            + "you sure of what you are doing? See also `//` (MergeSlabs), `>` "
            + "(SerialSlabs) and `concat_spectra`. If you're sure, use `force=True`"
        )

    # Get s1 units
    Iunit1 = s1.units[var]
    wunit1 = s1.get_waveunit()

    # Resample s2 on s1
    s2 = s2.resample(s1, inplace=False)

    # Add, change output unit if needed.
    w1, I1 = s1.get(var=var, Iunit=Iunit1, wunit=wunit1)
    w2, I2 = s2.get(var=var, Iunit=Iunit1, wunit=wunit1)

    name = s1.get_name() + "+" + s2.get_name()

    sub = Spectrum.from_array(w1, I1 + I2, var, wunit=wunit1, unit=Iunit1, name=name)
    #    warn("Conditions of the left spectrum were copied in the substraction.", Warning)
    return sub


def substract_spectra(s1, s2, var=None):
    """Return a new spectrum with ``s2`` substracted from ``s1``.
    Equivalent to::

        s1 - s2

    Parameters
    ----------

    s1, s2: Spectrum objects
        Spectrum you want to substract
    var: str
        quantity to manipulate: 'radiance', 'transmittance', ... If ``None``,
        get the unique spectral quantity of ``s1``, or the unique spectral
        quantity of ``s2``, or raises an error if there is any ambiguity

    Returns
    -------

    s: Spectrum
        Spectrum object with the same units and waveunits as ``s1``

    See Also
    --------

    :func:`~radis.spectrum.operations.add_spectra`

    """

    # Get variable
    if var is None:
        try:
            var = _get_unique_var(
                s2, var, inplace=False
            )  # unique variable of 2nd spectrum
        except KeyError:
            var = _get_unique_var(
                s1, var, inplace=False
            )  # if doesnt exist, unique variable of 1st spectrum
            # if it fails, let it fail
    # Make sure it is in both Spectra
    if var not in s1.get_vars():
        raise KeyError("Variable {0} not in Spectrum {1}".format(var, s1.get_name()))
    if var not in s2.get_vars():
        raise KeyError("Variable {0} not in Spectrum {1}".format(var, s1.get_name()))

    # Use same units
    Iunit1 = s1.units[var]
    wunit1 = s1.get_waveunit()

    # Resample s2 on s1
    s2 = s2.resample(s1, inplace=False)

    # Substract
    w1, I1 = s1.get(var=var, Iunit=Iunit1, wunit=wunit1)
    w2, I2 = s2.get(var=var, Iunit=Iunit1, wunit=wunit1)

    name = s1.get_name() + "-" + s2.get_name()

    sub = Spectrum.from_array(w1, I1 - I2, var, wunit=wunit1, unit=Iunit1, name=name)
    #    warn("Conditions of the left spectrum were copied in the substraction.", Warning)
    return sub


def concat_spectra(s1, s2, var=None):
    """Concatenate two spectra ``s1`` and ``s2`` side by side.

    Note: their spectral range should not overlap

    Returns
    -------

    s: Spectrum
        Spectrum object with the same units and waveunits as ``s1``

    Parameters
    ----------

    s1, s2: Spectrum objects
        Spectrum you want to concatenate
    var: str
        quantity to manipulate: 'radiance', 'transmittance', ... If ``None``,
        get the unique spectral quantity of ``s1``, or the unique spectral
        quantity of ``s2``, or raises an error if there is any ambiguity

    Notes
    -----

    .. warning::

        the output Spectrum has the sum of the spectral ranges of s1 and s2.
        It won't be evenly spaced. This means that you cannot apply a slit without
        side effects. Typically, you want to use this function for convolved
        quantities only, such as experimental spectra. Else, use
        :func:`~radis.los.slabs.MergeSlabs` with the options
        ``resample='full', out='transparent'``

    See Also
    --------

    :func:`~radis.spectrum.operations.add_spectra`,
    :func:`~radis.los.slabs.MergeSlabs`

    """

    # Get variable
    if var is None:
        try:
            var = _get_unique_var(
                s2, var, inplace=False
            )  # unique variable of 2nd spectrum
        except KeyError:
            var = _get_unique_var(
                s1, var, inplace=False
            )  # if doesnt exist, unique variable of 1st spectrum
            # if it fails, let it fail
    # Make sure it is in both Spectra
    if var not in s1.get_vars():
        raise KeyError("Variable {0} not in Spectrum {1}".format(var, s1.get_name()))
    if var not in s2.get_vars():
        raise KeyError("Variable {0} not in Spectrum {1}".format(var, s1.get_name()))

    if var in ["transmittance_noslit", "transmittance"]:
        warn(
            "It does not make much physical sense to sum transmittances. Are "
            + "you sure of what you are doing? See also // (MergeSlabs) and > "
            + "(SerialSlabs)"
        )

    # Use same units
    Iunit1 = s1.units[var]
    wunit1 = s1.get_waveunit()

    # Get the value, on the same wunit)
    w1, I1 = s1.get(var=var, copy=False)  # @dev: faster to just get the stored value.
    # it's copied in hstack() below anyway).
    w2, I2 = s2.get(var=var, Iunit=Iunit1, wunit=wunit1)

    if not (w1.max() < w2.min() or w2.max() > w1.min()):
        raise ValueError(
            "You cannot use concat_spectra for overlapping spectral ranges. "
            + "Got: {0:.2f}-{1:.2f} and {2:.2f}-{3:.2f} {4}. ".format(
                w1.min(), w1.max(), w2.min(), w2.max(), wunit1
            )
            + "Use MergeSlabs instead, with the correct `out=` parameter "
            + "for your case"
        )

    w_tot = hstack((w1, w2))
    I_tot = hstack((I1, I2))

    name = s1.get_name() + "&" + s2.get_name()  # use "&" instead of "+"

    concat = Spectrum.from_array(
        w_tot, I_tot, var, wunit=wunit1, unit=Iunit1, name=name
    )

    return concat


def offset(s, offset, unit, name=None, inplace=False):
    # type: (Spectrum, float, str, str, bool) -> Spectrum
    """Offset the spectrum by a wavelength or wavenumber

    Parameters
    ----------
    s: Spectrum
        Spectrum you want to modify
    offset: float
        Constant to add to all quantities in the Spectrum.
    unit: 'nm' or 'cm-1'
        unit for ``offset``.

    Other Parameters
    ----------------
    name: str
        name of output spectrum
    inplace: bool
        if ``True``, modifies ``s`` directly. Else, returns a copy.
        Default ``False``

    Returns
    -------
    s : Spectrum
        Spectrum object where cst is added to intensity of s['var']
        If ``inplace=True``, ``s`` has been modified directly.

    See Also
    --------
    call as a Spectrum method directly: :py:meth:`~radis.spectrum.spectrum.Spectrum.offset`
    """

    stored_waveunit = s.get_waveunit()

    # Convert offset to correct unit:
    if stored_waveunit == "cm-1":
        if unit == "nm":
            # Note @EP: here we're offsetting by a constant value in 'nm', which is
            # not a constant value in 'cm-1'. The offset is an array
            offset_q = -dnm_air2dcm(offset, s.get_wavelength())  # this is an array
        elif unit == "nm_vac":
            offset_q = -dnm2dcm(offset, s.get_wavelength())  # this is an array
        elif unit == "cm-1":
            offset_q = offset
        else:
            raise ValueError(unit)
    elif stored_waveunit == "nm":  # wavelength air
        if unit == "nm":
            offset_q = offset
        elif unit == "nm_vac":
            # Note @EP: strictly speaking, the offset should change a little bit
            # as nm_vac > nm depends on the wavelength. Neglected here. # TODO ?
            offset_q = offset
        elif unit == "cm-1":
            offset_q = -dcm2dnm_air(offset, s.get_wavenumber())  # this is an array
        else:
            raise ValueError(unit)
    elif stored_waveunit == "nm_vac":  # wavelength vacuum
        if unit == "nm":
            # Note @EP: strictly speaking, the offset should change a little bit
            # as nm > nm_vac depends on the wavelength. Neglected here. # TODO ?
            offset_q = offset
        elif unit == "nm_vac":
            offset_q = offset
        elif unit == "cm-1":
            offset_q = -dcm2dnm(offset, s.get_wavenumber())  # this is an array
        else:
            raise ValueError(unit)
    else:
        raise ValueError(stored_waveunit)

    if not inplace:
        s = s.copy()

    # Update all variables
    s._q["wavespace"] += offset_q
    # @dev: updates the Spectrum directly because of copy=False

    if name:
        s.name = name

    return s


def get_baseline(
    s, var="radiance", algorithm="als", wunit="default", Iunit="default", **kwargs
):
    """Calculate and returns a baseline

    Parameters
    ----------
    s: Spectrum
        Spectrum which needs a baseline
    var: str
        on which spectral quantity to read the baseline. Default ``'radiance'``.
        See :py:data:`~radis.spectrum.utils.SPECTRAL_QUANTITIES`
    algorithm: 'als', 'polynomial'
        Asymmetric least square or Polynomial fit
    **kwargs: dict
        additional parameters to send to the algorithm. By default,
        for 'polynomial':

            **kwargs = **{"deg": 1, "max_it":500}

        for 'als':

            **kwargs = {"asymmetry_param": 0.05,
                        "smoothness_param": 1e6}

    Returns
    -------
    baseline: Spectrum
        Spectrum object where intenisity is the baseline of s is computed by peakutils

    Examples
    --------

    .. minigallery:: radis.spectrum.spectrum.Spectrum.get_baseline

    See Also
    --------

    :py:func:`~radis.spectrum.operations.sub_baseline`,
    :py:meth:`~radis.spectrum.spectrum.Spectrum.get_baseline`
    """
    if wunit == "default":
        wunit = s.get_waveunit()
    if Iunit == "default":
        Iunit = s.units[var]

    w1, I1 = s.get(var=var, wunit=wunit, Iunit=Iunit)

    if algorithm == "polynomial":
        import peakutils

        polyargs = {"deg": 1, "max_it": 500}
        polyargs.update(kwargs)
        baseline = peakutils.baseline(I1, **polyargs)
    elif algorithm == "als":
        from radis.misc.signal import als_baseline

        alsargs = {"asymmetry_param": 0.05, "smoothness_param": 1e6}
        alsargs.update(kwargs)
        baseline = als_baseline(I1, **alsargs)

    baselineSpectrum = Spectrum.from_array(
        w1, baseline, var, wunit=wunit, unit=Iunit, name=s.get_name() + "_baseline"
    )
    return baselineSpectrum


# %% Tests

if __name__ == "__main__":
    from radis.test.spectrum.test_operations import _run_testcases

    _run_testcases(plot=True)
