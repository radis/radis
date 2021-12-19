# -*- coding: utf-8 -*-
"""

-------------------------------------------------------------------------------


"""

import warnings


def Unit(st, *args, **kwargs):
    """Radis evaluation of an unit, using :py:class:`~astropy.units.Unit`

    Changes compare to Astropy standards:

    - "µm" is accepted and converted to "um"
    - "/molecule" or "/molec" are removed, but do not raise an error
    - we do not raise a warning if multiple slashes, e.g. "mW/cm2/sr/nm"

    Examples
    --------

    ::

        from radis.phys.units import Unit as u
        a = 200 * u("mW/cm2/sr/nm")
        a += 0.1 * u("W/cm2/sr/nm")
    """

    import astropy.units as u

    try:
        st = st.replace("µ", "u")
        st = st.replace("/molecule", "")
        st = st.replace("/molec", "")
    except AttributeError:
        pass

    with warnings.catch_warnings():
        # Ignore some warnings : Cause all warnings to always be triggered.
        warnings.filterwarnings(
            "ignore", ".*multiple slashes.*", category=u.UnitsWarning
        )
        unit = u.Unit(st, *args, **kwargs)
    return unit


def conv2(quantity, fromunit, tounit):
    """Converts `quantity` from unit `fromunit` to unit `tounit`

    Parameters
    ----------
    quantity: array
        quantity to convert
    fromunit: str
        input unit
    tounit: str
        output unit

    Note
    ----

    1.

    The output is still non dimensional. We don't transform `quantity`
    into a dimensioned array because this may create a performance
    drop in computationaly-expensive task. Instead, we assume we know for
    sure the units in which some of our quantities will be created, and just
    want to let the users choose another output unit
    """

    import astropy.units as u

    try:
        a = quantity * Unit(fromunit)
        a = a.to(Unit(tounit))

    except u.UnitConversionError as err:
        suggestion = ""
        if fromunit == "" or tounit == "":
            suggestion += " Is one of the two arrays normalized ?"
        raise TypeError(
            f"Cannot convert quantity `{fromunit}` to the specified unit `{tounit}`. Please check the dimensions."
            + suggestion
        ) from err

    return a.value


def is_homogeneous(unit1, unit2):
    """Tells if unit1 and unit2 are homogeneous, using the Astropy library.

    Parameters
    ----------
    unit1, unit2: str
        units
    """

    import astropy.units as u

    try:
        1 * Unit(unit1) + 1 * Unit(unit2)
        return True
    except u.UnitConversionError:
        return False


# %% Emission density conversion functions (with changes from ~1/nm in ~1/cm-1)


def convert_emi2cm(j_nm, wavenum, Iunit0, Iunit):
    """Convert spectral emission density in wavelength base (typically
    ~mW/cm3/sr/nm) to spectral emission density in wavenumber base
    (~mW/cm3/sr/cm-1)

    Parameters
    ----------
    j_nm: array
        spectral emission density in Iunit0 unit (~ per wavelength) )
    wavenum: array (cm-1)
        wavenumber
    Iunit0: str
        unit (~ per wavelength) to convert from
    Iunit: str
        unit (~ per wavenumber) to convert to


    Notes
    -----
    Implementation:

    We use the variable substitution::

            dλ/dν = - 1/ν**2        [in SI]
            dλ/dν = - 10^7 / ν**2   [nm -> cm-1]

    We want   Lν.dν = Lλ.dλ
    so        Lν = - Lλ * 1e7 / ν^2


    Example
    -------
    Validation::

        w_nm, j_nm = s.get('emisscoeff', 'nm', 'mW/cm2/sr/nm')
        w_cm, j_cm = s.get('emisscoeff', 'cm', 'mW/cm2/sr/cm-1')
        print(trapz(y_nm, x_nm))
        print(trapz(y_cm, x_cm))

    Both integrals are to be the same
    """

    j_nm = conv2(j_nm, Iunit0, "mW/cm**3/sr/nm")

    # Convert ''mW/cm3/sr/nm' to 'mW/cm3/sr/cm-1'
    j_cm = j_nm * 1e7 / wavenum ** 2
    # note that we discard the - sign.
    # This means that one of the integrals `trapz(j_nm, wavelen)` or
    # `trapz(j_cm, wavenum)` will be < 0 as they are sorted in opposite order.
    # But both plots look alright

    # Convert to whatever was wanted
    j_cm = conv2(j_cm, "mW/cm3/sr/cm-1", Iunit)

    return j_cm


def convert_emi2nm(j_cm, wavenum, Iunit0, Iunit):
    """Convert spectral emission density in wavenumber base (typically
    ~mW/cm3/sr/cm-1) to spectral radiance in wavelength base (~mW/cm3/sr/nm)

    Parameters
    ----------
    j_cm: array
        spectral emission density in Iunit0 unit (~ per wavenumber)
    wavenum: array (cm-1)
        wavenumber
    Iunit0: str
        unit (~ per wavenumber) to convert from
    Iunit: str
        unit (~ per wavelength) to convert to


    Notes
    -----
    Implementation:

    We use the variable substitution::

            dλ/dν = - 1/ν**2        [in SI]
            dλ/dν = - 10^7 / ν**2   [nm -> cm-1]

    We want   Lν.dν = Lλ.dλ
    so        Lλ  = - Lν * 1e-7 * ν^2
    """

    j_cm = conv2(j_cm, Iunit0, "mW/cm3/sr/cm-1")

    # Convert 'mW/cm3/sr/cm-1' to 'mW/cm3/sr/nm'
    j_nm = j_cm * 1e-7 * wavenum ** 2
    # note that we discard the - sign.
    # This means that one of the integrals `trapz(L_nm, wavelen)` or
    # `trapz(L_cm, wavenum)` will be < 0 as they are sorted in opposite order.
    # But both plots look alright

    # Convert to whatever was wanted
    j_nm = conv2(j_nm, "mW/cm3/sr/nm", Iunit)

    return j_nm


# %% Radiance conversion functions (with changes from ~1/nm in ~1/cm-1)


def convert_rad2cm(l_nm, wavenum, Iunit0, Iunit):
    """Convert spectral radiance in wavelength base (~1/nm) to spectral
    radiance in wavenumber base (~1/cm-1)

    Parameters
    ----------
    l_nm: array
        spectral radiance in yunit0 unit (~ per wavelength)
    wavenum: array (cm-1)
        wavenumber
    Iunit0: str
        unit (~ per wavelength) to convert from
    Iunit: str
        unit (~ per wavenumber) to convert to


    Notes
    -----
    Implementation:

    We use the variable substitution:
            dλ/dν = - 1/ν**2        [in SI]
            dλ/dν = - 10^7 / ν**2   [nm -> cm-1]

    We want   Lν.dν = Lλ.dλ
    so        Lν = - Lλ * 1e7 / ν^2


    Example
    -------
    Validation::

        x_nm, y_nm = s.get('radiance_noslit', 'nm', 'mW/cm2/sr/nm')
        x_cm, y_cm = s.get('radiance_noslit', 'cm', 'mW/cm2/sr/cm-1')
        print(trapz(y_nm, x_nm))
        print(trapz(y_cm, x_cm))

    Both integrals are to be the same
    """

    l_nm = conv2(l_nm, Iunit0, "mW/cm2/sr/nm")

    # Convert ''mW/cm2/sr/nm' to 'mW/cm2/sr/cm-1'
    l_cm = l_nm * 1e7 / wavenum ** 2
    # note that we discard the - sign.
    # This means that one of the integrals `trapz(L_nm, wavelen)` or
    # `trapz(L_cm, wavenum)` will be < 0 as they are sorted in opposite order.
    # But both plots look alright

    # Convert to whatever was wanted
    l_cm = conv2(l_cm, "mW/cm2/sr/cm-1", Iunit)

    return l_cm


def convert_rad2nm(l_cm, wavenum, Iunit0, Iunit):
    """Convert spectral radiance in wavenumber base (~1/cm-1) to spectral
    radiance in wavelength base (~1/nm)

    Parameters
    ----------
    l_cm: array
        spectral radiance in yunit0 unit (~ per wavenumber)
    wavenum: array (cm-1)
        wavenumber
    Iunit0: str
        unit (~ per wavenumber) to convert from
    Iunit: str
        unit (~ per wavelength) to convert to

    Notes
    -----
    Implementation:

    We use a variable substitution::

            dλ/dν = - 1/ν**2        [in SI]
            dλ/dν = - 10^7 / ν**2   [nm -> cm-1]

    We want   Lν.dν = Lλ.dλ
    so        Lλ  = - Lν * 1e-7 * ν^2
    """

    l_cm = conv2(l_cm, Iunit0, "mW/cm2/sr/cm-1")

    # Convert 'mW/cm2/sr/cm-1' to 'mW/cm2/sr/nm'
    l_nm = l_cm * 1e-7 * wavenum ** 2
    # note that we discard the - sign.
    # This means that one of the integrals `trapz(L_nm, wavelen)` or
    # `trapz(L_cm, wavenum)` will be < 0 as they are sorted in opposite order.
    # But both plots look alright

    # Convert to whatever was wanted
    l_nm = conv2(l_nm, "mW/cm2/sr/nm", Iunit)

    return l_nm


def convert_universal(
    I,
    from_unit,
    to_unit,
    spec=None,
    per_nm_is_like="mW/sr/cm2/nm",
    per_cm_is_like="mW/sr/cm2/cm-1",
):
    """Return variable var in whatever unit, and converts to to_unit Also deal
    with cases where var is in ~1/nm (per_nm_is_like) or ~1/cm-1
    (per_cm_is_like)

    Parameters
    ----------
    var: str
        variable to get. Usually 'radiance' or 'radiance_noslit'
    to_unit: str
        unit to convert variable to

    Other Parameters
    ----------------
    spec: :class:`~radis.spectrum.spectrum.Spectrum` object
        needed to get wavenumber in case we need to do a change of variable
        within the integral

    Notes
    -----
    wavenumber is needed in case we convert from ~1/nm to ~1/cm-1 (requires
    a change of variable in the integral)
    """
    import astropy.units as u

    Iunit0 = from_unit
    Iunit = to_unit
    try:
        if is_homogeneous(Iunit0, per_nm_is_like) and is_homogeneous(
            Iunit, per_cm_is_like
        ):
            w_cm = spec.get_wavenumber()
            I = convert_rad2cm(I, w_cm, Iunit0, Iunit)
            # note that there may still be a new DimensionalityError
            # raise here if the input was non-sense.
        elif is_homogeneous(Iunit0, per_cm_is_like) and is_homogeneous(
            Iunit, per_nm_is_like
        ):
            w_cm = spec.get_wavenumber()
            I = convert_rad2nm(I, w_cm, Iunit0, Iunit)
            # note that there may still be a new DimensionalityError
            # raise here if the input was non-sense.
        else:  # general case: just convert
            I = conv2(I, Iunit0, Iunit)
    except u.UnitConversionError:
        raise TypeError(Iunit0, Iunit)

    return I


# %% Test
if __name__ == "__main__":
    from radis.test.test_phys import test_units

    print(("Test :", test_units()))
