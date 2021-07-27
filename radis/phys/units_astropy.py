# -*- coding: utf-8 -*-
"""

-------------------------------------------------------------------------------


"""


def convert_and_strip_units(quantity, output_unit=None, digit=10):
    """Strips units and return the numerical value.

    Parameters
    ----------
    quantity : int or float or list or None or `~astropy.units.quantity.Quantity`
        Numerical quantity. Pass it without including units if default units are intended.
    output_unit : `~astropy.units.core.UnitBase` or `~astropy.units.quantity.Quantity`
        The default units in which the quantity is to be converted before extracting the value.

    Other Parameters
    ----------------

    round: int
        round to that number of digit after conversion. Default ``10``. Set to
        ``None`` to disable.

    Returns
    -------
    int or float or None or ~numpy.ndarray
        The numerical value extracted from ``quantity``

    Examples
    --------

        >>> convert_and_strip_units(4.5 * u.um, u.nm)
        <<< 4500.0

        >>> convert_and_strip_units(200000 * (u.m)**-1, 1/u.cm)
        <<< 2000.0

    Raises
    ------
    TypeError
        Raised when ``quantity`` is a astropy.units quantity and ``output_unit`` is ``None``.
    """

    import astropy.units as u

    if isinstance(quantity, u.Quantity):
        from astropy.units.imperial import deg_F

        # TODO : make it possible to test if not adimensionned, before loading.
        # This would allow not to have to load Astropy.units on start-up?
        if output_unit in (u.deg_C, deg_F, u.K):
            quantity = quantity.to_value(output_unit, equivalencies=u.temperature())
        elif isinstance(output_unit, (u.UnitBase, u.Quantity)):
            quantity = quantity.to_value(output_unit, equivalencies=u.spectral())
        else:
            raise TypeError(
                "'output_unit' parameter is not a valid astropy unit: {0}".format(
                    output_unit
                )
            )

        if digit:
            quantity = round(quantity, digit)

    return quantity


if __name__ == "__main__":
    import pytest

    pytest.main(["../test/phys/test_units.py"])
