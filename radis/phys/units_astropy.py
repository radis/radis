import astropy.units as u


def convert_and_strip_units(quantity, output_unit=None):
    """
    Strips units and return the numerical value.

    Parameters
    ----------
    quantity : int or float or list or None or ~astropy.units.quantity.Quantity
        Numerical quantity. Pass it without including units if default units are intended.
    output_unit : ~astropy.units.core.UnitBase
        The default units in which the quantity is to be converted before extracting the value.

    Returns
    -------
    int or float or None or ~numpy.ndarray
        The numerical value extracted from ``quantity``

    Raises
    ------
    TypeError
        Raised when ``quantity`` is a astropy.units quantity and ``output_unit`` is ``None``.
    
    """
    if isinstance(quantity, u.Quantity):
        if isinstance(output_unit, u.UnitBase):
            if output_unit in (u.deg_C, u.imperial.deg_F, u.K):
                return quantity.to(output_unit, equivalencies = u.temperature()).value
            else:
                return quantity.to(output_unit).value
        else:
            raise TypeError("'output_unit' parameter is not a valid astropy unit")
    else:
        return quantity
