# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 12:55:28 2021

@author: erwan

Test the :py:func:`~radis.tools.gascomp.get_eq_mole_fraction` function.
"""

import pytest


@pytest.mark.fast
def test_get_eq_mole_fraction_without_cantera():
    """Test that ImportError is raised when Cantera is not installed."""
    try:
        import cantera  # noqa: F401

        pytest.skip("Cantera is installed, skipping import error test")
    except ImportError:
        pass

    from radis.tools.gascomp import get_eq_mole_fraction

    with pytest.raises(ImportError, match="Cantera is needed"):
        get_eq_mole_fraction("CO2:1", 3000, 1e5)


@pytest.mark.fast
def test_get_eq_mole_fraction_with_cantera():
    """Test equilibrium mole fraction calculation with Cantera."""
    ct = pytest.importorskip("cantera")  # noqa: F841

    from radis.misc.basics import all_in
    from radis.tools.gascomp import get_eq_mole_fraction

    gas = get_eq_mole_fraction("CO2:1", 3000, 1e5)

    # Check that expected species are present
    assert all_in(["CO", "CO2", "O", "O2"], gas.keys())

    # Check that mole fractions sum to ~1
    total = sum(gas.values())
    assert abs(total - 1.0) < 1e-6

    # CO2 should still be dominant at 3000K
    assert gas["CO2"] > 0.5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
