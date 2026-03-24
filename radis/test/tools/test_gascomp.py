# -*- coding: utf-8 -*-
"""Test the :py:func:`~radis.tools.gascomp.get_eq_mole_fraction` function."""

import pytest


@pytest.mark.fast
def test_get_eq_mole_fraction_without_cantera():
    """Test that ImportError is raised when Cantera is not installed."""
    try:
        import cantera  # noqa: F401

        pytest.skip("Cantera is installed")
    except ImportError:
        pass
    from radis.tools.gascomp import get_eq_mole_fraction

    with pytest.raises(ImportError, match="Cantera is needed"):
        get_eq_mole_fraction("CO2:1", 3000, 1e5)


@pytest.mark.fast
def test_get_eq_mole_fraction_with_cantera():
    """Test equilibrium mole fraction calculation with Cantera."""
    pytest.importorskip("cantera")
    from radis.misc.basics import all_in
    from radis.tools.gascomp import get_eq_mole_fraction

    gas = get_eq_mole_fraction("CO2:1", 3000, 1e5)
    assert all_in(["CO", "CO2", "O", "O2"], gas.keys())
    assert abs(sum(gas.values()) - 1.0) < 1e-6
    assert gas["CO2"] > 0.5
