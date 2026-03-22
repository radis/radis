# -*- coding: utf-8 -*-
"""
Tests for radis/tools/plot_tools.py

Coverage tests for ParamRange class.
The add_ruler function requires matplotlib GUI and is tested visually.
"""

import pytest


def test_param_range_init_defaults():
    """Test ParamRange initialization with default values."""
    from radis.tools.plot_tools import ParamRange

    pr = ParamRange()
    assert pr.valmin == 0
    assert pr.valmax == 1
    assert pr.valinit == 0  # defaults to valmin when not specified
    assert pr.val == 0
    assert pr.name is None
    assert pr.widget is None


def test_param_range_init_custom():
    """Test ParamRange initialization with custom values."""
    from radis.tools.plot_tools import ParamRange

    pr = ParamRange(valmin=100, valmax=500, valinit=300)
    assert pr.valmin == 100
    assert pr.valmax == 500
    assert pr.valinit == 300
    assert pr.val == 300
    assert pr.name is None
    assert pr.widget is None


def test_param_range_repr():
    """Test ParamRange string representation."""
    from radis.tools.plot_tools import ParamRange

    pr = ParamRange(valmin=0, valmax=100, valinit=50)
    repr_str = repr(pr)

    # Check that the repr contains the expected values
    assert "ParamRange" in repr_str
    assert "0" in repr_str
    assert "100" in repr_str
    assert "50" in repr_str


def test_param_range_valinit_defaults_to_valmin():
    """Test that valinit defaults to valmin when not provided."""
    from radis.tools.plot_tools import ParamRange

    pr = ParamRange(valmin=42, valmax=100)
    assert pr.valinit == 42
    assert pr.val == 42


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
