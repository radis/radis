# -*- coding: utf-8 -*-
"""Tests for radis/tools/plot_tools.py - Coverage tests for ParamRange class."""

import pytest


@pytest.mark.fast
def test_param_range():
    """Test ParamRange initialization and behavior."""
    from radis.tools.plot_tools import ParamRange

    # Test defaults
    pr = ParamRange()
    assert pr.valmin == 0 and pr.valmax == 1 and pr.valinit == 0 and pr.val == 0

    # Test custom values
    pr2 = ParamRange(valmin=100, valmax=500, valinit=300)
    assert pr2.valmin == 100 and pr2.valmax == 500 and pr2.val == 300

    # Test repr contains expected values
    repr_str = repr(pr2)
    assert "ParamRange" in repr_str and "100" in repr_str and "500" in repr_str

    # Test valinit defaults to valmin
    pr3 = ParamRange(valmin=42, valmax=100)
    assert pr3.valinit == 42 and pr3.val == 42
