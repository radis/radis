# -*- coding: utf-8 -*-
"""
Deprecated module - use radis.lbl.populations instead.

This module exists for backward compatibility only.
"""

import warnings

from . import populations

# Explicitly expose PopulationFactory for those using the new name via old path
PopulationFactory = populations.PopulationFactory


def __getattr__(name):
    """Lazy attribute access with deprecation warning for BaseFactory."""
    if name == "BaseFactory":
        warnings.warn(
            "BaseFactory has been renamed to PopulationFactory. "
            "Update your import to: from radis.lbl.populations import PopulationFactory",
            DeprecationWarning,
            stacklevel=2,
        )
        return populations.PopulationFactory

    # Forward other attributes from populations without warning
    # (they weren't renamed, just moved)
    if hasattr(populations, name):
        return getattr(populations, name)

    raise AttributeError(f"module 'radis.lbl.base' has no attribute '{name}'")


def __dir__():
    """Support for dir() and tab-completion."""
    return dir(populations) + ["BaseFactory"]
