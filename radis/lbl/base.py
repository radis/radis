# -*- coding: utf-8 -*-
"""
Deprecated module - use radis.lbl.populations instead.

This module exists for backward compatibility only.
"""

import warnings

from . import populations

# Explicitly expose PopulationFactory for those using the new name via old path
PopulationFactory = populations.PopulationFactory

__all__ = [name for name in dir(populations) if not name.startswith("_")]
for _legacy_name in ("PopulationFactory", "BaseFactory"):
    if _legacy_name not in __all__:
        __all__.append(_legacy_name)


def __getattr__(name):
    """Lazy attribute access with deprecation warning for BaseFactory."""
    if name == "BaseFactory":
        warnings.warn(
            "BaseFactory has been renamed to PopulationFactory. "
            "Update your import to: from radis.lbl.populations import PopulationFactory",
            DeprecationWarning,
            stacklevel=2,
        )
        # Cache to avoid repeated warnings on subsequent access
        globals()["BaseFactory"] = populations.PopulationFactory
        return populations.PopulationFactory

    # Forward other attributes from populations without warning
    # (they weren't renamed, just moved)
    if hasattr(populations, name):
        value = getattr(populations, name)
        globals()[name] = value
        return value

    raise AttributeError(f"module 'radis.lbl.base' has no attribute '{name}'")
