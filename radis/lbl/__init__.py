# -*- coding: utf-8 -*-
"""Core of the line-by-line calculations"""

import warnings

# prevent cyclic imports:
from . import bands, base, broadening, calc, factory, labels, loader, overp, populations
from .calc import calc_spectrum, spectrum_test
from .factory import SpectrumFactory
from .overp import LevelsList
from .populations import PopulationFactory

__all__ = [
    "LevelsList",
    "PopulationFactory",
    "SpectrumFactory",
    "calc_spectrum",
    "spectrum_test",
]


def __getattr__(name):
    """Lazy access with deprecation warning for BaseFactory."""
    if name == "BaseFactory":
        warnings.warn(
            "BaseFactory is deprecated, use PopulationFactory instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        # Cache to avoid duplicate warnings during from-import resolution.
        globals()["BaseFactory"] = PopulationFactory
        return PopulationFactory
    raise AttributeError(f"module 'radis.lbl' has no attribute '{name}'")
