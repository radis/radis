# -*- coding: utf-8 -*-
"""Core of the line-by-line calculations"""

# prevent cyclic imports:
from . import bands, broadening, calc, factory, labels, loader, overp, populations
from .calc import calc_spectrum, spectrum_test
from .factory import SpectrumFactory
from .overp import LevelsList
from .populations import PopulationFactory

# Backward compatibility alias (deprecated - use PopulationFactory instead)
BaseFactory = PopulationFactory

__all__ = [
    "LevelsList",
    "PopulationFactory",
    "SpectrumFactory",
    "calc_spectrum",
    "spectrum_test",
    "BaseFactory",  # deprecated alias
]
