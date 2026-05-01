# -*- coding: utf-8 -*-
"""Tools : database of spectra, line survey, interface with Cantera."""


from .database import SpecDatabase, load_spec, plot_spec, save
from .gascomp import get_eq_mole_fraction
from .slit import (
    convolve_with_slit,
    crop_slit,
    get_effective_FWHM,
    get_FWHM,
    plot_slit,
    recenter_slit,
)
from .uncertainty import (
    ConfidenceBandResult,
    SensitivityAnalyzer,
    UncertaintyModel,
    UncertaintyPropagator,
    decode_geisa_uncertainty,
    decode_hitran_uncertainty,
    parse_hitran_ierr,
)

__all__ = [
    "SpecDatabase",
    "load_spec",
    "plot_spec",
    "save",
    "get_eq_mole_fraction",
    "plot_slit",
    "get_effective_FWHM",
    "get_FWHM",
    "decode_hitran_uncertainty",
    "decode_geisa_uncertainty",
    "parse_hitran_ierr",
    "UncertaintyModel",
    "UncertaintyPropagator",
    "SensitivityAnalyzer",
    "ConfidenceBandResult",
]
