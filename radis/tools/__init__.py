# -*- coding: utf-8 -*-
"""Tools : database of spectra, line survey, interface with Cantera."""


from .database import SpecDatabase, load_spec, plot_spec, save
from .gascomp import get_eq_mole_fraction
from .read_wav_index import get_wavno_for_byte_range, get_wavno_lower_offset
from .slit import (
    convolve_with_slit,
    crop_slit,
    get_effective_FWHM,
    get_FWHM,
    plot_slit,
    recenter_slit,
)

__all__ = [
    "SpecDatabase",
    "load_spec",
    "plot_spec",
    "save",
    "get_eq_mole_fraction",
    "get_wavno_for_byte_range",
    "get_wavno_lower_offset",
    "plot_slit",
    "get_effective_FWHM",
    "get_FWHM",
]
