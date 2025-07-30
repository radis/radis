# -*- coding: utf-8 -*-
"""Tools : database of spectra, line survey, interface with Cantera."""


from .database import SpecDatabase, load_spec, plot_spec, save
from .gascomp import get_eq_mole_fraction
from .read_wav_index import get_wavno_lower_offset, offset_difference_from_lower_wavno
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
    "offset_difference_from_lower_wavno",
    "get_wavno_lower_offset",
    "plot_slit",
    "get_effective_FWHM",
    "get_FWHM",
]
