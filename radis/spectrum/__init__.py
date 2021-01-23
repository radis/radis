# -*- coding: utf-8 -*-
"""Spectrum-class module."""

from .compare import (
    get_diff,
    get_distance,
    get_ratio,
    get_residual,
    get_residual_integral,
    plot_diff,
)
from .models import calculated_spectrum, experimental_spectrum, transmittance_spectrum
from .operations import (
    PerfectAbsorber,
    Radiance,
    Radiance_noslit,
    Transmittance,
    Transmittance_noslit,
    get_baseline,
    sub_baseline,
)
from .spectrum import Spectrum, is_spectrum
