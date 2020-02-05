# -*- coding: utf-8 -*-
"""
Spectrum-class module 
"""

from .spectrum import Spectrum, is_spectrum
from .models import calculated_spectrum, experimental_spectrum, transmittance_spectrum
from .compare import (
    plot_diff,
    get_diff,
    get_distance,
    get_ratio,
    get_residual,
    get_residual_integral,
)
from .operations import (
    Transmittance,
    Radiance,
    Radiance_noslit,
    Transmittance_noslit,
    PerfectAbsorber,
    get_baseline,
    sub_baseline,
)
