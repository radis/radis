# -*- coding: utf-8 -*-
"""
Spectrum-class module 
"""

from .spectrum import (Spectrum, calculated_spectrum, experimental_spectrum,
                       transmittance_spectrum, is_spectrum)
from .compare import (plot_diff, get_diff, get_distance, get_ratio, get_residual)