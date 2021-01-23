# -*- coding: utf-8 -*-
"""Created on Tue May 26 11:52:15 2015.

Erwan Pannier EM2C, CentraleSupélec, 2015 CNRS UPR 288
"""


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
