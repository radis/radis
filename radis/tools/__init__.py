# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:52:15 2015

Erwan Pannier
EM2C, CentraleSupélec, 2015
CNRS UPR 288

"""

from __future__ import absolute_import, division, print_function, unicode_literals

from .slit import (plot_slit, get_effective_FWHM, get_FWHM, convolve_with_slit,
                   recenter_slit, crop_slit)
from .database import SpecDatabase, save, load_spec, plot_spec