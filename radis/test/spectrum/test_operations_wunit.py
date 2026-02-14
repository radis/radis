
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 14 2026

@author: Antigravity
"""

import astropy.units as u
import numpy as np
import pytest
from radis import Spectrum

def test_crop_wunit_explicit_new():
    """Test that crop works with wunit='nm_air' and 'nm_vac', which were added."""
    
    w = np.linspace(200, 300, 100) # nm 
    I = np.zeros_like(w)
    s = Spectrum.from_array(w, I, 'radiance_noslit', wunit='nm', unit='mW/cm2/sr/nm')
    
    # Test wunit='nm_air'
    s_air = s.crop(210, 290, wunit='nm_air', inplace=False)
    assert s_air.get_wavelength().min() >= 210
    assert s_air.get_wavelength().max() <= 290

    # Test wunit='nm_vac'
    # Note: conversion happens which slightly shifts wavelengths if stored as 'nm' (air) and requested as 'nm_vac'
    # but crop function handles this internally.
    s_vac = s.crop(210, 290, wunit='nm_vac', inplace=False)
    # Just ensure it runs without error as the reproduction script failed here previously
    assert len(s_vac.get_wavelength()) < len(s.get_wavelength())

    print("New wunit options test passed.")

if __name__ == "__main__":
    test_crop_wunit_explicit_new()
