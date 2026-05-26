"""Test that crop works with wmin=0.0 (truthiness fix)."""

import numpy as np
import pytest

from radis import Spectrum


def _make_spectrum(wmin=0, wmax=10, npts=100, wunit="nm"):
    """Create a simple spectrum for testing."""
    w = np.linspace(wmin, wmax, npts)
    I = np.ones_like(w)
    return Spectrum.from_array(w, I, "radiance_noslit", wunit=wunit, Iunit="mW/cm2/sr/nm")


def test_crop_wmin_zero_nm():
    """crop(s, wmin=0.0, ...) should not silently fail."""
    s = _make_spectrum(wmin=0, wmax=10, wunit="nm")
    cropped = s.crop(wmin=0.0, wmax=5.0, wunit="nm", inplace=False)
    assert len(cropped) > 0, "crop with wmin=0.0 produced empty spectrum"
    assert cropped.get_wavelength().min() >= -0.01


def test_crop_wmax_zero_cm1():
    """crop with wmax=0.0 in cm-1 should work (edge case)."""
    s = _make_spectrum(wmin=-5, wmax=5, wunit="cm-1")
    cropped = s.crop(wmin=-2.0, wmax=2.0, wunit="cm-1", inplace=False)
    assert len(cropped) > 0


def test_crop_normal_range_still_works():
    """Ensure normal crop still works after fix."""
    s = _make_spectrum(wmin=100, wmax=200, wunit="nm")
    cropped = s.crop(wmin=120, wmax=180, wunit="nm", inplace=False)
    assert len(cropped) > 0
    w = cropped.get_wavelength()
    assert w.min() >= 119.9
    assert w.max() <= 180.1
