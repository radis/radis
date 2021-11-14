# -*- coding: utf-8 -*-
"""
============================
Post-process using Specutils
============================

Find peaks or uncertainties using the :py:mod:`specutils` library. A Radis Spectrum
object can easily be converted to a ``specutils`` :py:class:`specutils.spectra.spectrum1d.Spectrum1D`
using :py:meth:`~radis.spectrum.spectrum.Spectrum.to_specutils`.

Below, we create a noisy spectrum based on a synthetic CO spectrum,
we convert it to :py:mod:`specutils`, add uncertainties by targeting a
noisy region, then determine the lines using :py:func:`~specutils.fitting.find_lines_threshold` :

"""

import astropy.units as u
import numpy as np

from radis import test_spectrum

""" We create a synthetic CO spectrum"""

s = (
    test_spectrum(molecule="CO", wavenum_min=2000, wavenum_max=2030)
    .apply_slit(1.5, "nm")
    .take("radiance")
)
s.trim()  # removes nans created by the slit convolution boundary effects
noise = np.random.normal(0.0, s.max().value * 0.03, len(s))
s_exp = s + noise

s_exp.plot()

#%% Determine the noise level by selecting a noisy region from the graph above.

spectrum = s_exp.to_specutils()

from specutils import SpectralRegion
from specutils.manipulation import noise_region_uncertainty

noise_region = SpectralRegion(2010.5 / u.cm, 2009.5 / u.cm)
spectrum = noise_region_uncertainty(spectrum, noise_region)


#%% Find lines :

from specutils.fitting import find_lines_threshold

lines = find_lines_threshold(spectrum, noise_factor=2)

print(lines)


s_exp.plot(lw=2, show_ruler=True)
import matplotlib.pyplot as plt

for line in lines.to_pandas().line_center.values:
    plt.axvline(line, color="r", zorder=-1)
s.plot(nfig="same")
plt.axvspan(noise_region.lower.value, noise_region.upper.value, color="b", alpha=0.1)


#%% Note: we can also create a RADIS spectrum object from the Specutils Spectrum1D :

from radis import Spectrum

s2 = Spectrum.from_specutils(spectrum)

s2.plot(Iunit="mW/cm2/sr/nm", wunit="nm")
s_exp.plot(Iunit="mW/cm2/sr/nm", wunit="nm", nfig="same")
assert s_exp == s2
