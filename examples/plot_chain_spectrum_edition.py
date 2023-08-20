# -*- coding: utf-8 -*-
"""
==============================================
Chain Editing and Lineshape Fitting a Spectrum
==============================================

RADIS includes a powerful chaining syntax to edit (crop/offset/multiply/etc) calculated
and experimental spectra. Some examples are given below.

We also do some lineshape fitting using the :py:mod:`specutils`
:py:func:`specutils.fitting.fit_lines` routine and some :py:mod:`astropy`
models, among :py:class:`astropy.modeling.functional_models.Gaussian1D`,
:py:class:`astropy.modeling.functional_models.Lorentz1D` or :py:class:`astropy.modeling.functional_models.Voigt1D`

See more loading and post-processing functions on the :ref:`Spectrum page <label_spectrum>`.

"""


from radis import Spectrum
from radis.test.utils import getTestFile

s = Spectrum.from_mat(
    getTestFile("trimmed_1857_VoigtCO_Minesi.mat"),
    "absorbance",
    wunit="cm-1",
    unit="",
    index=10,
)

# %%
# Plot default Spectrum :
s.plot()


# %%
# Now crop to the range we want to study
s.crop(2010.66, 2010.80).plot()


# %%
# It seems the baseline is slightly offset (negative). Fix it with algebraic
# operations
s += 0.003  # could have used an Astropy unit here too !
s.plot()


# %%
# Some operations, such as crop, happen "in-place" by default, i.e. they modify the Spectrum.
# If you don't want to modify the Spectrum make sure you specify ``inplace=False``
# For instance we compare the current Spectrum to a mock-up new experimental
# spectrum obtained by adding some random noise (15%) on top of the the previous one,
# with some offset and a coarser grid.
#
# Note the use of ``len(s)`` to get the number of spectral points
import matplotlib.pyplot as plt
import numpy as np

from radis.phys.units import Unit

s.normalize(inplace=False).plot(lw=1.5)
noise_array = 0.15 * (
    np.random.rand(len(s)) * Unit("")
)  # must be dimensioned to be multipled to spectra
noise_offset = np.random.rand(1) * 0.01
s2 = (-1 * (s.normalize(inplace=False) + noise_array)).offset(
    noise_offset, "cm-1", inplace=False
)
# Resample the spectrum on a coarser (1 every 3 points) grid
s2.resample(s2.get_wavenumber()[::3], inplace=True, energy_threshold=None)

s2.plot(nfig="same", lw=1.5)
plt.axhline(0)  # draw horizontal line

# %%
# Back to our original spectrum, we get the line positions using the
# :py:meth:`~radis.spectrum.spectrum.Spectrum.max` or :py:meth:`~radis.spectrum.spectrum.Spectrum.argmax` function
print("Maximum identified at ", s.argmax())


# %% Fit a Lorentzian
# For a more accurate measurement of the line position, we fit a Lorentzian lineshape
# using Astropy & Specutil models (with straightforward Radis ==> Specutil conversion)
# and compare the fitted line center.

from astropy import units as u
from astropy.modeling import models
from specutils.fitting import fit_lines

# Fit the spectrum and calculate the fitted flux values (``y_fit``)
g_init = models.Lorentz1D(amplitude=s.max(), x_0=s.argmax(), fwhm=0.2)
g_fit = fit_lines(s.to_specutils(), g_init)
w_fit = s.get_wavenumber() / u.cm
y_fit = g_fit(w_fit)

print("Fitted Line Center :  ", g_fit.x_0.value)

# Plot the original spectrum and the fitted.
import matplotlib.pyplot as plt

s.plot(lw=2)
plt.plot(w_fit, y_fit, label="Fit result")
plt.grid(True)
plt.legend()

#%%
# We can also mesure the area under the line :
import numpy as np

print("Absorbance of original line : ", s.get_integral("absorbance"))
print("Absorbance of fitted line   :", np.trapz(y_fit, w_fit))
#


#%%
# Finally note that the fitting routine can be achieved directly
# using the :py:meth:`~radis.spectrum.spectrum.Spectrum.fit_model` function:
from astropy.modeling import models

s.fit_model(models.Lorentz1D(), plot=True)

#%%
# We see that fitting a Voigt profile yields substantially better results
from astropy.modeling import models

s.fit_model(models.Voigt1D(), plot=True)
