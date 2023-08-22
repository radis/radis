# -*- coding: utf-8 -*-
"""
==================================
Fit Multiple Voigt Lineshapes
==================================

Direct-access functions to fit a sum of three Voigt lineshapes on an
experimental spectrum.

This uses the underlying fitting routines of :py:mod:`specutils`
:py:func:`specutils.fitting.fit_lines` routine and some :py:mod:`astropy`
models, among :py:class:`astropy.modeling.functional_models.Gaussian1D`,
:py:class:`astropy.modeling.functional_models.Lorentz1D` or :py:class:`astropy.modeling.functional_models.Voigt1D`

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
# Fix baseline & Fit 3 Voigts profiles :
from astropy.modeling import models

s += 0.002
gfit, y_err = s.fit_model(
    [models.Voigt1D() for _ in range(3)], confidence=0.9545, plot=True, verbose=True
)

for mod in gfit:
    print(mod)
    print(mod.area)
