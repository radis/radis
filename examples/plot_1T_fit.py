# -*- coding: utf-8 -*-
"""
.. _example_one_temperature_fit:

=================
1 temperature fit
=================

Quickly fit an experimental spectrum with a one-temperature model,
directly from :py:class:`~radis.lbl.factory.SpectrumFactory`,
with :py:meth:`~radis.lbl.factory.SpectrumFactory.fit_spectrum`

The method requires a fitting model. An example model is provided in :py:mod:`radis.tools.fitting` :
:py:func:`~radis.tools.fitting.LTEModel`. Other models can be used; such as in
the :ref:`multi-temperature fit example<example_multi_temperature_fit>`

More advanced tools for interactive fitting of multi-dimensional, multi-slabs
spectra can be found in :py:mod:`fitroom`.
Finally, the :ref:`GPU-accelerated example<example_real_time_gpu_spectra>` shows
how to obtain real-time interactive spectra.

"""

from radis import SpectrumFactory, load_spec

#%% Get Fitted Data
# Here we get an experimental spectrum from RADIS test cases. Use your own instead.
from radis.test.utils import getTestFile, setup_test_line_databases

setup_test_line_databases()

# fit range
wlmin = 4167
wlmax = 4180

s_exp = (
    load_spec(getTestFile("CO2_measured_spectrum_4-5um.spec"))
    .crop(wlmin, wlmax, "nm")
    .normalize()
    .sort()
    .offset(-0.2, "nm")
)

#%%
# Customize the :py:func:`~radis.tools.fitting.LTEModel` for our case: we add a slit
# (non fittable parameter) and normalize it

from radis.tools.fitting import LTEModel


def LTEModel_withslitnorm(factory, fit_parameters, fixed_parameters):
    s = LTEModel(factory, fit_parameters, fixed_parameters)
    # we could also have added a fittable parameter, such as an offset,
    # or made the slit width a fittable parameter.
    # ... any paramter in model_input will be fitted.
    # s.offset(model_input["offset"], 'nm')
    s.apply_slit(1.4, "nm")
    return s.take("radiance").normalize()


#%% Calculate
# using :py:meth:`~radis.lbl.factory.SpectrumFactory.fit_spectrum`

import astropy.units as u

sf = SpectrumFactory(
    wlmin * u.nm,
    wlmax * u.nm,
    wstep=0.001,  # cm-1
    pressure=1 * 1e-3,  # bar
    cutoff=1e-25,
    isotope="1,2",
    path_length=10,  # cm-1
    mole_fraction=1,
    truncation=1,  # cm-1
)
sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
sf.warnings["HighTemperatureWarning"] = "ignore"
sf.load_databank("HITRAN-CO2-TEST")

s_best, best = sf.fit_spectrum(
    s_exp.take("radiance"),
    model=LTEModel_withslitnorm,
    fit_parameters={
        "Tgas": 300,
        # "offset": 0
    },
    bounds={
        "Tgas": [300, 2000],
        # "offset": [-1, 1],
    },
    plot=True,
    solver_options={
        "maxiter": 15,  # 👈 increase to let the fit converge
        "ftol": 1e-15,
    },
    verbose=2,
)
# plot_diff(s_exp, s_best)  # , show_ruler=True)
