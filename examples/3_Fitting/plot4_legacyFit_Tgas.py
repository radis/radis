# -*- coding: utf-8 -*-
"""
.. _example_one_temperature_fit:

=================
Legacy #1: Temperature fit of CO2 spectrum
=================

Quickly fit an experimental spectrum with a one-temperature model,
directly from :py:class:`~radis.lbl.factory.SpectrumFactory`,
with :py:meth:`~radis.lbl.factory.SpectrumFactory.fit_legacy`

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
    # ... any parameter in model_input will be fitted.
    s.offset(0, "nm")  # or we could have used a fittable parameter below :
    # s.offset(model_input["offset"], 'nm')

    # Alternative: with a wavelength offset
    # WARNING: very sensitive parameter
    # copy_parameters = fit_parameters.copy()
    # offset = copy_parameters["offset"]
    # copy_parameters.pop("offset")
    # s = LTEModel(factory, copy_parameters, fixed_parameters)
    # s.offset(offset, 'nm')

    # Comment:
    # We could also have made the slit width a fittable parameter.
    # ... any parameter in model_input will be fitted.
    # Here we simply employ a fixed slit.
    s.apply_slit(1.4, "nm")
    return s.take("radiance").normalize()


#%% Calculate
# using :py:meth:`~radis.lbl.factory.SpectrumFactory.fit_legacy`

import astropy.units as u

sf = SpectrumFactory(
    wlmin * u.nm,
    wlmax * u.nm,
    molecule="CO2",
    wstep=0.001,  # cm-1
    pressure=1 * 1e-3,  # bar
    cutoff=1e-25,
    isotope="1,2",
    path_length=10,  # cm-1
    mole_fraction=1,
    truncation=1,  # cm-1
    verbose=0,
)
sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
sf.warnings["HighTemperatureWarning"] = "ignore"
sf.load_databank(
    "HITRAN-CO2-TEST"
)  # see 'fetch_databank' below for a more general application
# sf.fetch_databank("hitemp") #use "hitemp" or another database

s_best, best = sf.fit_legacy(
    s_exp.take("radiance"),
    model=LTEModel_withslitnorm,
    fit_parameters={
        "Tgas": 1450,
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
        # "gtol": 1e-10,
        # "eps":1e-5
    },
    verbose=2,
)
# plot_diff(s_exp, s_best)  # , show_ruler=True)
