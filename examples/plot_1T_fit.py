# -*- coding: utf-8 -*-
"""
=================
1 temperature fit
=================

Quickly fit an experimental spectrum with a one-temperature model,
directly from :py:class:`~radis.lbl.factory.SpectrumFactory`,
with :py:meth:`~radis.lbl.factory.SpectrumFactory.fit_spectrum`

The method requires a fitting model. An example model is provided in :py:mod:`radis.tools.fitting` :
:py:func:`~radis.tools.fitting.TrotModel`. Other models can be shown.

More advanced tools for interactive fitting of multi-dimensional, multi-slabs
spectra can be found in `Fitroom <https://github.com/radis/fitroom>`__
(access on request by asking on Slack)

"""

from radis import SpectrumFactory, load_spec, plot_diff

# %% Get Fitted Data
from radis.test.utils import getTestFile, setup_test_line_databases
from radis.tools.fitting import TrotModel  # Tvib12Tvib3TrotModel

setup_test_line_databases()

# fit range
wlmin = 4167
wlmax = 4180

s_exp = (
    load_spec(getTestFile("CO2_measured_spectrum_4-5um.spec"))
    .crop(wlmin, wlmax, "nm")
    .normalize()
    .sort()
    .offset(-0.5, "nm")
)


def TrotModel_norm(*args, **kwargs):
    s = TrotModel(*args, **kwargs)
    s.apply_slit(1.5, "nm")
    return s.take("radiance").normalize()


# %% Calculate

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
    broadening_max_width=1,  # cm-1
)
sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
sf.load_databank("HITRAN-CO2-TEST")

s_best, best = sf.fit_spectrum(
    s_exp.take("radiance"),
    model=TrotModel_norm,
    fit_parameters={
        "Trot": 300,
    },
    bounds={"Trot": [300, 2000]},
    plot=True,
    solver_options={
        "maxiter": 100,  # 👈 increase to let the fit converge
        "ftol": 1e-15,
    },
)
plot_diff(s_exp, s_best)  # , show_ruler=True)
