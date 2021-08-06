# -*- coding: utf-8 -*-
"""
=====================
Multi-temperature Fit
=====================

A method to fit an experimental spectrum directly from :py:class:`~radis.lbl.factory.SpectrumFactory`,
with :py:meth:`~radis.lbl.factory.SpectrumFactory.fit_spectrum`

The method requires a fitting model. An example model is provided in :py:mod:`radis.tools.fitting` :
:py:func:`~radis.tools.fitting.Tvib12Tvib3TrotModel`. Other models can be shown.

More advanced tools for interactive fitting of multi-dimensional, multi-slabs
spectra can be found in `Fitroom <https://github.com/radis/fitroom>`__
(access on request by asking on Slack)

Typical output is similar to the
`radis-examples Multi-temperature fit <https://github.com/radis/radis-examples#1-multi-temperature-fit>`__ :

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/docs/multi-temperature-fit.gif

"""

from os.path import join

from radis import Spectrum, SpectrumFactory

# %% Get Fitted Data
from radis.test.utils import getValidationCase, setup_test_line_databases
from radis.tools.fitting import Tvib12Tvib3TrotModel

setup_test_line_databases()
# Data from Dang, adapted by Klarenaar, digitized by us
s_exp = Spectrum.from_txt(
    getValidationCase(
        join("test_CO2_3Tvib_vs_klarenaar_data", "klarenaar_2017_digitized_data.csv")
    ),
    "transmittance_noslit",
    waveunit="cm-1",
    unit="",
    delimiter=",",
    name="Klarenaar 2017",
)

# %% Calculate

sf = SpectrumFactory(
    2284.2,
    2284.6,
    wstep=0.001,  # cm-1
    pressure=20 * 1e-3,  # bar
    db_use_cached=True,
    lvl_use_cached=True,
    cutoff=1e-25,
    isotope="1,2",
    path_length=10,  # cm-1
    mole_fraction=0.1 * 28.97 / 44.07,
    broadening_max_width=1,  # cm-1
    medium="vacuum",
    export_populations=None,  # 'vib',
)
sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
sf.load_databank("HITEMP-CO2-TEST")

s_best, best = sf.fit_spectrum(
    s_exp.take("transmittance_noslit"),
    model=Tvib12Tvib3TrotModel,
    fit_parameters={
        "T12": 517,
        "T3": 2641,
        "Trot": 491,
    },
    bounds={"T12": [300, 2000], "T3": [300, 5000], "Trot": [300, 2000]},
    plot=True,
    maxiter=30,  # 👈 increase to let the fit converge
)
