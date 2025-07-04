# -*- coding: utf-8 -*-
"""
.. _example_multi_temperature_fit:

=====================
Legacy #2: non-equilibrium CO2 (Tvib_12, Tvib_3, Trot)
=====================

A `3 temperature fitting example <https://github.com/radis/radis-examples/tree/master/multi-temperature-fit>`__ .
reproducing the validation case of Klarenaar 2017 [1]_, who calculated a transmittance
spectrum from the initial data of Dang 1973 [2]_, with a 1 rotational temperature +
3 vibrational temperature (Treanor distributions) model.

|CO2| Energies are calculated from Dunham developments in an uncoupled harmonic oscillator - rigid rotor model.
The example is based on one of
`RADIS validation cases <https://github.com/radis/radis/blob/master/radis/test/validation/test_CO2_3Tvib_vs_klarenaar.py>`_.
It makes use of the RADIS `Spectrum <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`_
class and the associated compare and load functions

Typical output is similar to the
`radis-examples Multi-temperature fit <https://github.com/radis/radis-examples#1-multi-temperature-fit>`__ :

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/docs/multi-temperature-fit.gif

The method requires a fitting model. An example model is provided in :py:mod:`radis.tools.fitting` :
:py:func:`~radis.tools.fitting.Tvib12Tvib3Trot_NonLTEModel`. Other models can be used, such
as in the :ref:`one-temperature fit example <example_one_temperature_fit>`

More advanced tools for interactive fitting of multi-dimensional, multi-slabs
spectra can be found in :py:mod:`fitroom`.
Finally, the :ref:`GPU-accelerated example<example_real_time_gpu_spectra>` shows
how to obtain real-time interactive spectra.

.. [1] Klarenaar et al 2017, "Time evolution of vibrational temperatures in a CO2 glow
       discharge measured with infrared absorption spectroscopy" doi/10.1088/1361-6595/aa902e

.. [2] Dang et al 1982, "Detailed vibrational population distributions in a CO2 laser
        discharge as measured with a tunable diode laser" doi/10.1007/BF00694640
"""

from os.path import join

from radis import Spectrum, SpectrumFactory

# %%
# Get Fitted Data
from radis.test.utils import getValidationCase, setup_test_line_databases
from radis.tools.fitting import Tvib12Tvib3Trot_NonLTEModel

setup_test_line_databases()
# Data from Dang, adapted by Klarenaar, digitized by us
s_exp = Spectrum.from_txt(
    getValidationCase(
        join("test_CO2_3Tvib_vs_klarenaar_data", "klarenaar_2017_digitized_data.csv")
    ),
    "transmittance_noslit",
    wunit="cm-1",
    unit="",
    delimiter=",",
    name="Klarenaar 2017",
)

# %%
# Calculate

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
    truncation=1,  # cm-1
    medium="vacuum",
    export_populations=None,  # 'vib',
    # parsum_mode="tabulation"
)
sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
sf.warnings["PerformanceWarning"] = "ignore"
sf.load_databank("HITEMP-CO2-TEST")

s_best, best = sf.fit_legacy(
    s_exp.take("transmittance_noslit"),
    model=Tvib12Tvib3Trot_NonLTEModel,
    fit_parameters={
        "T12": 517,
        "T3": 2641,
        "Trot": 491,
    },
    bounds={"T12": [300, 2000], "T3": [300, 5000], "Trot": [300, 2000]},
    fixed_parameters={"vib_distribution": "treanor"},
    plot=True,
    solver_options={
        "method": "TNC",
        "maxiter": 80,  # 👈 increase to let the fit converge
    },
)
