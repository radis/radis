# -*- coding: utf-8 -*-
"""
.. _example_real_time_gpu_spectra:

===============================================
Real-time GPU Accelerated Spectra (Interactive)
===============================================

Example using GPU sliders and GPU calculation with :py:meth:`~radis.lbl.SpectrumFactory.eq_spectrum_gpu_intereactive`

This method requires a GPU - Currently, only Nvidia GPU's are supported.
For more information on how to setup your system to run GPU-accelerated methods
using CUDA, check :ref:`GPU Spectrum Calculation on RADIS <label_radis_gpu>`

.. note::

    In some cases matplotlib immediately closes the window and returns; this is solved
    by running python in interactive mode as follows: ``python -i plot_gpu_widgets.py``.

    in the example below, the code runs on the GPU by default. In case no Nvidia GPU is
    detected, the code will instead be ran on CPU. This can be toggled manually by setting
    the ``backend`` keyword either to ``'gpu-cuda'`` or ``'cpu-cuda'``.
    The run time reported below is for CPU.

"""

from radis import SpectrumFactory
from radis.tools.plot_tools import ParamRange

sf = SpectrumFactory(
    2150,
    2450,  # cm-1
    molecule="CO2",
    isotope="1,2,3",
    wstep=0.002,
)

sf.fetch_databank("hitran")

s = sf.eq_spectrum_gpu_interactive(
    var="radiance",
    Tgas=ParamRange(300.0, 2500.0, 1500.0),  # K
    pressure=ParamRange(0.1, 2, 1),  # bar
    mole_fraction=ParamRange(0, 1, 0.8),
    path_length=ParamRange(0, 1, 0.2),  # cm
    slit_function=ParamRange(0, 1.5, 0.5),  # cm-1
    # device_id='nvidia',
    plotkwargs={"wunit": "nm"},  # "nfig": "same",
)
