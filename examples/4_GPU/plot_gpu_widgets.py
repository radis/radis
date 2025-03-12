# -*- coding: utf-8 -*-
"""
.. _example_real_time_gpu_spectra:

===============================================
Real-time GPU Accelerated Spectra (Interactive)
===============================================

Example using GPU sliders and GPU calculation with :py:meth:`~radis.lbl.SpectrumFactory.eq_spectrum_gpu_intereactive`

.. note::

    In some cases matplotlib immediately closes the window and returns; this is solved
    by running python in interactive mode as follows: ``python -i plot_gpu_widgets.py``.

    Systems with a dedicated GPU often have multiple devices available, because they also have an integrated GPU in the
    main processor. This can be selected by chosing a different device_id.
    Look at the device overview printed when running the GPU spectrum to see what options are available.

    s.gpu_exit() does not have to be called explicitly because it is called when the interactive window is closed.

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
