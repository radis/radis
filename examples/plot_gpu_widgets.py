# -*- coding: utf-8 -*-
"""
.. _example_real_time_gpu_spectra:

===============================================
Real-time GPU Accelerated Spectra (Interactive)
===============================================

Example using GPU sliders and GPU calculation with :py:meth:`~radis.lbl.SpectrumFactory.eq_spectrum_gpu`

This method requires CUDA compatible hardware to execute.
For more information on how to setup your system to run GPU-accelerated methods
using CUDA and Cython, check :ref:`GPU Spectrum Calculation on RADIS <label_radis_gpu>`

.. note::

    in the example below, the GPU code runs on CPU, using the parameter ``emulate=True``.
    In your environment, to run the GPU code with the full power of the GPU, remove this line
    or set  ``emulate=False`` (default)

"""

from radis import SpectrumFactory

# from radis.test.utils import getTestFile
from radis.tools.plot_tools import ParamRange

## Add an experimental file :
# my_file = getTestFile("CO2_measured_spectrum_4-5um.spec")  # for the example here
# s_exp = load_spec(my_file)
# s_exp.crop(4120, 4790).plot(Iunit="mW/cm2/sr/nm")
#
## This spectrum is significantly absorbed by atmospheric CO2
## so it will never match the synthetic spectrum.
## TODO: find different spectrum for this example.


sf = SpectrumFactory(
    2150,
    2450,  # cm-1
    molecule="CO2",
    isotope="1,2,3",
    wstep=0.002,
)

sf.fetch_databank("hitemp")

s = sf.eq_spectrum_gpu_interactive(
    var="radiance",
    Tgas=ParamRange(300.0, 2500.0, 1100.0),  # K
    pressure=ParamRange(0.1, 2, 1),  # bar
    mole_fraction=ParamRange(0, 1, 0.8),
    path_length=ParamRange(0, 1, 0.2),  # cm
    slit_FWHM=ParamRange(0, 1.5, 0.24),  # cm-1
    emulate=True,  # if True, runs CPU code on GPU. Set to False or remove to run on the GPU
    plotkwargs={"wunit": "nm"},  # "nfig": "same",
)

print(s)
