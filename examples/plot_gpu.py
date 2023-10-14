# -*- coding: utf-8 -*-
"""
=======================
GPU Accelerated Spectra
=======================

Example using GPU calculation with :py:meth:`~radis.lbl.SpectrumFactory.eq_spectrum_gpu`

This method requires a GPU - Currently, only Nvidia GPU's are supported.
For more information on how to setup your system to run GPU-accelerated methods
using CUDA, check :ref:`GPU Spectrum Calculation on RADIS <label_radis_gpu>`

.. note::

    in the example below, the code runs on the GPU by default. In case no Nvidia GPU is
    detected, the code will instead be ran on CPU. This can be toggled manually by setting
    the ``backend`` keyword either to ``'gpu-cuda'`` or ``'cpu-cuda'``.
    The run time reported below is for CPU.

"""

from radis import SpectrumFactory, plot_diff

sf = SpectrumFactory(
    2150,
    2450,  # cm-1
    molecule="CO2",
    isotope="1,2,3",
    wstep=0.002,
)

sf.fetch_databank("hitemp")

T = 1500.0  # K
p = 1.0  # bar
x = 0.8
l = 0.2  # cm
w_slit = 0.5  # cm-1

# s_cpu = sf.eq_spectrum(
    # name="CPU",
    # Tgas=T,
    # pressure=p,
    # mole_fraction=x,
    # path_length=l,
# )
# s_cpu.apply_slit(w_slit, unit="cm-1")

s_gpu = sf.eq_spectrum_gpu(
    name="GPU",
    Tgas=T,
    pressure=p,
    mole_fraction=x,
    path_length=l,
    backend="gpu-cuda",
)
#s_gpu.apply_slit(w_slit, unit="cm-1")

s_gpu.plot()
#plot_diff(s_cpu, s_gpu, var="radiance", wunit="nm", method="diff")
