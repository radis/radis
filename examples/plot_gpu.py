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

from radis import SpectrumFactory

sf = SpectrumFactory(
    2150,
    2450,  # cm-1
    molecule="CO2",
    isotope="1,2,3",
    wstep=0.002,
)

sf.fetch_databank("hitemp")

s = sf.eq_spectrum_gpu(
    Tgas=1100.0,  # K
    pressure=1,  # bar
    mole_fraction=0.8,
    path_length=0.2,  # cm
)
print(s)

s.apply_slit(0.5)  # cm-1
s.plot("radiance", show=True)
