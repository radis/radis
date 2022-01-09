# -*- coding: utf-8 -*-
"""
=======================
GPU Accelerated Spectra
=======================

Example using GPU calculation with :py:meth:`~radis.lbl.SpectrumFactory.eq_spectrum_gpu`

.. note::

    in the example below, the GPU code runs on CPU, using the parameter ``emulate=True``.
    In your environment, to run the GPU code with the full power of the GPU, remove this line
    or set  ``emulate=False`` (default)

"""

from radis import SpectrumFactory

sf = SpectrumFactory(
    2100,
    2400,  # cm-1
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
    emulate=True,  # if True, runs CPU code on GPU. Set to False or remove to run on the GPU
)
print(s)

s.apply_slit(0.5)  # cm-1
s.plot("radiance", show=True)
