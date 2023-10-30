# -*- coding: utf-8 -*-
"""
=======================
GPU Accelerated Spectra (recalc_gpu() demo)
=======================

Example using GPU calculation with :py:meth:`~radis.spectrum.spectrum.Spectrum.recalc_gpu`

After producing a spectrum object with sf.eq_spectrum_gpu(), new spectra can be produced quickly
with spectrum.recalc_gpu().

.. note::

    make sure you pass ``exit_gpu=False`` when producing the spectrum object, otherwise
    it will destroy the GPU context which is needed for spectrum.recalc_gpu().
    Also be sure to call gpu_exit() at the end.

"""

from radis import SpectrumFactory
from radis.gpu.gpu import gpu_exit

sf = SpectrumFactory(
    2150,
    2450,  # cm-1
    molecule="CO2",
    isotope="1,2,3",
    wstep=0.002,
)

sf.fetch_databank("hitemp")

T_list = [1000.0, 1250.0, 1500.0, 1750.0, 2000.0]

s = sf.eq_spectrum_gpu(
    Tgas=T_list[0],  # K
    pressure=1,  # bar
    mole_fraction=0.8,
    path_length=0.2,  # cm
    exit_gpu=False,
)
s.apply_slit(0.5, unit="cm-1")  # cm-1
print("Plot0 finished in {:6.1f} ms".format(s.conditions["calculation_time"] * 1e3))
s.plot("absorbance", wunit="nm", show=False)

for i, T in enumerate(T_list[1:]):
    s.recalc_gpu(Tgas=T)
    print(
        "Plot{:d} finished in {:6.1f} ms".format(
            i + 1, s.conditions["calculation_time"] * 1e3
        )
    )
    show = True if T == T_list[-1] else False
    s.plot("absorbance", wunit="nm", show=show, nfig="same")

gpu_exit()
