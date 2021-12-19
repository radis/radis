# -*- coding: utf-8 -*-
"""
Example using GPU sliders and GPU calculation with :py:meth:`~radis.lbl.SpectrumFactory.eq_spectrum_gpu`

Use ``emulate=True`` to run the GPU code on CPU, and ``emulate=False`` (default)
to run it with the full power of the GPU

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
    emulate=False,  # runs on GPU
)
print(s)

s.apply_slit(0.5)  # cm-1
s.plot("radiance", show=True)
