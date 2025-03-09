# -*- coding: utf-8 -*-
"""
=======================
GPU Accelerated Spectra (recalc_gpu() demo)
=======================

Example using GPU calculation with :py:meth:`~radis.spectrum.spectrum.Spectrum.recalc_gpu`

After producing a spectrum object with sf.eq_spectrum_gpu(), new spectra can be produced quickly
with spectrum.recalc_gpu().

.. note::

    Systems with a dedicated GPU often have multiple devices available, because they also have an integrated GPU in the
    main processor. This can be selected by chosing a different device_id.
    Look at the device overview printed when running the GPU spectrum to see what options are available.

    Make sure to call gpu_exit() at the end to release all GPU resources.
    This can't be done automatically because Radis doesn't know how long you want to keep using the GPU.
    If you want to immediately exit the GPU calculations, the keyword exit_gpu=True can be passed to
    sf.eq_spectrum_gpu(), but this is uncommon because it doesn't leverage the full power of GPU calculations.

"""

from radis import SpectrumFactory

sf = SpectrumFactory(
    2150,
    2450,  # cm-1
    # molecule="CO2",
    molecule="CO",
    isotope="1,2,3",
    wstep=0.002,
)


sf.fetch_databank("hitran")
##sf.fetch_databank("exomol")

T_list = [1000.0, 1250.0, 1500.0, 1750.0, 2000.0]

s = sf.eq_spectrum_gpu(
    Tgas=T_list[0],  # K
    pressure=1,  # bar
    mole_fraction=0.8,
    path_length=0.2,  # cm
    # device_id='intel',
    # device_id='nvidia',
)
s.apply_slit(0.5, unit="cm-1")  # cm-1
print("Plot0 finished in {:6.1f} ms".format(s.conditions["calculation_time"] * 1e3))
s.plot("radiance", wunit="nm", show=False)

for i, T in enumerate(T_list[1:]):
    s.recalc_gpu(Tgas=T)
    print(
        "Plot{:d} finished in {:6.1f} ms".format(
            i + 1, s.conditions["calculation_time"] * 1e3
        )
    )
    show = True if T == T_list[-1] else False
    s.plot("radiance", wunit="nm", show=show, nfig="same")

s.exit_gpu()
