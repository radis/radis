# -*- coding: utf-8 -*-
"""
=======================
GPU Accelerated Spectra
=======================

Example using GPU calculation with :py:meth:`~radis.lbl.SpectrumFactory.eq_spectrum_gpu`

.. note::

    Systems with a dedicated GPU often have multiple devices available, because they also have an integrated GPU in the 
    main processor. This can be selected by chosing a different device_id.
    Look at the device overview printed when running the GPU spectrum to see what options are available.

    Make sure to call gpu_exit() at the end to release all GPU resources.
    This can't be done automatically because Radis doesn't know how long you want to keep using the GPU.
    If you want to immediately exit the GPU calculations, the keyword exit_gpu=True can be passed to
    sf.eq_spectrum_gpu(), but this is uncommon because it doesn't leverage the full power of GPU calculations.

"""

from radis import SpectrumFactory, plot_diff

sf = SpectrumFactory(
    2150,
    2450,  # cm-1
    molecule="CO2",
    isotope="1",
    wstep=0.002,
)

sf.fetch_databank(
    source="hitran"
)  # use hitemp or exomol for accuracy at high tempertatures

T = 1500.0  # K
p = 1.0  # bar
x = 0.8
l = 0.2  # cm
w_slit = 0.5  # cm-1

s_cpu = sf.eq_spectrum(
    name="CPU",
    Tgas=T,
    pressure=p,
    mole_fraction=x,
    path_length=l,
)
s_cpu.apply_slit(w_slit, unit="cm-1")

s_gpu = sf.eq_spectrum_gpu(
    name="GPU",
    Tgas=T,
    pressure=p,
    mole_fraction=x,
    path_length=l,
    # device_id='nvidia',
    exit_gpu=True,
)
s_gpu.apply_slit(w_slit, unit="cm-1")

plot_diff(s_cpu, s_gpu, var="emissivity", wunit="nm", method="diff")

# s_gpu.exit_gpu()
