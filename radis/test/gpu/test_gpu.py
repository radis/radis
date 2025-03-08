# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 13:34:42 2020

@author: pankaj

------------------------------------------------------------------------

"""
#%%
import warnings

import pytest
from numpy import allclose

import radis
from radis import SpectrumFactory, get_residual
from radis.misc.printer import printm
from radis.misc.warning import NoGPUWarning
from radis.test.utils import getTestFile


@pytest.mark.fast
def test_eq_spectrum_gpu(device_id=0, verbose=False, plot=False, *args, **kwargs):
    """Compare Spectrum calculated in the GPU code
    :py:func:`radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu` to Spectrum
    calculated with the CPU code :py:func:`radis.lbl.factory.SpectrumFactory.eq_spectrum`
    """

    T = 1000
    p = 0.1
    wstep = 0.001
    wmin = 2284.0  # cm-1
    wmax = 2285.0  # cm-1

    # Compare GPU & CPU without sparse-LDM (not implemented in GPU yet)
    radis.config["SPARSE_WAVERANGE"] = False

    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=0.01,  # until self and air broadening is implemented
        path_length=1,  # doesnt change anything
        wstep=wstep,
        pressure=p,
        isotope="1,2",
        warnings={
            "MissingSelfBroadeningWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "HighTemperatureWarning": "ignore",
            "GaussianBroadeningWarning": "ignore",
        },
    )
    sf.params.broadening_method = "fft"
    sf.load_databank(
        path=getTestFile("cdsd_hitemp_09_fragment.txt"),
        format="cdsd-hitemp",
        parfuncfmt="hapi",
    )

    s_cpu = sf.eq_spectrum(Tgas=T, name="CPU")
    s_cpu.name += f" [{s_cpu.c['calculation_time']:.2f}s]"
    s_gpu = sf.eq_spectrum_gpu(
        Tgas=T,
        name="GPU",
        exit_gpu=True,
    )
    s_gpu.name += f"[{s_gpu.c['calculation_time']:.2f}s]"
    s_cpu.crop(wmin=2284.2, wmax=2284.8)  # remove edge lines
    s_gpu.crop(wmin=2284.2, wmax=2284.8)
    if plot:
        s_cpu.compare_with(s_gpu, spectra_only=True, plot=plot)
    assert get_residual(s_cpu, s_gpu, "abscoeff") < 1.4e-5
    assert get_residual(s_cpu, s_gpu, "radiance_noslit") < 7.3e-6
    assert get_residual(s_cpu, s_gpu, "transmittance_noslit") < 1.4e-5

    if verbose >= 2:
        print(s_gpu)
        print("\n" * 2)
        print(s_cpu)
        print("\n" * 2)
        s_gpu.print_perf_profile()
        print("\n" * 2)
        s_cpu.print_perf_profile()


@pytest.mark.needs_cuda
def test_eq_spectrum_gpu_cuda(plot=False, *args, **kwargs):
    """Compare Spectrum calculated in the GPU code
    :py:func:`radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu` to Spectrum
    calculated with the CPU code :py:func:`radis.lbl.factory.SpectrumFactory.eq_spectrum`,
    specifically using an Nvidia GPU.
    """
    # Ensure that GPU is not deactivated (which triggers a NoGPUWarning)
    with warnings.catch_warnings():
        warnings.simplefilter("error", category=NoGPUWarning)
        test_eq_spectrum_gpu(device_id="nvidia", plot=plot, *args, **kwargs)


@pytest.mark.fast
def test_multiple_gpu_calls(plot=False, hard_test=True):
    from radis import SpectrumFactory

    fixed_conditions = {
        "wmin": 2000,
        "wmax": 2010,
        "pressure": 0.1,
        "wstep": 0.001,
        "mole_fraction": 0.01,
    }
    if plot:
        import matplotlib.pyplot as plt

    sf = SpectrumFactory(**fixed_conditions, broadening_method="fft", molecule="CO")
    sf.fetch_databank("hitran")

    s_gpu = sf.eq_spectrum_gpu(
        Tgas=300.0,
        path_length=1.0,
        backend="vulkan",
        diluent={"air": 0.99},  # K  # runs on GPU
    )
    wl, A1 = s_gpu.get("absorbance")
    A2 = s_gpu.recalc_gpu("absorbance", path_length=3.0)
    s_gpu.exit_gpu()

    if plot:
        plt.plot(wl, 3 * A1, "k", lw=3)
        plt.plot(wl, A2, "r", lw=1)

    assert allclose(3 * A1, A2, rtol=1e-4)
    assert allclose(3 * A1, A2, atol=1e-6)

    if not hard_test:
        if plot:
            plt.show()
        return

    # The second part of the tests spawns a new gpuApp.
    # this is only possible if all cleaning up went well the first time.

    s_gpu3 = sf.eq_spectrum_gpu(
        Tgas=300.0,
        path_length=2.0,
        backend="vulkan",
        diluent={"air": 0.99},  # K  # runs on GPU
        exit_gpu=True,
    )
    wl, A3 = s_gpu3.get("absorbance")

    if plot:
        import matplotlib.pyplot as plt

        plt.plot(wl, 1.5 * A3, "m+")

    assert allclose(2 * A1, A3, rtol=1e-4)

    if plot:
        plt.show()


# --------------------------
if __name__ == "__main__":

    # test_eq_spectrum_gpu(plot=True, verbose=2)
    # test_eq_spectrum_gpu_cuda(plot=True)
    # test_multiple_gpu_calls(plot=True, hard_test=True)

    printm("Testing GPU spectrum calculation:", pytest.main(["test_gpu.py"]))
