# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 13:34:42 2020

@author: pankaj

------------------------------------------------------------------------

"""
#%%
import warnings

import pytest

import radis
from radis import SpectrumFactory, get_residual
from radis.misc.printer import printm
from radis.misc.warning import NoGPUWarning
from radis.test.utils import getTestFile


@pytest.mark.fast
def test_eq_spectrum_emulated_gpu(
    backend="cpu-cuda", verbose=False, plot=False, *args, **kwargs
):
    """Compare Spectrum calculated in the emulated-GPU code
    :py:func:`radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu` to Spectrum
    calculated with the CPU code :py:func:`radis.lbl.factory.SpectrumFactory.eq_spectrum`
    """

    print("Backend: ", backend)
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
        backend=backend,
        name="GPU (emulate)" if backend == "cpu-cuda" else "GPU",
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
def test_eq_spectrum_gpu(plot=False, *args, **kwargs):
    """Compare Spectrum calculated in the GPU code
    :py:func:`radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu` to Spectrum
    calculated with the CPU code :py:func:`radis.lbl.factory.SpectrumFactory.eq_spectrum`
    """
    # Ensure that GPU is not deactivated (which triggers a NoGPUWarning)
    with warnings.catch_warnings():
        warnings.simplefilter("error", category=NoGPUWarning)
        test_eq_spectrum_emulated_gpu(backend="gpu-cuda", plot=plot, *args, **kwargs)


@pytest.mark.fast
def test_multiple_gpu_calls():
    from radis import SpectrumFactory

    fixed_conditions = {
        "path_length": 1,
        "wmin": 2000,
        "wmax": 2010,
        "pressure": 0.1,
        "wstep": 0.001,
        "mole_fraction": 0.01,
    }

    sf = SpectrumFactory(**fixed_conditions, broadening_method="fft", molecule="CO")
    sf.fetch_databank("hitran")

    s1_gpu = sf.eq_spectrum_gpu(
        Tgas=300, backend="gpu-cuda", diluent={"air": 0.99}  # K  # runs on GPU
    )
    s2_gpu = sf.eq_spectrum_gpu(
        Tgas=300, backend="gpu-cuda", diluent={"air": 0.99}  # K  # runs on GPU
    )

    assert abs(s1_gpu.get_power() - s2_gpu.get_power()) / s1_gpu.get_power() < 1e-5
    assert s1_gpu.get_power() > 0


# --------------------------
if __name__ == "__main__":

    # test_eq_spectrum_gpu(plot=True)
    test_eq_spectrum_emulated_gpu(plot=True, verbose=2)

    printm("Testing GPU spectrum calculation:", pytest.main(["test_gpu.py"]))
