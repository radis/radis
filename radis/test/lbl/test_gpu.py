# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 13:34:42 2020

@author: pankaj

------------------------------------------------------------------------

"""
import pytest

from radis import SpectrumFactory, get_residual
from radis.misc.printer import printm
from radis.test.utils import getTestFile


def test_eq_spectrum_emulated_gpu(emulate=True, plot=False, *args, **kwargs):

    print("Emulate: ", emulate)
    T = 1000
    p = 0.1
    wstep = 0.001
    wmin = 2284.0  # cm-1
    wmax = 2285.0  # cm-1

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
    sf._broadening_method = "fft"
    sf.load_databank(
        path=getTestFile("cdsd_hitemp_09_fragment.txt"),
        format="cdsd-4000",
        parfuncfmt="hapi",
    )

    s_cpu = sf.eq_spectrum(Tgas=T, name="CPU")
    s_gpu = sf.eq_spectrum_gpu(
        Tgas=T, emulate=emulate, name="GPU (emulate)" if emulate else "GPU"
    )
    s_cpu.crop(wmin=2284.2, wmax=2284.8)  # remove edge lines
    s_gpu.crop(wmin=2284.2, wmax=2284.8)
    if plot:
        s_cpu.compare_with(s_gpu, spectra_only=True, plot=plot)
    assert get_residual(s_cpu, s_gpu, "abscoeff") < 1.4e-5
    assert get_residual(s_cpu, s_gpu, "radiance_noslit") < 7.3e-6
    assert get_residual(s_cpu, s_gpu, "transmittance_noslit") < 1.4e-5


@pytest.mark.needs_cuda
def test_eq_spectrum_gpu(plot, *args, **kwargs):
    test_eq_spectrum_emulated_gpu(emulate=False, plot=plot, *args, **kwargs)


# --------------------------
if __name__ == "__main__":

    printm("Testing GPU spectrum calculation:", pytest.main(["test_gpu.py"]))
