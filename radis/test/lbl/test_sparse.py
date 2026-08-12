# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 17:34:58 2025

@author: Nicolas Minesi
"""

import numpy as np

import radis
from radis import SpectrumFactory, get_residual, plot_diff

radis.config["DEBUG_MODE"] = True
sf_args = {
    # "wavenum_min": 3800,
    # "wavenum_max": 4500,
    "wavenum_min": 500,
    "wavenum_max": 10000,
    "molecule": "CO",
    # "wavenum_min": 2008,
    # "wavenum_max": 2009,
    ####
    "path_length": 0.1,
    "mole_fraction": 0.2,
    "isotope": 1,
    "pressure": 1e-5,
    "wstep": 0.004,
    "save_memory": False,
    "verbose": False,
}


def test_sparse_vs_regular_methods(plot=True):
    # save original config
    init_config = radis.config["MULTI_SPARSE_GRID"]
    sf = SpectrumFactory(**sf_args)

    for broadening_method in ["voigt_poly"]:
        # for broadening_method in ["voigt", 'convolve']: #the convolve method is simply too long

        # Compute Multi-grid
        radis.config["MULTI_SPARSE_GRID"] = True
        sf = SpectrumFactory(
            **sf_args, optimization=None, broadening_method=broadening_method
        )
        sf.fetch_databank("hitran")
        s_multi = sf.eq_spectrum(700)
        s_multi.name = f"Multi grid : {s_multi.c['calculation_time']:.1f}s"

        for optim in [None, "simple", "min-RMS"]:
            # for optim in ["simple", "min-RMS"]:
            # NOW compute a Spectrum
            if optim == None and broadening_method == "convolve":
                pass
            radis.config["MULTI_SPARSE_GRID"] = False
            sf = SpectrumFactory(
                **sf_args, optimization=optim, broadening_method=broadening_method
            )
            sf.fetch_databank("hitran")
            s_single = sf.eq_spectrum(700)
            s_single.name = f"1 grid : {s_single.c['calculation_time']:.1f}s"

            plot_diff(
                s_single,
                s_multi,
                "abscoeff",
                yscale="log",
                method="diff",
                title=f"Optimization x broad. method: {optim} x {broadening_method}",
            )

            residual = get_residual(s_single, s_multi, "abscoeff")
            print(f"******** {optim}: {residual:.2e} ********")
            assert np.isclose(residual, 0, atol=2e-10)

    radis.config["MULTI_SPARSE_GRID"] = init_config


def test_sparse_with_optimization_raises_error():
    """Test that MULTI_SPARSE_GRID with optimization raises NotImplementedError."""

    import pytest

    radis.config["MULTI_SPARSE_GRID"] = True

    for optim in ["simple", "min-RMS"]:
        with pytest.raises(ValueError):
            sf = SpectrumFactory(
                **sf_args,
                optimization=optim,
            )

    for broadening_method in ["voigt", "convolve"]:
        with pytest.raises(ValueError):
            sf = SpectrumFactory(
                **sf_args,
                optimization=optim,
            )
            sf.fetch_databank("hitran")  # to avoid flake8 error in linting


if __name__ == "__main__":
    # save original config
    init_config0 = radis.config["MULTI_SPARSE_GRID"]

    # sparse_vs_regular(optim='simple', broadening_method='fft')
    test_sparse_vs_regular_methods(plot=True)
    test_sparse_with_optimization_raises_error()
    # put back initial config
    radis.config["MULTI_SPARSE_GRID"] = init_config0
