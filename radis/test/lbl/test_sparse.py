# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 17:34:58 2025

@author: Nicolas Minesi
"""

import numpy as np
import pytest

import radis
from radis import SpectrumFactory, get_residual, plot_diff

radis.config["DEBUG_MODE"] = True
sf_args = {
    # "wavenum_min": 3800,
    # "wavenum_max": 4500,
    # "wavenum_min": 500,
    # "wavenum_max": 10000,
    "molecule": "CO",
    "wavenum_min": 2000,
    "wavenum_max": 2300,
    ####
    "path_length": 0.1,
    "mole_fraction": 0.2,
    "isotope": 1,
    "pressure": 1e-5,
    "wstep": 0.004,
    "save_memory": False,
    "verbose": False,
}


@pytest.mark.fast
def test_sparse_vs_regular_methods(plot=False):
    # save original config
    init_config = radis.config["SPARSE_WAVERANGE"]
    sf = SpectrumFactory(**sf_args)

    for broadening_method in ["voigt_poly"]:
        # for broadening_method in ["voigt", 'convolve']: #the convolve method is simply too long

        # Compute Multi-grid
        radis.config["SPARSE_WAVERANGE"] = "multi_sparse_grid"
        sf = SpectrumFactory(
            **sf_args, optimization=None, broadening_method=broadening_method
        )
        sf.fetch_databank("hitran")
        s_multi = sf.eq_spectrum(700)
        s_multi.name = f"Multi grid : {s_multi.c['calculation_time']:.1f}s"

        for type_sparse in ["False", "simple"]:
            radis.config["SPARSE_WAVERANGE"] = type_sparse

            for optim in [None, "simple", "min-RMS"]:
                # for optim in ["simple", "min-RMS"]:
                # NOW compute a Spectrum
                sf = SpectrumFactory(
                    **sf_args, optimization=optim, broadening_method=broadening_method
                )
                # sf.fetch_databank("hitran")
                sf.load_databank("HITRAN-CO-TEST")
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
                assert np.isclose(residual, 0, atol=1e-8)  # 10)

    radis.config["SPARSE_WAVERANGE"] = init_config


@pytest.mark.fast
def test_sparse_with_optimization_raises_error():
    """Test that multi_sparse_grid with non-valid optimization or broadening_method raises ValueError."""

    import pytest

    init_config = radis.config["SPARSE_WAVERANGE"]
    radis.config["SPARSE_WAVERANGE"] = "multi_sparse_grid"

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

    radis.config["SPARSE_WAVERANGE"] = init_config


@pytest.mark.fast
def test_sparse_waverange_single_mode():
    """A change in the config should change the attribute sparse_waverange."""

    init_config = radis.config["SPARSE_WAVERANGE"]

    try:
        for value, expected_sparse_waverange in [
            (False, False),
            ("auto", "simple"),  # Legacy value - same behavior
            ("True", "simple"),  # Legacy value - closest bahvior
            ("simple", "simple"),
            ("multi_sparse_grid", "multi_sparse_grid"),
        ]:
            radis.config["SPARSE_WAVERANGE"] = value
            sf = SpectrumFactory(**sf_args)
            print(value, expected_sparse_waverange)
            assert sf.sparse_waverange == expected_sparse_waverange
    finally:
        radis.config["SPARSE_WAVERANGE"] = init_config


if __name__ == "__main__":
    # save original config
    init_config0 = radis.config["SPARSE_WAVERANGE"]

    # test_sparse_waverange_single_mode()
    test_sparse_vs_regular_methods(plot=True)
    # test_sparse_with_optimization_raises_error()

    # put back initial config
    radis.config["SPARSE_WAVERANGE"] = init_config0
