# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 17:34:58 2025

@author: Nicolas Minesi
"""

import numpy as np

import radis

radis.config["DEBUG_MODE"] = True


def test_sparse_vs_regular():
    from radis import SpectrumFactory

    sf = SpectrumFactory(
        # wavenum_min=3800,
        # wavenum_max=4500,
        # wavenum_min=500,
        # wavenum_max=10000,
        molecule="CO",
        wavenum_min=2008,
        wavenum_max=2009,
        Tgas=300,
        path_length=0.1,
        mole_fraction=0.2,
        isotope=1,
        pressure=1e-5,
        wstep=0.004,
        databank="hitran",  # or 'hitemp'
        optimization=None,
        # also measure interpolation time
        return_factory=True,
        save_memory=False,
    )

    #%%

    # #Test:

    # %%
    # keep only few lines
    # sf.df0= sf.df0.iloc[150:153]
    # keep only so many lines
    # sf.df0 = sf.df0.iloc[np.arange(len(sf.df0))[::10]]

    #%%

    # NOW compute a Spectrum
    import radis

    radis.config["MULTI_SPARSE_GRID"] = False
    s_single = sf.eq_spectrum(700)
    s_single.plot("abscoeff", yscale="log", lw=4, nfig=10)

    # %%

    # NOW compute a Spectrum
    import radis

    radis.config["MULTI_SPARSE_GRID"] = True
    s_multi = sf.eq_spectrum(700)

    s_multi.plot("abscoeff", yscale="log", nfig=10, lw=2)

    # %% Plot the graphs

    # import matplotlib.pyplot as plt
    # from radis.lbl.factory import _generate_wavenumber_range_sparse
    # wstep_calc_narrow = self.params.wstep
    # truncation = self.params.truncation
    # neighbour_lines = self.params.neighbour_lines

    # plt.figure()
    # plt.xlabel("Wavenumber (cm-1)")
    # plt.ylabel("Wstep")
    # plt.yscale("log")
    # wavenumber_arrays = []

    # for i, (wstep, lineshape_half_width) in enumerate(zip(self._wstep_multigrid,
    #                                                       self._truncation_multigrid)):
    #     print("wstep", wstep)
    #     print("lineshape_half_width", lineshape_half_width)
    #     (
    #         wavenumber,
    #         wavenumber_calc,
    #         woutrange,
    #         ix_ranges,
    #     ) = _generate_wavenumber_range_sparse(
    #         self.input.wavenum_min,
    #         self.input.wavenum_max,
    #         wstep,
    #         neighbour_lines,
    #         self.df1.wav,
    #         lineshape_half_width,
    #     )
    #     wavenumber_arrays.append(wavenumber)

    #     from publib import keep_color
    #     for wavenumber_group in wavenumber: # iterate over all wavenumber groups :
    #         plt.plot(wavenumber_group, np.ones_like(wavenumber_group)*wstep, "-o", ms=3)
    #         if wavenumber_group[-1] != wavenumber[-1][-1]: # not the last element :
    #             keep_color()
    #     # Plot line centers :
    #     for line in self.df1.wav:
    #         plt.axvline(line,color='k', alpha=0.05, zorder=-1)
    # plt.title("3 grids & Line centers\n")

    # %% plot diff
    s_multi.name = f"Multi grid : {s_multi.c['calculation_time']:.1f}s"
    # s.name = f"{s.c['calculation_time']:.1f}s"
    s_single.name = f"1 grid : {s_single.c['calculation_time']:.1f}s"

    from radis import get_residual, plot_diff

    plot_diff(s_single, s_multi, "abscoeff", yscale="log", method="diff")

    residual = get_residual(s_single, s_multi, "abscoeff")
    assert np.isclose(residual, 0, atol=1e-8)


if __name__ == "__main__":
    test_sparse_vs_regular()
