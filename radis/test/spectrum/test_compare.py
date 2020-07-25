# -*- coding: utf-8 -*-
"""

Examples
--------

Run all tests::

    pytest      # (in command line, in project folder)

Run only fast tests (i.e: tests that a  'fast' label)::

    pytest -m fast

-------------------------------------------------------------------------------


"""

from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from radis.test.utils import getTestFile
from radis.tools.database import load_spec
from radis.spectrum.compare import get_distance, plot_diff
from radis import calc_spectrum, plot_diff, Radiance_noslit, get_residual

# Test routines


def test_compare_methods(verbose=True, plot=True, close_plots=True, *args, **kwargs):
    """ Just run all Spectrum compare methods to check they work"""

    if plot and close_plots:
        import matplotlib.pyplot as plt

        plt.close("all")

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)

    # limits to a single line, because get_distance()
    s.resample(np.linspace(2193, 2193.8, 100))
    # is very computationaly heavy
    s.update("radiance_noslit")
    s_noabsorption = s.copy()
    s.name = "solve RTE"
    s_noabsorption.name = "optically thin"

    # rescale, normal
    #    s.rescale_mole_fraction(10)
    s.rescale_path_length(10)

    # rescale, optically thin mode
    s_noabsorption.conditions["self_absorption"] = False
    #    s_noabsorption.rescale_mole_fraction(10)
    s_noabsorption.rescale_path_length(10)

    # Compare
    # should be added in an example with experimental fit of bandhead
    get_distance(s, s_noabsorption, "radiance_noslit")
    title = "CO x={0:.2f}, L={1:.2f}cm".format(
        s.conditions["mole_fraction"], s.conditions["path_length"]
    )
    if plot:
        plot_diff(s, s_noabsorption, method="diff", diff_window=1, title=title)
        plot_diff(s, s_noabsorption, method="ratio", normalize=True, title=title)


# %%


def test_plot_compare_with_nan(
    plot=True, verbose=True, close_plots=True, *args, **kwargs
):

    if plot and close_plots:
        import matplotlib.pyplot as plt

        plt.close("all")

    s = calc_spectrum(
        1900,
        2300,  # cm-1
        molecule="CO",
        isotope="1,2,3",
        pressure=1.01325,  # bar
        Tgas=700,  # K
        mole_fraction=0.1,
        path_length=1,  # cm
    )

    s = Radiance_noslit(s)
    s._q["radiance_noslit"][0] = np.nan

    # Test get_residual methods when there are nans in the spectra
    diff1 = get_residual(
        s,
        s * 1.2,
        var="radiance_noslit",
        ignore_nan=True,
        normalize=True,
        normalize_how="mean",
    )
    diff2 = get_residual(
        s,
        s * 1.2,
        var="radiance_noslit",
        ignore_nan=True,
        normalize=True,
        normalize_how="area",
    )
    # TODO Write similar testcase for normalize_how="mean"

    # Test Plot function when there are Nans in the spectrum:
    if plot:
        s.plot(normalize=True)
        s.plot(normalize=(2200, 2250))

    # Test plot_diff function when there are Nans in the spectra:
    if plot:
        plot_diff(s, s * 1.2, normalize=True)
        plot_diff(s, s * 1.2, "radiance_noslit", normalize=(2000, 2100))


def _run_testcases(plot=True, verbose=True, warnings=True, *args, **kwargs):
    """ Test procedures
    """

    # Test all Spectrum compare methods
    # ----------------------------------
    test_compare_methods(verbose=verbose, plot=plot, *args, **kwargs)
    test_plot_compare_with_nan(verbose=verbose, plot=True, *args, **kwargs)
    return True


if __name__ == "__main__":
    print("Test_compare.py: ", _run_testcases())
