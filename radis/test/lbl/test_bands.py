# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 08:56:12 2017

@author: erwan

Examples
--------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)

>>> pytest -m fast


----------


"""

import matplotlib.pyplot as plt
import pytest

from radis import get_residual
from radis.lbl import SpectrumFactory
from radis.los import MergeSlabs
from radis.misc.printer import printm


@pytest.mark.needs_connection
def test_plot_all_CO2_bandheads(verbose=True, plot=False, *args, **kwargs):
    """In this test we use the :meth:`~radis.lbl.bands.BandFactory.non_eq_bands`
    method to calculate separately all vibrational bands of CO2, and compare
    them with the final Spectrum.

    """

    # Note: only with iso1 at the moment
    # ('hitran_co2_626_bandhead_4165_4200nm' has only iso=1)

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    Tgas = 1000

    sf = SpectrumFactory(
        wavelength_min=4165,
        wavelength_max=4200,
        mole_fraction=1,
        path_length=0.3,
        cutoff=1e-23,
        molecule="CO2",
        truncation=5,
        neighbour_lines=5,
        isotope="1",
        optimization=None,
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.warnings["HighTemperatureWarning"] = "ignore"
    sf.load_databank("HITRAN-CO2-TEST")

    s_tot = sf.non_eq_spectrum(Tvib=Tgas, Trot=Tgas)
    s_bands = sf.non_eq_bands(Tvib=Tgas, Trot=Tgas)

    if verbose:
        printm("{0} bands in spectrum".format(len(s_bands)))

    assert len(s_bands) == 3

    # Ensure that recombining gives the same
    s_merged = MergeSlabs(*list(s_bands.values()))
    assert s_tot.compare_with(s_merged, "radiance_noslit", plot=False)

    if verbose:
        printm("Recombining bands give the same Spectrum")

    # %%
    if plot:
        s_tot.apply_slit(1, "nm")
        s_tot.name = "Full spectrum"
        s_tot.plot(wunit="nm", lw=3)
        for band, s in s_bands.items():
            s.plot(wunit="nm", nfig="same")
        plt.legend(loc="upper left")
        plt.ylim(ymax=0.25)

    # %% Compare with equilibrium bands now
    s_bands_eq = sf.eq_bands(Tgas)
    s_merged_eq = MergeSlabs(*list(s_bands_eq.values()))

    assert get_residual(s_tot, s_merged_eq, "radiance_noslit") < 1.5e-5

    return True


def run_testcases(verbose=True, plot=False, warnings=True, *args, **kwargs):

    test_plot_all_CO2_bandheads(plot=plot)

    return True


if __name__ == "__main__":

    printm("test_overp.py:", run_testcases(plot=True))
