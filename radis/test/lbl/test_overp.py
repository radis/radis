# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 08:56:12 2017

@author: erwan

Use
------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)

>>> pytest -m fast


----------


"""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from radis.lbl import SpectrumFactory
from radis.lbl.overp import LevelsList
from radis.misc.printer import printm
from radis.spectrum import plot_diff
from radis.test.utils import setup_test_line_databases


@pytest.mark.needs_config_file
@pytest.mark.needs_db_CDSD_HITEMP_PCN
# @pytest.mark.needs_connection
def test_direct_overpopulation_vs_recombined_bands(
    verbose=True, plot=False, warnings=True, rtol=0.05, *args, **kwargs
):
    """Compare a non-equilibrium spectrum calculated directly with overpopulations,
    or by recombining pre-calculated vibrational bands.

    The later allows for almost instantaneous changes of the overpopulation factors,
    (mostly useful in fitting algorithms), but is only valid for optically thin emission spectra

    Expected output:

    when x_CO2 = 1e-3, radiance in both cases match
    when x_CO2 = 1, they dont

    """

    #    Notes
    #    -----
    #
    #    On the old NeQ package the test used [HITEMP-2010]_
    #
    #    Starting from RADIS 1.0.1, the test is run on [HITRAN-2016]_, which
    #    is not valid for these temperatures but can be more conveniently
    #    downloaded automatically and thus executed everytime with `Travis CI <https://travis-ci.com/radis/radis>`_
    #

    # Note: only with iso1 at the moment

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    # Generate factory
    iso = 1
    sf = SpectrumFactory(
        wavelength_min=4220,
        wavelength_max=4280,
        mole_fraction=1e-3,
        path_length=10,
        cutoff=1e-25,
        molecule="CO2",
        isotope=iso,
        wstep=0.01,
        truncation=2.5,
        medium="air",
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.load_databank("CDSD-HITEMP-PCN")
    #        sf.fetch_databank()   # uses HITRAN: not really valid at this temperature, but runs on all machines without install

    # Generate bands to recombine
    parfunc = sf.parsum_calc["CO2"][iso]["X"]
    Tref = 1500
    # , return_lines=False)
    s_bands = sf.non_eq_bands(Tvib=Tref, Trot=Tref)
    lvlist = LevelsList(parfunc, s_bands, sf.params.levelsfmt)

    # Compare ab initio and recombined from bands at M = 1e-3
    s_recombined = lvlist.non_eq_spectrum(
        Tvib=Tref, Trot=Tref, overpopulation={"(4,1,3)": 3}
    )
    sref = sf.non_eq_spectrum(Tvib=Tref, Trot=Tref, overpopulation={"(4,1,3)": 1})

    if verbose:
        printm(
            "Testing x_CO2 = 1e-3: ab initio ~ recombined bands (<{0:.1f}%):\t".format(
                rtol * 100
            )
        )

    if plot:
        plot_diff(
            sref,
            s_recombined,
            var="radiance_noslit",
            label1="ab initio",
            label2="recombined bands",
            title="x_CO2 = 1e-3",
        )

    assert np.allclose(
        s_recombined.get_radiance_noslit(), sref.get_radiance_noslit(), rtol=rtol
    )

    # Rescale and try again for x_CO2 = 1
    s_recombined.rescale_mole_fraction(1)
    sref.rescale_mole_fraction(1)

    if plot:
        plot_diff(
            sref,
            s_recombined,
            var="radiance_noslit",
            label1="ab initio",
            label2="recombined bands",
            title="x_CO2 = 1",
        )

    if verbose:
        printm(
            "Testing x_CO2 = 1: ab initio ~ recombined bands (<{0:.1f}%):\t{1}".format(
                rtol * 100,
                np.allclose(
                    s_recombined.get_radiance_noslit(),
                    sref.get_radiance_noslit(),
                    rtol=rtol,
                ),
            )
        )

    with pytest.raises(AssertionError):
        assert np.allclose(
            s_recombined.get_radiance_noslit(), sref.get_radiance_noslit(), rtol=rtol
        )

    return True


def test_3Tvib_vs_1Tvib(verbose=True, plot=False, warnings=True, *args, **kwargs):
    """Compare 3-vibrational Temperature algorithm with 1 vibrational temperature
    algorithm, at equilibrium. Expect same output.
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()

    T = 1500

    iso = [1]
    sf = SpectrumFactory(
        wavenum_min=2380,
        wavenum_max=2400,
        pressure=20 * 1e-3,
        cutoff=1e-25,
        isotope=iso,  # ,2',
        path_length=10,
        mole_fraction=0.1,
        truncation=0.5,
        medium="vacuum",
        wstep=0.001,
        verbose=verbose,
    )
    sf.warnings.update(
        {
            "MissingSelfBroadeningWarning": "ignore",
            "VoigtBroadeningWarning": "ignore",
            "HighTemperatureWarning": "ignore",
        }
    )
    sf.load_databank("HITRAN-CO2-TEST", load_energies=True)

    # Compare energies
    for I in iso:
        energies = sf.get_energy_levels("CO2", I, "X")
        assert (energies.Evib == energies.Evib1 + energies.Evib2 + energies.Evib3).all()
    if verbose:
        printm("Tested Evib == Evib1+Evib2+Evib3: OK")

    s3 = sf.non_eq_spectrum((T, T, T), T)
    s1 = sf.non_eq_spectrum(T, T)
    s3.name = "Tvib=({0},{0},{0}) K, Trot={0} K".format(T)
    s1.name = "Tvib={0} K, Trot={0} K".format(T)

    if plot:
        s1.plot(
            "transmittance_noslit",
            wunit="cm-1",
            color="k",
            lw=3,
            label="Tvib = {0}K".format(T),
        )

        s3.plot(
            "transmittance_noslit",
            wunit="cm-1",
            color="r",
            lw=3,
            nfig="same",
            label="Tvib1 = Tvib2 = Tvib3 = {0}K".format(T),
        )

    assert s1.compare_with(s3, spectra_only=True, verbose=verbose, plot=plot)
    if verbose:
        printm("Tested Spectra 3Tvib(T1=T2=T3=T) and 1Tvib(T) are the same: OK")

    return True


def run_testcases(verbose=True, plot=False, warnings=True, *args, **kwargs):

    test_direct_overpopulation_vs_recombined_bands(plot=plot)
    test_3Tvib_vs_1Tvib(plot=plot)

    return True


if __name__ == "__main__":
    printm("test_overp.py:", run_testcases(plot=True))
