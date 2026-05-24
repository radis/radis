# -*- coding: utf-8 -*-
"""

Examples
--------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests

>>> pytest -m fast

"""

from os.path import basename

import matplotlib.pyplot as plt
import numpy as np
import pytest

from radis.misc.printer import printm
from radis.test.utils import setup_test_line_databases

fig_prefix = basename(__file__) + ": "

from radis.misc.utils import NotInstalled, not_installed_vaex_args

try:
    import vaex
except ImportError:
    vaex = NotInstalled(*not_installed_vaex_args)


# %% Test routines
# @pytest.mark.fast #Removed from "fast" tag because it seems this test needs to be performed after other ones (minouHub)
def test_populations(verbose=True, *args, **kwargs):
    """Test that vib and rovib populations are calculated correctly"""

    from radis.lbl import SpectrumFactory
    from radis.misc.basics import all_in

    export = ["vib", "rovib"]

    sf = SpectrumFactory(
        2000,
        2300,
        export_populations=export,
        cutoff=1e-25,
        isotope="1",
    )
    sf.warnings.update(
        {"MissingSelfBroadeningWarning": "ignore", "VoigtBroadeningWarning": "ignore"}
    )
    sf.load_databank("HITRAN-CO-TEST")
    sf.misc.export_rovib_fraction = True
    # we test that "tabulation" and "export_population" are incompatible

    if isinstance(vaex, NotInstalled):
        # Vaex not available, just-in-time partition functions are only implemented in Vaex as of Radis 0.15"
        # SKIPPING TEST
        pass
    else:

        sf.params.parsum_mode = "tabulation"
        with pytest.raises(ValueError) as err:
            s = sf.non_eq_spectrum(2000, 2000)
        assert (
            str(err.value)
            == "Cannot update populations of individual levels with `tabulation` mode. Choose `update_populations=False` or `mode='full summation'`"
        )

    sf.params.parsum_mode = "full summation"  # won't be default at some point
    s = sf.non_eq_spectrum(2000, 2000)

    pops = sf.get_populations(export)
    if not all_in(["rovib", "vib"], list(pops["CO"][1]["X"].keys())):
        raise AssertionError(
            "vib and rovib levels should be defined after non_eq_spectrum calculation!"
        )
    if not "nvib" in list(pops["CO"][1]["X"]["vib"].keys()):
        raise AssertionError(
            "Vib populations should be defined after non_eq_spectrum calculation!"
        )

    s = sf.eq_spectrum(300)

    pops = sf.get_populations(export)
    if "nvib" in list(pops["CO"][1]["X"]["vib"].keys()):
        raise AssertionError(
            "Vib levels should not be defined anymore after eq_spectrum calculation!"
        )

    # Any of these is True and something went wrong
    s = sf.non_eq_spectrum(2000, 2000)
    s2 = sf.non_eq_spectrum(300, 300)

    # printm(all(s2.get_vib_levels(isotope=1) == s.get_vib_levels(isotope=1)))
    assert not (s2.get_vib_levels() is s.get_vib_levels())
    assert not (s2.get_rovib_levels() == s.get_rovib_levels()).all().all()
    assert not (s2.get_rovib_levels() is s.get_rovib_levels())


@pytest.mark.fast
def test_rescaling_path_length(
    debug=False, plot=False, verbose=True, warnings=True, *args, **kwargs
):
    """Test rescaling functions"""

    if plot:  # Make sure matplotlib is interactive so that test are not stuck
        plt.ion()

    from radis.lbl import SpectrumFactory

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    Tgas = 1500
    sf = SpectrumFactory(
        wavelength_min=4400,
        wavelength_max=4800,
        mole_fraction=0.01,
        #                         path_length=0.1,
        cutoff=1e-25,
        wstep=0.005,
        isotope=[1],
        self_absorption=True,
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    #        sf.warnings['NegativeEnergiesWarning'] = 'ignore'
    sf.load_databank("HITRAN-CO-TEST")
    s1 = sf.non_eq_spectrum(Tgas, Tgas, path_length=0.01)
    s2 = sf.non_eq_spectrum(Tgas, Tgas, path_length=3)
    s1.rescale_path_length(3)

    if plot:
        fig = plt.figure(fig_prefix + "Rescaling path length")
        s2.plot("radiance_noslit", nfig=fig.number, lw=3, label="L=3m")
        s1.plot(
            "radiance_noslit",
            nfig=fig.number,
            color="r",
            label="L=0.01m, rescaled to 3m",
        )
        plt.title("Non optically thin rescaling")
        plt.legend()
        plt.tight_layout()

    if verbose:
        printm("Test rescaling:")
        printm(f"... Difference: {abs(s1.get_power() / s2.get_power() - 1) * 100:.2f}%")

    assert np.isclose(s2.get_power(), s1.get_power(), 2e-3)


@pytest.mark.fast
def test_rescaling_mole_fraction(
    debug=False, plot=False, verbose=True, warnings=True, *args, **kwargs
):
    """Test rescaling functions"""

    from radis.lbl import SpectrumFactory

    if plot:  # Make sure matplotlib is interactive so that test are not stuck
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    Tgas = 1500
    sf = SpectrumFactory(
        wavelength_min=4400,
        wavelength_max=4800,
        #                     mole_fraction=1,
        path_length=0.1,
        mole_fraction=0.01,
        cutoff=1e-25,
        wstep=0.005,
        isotope=[1],
        self_absorption=True,
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.load_databank("HITRAN-CO-TEST")
    error = []
    N = [1e-3, 1e-2, 1e-1, 0.3, 0.6, 1]  # first is ref
    for Ni in N:
        s1 = sf.non_eq_spectrum(Tgas, Tgas, mole_fraction=N[0])
        sN = sf.non_eq_spectrum(Tgas, Tgas, mole_fraction=Ni)
        s1.rescale_mole_fraction(Ni)
        error.append(sN.get_power() / s1.get_power())

    if plot:
        plt.figure(fig_prefix + "Rescaling mole fractions")
        plt.plot(N, error, "-ok")
        plt.scatter(
            N[0],
            error[0],
            s=200,
            facecolors="none",
            edgecolors="r",
            label="reference",
        )
        plt.xlabel("Mole fraction")
        plt.ylabel("scaled energy / ab initio energy")
        plt.xscale("log")
        plt.legend()
        plt.title("Effect of scaling mole fraction w/o lineshape update")
        plt.tight_layout()

    # less than 1% error when rescaling from 1e-3 to 0.6
    assert abs(error[-2] - 1) < 0.01


@pytest.mark.fast
def test_medium(plot=False, verbose=True, debug=False, warnings=True, *args, **kwargs):
    """Test effect of propagating medium"""

    from radis.lbl.factory import SpectrumFactory

    if plot:  # Make sure matplotlib is interactive so that test are not stuck
        plt.ion()

    T = 300

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    pl = SpectrumFactory(
        wavenum_min=2171.5,
        wavenum_max=2174,
        mole_fraction=0.01,
        medium="vacuum",
        isotope="1,2",
    )
    pl.warnings["MissingSelfBroadeningWarning"] = "ignore"
    pl.load_databank("HITRAN-CO-TEST")
    s = pl.non_eq_spectrum(Tvib=T, Trot=T)  # , Ttrans=300)

    pla = SpectrumFactory(
        wavenum_min=2171.5,
        wavenum_max=2174,
        mole_fraction=0.01,
        medium="air",
        isotope="1,2",
    )
    pla.load_databank("HITRAN-CO-TEST")
    s_air = pla.non_eq_spectrum(Tvib=T, Trot=T)  # , Ttrans=300)

    if plot:
        plt.figure(fig_prefix + "Propagating medium conversions")
        s.plot(wunit="nm_vac", nfig="same", lw=3, label="vacuum")
        s.plot(wunit="nm", nfig="same", label="air")

    assert np.allclose(s.get_wavenumber(), s_air.get_wavenumber())
    assert np.allclose(
        s.get_wavelength(medium="vacuum"), s_air.get_wavelength(medium="vacuum")
    )
    assert np.allclose(
        s.get_wavelength(medium="air"), s_air.get_wavelength(medium="air")
    )
    assert all(s.get_wavelength(medium="vacuum") > s.get_wavelength(medium="air"))


@pytest.mark.fast
def test_plot_wunit_nm_medium_consistency(verbose=True, *args, **kwargs):
    """Test that plot(wunit='nm') respects the medium from calculation
    and that get_wavenumber_range accepts 'nm_vac'/'nm_air'.

    Regression test for https://github.com/radis/radis/issues/707
    """
    from radis import Spectrum
    from radis.lbl.base import get_wavenumber_range

    # --- Test 1: get_wavenumber_range accepts 'nm_vac' and 'nm_air' ---
    wmin_nm = 3403.0
    wmax_nm = 3405.0

    wn_min_vac, wn_max_vac = get_wavenumber_range(
        wmin=wmin_nm, wmax=wmax_nm, wunit="nm_vac"
    )
    wn_min_vac2, wn_max_vac2 = get_wavenumber_range(
        wmin=wmin_nm, wmax=wmax_nm, wunit="nm", medium="vacuum"
    )
    assert np.isclose(wn_min_vac, wn_min_vac2)
    assert np.isclose(wn_max_vac, wn_max_vac2)

    wn_min_air, wn_max_air = get_wavenumber_range(
        wmin=wmin_nm, wmax=wmax_nm, wunit="nm_air"
    )
    wn_min_air2, wn_max_air2 = get_wavenumber_range(
        wmin=wmin_nm, wmax=wmax_nm, wunit="nm", medium="air"
    )
    assert np.isclose(wn_min_air, wn_min_air2)
    assert np.isclose(wn_max_air, wn_max_air2)

    assert not np.isclose(
        wn_min_vac, wn_min_air
    ), "nm_vac and nm_air should produce different wavenumber conversions"

    # --- Test 2: get(wunit='nm') respects medium from conditions ---
    w_cm = np.linspace(2937.0, 2938.0, 100)
    I = np.random.uniform(0.8, 1.0, len(w_cm))

    s_vac = Spectrum.from_array(
        w_cm,
        I,
        "transmittance_noslit",
        wunit="cm-1",
        Iunit="",
        conditions={"medium": "vacuum"},
    )
    s_air = Spectrum.from_array(
        w_cm,
        I,
        "transmittance_noslit",
        wunit="cm-1",
        Iunit="",
        conditions={"medium": "air"},
    )

    w_nm_from_vac, _ = s_vac.get("transmittance_noslit", wunit="nm")
    w_nm_from_air, _ = s_air.get("transmittance_noslit", wunit="nm")
    w_nm_vac_explicit, _ = s_vac.get("transmittance_noslit", wunit="nm_vac")

    # wunit='nm' on a vacuum spectrum should return vacuum wavelengths
    assert np.allclose(
        w_nm_from_vac, w_nm_vac_explicit
    ), "get(wunit='nm') on a vacuum spectrum should return vacuum wavelengths"

    # wunit='nm' on air vs vacuum spectrum should give different values
    assert not np.allclose(
        w_nm_from_vac, w_nm_from_air
    ), "get(wunit='nm') should give different results for air vs vacuum spectra"

    # Vacuum wavelengths > air wavelengths
    assert all(
        w_nm_from_vac > w_nm_from_air
    ), "Vacuum wavelengths should be larger than air wavelengths"

    if verbose:
        printm("test_plot_wunit_nm_medium_consistency: OK")


def _run_testcases(
    plot=False, verbose=True, debug=False, warnings=True, *args, **kwargs
):
    """Test procedures

    Parameters
    ----------

    debug: boolean
        swamps the console namespace with local variables. Default ``False``

    """

    # Test populations
    # ----------
    test_populations(verbose=verbose, *args, **kwargs)

    # Test conversion of intensity cm-1 works
    # -------------
    #    test_intensity_conversion(debug=debug, verbose=verbose, *args, **kwargs)   # Moved in RADIS

    # Test updating / rescaling functions (no self absorption)
    # ---------
    #    test_rescaling_function(debug=debug, *args, **kwargs)     # Moved in RADIS
    test_rescaling_path_length(
        plot=plot, verbose=verbose, debug=debug, warnings=warnings, *args, **kwargs
    )
    test_rescaling_mole_fraction(
        plot=plot, verbose=verbose, debug=debug, warnings=warnings, *args, **kwargs
    )

    # Test propagating medium
    test_medium(
        plot=plot, verbose=verbose, debug=debug, warnings=warnings, *args, **kwargs
    )


if __name__ == "__main__":
    printm("Test spectrum: ", _run_testcases(debug=False, plot=True))
