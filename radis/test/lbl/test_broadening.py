# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 17:14:01 2018

@author: erwan

Test broadening against HAPI and tabulated data

We're looking at CO(0->1) line 'R1' at 2150.86 cm-1

"""

from os.path import dirname, join

import matplotlib.pyplot as plt
import pytest
from numpy import isclose

from radis import get_residual, get_residual_integral, plot_diff
from radis.lbl.factory import SpectrumFactory
from radis.misc.printer import printm
from radis.spectrum.spectrum import Spectrum
from radis.test.utils import setup_test_line_databases


@pytest.mark.fast
@pytest.mark.needs_connection
def test_broadening_vs_hapi(rtol=1e-2, verbose=True, plot=False, *args, **kwargs):
    """
    Test broadening against HAPI and tabulated data

    We're looking at CO(0->1) line 'R1' at 2150.86 cm-1
    """
    from hapi import absorptionCoefficient_Voigt, db_begin, fetch, tableList

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    # Conditions
    T = 3000
    p = 0.0001
    wstep = 0.001
    wmin = 2150  # cm-1
    wmax = 2152  # cm-1
    truncation = 5  # cm-1
    neighbour_lines = 5  # cm-1

    # %% HITRAN calculation
    # -----------

    # Generate HAPI database locally
    hapi_data_path = join(dirname(__file__), __file__.replace(".py", "_HAPIdata"))

    db_begin(hapi_data_path)
    if not "CO" in tableList():  # only if data not downloaded already
        fetch("CO", 5, 1, wmin - truncation / 2, wmax + truncation / 2)
        # HAPI doesnt correct for side effects

    # Calculate with HAPI
    nu, coef = absorptionCoefficient_Voigt(
        SourceTables="CO",
        Environment={
            "T": T,
            "p": p / 1.01325,
        },  # K  # atm
        WavenumberStep=wstep,
        HITRAN_units=False,
    )

    s_hapi = Spectrum.from_array(
        nu, coef, "abscoeff", "cm-1", "cm-1", conditions={"Tgas": T}, name="HAPI"
    )

    # %% Calculate with RADIS
    # ----------
    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=1,
        path_length=1,  # doesnt change anything
        wstep=wstep,
        pressure=p,
        truncation=truncation,
        neighbour_lines=neighbour_lines,
        isotope=[1],
        warnings={
            "MissingSelfBroadeningWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "HighTemperatureWarning": "ignore",
            "GaussianBroadeningWarning": "ignore",
        },
    )
    sf.load_databank(
        path=join(hapi_data_path, "CO.data"), format="hitran", parfuncfmt="hapi"
    )
    #    s = pl.non_eq_spectrum(Tvib=T, Trot=T, Ttrans=T)
    s = sf.eq_spectrum(Tgas=T, name="RADIS")

    if plot:  # plot broadening of line of largest linestrength
        sf.plot_broadening(i=sf.df1.S.idxmax())

    # Plot and compare
    res = abs(get_residual_integral(s, s_hapi, "abscoeff"))
    if plot:
        plot_diff(
            s,
            s_hapi,
            var="abscoeff",
            title="{0} bar, {1} K (residual {2:.2g}%)".format(p, T, res * 100),
            show_points=False,
        )
        plt.xlim((wmin, wmax))
    if verbose:
        printm("residual:", res)
    assert res < rtol


@pytest.mark.fast
def test_broadening_methods_different_conditions(
    verbose=True, plot=False, *args, **kwargs
):
    """
    Test direct Voigt broadening vs convolution of Gaussian x Lorentzian
    for different spectral grid resolution

    Notes
    -----

    Reference broadening calculated manually with the HWHM formula of
    `HITRAN.org <https://hitran.org/docs/definitions-and-units/>`_
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    # Conditions
    wstep = 0.005
    wmin = 2150.4  # cm-1
    wmax = 2151.4  # cm-1
    truncation = 1  # cm-1

    for (T, p, fwhm_lorentz, fwhm_gauss) in [
        # K, bar, expected FWHM for Lotentz, gauss (cm-1)
        (3000, 1, 0.02849411, 0.01594728),
        (300, 1, 0.16023415, 0.00504297),
        (3000, 0.01, 0.00028494, 0.01594728),
    ]:

        # %% Calculate with RADIS
        # ----------
        sf = SpectrumFactory(
            wavenum_min=wmin,
            wavenum_max=wmax,
            mole_fraction=1,
            path_length=1,  # doesnt change anything
            wstep=wstep,
            pressure=p,
            truncation=truncation,
            isotope="1",
            verbose=False,
            warnings={
                "MissingSelfBroadeningWarning": "ignore",
                "NegativeEnergiesWarning": "ignore",
                "HighTemperatureWarning": "ignore",
                "OutOfRangeLinesWarning": "ignore",
                "GaussianBroadeningWarning": "ignore",
                "CollisionalBroadeningWarning": "ignore",
            },
        )
        sf.load_databank("HITRAN-CO-TEST")
        # Manually filter line database, keep one line only:
        sf.df0.drop(sf.df0[sf.df0.vu != 1].index, inplace=True)
        assert isclose(sf.df0.wav, 2150.856008)

        # Calculate spectra (different broadening methods)
        sf.params.broadening_method = "voigt"
        s_voigt = sf.eq_spectrum(Tgas=T, name="direct")

        # assert broadening FWHM are correct
        assert isclose(2 * float(sf.df1.hwhm_gauss), fwhm_gauss)
        assert isclose(2 * float(sf.df1.hwhm_lorentz), fwhm_lorentz)

        sf.params.broadening_method = "convolve"
        s_convolve = sf.eq_spectrum(Tgas=T, name="convolve")

        # assert broadening FWHM are correct
        assert isclose(2 * float(sf.df1.hwhm_gauss), fwhm_gauss)
        assert isclose(2 * float(sf.df1.hwhm_lorentz), fwhm_lorentz)

        res = get_residual(s_voigt, s_convolve, "abscoeff")

        if verbose:
            print(
                "{0} K, {1} bar: FWHM lorentz = {2:.3f} cm-1, FWHM gauss = {3:.3f} cm-1".format(
                    T, p, 2 * float(sf.df1.hwhm_lorentz), 2 * float(sf.df1.hwhm_gauss)
                )
            )

        if plot:
            plot_diff(
                s_voigt,
                s_convolve,
                "abscoeff",
                title=r"T {0} K, p {1} bar: w$_\mathrm{{L}}$ {2:.3f}, w$_\mathrm{{G}}$ {3:.3f} cm$^{{-1}}$".format(
                    T, p, 2 * float(sf.df1.hwhm_lorentz), float(sf.df1.hwhm_gauss)
                ),
            )

        # assert all broadening methods match
        assert res < 2e-4


@pytest.mark.fast
def test_broadening_methods_different_wstep(verbose=True, plot=False, *args, **kwargs):
    """
    Test direct Voigt broadening vs convolution of Gaussian x Lorentzian
    for different spectral grid resolution
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    # Conditions
    T = 3000
    p = 1
    wmin = 2150  # cm-1
    wmax = 2152  # cm-1
    truncation = 5  # cm-1
    neighbour_lines = 5  # cm-1

    for i, wstep in enumerate([0.01, 0.1, 0.5]):

        # %% Calculate with RADIS
        # ----------
        sf = SpectrumFactory(
            wavenum_min=wmin,
            wavenum_max=wmax,
            mole_fraction=1,
            path_length=1,  # doesnt change anything
            wstep=wstep,
            pressure=p,
            truncation=truncation,
            neighbour_lines=neighbour_lines,
            isotope="1",
            optimization=None,
            verbose=False,
            warnings={
                "MissingSelfBroadeningWarning": "ignore",
                "NegativeEnergiesWarning": "ignore",
                "HighTemperatureWarning": "ignore",
                "GaussianBroadeningWarning": "ignore",
                "AccuracyError": "ignore",
                "AccuracyWarning": "ignore",
            },
        )  # 0.2)
        sf.load_databank("HITRAN-CO-TEST")
        #    s = pl.non_eq_spectrum(Tvib=T, Trot=T, Ttrans=T)
        sf.params.broadening_method = "voigt"
        s_voigt = sf.eq_spectrum(Tgas=T, name="direct")

        sf.params.broadening_method = "convolve"
        s_convolve = sf.eq_spectrum(Tgas=T, name="convolve")

        res = get_residual(s_voigt, s_convolve, "abscoeff")

        if verbose:
            print("Residual:", res)

        # plot the last one
        if plot:
            plot_diff(
                s_voigt,
                s_convolve,
                "abscoeff",
                nfig="test_voigt_broadening_methods" + str(i),
                title="P {0} bar, T {1} K, wstep {2} cm-1".format(p, T, wstep),
            )

        assert res < 2e-4


@pytest.mark.fast
def test_broadening_LDM(verbose=True, plot=False, *args, **kwargs):
    """
    Test use of lineshape template for broadening calculation.

    Ensures that results are the same with and without LDM.
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    # Conditions
    T = 3000
    p = 1
    wstep = 0.002
    wmin = 2150  # cm-1
    wmax = 2152  # cm-1
    truncation = 5  # cm-1
    neighbour_lines = 5  # cm-1

    # %% Calculate with RADIS
    # ----------
    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=1,
        path_length=1,  # doesnt change anything
        wstep=wstep,
        pressure=p,
        truncation=truncation,
        neighbour_lines=neighbour_lines,
        optimization=None,
        isotope="1",
        verbose=False,
        warnings={
            "MissingSelfBroadeningWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "HighTemperatureWarning": "ignore",
            "GaussianBroadeningWarning": "ignore",
        },
    )
    sf.load_databank("HITRAN-CO-TEST")

    # Reference: calculate without LDM
    assert sf.params["optimization"] is None
    s_ref = sf.eq_spectrum(Tgas=T)
    s_ref.name = "Reference ({0:.2f}s)".format(s_ref.conditions["calculation_time"])

    # LDM:
    sf.params["optimization"] = "simple"
    sf.params.broadening_method = "convolve"
    s_ldm = sf.eq_spectrum(Tgas=T)
    s_ldm.name = "LDM ({0:.2f}s)".format(s_ldm.conditions["calculation_time"])
    # LDM Voigt with Whiting approximation:
    sf.params.broadening_method = "voigt"
    s_ldm_voigt = sf.eq_spectrum(Tgas=T)
    s_ldm_voigt.name = "LDM Whiting ({0:.2f}s)".format(
        s_ldm_voigt.conditions["calculation_time"]
    )

    # Compare
    res = get_residual(s_ref, s_ldm, "abscoeff")
    res_voigt = get_residual(s_ldm, s_ldm_voigt, "abscoeff")

    if verbose:
        print("Residual:", res)

    # plot the last one
    if plot:
        plot_diff(s_ref, s_ldm, "abscoeff")
        plt.legend()
        plot_diff(s_ldm, s_ldm_voigt, "abscoeff")

    assert res < 1.2e-5
    assert res_voigt < 1e-5


@pytest.mark.fast
def test_broadening_LDM_FT(verbose=True, plot=False, *args, **kwargs):
    """
    Test use of LDM with and without Fourier Transform

    Ensures that results are the same, and compare calculation times.
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    # Conditions
    T = 3000
    p = 1
    wstep = 0.002
    wmin = 2000  # cm-1
    wmax = 2300  # cm-1
    truncation = 10  # cm-1
    neighbour_lines = 10  # cm-1

    # %% Calculate with RADIS
    # ----------
    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=1,
        path_length=1,  # doesnt change anything
        wstep=wstep,
        pressure=p,
        truncation=truncation,
        neighbour_lines=neighbour_lines,
        isotope="1",
        verbose=verbose,
        optimization="simple",
        warnings={
            "MissingSelfBroadeningWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "HighTemperatureWarning": "ignore",
            "GaussianBroadeningWarning": "ignore",
        },
    )
    sf.load_databank("HITRAN-CO-TEST")

    # LDM, real space
    if verbose:
        print("\nConvolve version \n")
    sf._broadening_method = "convolve"
    s_ldm = sf.eq_spectrum(Tgas=T)
    s_ldm.name = "LDM ({0:.2f}s)".format(s_ldm.conditions["calculation_time"])

    # LDM , with Fourier
    if verbose:
        print("\nFFT version \n")
    sf.params.broadening_method = "fft"
    s_ldm_fft = sf.eq_spectrum(Tgas=T)
    s_ldm_fft.name = "LDM FFT ({0:.2f}s)".format(
        s_ldm_fft.conditions["calculation_time"]
    )

    # Compare
    res = get_residual(s_ldm, s_ldm_fft, "abscoeff")

    if verbose:
        print("Residual:", res)

    # plot
    if plot:
        plot_diff(s_ldm, s_ldm_fft, "abscoeff")
        plt.legend()

    assert res < 5e-6


@pytest.mark.fast
def test_broadening_LDM_noneq(verbose=True, plot=False, *args, **kwargs):
    """
    Test Noneq version of LDM and makes sure it gives the same results as the eq
    one when used with Tvib=Trot

    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    setup_test_line_databases()  # add HITRAN-CO2-TEST in ~/radis.json if not there

    # Conditions
    p = 1
    wstep = 0.002
    wmin = 2380  # cm-1
    wmax = 2400  # cm-1
    truncation = 5  # cm-1
    neighbour_lines = 5  # cm-1

    # %% Calculate with RADIS
    # ----------
    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=1,
        path_length=1,  # doesnt change anything
        wstep=wstep,
        pressure=p,
        truncation=truncation,
        neighbour_lines=neighbour_lines,
        isotope="1",
        verbose=3,
        warnings={
            "MissingSelfBroadeningWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "HighTemperatureWarning": "ignore",
            "GaussianBroadeningWarning": "ignore",
        },
    )
    sf.load_databank("HITRAN-CO2-TEST", load_columns="noneq")

    # LDM:
    sf.params["optimization"] = "simple"
    s_ldm_eq = sf.eq_spectrum(Tgas=3000)
    s_ldm_eq.name = "LDM eq ({0:.2f}s)".format(s_ldm_eq.conditions["calculation_time"])

    s_ldm_noneq = sf.non_eq_spectrum(Tvib=3000, Trot=3000)
    s_ldm_noneq.name = "LDM noneq ({0:.2f}s)".format(
        s_ldm_noneq.conditions["calculation_time"]
    )

    # Compare
    res = get_residual(s_ldm_eq, s_ldm_noneq, "radiance_noslit")

    if verbose:
        print("Residual:", res)

    # plot
    if plot:
        plot_diff(s_ldm_eq, s_ldm_noneq)

    assert res <= 1e-4


def test_truncations_and_neighbour_lines(*args, **kwargs):
    """Test new truncations introduced in https://github.com/radis/radis/issues/340

    So far: test that all functions work
    More assertion could be added.
    """
    # TODO  . Check on databases with one line that truncatino is achieved, etc.

    #%% Test truncation > lines

    conditions = {
        "wavenum_min": 2100,
        "wavenum_max": 2300,
        "molecule": "CO",
        "isotope": "1,2,3",
        "pressure": 20,  # bar
        "Tgas": 700,
        "mole_fraction": 0.1,
        "path_length": 1,
        "databank": "HITRAN-CO-TEST",
    }

    # No Neighbourling lines, only truncation
    from radis import calc_spectrum

    s_lbl_voigt_trunc10 = calc_spectrum(
        **conditions,
        optimization=None,
        truncation=10,
        neighbour_lines=0,
    )
    assert s_lbl_voigt_trunc10.conditions["truncation"] == 10

    #%% Test no truncation

    from radis import calc_spectrum

    s_lbl_voigt_notrunc = calc_spectrum(
        **conditions,
        optimization=None,
        truncation=None,
    )
    assert s_lbl_voigt_notrunc.conditions["truncation"] is None

    # Effect of neighbour lines on both edges of the spectrum:
    assert (
        s_lbl_voigt_notrunc.get("abscoeff")[1][0]
        > s_lbl_voigt_trunc10.get("abscoeff")[1][0]
    )
    assert (
        s_lbl_voigt_notrunc.get("abscoeff")[1][-1]
        > s_lbl_voigt_trunc10.get("abscoeff")[1][-1]
    )

    #%%

    # Check incompatibilities correctly raise errors:

    import pytest

    with pytest.raises(NotImplementedError) as err:

        calc_spectrum(
            **conditions,
            optimization=None,
            broadening_method="fft",  # = "fft"
            truncation=10,  # truncation != 0
        )
    assert "Lines cannot be truncated with `broadening_method='fft'`" in str(err.value)

    #%% same with optimization=simple (DIT)
    with pytest.raises(NotImplementedError) as err:
        calc_spectrum(
            **conditions,
            optimization="simple",
            broadening_method="fft",  # = "fft"
            truncation=10,  # truncation != 0
        )
    assert "Lines cannot be truncated with `broadening_method='fft'`" in str(err.value)

    #%% Test there is DeprecationError raised
    with pytest.raises(DeprecationWarning) as err:
        calc_spectrum(
            **conditions,
            broadening_max_width=20,  # = "fft"
        )
    assert (
        "To keep the current behavior, replace `broadening_max_width=20` with `truncation=10.0"
        in str(err.value)
    )

    #%% Truncation, no neighbour lines : should be ~ same result that without DIT Optimization
    s_dit_voigt_trunc10 = calc_spectrum(
        **conditions,
        optimization="simple",
        broadening_method="voigt",  # = "fft"
        truncation=10,  # truncation != 0
        neighbour_lines=0,
    )

    assert get_residual(s_lbl_voigt_trunc10, s_dit_voigt_trunc10, "abscoeff") < 2e-5

    # %% LBL Optimization; with neighbour lines
    # Also check it works with floating point division

    s_lbl_voigt_trunc10 = calc_spectrum(
        **conditions,
        optimization=None,
        broadening_method="voigt",
        wstep=0.31,  # floating point division can create issues. We test here.
        truncation=300,
        neighbour_lines=15,
    )
    assert s_lbl_voigt_trunc10.conditions["truncation"] == 300


@pytest.mark.fast
def test_broadening_warnings(*args, **kwargs):
    """Test AccuracyWarning and AccuracyErrors are properly triggered.

    Inspired by https://github.com/radis/radis/issues/186

    Test :py:meth:`~radis.lbl.broadening.BroadenFactory._check_accuracy`

    Examples
    --------

    ::
            AccuracyError: Some lines are too narrow (FWHM ~ 0.0011 cm⁻¹) for
            the current spectral grid (wstep=0.01). Please reduce wstep to (at least)
            below 0.00055 cm⁻¹ or (suggested) 0.00022 cm⁻¹
    """
    import astropy.units as u

    from radis.misc.warning import AccuracyError, AccuracyWarning

    conditions = {
        # "wavenum_min": 667.58 / u.cm,
        # "wavenum_max": 667.61 / u.cm,
        "wavelength_min": 4170 * u.nm,
        "wavelength_max": 4180 * u.nm,
        "molecule": "CO2",
        "isotope": "1",
        "neighbour_lines": 0,
        "verbose": False,
    }

    # Try with low resolution, expect error :
    with pytest.raises(AccuracyError):
        sf = SpectrumFactory(**conditions, wstep=0.02)
        # sf.fetch_databank("hitran")
        sf.load_databank("HITRAN-CO2-TEST")

        sf.eq_spectrum(
            Tgas=234.5,
            pressure=0.00646122 * u.bar,
            mole_fraction=400e-6,
        )

    with pytest.warns(AccuracyWarning):
        sf = SpectrumFactory(**conditions, wstep=0.002)
        # sf.fetch_databank("hitran")
        sf.load_databank("HITRAN-CO2-TEST")

        sf.eq_spectrum(
            Tgas=234.5,
            pressure=0.00646122 * u.bar,
            mole_fraction=400e-6,
        )


# @pytest.mark.needs_config_file
# @pytest.mark.needs_db_HITEMP_CO2_DUNHAM
@pytest.mark.needs_connection
def test_abscoeff_continuum(
    plot=False, verbose=2, warnings=True, threshold=0.1, *args, **kwargs
):
    """
    Test calculation with pseudo-continuum

    Assert results on abscoeff dont change


    Notes
    -----

    Uses HITRAN so it can deployed and tested on `Travis CI <https://travis-ci.com/radis/radis>`_, but we should switch
    to HITEMP if some HITEMP files can be downloaded automatically at the
    execution time.

    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    if verbose:
        printm(">>> test_abscoeff_continuum")

    sf = SpectrumFactory(
        wavelength_min=4200,
        wavelength_max=4500,
        cutoff=1e-23,
        molecule="CO2",
        isotope="1,2",
        truncation=5,
        path_length=0.1,
        mole_fraction=1e-3,
        medium="vacuum",
        optimization=None,
        verbose=verbose,
    )
    sf.warnings.update(
        {
            "MissingSelfBroadeningWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "LinestrengthCutoffWarning": "ignore",
            "HighTemperatureWarning": "ignore",
        }
    )
    sf.fetch_databank(
        "hitran"
    )  # uses HITRAN: not really valid at this temperature, but runs on all machines without install
    #        sf.load_databank('HITEMP-CO2-DUNHAM')       # to take a real advantage of abscoeff continuum, should calculate with HITEMP
    sf._export_continuum = True  # activate it

    # Calculate one without pseudo-continuum
    sf.params.pseudo_continuum_threshold = 0
    s1 = sf.eq_spectrum(Tgas=2000)
    s1.name = "All lines resolved ({0}) ({1:.1f}s)".format(
        s1.conditions["lines_calculated"], s1.conditions["calculation_time"]
    )
    assert s1.conditions["pseudo_continuum_threshold"] == 0

    # Calculate one with pseudo-continuum
    sf.params.pseudo_continuum_threshold = threshold
    s2 = sf.eq_spectrum(Tgas=2000)
    s2.name = "Semi-continuum + {0} lines ({1:.1f}s)".format(
        s2.conditions["lines_calculated"], s2.conditions["calculation_time"]
    )
    assert s2.conditions["pseudo_continuum_threshold"] == threshold
    assert "abscoeff_continuum" in s2.get_vars()

    # Plot
    if plot:
        plot_diff(
            s1,
            s2,
            "radiance_noslit",
            Iunit="µW/cm2/sr/nm",
            nfig="test_abscoeff_continuum: diff",
        )

        s2.plot(
            "abscoeff",
            label="Full spectrum",
            nfig="test_abscoeff_continuum: show continuum",
            force=True,
        )
        s2.plot(
            "abscoeff_continuum",
            nfig="same",
            label="Pseudo-continuum ({0} lines in continuum))".format(
                s2.conditions["lines_in_continuum"]
            ),
            force=True,
        )
        plt.legend()

    # Compare
    res = get_residual(s1, s2, "abscoeff")
    if verbose:
        printm("residual:", res)

    globals().update(locals())

    assert res < 1.32e-6


# @pytest.mark.needs_config_file
# @pytest.mark.needs_db_HITEMP_CO2_DUNHAM
@pytest.mark.needs_connection
def test_noneq_continuum(plot=False, verbose=2, warnings=True, *args, **kwargs):
    """
    Test calculation with pseudo-continuum under nonequilibrium

    Assert results on emisscoeff dont change


    Notes
    -----

    Uses HITRAN so it can deployed and tested on `Travis CI <https://travis-ci.com/radis/radis>`_, but we should switch
    to HITEMP if some HITEMP files can be downloaded automatically at the
    execution time.

    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    if verbose:
        printm(">>> test_noneq_continuum")

    sf = SpectrumFactory(
        wavelength_min=4200,
        wavelength_max=4500,
        cutoff=1e-23,
        molecule="CO2",
        isotope="1,2",
        truncation=5,
        neighbour_lines=10,
        path_length=0.1,
        mole_fraction=1e-3,
        medium="vacuum",
        optimization=None,
        verbose=verbose,
    )
    sf.warnings.update(
        {
            "MissingSelfBroadeningWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "LinestrengthCutoffWarning": "ignore",
            "HighTemperatureWarning": "ignore",
        }
    )
    sf.fetch_databank(
        "hitran", load_columns="noneq"
    )  # uses HITRAN: not really valid at this temperature, but runs on all machines without install
    sf._export_continuum = True  # activate it

    # Calculate one without pseudo-continuum
    sf.params.pseudo_continuum_threshold = 0
    s1 = sf.non_eq_spectrum(Tvib=2000, Trot=1000)
    s1.name = "All lines resolved ({0}) ({1:.1f}s)".format(
        s1.conditions["lines_calculated"], s1.conditions["calculation_time"]
    )
    assert s1.conditions["pseudo_continuum_threshold"] == 0

    # Calculate one with pseudo-continuum
    sf.params.pseudo_continuum_threshold = 0.05
    s2 = sf.non_eq_spectrum(Tvib=2000, Trot=1000)
    s2.name = "Semi-continuum + {0} lines ({1:.1f}s)".format(
        s2.conditions["lines_calculated"], s2.conditions["calculation_time"]
    )
    assert s2.conditions["pseudo_continuum_threshold"] == 0.05
    assert "abscoeff_continuum" in s2.get_vars()
    assert "emisscoeff_continuum" in s2.get_vars()

    # Plot
    if plot:
        plot_diff(
            s1,
            s2,
            "radiance_noslit",
            Iunit="µW/cm2/sr/nm",
            nfig="test_noneq_continuum: diff",
        )

        plt.figure("test_noneq_continuum: show continuum").clear()
        s2.plot(
            "emisscoeff", label=s2.name, nfig="test_noneq_continuum: show continuum"
        )
        s2.plot(
            "emisscoeff_continuum",
            nfig="same",
            label="Pseudo-continuum (aggreg. {0:g} lines)".format(
                s2.conditions["lines_in_continuum"]
            ),
            force=True,
        )

    # Compare
    res = get_residual(s1, s2, "abscoeff") + get_residual(s1, s2, "emisscoeff")
    if verbose:
        printm("residual:", res)

    assert res < 5.2e-6


def _run_testcases(plot=False, verbose=True, *args, **kwargs):

    # Test broadening
    test_broadening_vs_hapi(plot=plot, verbose=verbose, *args, **kwargs)
    test_broadening_methods_different_conditions(
        plot=plot, verbose=verbose, *args, **kwargs
    )
    test_broadening_methods_different_wstep(plot=plot, verbose=verbose, *args, **kwargs)
    test_broadening_LDM(plot=plot, verbose=verbose, *args, **kwargs)
    test_broadening_LDM_FT(plot=plot, verbose=3, *args, **kwargs)
    test_broadening_LDM_noneq(plot=plot, verbose=verbose, *args, **kwargs)
    test_truncations_and_neighbour_lines(*args, **kwargs)

    # Test warnings
    test_broadening_warnings(*args, **kwargs)

    # Test pseudo-continuum
    test_abscoeff_continuum(plot=plot, verbose=verbose, *args, **kwargs)
    test_noneq_continuum(plot=plot, verbose=verbose, *args, **kwargs)

    return True


if __name__ == "__main__":
    printm("test_broadening: ", _run_testcases(plot=True, verbose=True, debug=False))
