# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 09:59:55 2017

@author: erwan

Examples
--------

Run all tests::

    pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)::

    pytest -m fast

------------------------------------------------------------------------

"""

from os.path import basename

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pytest

import radis
from radis import config
from radis.lbl.factory import SpectrumFactory
from radis.misc.printer import printm
from radis.test.utils import setup_test_line_databases

fig_prefix = basename(__file__) + ": "

from radis.misc.utils import NotInstalled, not_installed_vaex_args

try:
    import vaex
except ImportError:
    vaex = NotInstalled(*not_installed_vaex_args)

# %% ======================================================================
# TEST
# ---------------------


@pytest.mark.needs_config_file
@pytest.mark.needs_db_CDSD_HITEMP
def test_spec_generation(
    plot=True,
    verbose=2,
    warnings=True,
    update_reference_spectrum=False,
    *args,
    **kwargs,
):
    """Test spectrum generation
    Can be used as a base to generate spectra in your codes

    Non-regression test: compare with past version (see conditions below)

    Compare results from a reference case to results calculated on 30/12/2017
    This is not a validation case (30/12/2017 results are not a physically validated
    case), but it makes sure results dont change over time

    Conditions (30/12/2017)::

        Physical Conditions
        ----------------------------------------
           Tgas                 300 K
           Trot                 300 K
           Tvib                 300 K
           pressure             1.01325 bar
           isotope              1,2
           mole_fraction        1
           molecule             CO2
           path_length          1 cm
           wavelength_max       4400.0 nm
           wavelength_min       4150.0 nm
           wavenum_max          2409.6385542168673 cm-1
           wavenum_min          2272.7272727272725 cm-1
        Computation Parameters
        ----------------------------------------
           Tref                 296 K
           broadening_max_width  10 cm-1
           cutoff               1e-25 cm-1/(#.cm-2)
           db_assumed_sorted    True
           db_use_cached        True
           dbformat             cdsd
           dbpath               # USER-DEPENDANT: CDSD-HITEMP
           fillmissinglevelswithzero  False
           levelsfmt            cdsd
           levelspath           # USER-DEPENDANT: CDSD-4000
           medium               vacuum
           parfuncfmt           cdsd
           parfuncpath          # USER-DEPENDANT: CDSD-4000
           rot_distribution     boltzmann
           self_absorption      True
           vib_distribution     boltzmann
           wavenum_max_calc     2414.6385542168673 cm-1
           wavenum_min_calc     2267.7272727272725 cm-1
           waveunit             cm-1
           wstep                0.01 cm-1
        ----------------------------------------

    Notes
    -----

    Performance test. How long it took to calculate this Spectrum?
    Test with cutoff 1e-25, broadening_max_width=10

    - 0.9.15: >>> 33s

    - 0.9.16*: (replaced groupby().apply() with iteration over indexes) >>> 32s
            [but large impact expected on big files]

    - 0.9.16*: (upgraded cache files to h5) >>> 25s

    - 0.9.16*: (also added h5 cache file for levels) >>> 21s

    - 0.9.16*: (with Whiting slit voigt function) >>> 5.8s

    Test with cutoff 1e-27, broadening_max_width=50 :
    ("Spectrum calculated in ... ", including database loading time)

    - 0.9.16*: (same code as last) >>> 12.5s including 7.6s of broadening

    - 0.9.16**: (with pseudo_continuum_threshold=0.01) >>> 7.8s including 2.3s of broadening

    - 0.9.18 (normal code, no pseudo continuum). >>> ?

    - 0.9.21 (normal code) >>> 13.7s, including 8.7s of broadening
             (with pseudo_continuum_threshold=0.01) >>> 4.3s, including 2.6s of broadening

    - 0.9.21*              >>> 14.0s  (added the manual lineshape normalization instead of
                                       Whitings's polynomial)

    - 0.9.22 (normal code) >>> 11.3s   (without energy level lookup, for eq. calculations)
             (with pseudo_continuum_threshold=0.01) >>> 5.9s

    - 0.9.23 (normal code) >>> 7.2s   (added jit in Voigt broadening)
                           >>> 7.1s   (chunksize = None)  (and much faster for more lines)
             (with pseudo_continuum_threshold=0.01) >>> 4.9s

    RADIS:

    - 0.9.19 (normal code) >>> 6.3 s

    - 0.9.20 (normal code) >>> 6.3 s
             (with pseudo_continuum_threshold=0.01) >>> ???
             (with LDM) >>> 2.3 s

    - 0.9.26 (normal code) >>> 7.6 s
             (with pseudo_continuum_threshold=0.01) >>> 2.73s
             (with LDM) >>> 0.25 s

    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    from time import time

    t0 = time()
    if verbose:
        printm(">>> _test_spec_generation")

    # This is how you get a spectrum (see calc.py for front-end functions
    # that do just that)
    sf = SpectrumFactory(
        wavelength_min=4150,
        wavelength_max=4400,
        cutoff=1e-27,
        isotope="1,2",
        truncation=25,
        optimization=None,
        # optimization="min-RMS",
        # pseudo_continuum_threshold=0.01,
        medium="vacuum",
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.load_databank("HITEMP-CO2-DUNHAM")
    s = sf.eq_spectrum(Tgas=300, name="test_spec_generation")
    if verbose:
        printm(f">>> _test_spec_generation: Spectrum calculated in {time() - t0:.2f}s")

    if plot:
        plt.figure(fig_prefix + "Reference spectrum CDSD-HITEMP (radiance)")
        # Iunit is arbitrary. Use whatever makes sense
        s.plot("radiance_noslit", Iunit="µW/cm2/sr/nm", nfig="same")
    s.rescale_path_length(0.01)

    # Here we get some extra informations:
    if plot:
        sf.plot_broadening(i=0)  # show broadening of one line
        plt.xlim((2267.20, 2268.30))

    # Compare with hardcoded results
    # ... code previously used to export hardcoded results:
    # ... and header contains all input conditions:
    #        np.savetxt('output.txt', np.vstack(s.get('abscoeff', wunit='nm')).T[::10])
    #        print(s)
    # ................
    from radis import get_version
    from radis.test.utils import getTestFile

    wref, Iref = np.loadtxt(getTestFile("CO2abscoeff_300K_4150_4400nm.txt")).T
    match_reference = np.allclose(s.get("abscoeff", wunit="nm")[1][::10], Iref)
    if not match_reference:
        # give some more information before raising error
        printm(
            f"Error: {np.mean(abs(s.get('abscoeff', wunit='nm')[1][::10] / Iref - 1)) * 100:.2f}%"
        )
        # Store the faulty spectrum
        s.store(
            f"test_factory_failed_{radis.get_version()}.spec",
            if_exists_then="replace",
        )

    # Use "update_reference_spectrum=True" to update reference case :
    if update_reference_spectrum:
        wsave, Isave = s.get("abscoeff", wunit="nm")
        import io
        from contextlib import redirect_stdout

        with io.StringIO() as buf, redirect_stdout(buf):
            print(s)
            s_details = buf.getvalue()
        np.savetxt(
            getTestFile("CO2abscoeff_300K_4150_4400nm.txt"),
            np.vstack((wsave[::10], Isave[::10])).T,
            header=f"RADIS {get_version(add_git_number=False)}\n\n{s_details}",
        )

    # Plot comparison
    if plot:
        plt.figure(fig_prefix + "Reference spectrum (abscoeff)")
        # , show_points=True)  # show_points to have an
        s.plot(
            "abscoeff",
            wunit="nm",
            nfig="same",
            lw=3,
            label="RADIS, this version",
        )
        # idea of the resolution
        plt.plot(wref, Iref, "or", ms=3, label="version RADIS 0.9.26 (13/12/20)")
        plt.legend()
        plt.title(f"All close: {match_reference}")
        plt.tight_layout()

    # Another example, at higher temperature.
    # Removed because no test is associated with it and it takes time for
    # nothing
    #        s2 = sf.non_eq_spectrum(Tvib=1000, Trot=300)
    #        if plot: s2.plot('abscoeff', wunit='nm')

    if verbose:
        printm(
            f"Spectrum calculation (no database loading) took {s.conditions['calculation_time']:.1f}s\n"
        )
        printm(f"_test_spec_generation finished in {time() - t0:.1f}s\n")

    assert match_reference


# Test power integral


@pytest.mark.fast
def test_power_integral(verbose=True, warnings=True, *args, **kwargs):
    """Test direct calculation of power integral from Einstein coefficients
    matches integration of broadened spectrum in the optically thin case

    We compare:

    - direct calculation of power integral with equilibrium code
        :meth:`~radis.lbl.SpectrumFactory.optically_thin_power` (T)
    - direct calculation of power integral with non equilibrium code
        :meth:`~radis.lbl.SpectrumFactory.optically_thin_power` (T,T)
    - numerical integration of non equilibrium spectrum under optically thin conditions:
        :meth:`~radis.spectrum.spectrum.Spectrum.get_power`

    Test passes if error < 0.5%

    """

    if verbose:
        printm(">>> _test_power_integral")

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    sf = SpectrumFactory(
        wavelength_min=4300,
        wavelength_max=4666,
        wstep=0.001,
        cutoff=1e-30,
        path_length=10,
        mole_fraction=400e-6,
        isotope=[1],
        truncation=5,
        verbose=verbose,
    )
    sf.warnings.update(
        {
            "MissingSelfBroadeningWarning": "ignore",
            "OutOfRangeLinesWarning": "ignore",
            "HighTemperatureWarning": "ignore",
        }
    )
    sf.load_databank("HITRAN-CO-TEST", db_use_cached=True)
    unit = "µW/sr/cm2"
    T = 600

    # Calculate:

    # ... direct calculation of power integral with equilibrium code
    Peq = sf.optically_thin_power(Tgas=T, unit=unit)

    # ... direct calculation of power integral with non equilibrium code
    Pneq = sf.optically_thin_power(Tvib=T, Trot=T, unit=unit)

    # ... numerical integration of non equilibrium spectrum under optically thin
    # ... conditions
    sf.input.self_absorption = False
    s = sf.non_eq_spectrum(T, T)

    assert s.conditions["self_absorption"] == False

    # Compare
    err = abs(Peq - s.get_power(unit=unit)) / Peq
    if verbose:
        printm(f"Emission integral:\t{Peq:.4g}", unit)
        printm(f"Emission (noneq code):\t{Pneq:.4g}", unit)
        printm(f"Integrated spectrum:\t{s.get_power(unit=unit):.4g}", unit)
        printm(f"Error: {err * 100:.2f}%")

    assert err < 0.005


@pytest.mark.fast
def test_media_line_shift(plot=False, verbose=True, warnings=True, *args, **kwargs):
    """See wavelength difference in air and vacuum"""

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    if verbose:
        printm(">>> _test_media_line_shift")

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    sf = SpectrumFactory(
        wavelength_min=4500,
        wavelength_max=4600,
        wstep=0.001,
        cutoff=1e-30,
        path_length=0.1,
        mole_fraction=400e-6,
        isotope=[1],
        medium="vacuum",
        truncation=5,
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["GaussianBroadeningWarning"] = "ignore"
    sf.load_databank("HITRAN-CO-TEST")

    # Calculate a spectrum
    s = sf.eq_spectrum(2000)

    # Compare
    if plot:
        fig = plt.figure(fig_prefix + "Propagating media line shift")
        s.plot("radiance_noslit", wunit="nm_vac", nfig=fig.number, lw=2, label="Vacuum")
        plt.title("CO spectrum (2000 K)")
        s.plot(
            "radiance_noslit",
            wunit="nm",
            nfig=fig.number,
            lw=2,
            color="r",
            label="Air",
        )

    # ... there should be about ~1.25 nm shift at 4.5 µm:
    assert np.isclose(
        s.get("radiance_noslit", wunit="nm_vac")[0][0]
        - s.get("radiance_noslit", wunit="nm")[0][0],
        1.2540436086346745,
    )


@pytest.mark.fast
@pytest.mark.parametrize(
    ("input_wavelengths", "expected_wavelengths_nm"),
    [
        [(4300 * u.nm, 4.5 * u.um), (4300, 4500)],
        [(4500, 5000), (4500, 5000)],
    ],
)
def test_wavelength_units_conversion(
    input_wavelengths, expected_wavelengths_nm, verbose=True, *args, **kwargs
):
    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    wlmin, wlmax = input_wavelengths
    expected_wlmin, expected_wlmax = expected_wavelengths_nm
    sf = SpectrumFactory(
        wavelength_min=wlmin,
        wavelength_max=wlmax,
        wstep=0.01,
        cutoff=1e-30,
        pressure=1,
        path_length=1,
        mole_fraction=1,
        isotope=[1],
        verbose=verbose,
    )
    sf.load_databank("HITRAN-CO-TEST")
    s = sf.eq_spectrum(Tgas=300)
    assert np.isclose(s.get_wavelength().min(), expected_wlmin)
    assert np.isclose(s.get_wavelength().max(), expected_wlmax)


@pytest.mark.fast
@pytest.mark.parametrize(
    ("input_wavenumbers", "expected_wavenumbers_cm1"),
    [
        [(2000 * 1 / u.cm, 230000 * 1 / u.m), (2000, 2300)],
    ],
)
def test_wavenumber_units_conversion(
    input_wavenumbers, expected_wavenumbers_cm1, verbose=True, *args, **kwargs
):
    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    wmin, wmax = input_wavenumbers
    expected_wmin, expected_wmax = expected_wavenumbers_cm1
    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        wstep=0.01,
        cutoff=1e-30,
        pressure=1,
        path_length=1,
        mole_fraction=1,
        isotope=[1],
        verbose=verbose,
    )
    sf.load_databank("HITRAN-CO-TEST")
    s = sf.eq_spectrum(Tgas=300)
    assert np.isclose(s.get_wavenumber().min(), expected_wmin)
    assert np.isclose(s.get_wavenumber().max(), expected_wmax)


@pytest.mark.fast
@pytest.mark.parametrize(
    ("input_pressure, expected_pressure_bar"),
    [
        (1, 1),
        (1 * u.mbar, 1e-3),
        (10 * u.bar, 10),
    ],
)
def test_pressure_units_conversion(
    input_pressure, expected_pressure_bar, verbose=True, *args, **kwargs
):
    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    sf = SpectrumFactory(
        wavelength_min=4300,
        wavelength_max=4500,
        wstep=0.01,
        cutoff=1e-30,
        pressure=input_pressure,
        path_length=1,
        mole_fraction=1,
        isotope=[1],
        verbose=verbose,
        warnings={"AccuracyError": "ignore", "AccuracyWarning": "ignore"},
    )
    sf.load_databank("HITRAN-CO-TEST")
    s = sf.eq_spectrum(Tgas=300)
    assert np.isclose(s.conditions["pressure"], expected_pressure_bar)


@pytest.mark.fast
@pytest.mark.parametrize(
    ("input_pathlength, expected_pathlength_cm"),
    [
        (1, 1),
        (1 * u.m, 100),
        (1 * u.cm, 1),
    ],
)
def test_pathlength_units_conversion(
    input_pathlength, expected_pathlength_cm, verbose=True, *args, **kwargs
):
    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    sf = SpectrumFactory(
        wavelength_min=4300,
        wavelength_max=4500,
        wstep=0.01,
        cutoff=1e-30,
        pressure=1,
        path_length=input_pathlength,
        mole_fraction=1,
        isotope=[1],
        verbose=verbose,
    )
    sf.load_databank("HITRAN-CO-TEST")
    s = sf.eq_spectrum(Tgas=300)
    assert np.isclose(s.conditions["path_length"], expected_pathlength_cm)


@pytest.mark.fast
@pytest.mark.parametrize(
    ("input_temperature, expected_temperature_K"),
    [
        (300, 300),
        (300 * u.K, 300),
        (300 * u.deg_C, 573.15),
    ],
)
def test_temperature_units_conversion(
    input_temperature, expected_temperature_K, verbose=True, *args, **kwargs
):
    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    sf = SpectrumFactory(
        wavelength_min=4300,
        wavelength_max=4500,
        wstep=0.001,
        cutoff=1e-30,
        pressure=1,
        mole_fraction=1,
        isotope=[1],
        Tref=300 * u.K,
        verbose=verbose,
    )
    sf.load_databank("HITRAN-CO-TEST")
    s = sf.eq_spectrum(
        Tgas=input_temperature, pressure=20 * u.mbar, path_length=1 * u.mm
    )
    assert np.isclose(s.conditions["Tgas"], expected_temperature_K)
    assert np.isclose(s.conditions["path_length"], 0.1)  # cm
    assert np.isclose(s.conditions["pressure"], 0.02)


@pytest.mark.fast
def test_wstep_auto_method_sf(verbose=True, plot=False, *args, **kwargs):
    """Test to check that on computing several spectrum from the same Spectrum
    Factory object we get the different wstep for each case using auto method"""

    import radis
    from radis.misc.basics import round_off

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    sf = SpectrumFactory(
        wavelength_min=4400,
        wavelength_max=4800,
        mole_fraction=0.01,
        cutoff=1e-25,
        wstep="auto",
        isotope=[1],
        db_use_cached=True,
        self_absorption=True,
        verbose=verbose,
    )

    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.warnings["HighTemperatureWarning"] = "ignore"

    sf.load_databank("HITRAN-CO-TEST")

    s1 = sf.eq_spectrum(300, pressure=1)
    assert sf._wstep == "auto"

    wstep_calculated = s1.get_conditions()["wstep"]

    # Checking computed wstep and expected wstep are equal
    assert wstep_calculated == round_off(
        sf.min_width / radis.config["GRIDPOINTS_PER_LINEWIDTH_WARN_THRESHOLD"]
    )

    s2 = sf.eq_spectrum(300, pressure=0.2)
    assert sf._wstep == "auto"

    s3 = sf.eq_spectrum(300, pressure=0.001)
    assert sf._wstep == "auto"

    assert (
        s1.get_conditions()["wstep"]
        != s2.get_conditions()["wstep"]
        != s3.get_conditions()["wstep"]
    )


# @pytest.mark.fast #not fast, Nicolas Minesi 08/04/2024
def test_all_spectrum_using_wstep_auto(verbose=True, plot=False, *args, **kwargs):
    """Checks all methods to calculate Spectrum works with "auto" mode"""
    Tgas = 1000

    sf = SpectrumFactory(
        wavelength_min=4165,
        wavelength_max=4200,
        mole_fraction=1,
        path_length=0.3,
        cutoff=1e-23,
        molecule="CO2",
        isotope=1,
        wstep="auto",
        optimization=None,
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.warnings["HighTemperatureWarning"] = "ignore"
    sf.load_databank("HITRAN-CO2-TEST")

    sf.eq_spectrum(Tgas)
    wstep_1 = sf.params.wstep

    sf.non_eq_spectrum(Tvib=Tgas, Trot=Tgas)
    wstep_2 = sf.params.wstep

    sf.eq_bands(Tgas)
    wstep_3 = sf.params.wstep

    sf.non_eq_bands(Tvib=Tgas, Trot=Tgas)
    wstep_4 = sf.params.wstep

    assert wstep_1 == wstep_2 == wstep_3 == wstep_4
    assert sf._wstep == "auto"


@pytest.mark.skipif(isinstance(vaex, NotInstalled), reason="Vaex not available")
def test_vaex_and_pandas_spectrum():
    """Compares spectrum calculated using vaex and pandas are same.
    CONVOLUTED_QUANTITIES and dataframe df is compared to ensure the spectrum is same.

    Here are things that we ensure are same for both vaex and pandas that
    -Spectra calculated with Vaex & Pandas are the same
    -Spectra calculated with Vaex & Pandas using Database cutoff are the same
    Added in https://github.com/radis/radis/pull/580"""

    from radis import calc_spectrum

    # computing spectrum in vaex dataframe format
    config["DATAFRAME_ENGINE"] = "vaex"
    s_vaex, factory_s_vaex = calc_spectrum(
        1900,
        2300,  # cm-1
        molecule="CO",
        isotope="1,2,3",
        pressure=1.01325,  # bar
        Tgas=700,  # K
        # mole_fraction={'CO': 0.5, 'CO2': 0.3},
        mole_fraction=0.1,
        path_length=1,  # cm
        databank="hitran",  # or 'hitemp', 'geisa', 'exomol'
        # use_cached=False,
        return_factory=True,
    )

    # computing spectrum in pandas dataframe format
    config["DATAFRAME_ENGINE"] = "pandas"
    s_pandas, factory_s_pandas = calc_spectrum(
        1900,
        2300,  # cm-1
        molecule="CO",
        isotope="1,2,3",
        pressure=1.01325,  # bar
        Tgas=700,  # K
        mole_fraction=0.1,
        path_length=1,  # cm
        databank="hitran",  # or 'hitemp', 'geisa', 'exomol'
        # use_cached=False,
        return_factory=True,
    )

    assert np.allclose(
        s_vaex.get("absorbance"), s_pandas.get("absorbance"), equal_nan=True
    )

    for column in factory_s_pandas.df1.columns:
        assert np.all(
            factory_s_vaex.df1[column].to_numpy() == factory_s_pandas.df1[column]
        )

    # %% Additional test with other options (cutoff, ...)
    from radis import SpectrumFactory

    config["DATAFRAME_ENGINE"] = "vaex"
    sf_vaex = SpectrumFactory(
        wavelength_min=4200,
        wavelength_max=4500,
        cutoff=1e-23,
        molecule="CO2",
    )
    sf_vaex.fetch_databank("hitran")
    s_vaex = sf_vaex.eq_spectrum(Tgas=2000)  # failing on the 08/03/2023 - minouHub

    config["DATAFRAME_ENGINE"] = "pandas"
    sf_pd = SpectrumFactory(
        wavelength_min=4200,
        wavelength_max=4500,
        cutoff=1e-23,
        molecule="CO2",
    )
    sf_pd.fetch_databank("hitran")
    s_pd = sf_pd.eq_spectrum(Tgas=2000)
    for column in sf_pd.df1.columns:
        assert np.all(sf_vaex.df1[column].to_numpy() == sf_pd.df1[column])

    assert np.allclose(s_vaex.get("absorbance"), s_pd.get("absorbance"), equal_nan=True)


# %%
@pytest.mark.skipif(isinstance(vaex, NotInstalled), reason="Vaex not available")
def test_vaex_and_pandas_spectrum_noneq():
    """Compares the Spectrum calculated under non-equilibrium conditions .
    Comparison is made between CONVOLUTED_QUANTITIES and dataframe df.

    Here are things that we ensure are same for both vaex and pandas that
    -Spectra calculated with Vaex & Pandas are the same
    -Spectra calculated with Vaex & Pandas using Database cutoff are the same
    Added in https://github.com/radis/radis/pull/580"""

    from radis import calc_spectrum

    conditions = {
        "wmin": 2300,
        "wmax": 2320,
        "molecule": "CO",
        "isotope": 1,
        "pressure": 1.01325,
        "mole_fraction": 0.1,
        "wstep": "auto",
        "path_length": 1,
        "databank": "hitran",
        "verbose": 3,
        "return_factory": True,
    }
    # Calculating spectrum using vae
    config["DATAFRAME_TYPE"] = "vaex"
    s, factory_s = calc_spectrum(
        **conditions,
        Tgas=700,  # K
        Tvib=710,
        Trot=710,
    )

    # Calculating spectrum using pandas
    config["DATAFRAME_ENGINE"] = "pandas"
    s1, factory_s1 = calc_spectrum(
        **conditions,
        Tgas=700,  # K
        Tvib=710,
        Trot=710,
    )

    import numpy as np

    # Comparing different quantities
    assert np.allclose(s.get("absorbance"), s1.get("absorbance"), equal_nan=True)

    # Comparing different columns of dataframe df
    for column in factory_s1.df1.columns:
        assert np.all(factory_s1.df1[column] == factory_s.df1[column].to_numpy())


# --------------------------
if __name__ == "__main__":

    printm("Testing factory:", pytest.main(["test_factory.py", "--pdb"]))
#    printm('Testing factory:', pytest.main(['test_factory.py', '-k', 'test_wavenumber_units_conversion']))
