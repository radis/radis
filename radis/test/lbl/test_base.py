# -*- coding: utf-8 -*-
"""
Created on Mon May  7 17:34:52 2018

@author: erwan
"""

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pytest

import radis
from radis import get_residual, sPlanck
from radis.lbl import SpectrumFactory
from radis.lbl.base import get_waverange
from radis.misc.printer import printm
from radis.misc.progress_bar import ProgressBar
from radis.misc.utils import Default
from radis.test.utils import setup_test_line_databases


def test_linestrength_calculations(*args, **kwargs):
    """Compare and validate linestrength calculations from
    :py:func:`radis.lbl.base.linestrength_from_Einstein`,
    :py:meth:`radis.lbl.base.BaseFactory.calc_linestrength_eq`,
    :py:meth:`radis.lbl.base.BaseFactory.calc_linestrength_non_eq`,

    """

    # Test linestrength calculations
    from radis.io.hitran import hit2df
    from radis.test.utils import getTestFile

    df = hit2df(getTestFile("hitran_co_3iso_2000_2300cm.par"))

    from radis.db.classes import get_molecule_identifier
    from radis.db.molparam import MolParams
    from radis.lbl.base import linestrength_from_Einstein
    from radis.levels.partfunc import PartFuncTIPS

    molpar = MolParams()
    M = get_molecule_identifier("CO")

    Q_arr = df.iso.map({iso: PartFuncTIPS(M, iso).at(296) for iso in df.iso.unique()})
    Ia_arr = df.iso.map(
        {iso: molpar.get(M, iso, "abundance") for iso in df.iso.unique()}
    )

    # Reference linestrengths recomputed from Einstein coefficients match
    # within 0.1%
    S = linestrength_from_Einstein(df.A, df.gp, df.El, Ia_arr, df.wav, Q_arr, 296)
    assert np.allclose(df.int, S, rtol=1e-3, atol=0)

    #%% Now, at 1000 K
    Q_arr = df.iso.map({iso: PartFuncTIPS(M, iso).at(1000) for iso in df.iso.unique()})
    S_1000 = linestrength_from_Einstein(df.A, df.gp, df.El, Ia_arr, df.wav, Q_arr, 1000)

    from radis.lbl.factory import SpectrumFactory

    sf = SpectrumFactory(2000, 3000)
    sf.load_databank(
        path=getTestFile("hitran_co_3iso_2000_2300cm.par"),
        format="hitran",
        parfuncfmt="hapi",
    )

    # TODO : write an example of all the calculation steps in SpectrumFactory
    # sf._reinitialize()  # creates scaled dataframe df1 from df0  # TODO: make a public function.
    # sf.calc_linestrength_eq(1000)

    sf.eq_spectrum(1000)
    # Linestrengths computed from Einstein coefficients by linestrength_from_Einstein
    # match the one from the SpectrumFactory within 0.1%
    assert np.allclose(sf.df1.S, S_1000, rtol=1e-3, atol=0)

    #%% Using non-eq calculations
    S_eq = sf.df1.S.copy()

    sf.non_eq_spectrum(1000, 1000)
    # Linestrengths computed with equilibrium routine (scaling reference) match
    # the ones from the nonequilibrium routine (from Einstein coefficients)
    # within 0.1%
    assert np.allclose(sf.df1.S, S_eq, rtol=1e-3, atol=0)


@pytest.mark.fast
def test_export_populations(plot=True, verbose=True, warnings=True, *args, **kwargs):
    """Check populations calculated in the nonequilibrium module are exported in Spectrum.

    Compare with hardcoded values

    Examples
    --------
    ::
        sf = SpectrumFactory(export_populations="rovib")
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    if verbose:
        printm(">>> test_export_populations")

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
        export_populations="rovib",
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.load_databank("HITRAN-CO-TEST")

    # Populations cannot be calculated at equilibrium (no access to energy levels)
    s = sf.eq_spectrum(300)
    with pytest.raises(ValueError):  # .. we expect this error here!
        pops = s.get_populations("CO", isotope=1, electronic_state="X")

    # Calculate populations using the non-equilibrium module:
    s = sf.non_eq_spectrum(300, 300)
    pops = s.get_populations("CO", isotope=1, electronic_state="X")
    if not "n" in list(pops["rovib"].keys()):
        raise ValueError("Populations not calculated")
    if plot:
        s.plot_populations()

    # Test hardcoded populations
    assert np.isclose(pops["rovib"]["n"].iloc[0], 0.0091853446840826653)
    assert np.isclose(pops["rovib"]["n"].iloc[1], 0.027052543988733215)
    assert np.isclose(pops["rovib"]["n"].iloc[2], 0.04345502115897712)

    # Compare with factory
    # Plot populations:
    with pytest.raises(ValueError):
        sf.plot_populations()  # no isotope given: error expected
    if plot:
        sf.plot_populations("rovib", isotope=1)
        plt.close()  # no need to keep it open, we just tested the function


@pytest.mark.fast
def test_export_rovib_fractions(
    plot=True, verbose=True, warnings=True, *args, **kwargs
):
    """Compare rovib fraction (nu_vib, nl_vib, nu_rot, nl_rot) are calculated
    in the nonequilibrium module

    Examples
    --------
    ::
        sf = SpectrumFactory(export_lines=True, ...)
        sf.misc.export_rovib_fraction = True  # required from 0.9.30
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    if verbose:
        printm(">>> test_export_rovib_fractions")

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
        export_lines=True,
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.load_databank("HITRAN-CO-TEST")
    sf.misc.export_rovib_fraction = True  # required from 0.9.30

    # Calculate populations using the non-equilibrium module:
    s = sf.non_eq_spectrum(300, 300)
    # Test calculated quantities are there
    # assert "Qref" in sf.df1.attrs
    assert "Qvib" in sf.df1.attrs
    assert hasattr(sf.df1, "Qrotu")
    assert hasattr(sf.df1, "Qrotl")

    assert hasattr(s.lines, "nu_vib")
    assert hasattr(s.lines, "nl_vib")
    assert hasattr(s.lines, "nu_rot")
    assert hasattr(s.lines, "nl_rot")

    return True


@pytest.mark.fast
def test_populations_CO2_hamiltonian(
    plot=True, verbose=True, warnings=True, *args, **kwargs
):
    """Calculate nonequilibrium modes with the CO2 Hamiltonian

    ..warning::

        as we only use a reduced set of the CO2 effective Hamiltonian (< 3000 cm-1),
        many levels of the Line Database will not appear in the Levels Database.
        We will need to either filter the Line Database beforehand, or run it
        a first time and remove all levels not found.

        This database is obviously not to be used in a Production code!

    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    if verbose:
        printm(">>> _test_populations_CO2_hamiltonian")

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    sf = SpectrumFactory(
        wavenum_min=2283.7,
        wavenum_max=2285.1,
        wstep=0.001,
        cutoff=1e-30,
        path_length=0.1,
        mole_fraction=400e-6,
        isotope=[1],
        medium="vacuum",
        truncation=5,
        export_populations="rovib",
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.load_databank("HITEMP-CO2-HAMIL-TEST", load_energies=True, load_columns="noneq")

    # First run a calculation at equilibrium
    s = sf.eq_spectrum(300)
    s.name = "equilibrium"

    # Now generate vibrational energies for a 2-T model
    # ... Note that this is arbitrary. Lookup Pannier & Dubuet 2020 for more.
    levels = sf.parsum_calc["CO2"][1]["X"].df
    levels["Evib"] = levels.Evib1 + levels.Evib2 + levels.Evib3

    # Calculate populations using the non-equilibrium module:
    # This will crash the first time (see error in docstrings of the function).
    with pytest.raises(AssertionError):
        sf.non_eq_spectrum(300, 300)

    sf.df0.dropna(inplace=True)

    # Retry:
    s_noneq = sf.non_eq_spectrum(300, 300)
    s_noneq.name = "nonequilibrium"

    # Tests:

    #    s.compare_with(s_noneq, plot=plot)
    assert s_noneq.get_power() > 0
    # TODO: implement actual Assertions (right now we're just checking that
    # functions are properly parsed)

    # TODO. Fix below (and move in dedicated test):
    #    s.line_survey()

    return True


# @pytest.mark.needs_connection
def test_optically_thick_limit_1iso(verbose=True, plot=True, *args, **kwargs):
    """Test that we find Planck in the optically thick limit

    In particular, this test will fail if :

    - linestrength are not properly calculated

    - at noneq, linestrength and emission integrals are mixed up

    The test should be run for 1 and several isotopes, because different
    calculations paths are used internally, and this can lead to different
    errors.

    Also, this test is used to run with DEBUG_MODE = True, which will
    check that isotopes and molecule ids are what we expect in all the
    groupby() loops that make the production code very fast.

    Notes
    -----

    switched from large band calculation with [HITRAN-2016]_ to a calculation with
    the embedded [HITEMP-2010]_ fragment (shorter range, but no need to download files)

    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    # Force DEBUG_MODE
    DEBUG_MODE = radis.config["DEBUG_MODE"]
    radis.config["DEBUG_MODE"] = True

    try:

        wavenum_min = 2284.2
        wavenum_max = 2284.6

        P = 0.017  # bar
        wstep = 0.001  # cm-1

        Tgas = 1200

        # %% Generate some CO2 emission spectra
        # --------------
        sf = SpectrumFactory(
            wavenum_min=wavenum_min,
            wavenum_max=wavenum_max,
            molecule="CO2",
            mole_fraction=1,
            path_length=0.05,
            cutoff=1e-25,
            truncation=0.5,
            export_populations=False,  #'vib',
            export_lines=False,
            isotope=1,
            wstep=wstep,
            pseudo_continuum_threshold=0,
            pressure=P,
            verbose=False,
        )
        sf.warnings["NegativeEnergiesWarning"] = "ignore"
        sf.load_databank("HITEMP-CO2-TEST", load_columns="noneq")
        pb = ProgressBar(3, active=verbose)
        s_eq = sf.eq_spectrum(Tgas=Tgas, mole_fraction=1, name="Equilibrium")
        pb.update(1)
        s_2T = sf.non_eq_spectrum(
            Tvib=Tgas, Trot=Tgas, mole_fraction=1, name="Noneq (2T)"
        )
        pb.update(2)
        s_4T = sf.non_eq_spectrum(
            Tvib=(Tgas, Tgas, Tgas), Trot=Tgas, mole_fraction=1, name="Noneq (4T)"
        )
        pb.update(3)
        s_plck = sPlanck(
            wavelength_min=2000,  # =wavelength_min,
            wavelength_max=5000,  # =wavelength_max - wstep,   # there is a border effect on last point
            T=Tgas,
        )
        pb.done()

        # %% Post process:
        # MAke optically thick, and compare with Planck

        if plot:
            s_plck.plot(wunit="nm", Iunit="mW/cm2/sr/nm", lw=2)
        for s in [s_eq, s_2T, s_4T]:

            s.rescale_path_length(1e6)

            if plot:
                s.plot(wunit="nm", nfig="same", lw=4)

            if verbose:
                printm(
                    "Residual between opt. thick CO2 spectrum ({0}) and Planck: {1:.2g}".format(
                        s.name,
                        get_residual(s, s_plck, "radiance_noslit", ignore_nan=True),
                    )
                )

            #            assert get_residual(s, s_plck, 'radiance_noslit', ignore_nan=True) < 1e-3
            assert get_residual(s, s_plck, "radiance_noslit", ignore_nan=True) < 0.9e-4
        if plot:
            plt.legend()

        if verbose:
            printm("Tested optically thick limit is Planck (1 isotope): OK")

    finally:
        # Reset DEBUG_MODE
        radis.config["DEBUG_MODE"] = DEBUG_MODE


def test_optically_thick_limit_2iso(verbose=True, plot=True, *args, **kwargs):
    """Test that we find Planck in the optically thick limit

    In particular, this test will fail if :

    - linestrength are not properly calculated

    - at noneq, linestrength and emission integrals are mixed up

    The test should be run for 1 and several isotopes, because different
    calculations paths are used internally, and this can lead to different
    errors.

    Also, this test is used to run with DEBUG_MODE = True, which will
    check that isotopes and molecule ids are what we expect in all the
    groupby() loops that make the production code very fast.

    Notes
    -----

    switched from large band calculation with [HITRAN-2016]_ to a calculation with
    the embedded [HITEMP-2010]_ fragment (shorter range, but no need to download files)

    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    # Force DEBUG_MODE
    DEBUG_MODE = radis.config["DEBUG_MODE"]
    radis.config["DEBUG_MODE"] = True

    try:

        wavenum_min = 2284.2
        wavenum_max = 2284.6

        P = 0.017  # bar
        wstep = 0.001  # cm-1

        Tgas = 1200

        # %% Generate some CO2 emission spectra
        # --------------
        sf = SpectrumFactory(
            wavenum_min=wavenum_min,
            wavenum_max=wavenum_max,
            mole_fraction=1,
            path_length=0.05,
            cutoff=1e-25,
            truncation=0.5,
            export_populations=False,  #'vib',
            export_lines=False,
            molecule="CO2",
            isotope=[1, 2],
            wstep=wstep,
            pseudo_continuum_threshold=0,
            pressure=P,
            verbose=False,
        )
        sf.load_databank("HITEMP-CO2-TEST", load_columns="noneq")
        sf.warnings["NegativeEnergiesWarning"] = "ignore"
        #        sf.warnings['MissingSelfBroadeningWarning'] = 'ignore'
        pb = ProgressBar(3, active=verbose)
        s_eq = sf.eq_spectrum(Tgas=Tgas, mole_fraction=1, name="Equilibrium")
        pb.update(1)
        s_2T = sf.non_eq_spectrum(
            Tvib=Tgas, Trot=Tgas, mole_fraction=1, name="Noneq (2T)"
        )
        pb.update(2)
        s_4T = sf.non_eq_spectrum(
            Tvib=(Tgas, Tgas, Tgas), Trot=Tgas, mole_fraction=1, name="Noneq (4T)"
        )
        pb.update(3)
        s_plck = sPlanck(
            wavelength_min=2000,  # =wavelength_min,
            wavelength_max=5000,  # =wavelength_max - wstep,   # there is a border effect on last point
            T=Tgas,
        )
        pb.done()

        # %% Post process:
        # MAke optically thick, and compare with Planck

        for s in [s_eq, s_2T, s_4T]:

            s.rescale_path_length(1e6)

            if plot:

                nfig = "test_opt_thick_limit_2iso {0}".format(s.name)
                plt.figure(nfig).clear()
                s.plot(wunit="nm", nfig=nfig, lw=4)
                s_plck.plot(wunit="nm", nfig=nfig, Iunit="mW/cm2/sr/nm", lw=2)
                plt.legend()

            if verbose:
                printm(
                    "Residual between opt. thick CO2 spectrum ({0}) and Planck: {1:.2g}".format(
                        s.name,
                        get_residual(s, s_plck, "radiance_noslit", ignore_nan=True),
                    )
                )

            #            assert get_residual(s, s_plck, 'radiance_noslit', ignore_nan=True) < 1e-3
            assert get_residual(s, s_plck, "radiance_noslit", ignore_nan=True) < 0.9e-4

        if verbose:
            printm("Tested optically thick limit is Planck (2 isotopes): OK")

    finally:
        # Reset DEBUG_MODE
        radis.config["DEBUG_MODE"] = DEBUG_MODE


def test_get_waverange(*args, **kwargs):

    assert get_waverange(1 * u.um, 10 * u.um, medium="vacuum") == (1000, 10000)

    # 'wunit' is none
    # ...'wmin/wmax' is none
    # ...... wavenumber is passed > should return the same wavenumbers
    assert get_waverange(wavenum_min=10, wavenum_max=20, wunit=Default("cm-1")) == (
        10,
        20,
    )

    # ...... wavelength is passed > should convert and return the wavenumbers
    assert np.isclose(
        get_waverange(
            wavelength_min=1, wavelength_max=2, medium="vacuum", wunit=Default("cm-1")
        ),
        (5000000.0, 10000000.0),
    ).all()

    # ....... passed both wavenumber and wavelength > should throw error
    with pytest.raises(ValueError):
        get_waverange(
            wavenum_min=1000,
            wavenum_max=2000,
            wavelength_min=1,
            wavelength_max=2,
            medium="vacuum",
            wunit=Default("cm-1"),
        )

    # ...... passed neither wavenumber nor wavlength > should throw error
    with pytest.raises(ValueError):
        get_waverange(wunit=Default("cm-1"))

    # ... 'wmin/wmax' are passed as values without accompanying units
    # ...... wavenumber is passed > should throw error due to multiple input wavespace values
    with pytest.raises(ValueError):
        get_waverange(
            wavenum_min=1, wavenum_max=2, wmin=10, wmax=20, wunit=Default("cm-1")
        )

    # ...... wavelength is passed > should throw error due to multiple input wavespace values
    with pytest.raises(ValueError):
        get_waverange(wmin=10, wmax=20, wavelength_min=1, wavelength_max=2)

    # ...... passed both wavenumber and wavelength > should throw error due to multiple input wavespace values
    with pytest.raises(ValueError):
        get_waverange(
            wmin=1,
            wmax=2,
            wavenum_min=10,
            wavenum_max=20,
            wavelength_min=100,
            wavelength_max=200,
            wunit=Default("cm-1"),
        )

    # ...... passed neither wavenumber nor wavelength > should return wavenumber after converting wmin/wmax assuming default units
    assert get_waverange(wmin=10, wmax=20, wunit=Default("cm-1")) == (10.0, 20.0)

    # ... 'wmin/wmax' are passed as values with accompanying units
    # ...... wavenumber is passed > should throw error due to multiple input wavespace values
    with pytest.raises(ValueError):
        get_waverange(
            wavenum_min=10,
            wavenum_max=20,
            wmin=100 * (1 / u.cm),
            wmax=200 * (1 / u.cm),
            wunit=Default("cm-1"),
        )

    # ...... wavelength is passed > should throw error due to multiple input wavespace values
    with pytest.raises(ValueError):
        get_waverange(
            wmin=100 * (1 / u.cm),
            wmax=200 * (1 / u.cm),
            wavelength_min=10,
            wavelength_max=20,
            wunit=Default("cm-1"),
        )

    # ...... passed both wavenumber and wavelength > should throw error due to multiple input wavespace values
    with pytest.raises(ValueError):
        get_waverange(
            wavenum_min=10,
            wavenum_max=20,
            wmin=100 * (1 / u.cm),
            wmax=200 * (1 / u.cm),
            wavelength_min=10,
            wavelength_max=20,
            wunit=Default("cm-1"),
        )

    # ...... passed neither wavenumber nor wavelength > should return wavenumber after converting wmin/wmax from accompanying unit
    assert get_waverange(
        wmin=100 * (1 / u.cm), wmax=200 * (1 / u.cm), wunit=Default("cm-1")
    ) == (100.0, 200.0)

    assert np.isclose(
        get_waverange(
            wmin=1 * u.cm, wmax=2 * u.cm, medium="vacuum", wunit=Default("cm-1")
        ),
        (0.5, 1.0),
    ).all()

    # 'wunit' is not none
    # ... 'wmin/wmax' is none
    # ...... wavenumber is passed > should throw error as wunit can only be passed with wmin/wmax
    with pytest.raises(ValueError):
        get_waverange(wavenum_min=1, wavenum_max=2, wunit="cm")

    # ...... wavelength is passed > should throw error as wunit can only be passed with wmin/wmax
    with pytest.raises(ValueError):
        get_waverange(wavelength_min=1, wavelength_max=2, wunit="cm")

    # ...... passed both wavenumber and wavelength > should throw error as wunit can only be passed with wmin/wmax
    with pytest.raises(ValueError):
        get_waverange(
            wavenum_min=1,
            wavenum_max=2,
            wavelength_min=10,
            wavelength_max=20,
            wunit="cm",
        )

    # ...... passed neither wavenumber nor wavlength > should throw error
    with pytest.raises(ValueError):
        get_waverange(wunit="cm")

    # ... 'wmin/wmax' are passed as values without accompanying units
    # ...... wavenumber is passed > should throw error as only one set of wavespace parameters can be passed
    with pytest.raises(ValueError):
        get_waverange(wmin=1, wmax=2, wavenum_min=10, wavenum_max=20, wunit="cm")

    # ...... wavelength is passed > should throw error as only one set of wavespace parameters can be passed
    with pytest.raises(ValueError):
        get_waverange(wmin=1, wmax=2, wavelength_min=10, wavelength_max=20, wunit="cm")

    # ...... passed both wavenumber and wavelength > should throw error as only one set of wavespace parameters can be passed
    with pytest.raises(ValueError):
        get_waverange(
            wmin=1,
            wmax=2,
            wavenum_min=10,
            wavenum_max=20,
            wavelength_min=100,
            wavelength_max=200,
            wunit="cm",
        )

    # ...... passed neither wavenumber nor wavlength > should return wavenumber after converting wmin/wmax with given wunit
    assert np.isclose(
        get_waverange(wmin=1, wmax=2, wunit="cm"),
        (0.4998637271242577, 0.9997274542485038),
    ).all()

    assert get_waverange(wmin=1, wmax=2, wunit="cm-1") == (1.0, 2.0)

    # ... 'wmin/wmax' are passed as values with accompanying units
    # ...... wavenumber is passed > should throw error as only one set of wavespace parameters can be passed
    with pytest.raises(ValueError):
        get_waverange(
            wmin=1 * u.cm, wmax=2 * u.cm, wavenum_min=10, wavenum_max=20, wunit="cm"
        )

    # ...... wavelength is passed > should throw error as only one set of wavespace parameters can be passed
    with pytest.raises(ValueError):
        get_waverange(
            wmin=1 * u.cm,
            wmax=2 * u.cm,
            wavelength_min=10,
            wavelength_max=20,
            wunit="cm",
        )

    # ...... passed both wavenumber and wavelength > should throw error as only one set of wavespace parameters can be passed
    with pytest.raises(ValueError):
        get_waverange(
            wmin=1 * u.cm,
            wmax=2 * u.cm,
            wavenum_min=1,
            wavenum_max=2,
            wavelength_min=10,
            wavelength_max=20,
            wunit="cm",
        )

    # ...... passed neither wavenumber nor wavlength > should return wavenumber after converting wmin/wmax with given wunit
    assert np.isclose(
        get_waverange(wmin=1 * u.cm, wmax=2 * u.cm, wunit="cm"),
        (0.4998637271242577, 0.9997274542485038),
    ).all()

    assert get_waverange(wmin=1 * (1 / u.cm), wmax=2 * (1 / u.cm), wunit="cm-1") == (
        1.0,
        2.0,
    )
    assert np.isclose(
        get_waverange(wmin=1 * u.cm, wmax=2 * u.cm, wunit="cm"),
        (0.4998637271242577, 0.9997274542485038),
    ).all()
    assert np.isclose(
        get_waverange(
            wavelength_min=1 * u.cm, wavelength_max=2 * u.cm, wunit=Default("cm-1")
        ),
        (0.4998637271242577, 0.9997274542485038),
    ).all()

    # ... passed wmin/wmax with units different from wunit > should throw error due to conflicting units
    with pytest.raises(ValueError):
        get_waverange(wmin=1 * u.cm, wmax=2 * u.cm, wunit="cm-1")


def _run_testcases(verbose=True, plot=True):

    test_linestrength_calculations()
    test_export_populations(plot=plot, verbose=verbose)
    test_export_rovib_fractions(plot=plot, verbose=verbose)
    test_populations_CO2_hamiltonian(plot=plot, verbose=verbose)
    test_optically_thick_limit_1iso(plot=plot, verbose=verbose)
    test_optically_thick_limit_2iso(plot=plot, verbose=verbose)
    test_get_waverange()


if __name__ == "__main__":
    _run_testcases()
