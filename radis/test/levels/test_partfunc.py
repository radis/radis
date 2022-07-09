# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 10:16:13 2017

@author: erwan

Examples
--------

Run all tests::

    pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that a 'fast' label)::

    pytest -m fast


-------------------------------------------------------------------------------

"""

import os
from os.path import basename, exists, getmtime

import numpy as np
import pytest
from numpy import exp

from radis import SpectrumFactory
from radis.db.molecules import Molecules
from radis.levels.partfunc import PartFunc_Dunham, PartFuncTIPS
from radis.levels.partfunc_cdsd import PartFuncCO2_CDSDcalc, PartFuncCO2_CDSDtab
from radis.misc.printer import printm
from radis.misc.warning import DeprecatedFileWarning
from radis.phys.constants import hc_k
from radis.test.utils import getTestFile, setup_test_line_databases

fig_prefix = basename(__file__) + ": "


# %% Test

# never add @pytest.mark.fast so we don't delete cached files for 'fast' tests
def test_delete_all_cached_energies(verbose=True, warnings=True, *args, **kwargs):
    """Doesnt really test anything, but cleans all cached energy levels"""

    for molecule, isotopes in Molecules.items():
        for isotope, states in isotopes.items():
            for state, ElecState in states.items():
                # Delete cache file for Dunham expansions
                energies = PartFunc_Dunham(ElecState, use_cached="return")
                cachefile = energies.cachefile
                if exists(cachefile):
                    os.remove(cachefile)
                    print(
                        (
                            "Cleaned cached energies for {0}".format(
                                ElecState.get_fullname()
                            )
                        )
                    )


@pytest.mark.fast
def test_cache_file_generation_and_update(verbose=True, *args, **kwargs):
    """Test that cache file process works correctly, using CO as an example.

    Expected behavior:

    - generated when 'regen'
    - updated when using 'True' but input parameters have changed
    - an error is raised when 'force' is used but input parameters have changed

    """
    # TODO: move it as a test of cache_files.py in RADIS

    S = Molecules["CO"][1]["X"]

    from os.path import exists, getmtime

    if verbose:
        printm("Testing behavior for Energy cache file with use_cached=")

    # Create a first cache file (or use existing one)
    db = PartFunc_Dunham(
        S, vmax=12, vmax_morse=48, Jmax=300, use_cached=True, verbose=verbose
    )

    # .. test file was created
    cachefile = db.cachefile
    assert exists(cachefile)
    last_modif = getmtime(cachefile)
    if verbose:
        printm("... True: energy cache file has been used/created (vmax=12)")

    # Force recalculating (using same inputs, but use_cached='regen')
    db = PartFunc_Dunham(
        S, vmax=12, vmax_morse=48, Jmax=300, use_cached="regen", verbose=verbose
    )

    # ... test that cache file was updated
    cachefile = db.cachefile
    assert getmtime(cachefile) > last_modif
    last_modif = getmtime(cachefile)
    if verbose:
        print("... 'regen': energy cache file has been created again (vmax=12)")

    # Recompute with different parameters

    # ... if using 'force', ensures that an error is raised
    with pytest.raises(DeprecatedFileWarning):  # DeprecatedFileError is expected
        db = PartFunc_Dunham(
            S, vmax=11, vmax_morse=48, Jmax=300, use_cached="force", verbose=verbose
        )
    if verbose:
        printm(
            "... 'force': an error was correctly triggered for new input parameters (vmax=11)"
        )

    # ... if using 'True', ensures that cache file is updated
    db = PartFunc_Dunham(
        S, vmax=11, vmax_morse=48, Jmax=300, use_cached=True, verbose=verbose
    )
    assert getmtime(cachefile) > last_modif
    if verbose:
        printm("... True: cache file was updated for new input parameters (vmax=11)")


@pytest.mark.needs_config_file
@pytest.mark.fast
@pytest.mark.needs_db_CDSD_HITEMP_PC
def test_CDSD_calc_vs_tab(verbose=True, warnings=True, *args, **kwargs):
    """Test 1: compare calculated PartFunc to the tabulated one"""

    from radis.misc.config import getDatabankEntries

    iso = 1
    database = "CDSD-HITEMP-PC"

    # Compare tab to hardcoded
    parfunc = getDatabankEntries(database)["parfunc"]
    Qf = PartFuncCO2_CDSDtab(iso, parfunc)
    assert np.isclose(Qf.at(300), 291.0447781984652, rtol=0.001)
    assert np.isclose(Qf.at(3000), 114689.88454184022, rtol=0.001)

    if verbose:
        printm("Tested CDSD tabulated is correct: OK")

    # Compare tab to calculated
    energies = getDatabankEntries(database)["levels"]
    levelsfmt = getDatabankEntries(database)["levelsfmt"]
    Qfc = PartFuncCO2_CDSDcalc(
        energies[iso], levelsfmt=levelsfmt, isotope=iso, use_cached=True
    )
    assert np.isclose(Qf.at(300), Qfc.at(300), rtol=0.001)
    assert np.isclose(Qf.at(3000), Qfc.at(3000), rtol=0.001)

    if verbose:
        printm("Tested CDSD Q_calc vs Q_tab give same output: OK")

    return True


@pytest.mark.needs_config_file
@pytest.mark.fast
def test_reduced_CDSD_calc_vs_tab(verbose=True, warnings=True, *args, **kwargs):
    """Test 1: compare calculated PartFunc to the tabulated one

    Version where we use the reduced set of CO2 levels (< 3000 cm-1)"""
    from radis.misc.config import getDatabankEntries

    iso = 1
    database = "HITEMP-CO2-HAMIL-TEST"

    # Compare tab to hardcoded
    parfunc = getDatabankEntries(database)["parfunc"]
    Qf = PartFuncCO2_CDSDtab(iso, parfunc)
    assert np.isclose(Qf.at(300), 291.0447781984652, rtol=0.001)
    assert np.isclose(Qf.at(3000), 114689.88454184022, rtol=0.001)

    if verbose:
        printm("Tested CDSD tabulated is correct: OK")

    # Compare tab to calculated
    energies = getDatabankEntries(database)["levels"]
    levelsfmt = getDatabankEntries(database)["levelsfmt"]
    Qfc = PartFuncCO2_CDSDcalc(
        energies[iso],
        levelsfmt=levelsfmt,
        isotope=iso,
        use_cached=True,
        verbose=verbose,
    )
    assert np.isclose(
        Qf.at(300), Qfc.at(300), rtol=0.01
    )  # reduced rtol to accommodate for the reduced set of levels in the test set
    # assert np.isclose(Qf.at(3000), Qfc.at(3000), rtol=0.001)  # of course doesnt work with the reduced set of levels in the test set

    if verbose:
        printm("Tested CDSD Q_calc vs Q_tab give same output: OK")


@pytest.mark.fast
def test_calculatedQ_match_HAPI_CO(
    vmax=11, jmax=300, plot=False, verbose=True, *args, **kwargs
):
    """Tested that Q ab_initio (Dunham) match HAPI for CO at different temperatures"""

    vmax = 11
    vmax_morse = 48
    jmax = 300
    iso = 1

    S = Molecules["CO"][iso]["X"]

    # Dont use cached: force recalculating
    db = PartFunc_Dunham(
        S, vmax=vmax, vmax_morse=vmax_morse, Jmax=jmax, use_cached=False
    )  # , ZPE=1081.584383)
    assert not db.use_cached

    #    if plot: db.plot_states()

    hapi = PartFuncTIPS(M=5, I=1)  # CO  # isotope

    us = []
    hap = []
    T = np.linspace(300, 3000)
    for Ti in T:
        us.append(db.at(Ti))
        hap.append(hapi.at(Ti))

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        import matplotlib.pyplot as plt

        plt.ion()

        plt.figure(fig_prefix + "Partition function Dunham vs Precomputed")
        plt.plot(T, us, "ok", label="NeQ")
        plt.plot(T, hap, "or", label="HAPI")
        plt.legend()
        plt.xlabel("Temperature (K)")
        plt.ylabel("Partition function")
        plt.title(
            "Ab-initio partition function calculations\n"
            + "(vmax:{0},jmax:{1})".format(vmax, jmax)
        )
        plt.tight_layout()

    # Compare Calculated  vs HAPI
    assert np.allclose(us, hap, rtol=0.02)

    if verbose:
        printm("Tested Q_CO ab_initio (Dunham) matches HAPI for 300 - 3000 K: OK")

    return True


@pytest.mark.fast
def test_calculatedQ_match_HAPI(plot=False, verbose=True, *args, **kwargs):
    """Tested that Q ab_initio (Dunham) match HAPI for different molecules
    and isotopes"""

    # molecule, isotope, temperature, absolute tolerance
    for molecule, iso, T, atol, rtol in [
        ("CO2", 1, 300, 0.9, 0.02),
        ("CO2", 1, 1000, 5.0, 0.02),
        ("CO2", 1, 3000, 2150, 0.02),
        ("CO2", 2, 300, 1.8, 0.02),
        ("CO2", 2, 1000, 25, 0.02),
        ("CO2", 2, 3000, 4962, 0.02),
        ("CO", 1, 300, 0.16, 0.02),
        ("CO", 1, 1000, 1.9, 0.02),
        ("CO", 1, 3000, 27, 0.02),
        ("CO", 2, 300, 0.31, 0.02),
        ("CO", 2, 1000, 4.1, 0.02),
        ("CO", 2, 3000, 56.9, 0.02),
        ("CO", 3, 300, 0.16, 0.02),
        ("CO", 3, 1000, 2.1, 0.02),
        ("CO", 3, 3000, 28.6, 0.02),
    ]:

        S = Molecules[molecule][iso]["X"]

        # Dont use cached: force recalculating
        db = PartFunc_Dunham(S)

        from radis.db.classes import get_molecule_identifier

        hapi = PartFuncTIPS(M=get_molecule_identifier(molecule), I=iso)

        Q_radis = db.at(T)
        Q_hapi = hapi.at(T)

        if verbose:
            print(
                "Q({0},iso={1},{2}K)\tRADIS: {3:.2f}\tHAPI={4:.2f}".format(
                    molecule, iso, T, Q_radis, Q_hapi
                )
            )

        try:
            assert np.isclose(Q_radis, Q_hapi, atol=atol)
            assert np.isclose(Q_radis, Q_hapi, atol=atol)
        except AssertionError:
            raise AssertionError(
                "Partition function for {0}, iso={1} ".format(molecule, iso)
                + "at {0}K doesnt match in RADIS ".format(T)
                + "({0:.4f}) and HAPI ({1:.4f})".format(Q_radis, Q_hapi)
            )

    if verbose:
        printm("Tested Q ab_initio (Dunham) matches HAPI: OK")

    return True


@pytest.mark.needs_config_file
@pytest.mark.fast
@pytest.mark.needs_db_CDSD_HITEMP_PC
def test_CDSD_calc_vs_ref(warnings=True, verbose=True, *args, **kwargs):
    """Test partition functions calculated with CDSD energy levels against
    hardcoded values"""

    from radis.misc.config import getDatabankEntries

    iso = 1

    energies = getDatabankEntries("CDSD-HITEMP-PC")["levels"]
    levelsfmt = getDatabankEntries("CDSD-HITEMP-PC")["levelsfmt"]
    Qf = PartFuncCO2_CDSDcalc(
        energy_levels=energies[iso],
        isotope=iso,
        use_cached=True,
        levelsfmt=levelsfmt,
    )
    assert np.isclose(Qf.at(300), 291.0447781984652, rtol=0.001)
    assert np.isclose(Qf.at(3000), 114689.88454184022, rtol=0.001)
    assert np.isclose(Qf.at_noneq(300, 300), 291.0447781984652, rtol=0.001)
    assert np.isclose(Qf.at_noneq(3000, 3000), 114689.88454184022, rtol=0.001)
    assert np.isclose(
        Qf.at_noneq(3000, 3000, overpopulation={"(0,1)": 3}, returnQvibQrot=True)[0],
        120053.34252537244,
        rtol=0.001,
    )

    if verbose:
        printm("Tested Q_CDSD values are correct : OK")


@pytest.mark.fast
def test_reduced_CDSD_calc_noneq(verbose=True, warnings=True, *args, **kwargs):
    """Compare calculated partition function at equilibrium and nonequilibrium
    using the CDSD-format

    Examples
    --------

        assert Qfc.at(300) != Qfc.at_noneq_3Tvib((300, 300, 300), 300)

    After redefining Evib so that Evib + Erot = E:

        assert np.isclose(Qfc.at(300), Qfc.at_noneq(300, 300), rtol=0.001)
        assert np.isclose(Qfc.at(300), Qfc.at_noneq_3Tvib((300, 300, 300), 300), rtol=0.001)
    """

    from radis.misc.config import getDatabankEntries
    from radis.test.utils import define_Evib_as_sum_of_Evibi

    iso = 1
    database = "HITEMP-CO2-HAMIL-TEST"

    # Compare tab to calculated
    energies = getDatabankEntries(database)["levels"]
    levelsfmt = getDatabankEntries(database)["levelsfmt"]
    Qfc = PartFuncCO2_CDSDcalc(
        energies[iso],
        levelsfmt=levelsfmt,
        isotope=iso,
        use_cached=True,
        verbose=verbose,
    )

    # Note that partition functions are not equal because some of the energy
    # goes in coupling terms in the Hamiltonian formulation
    assert Qfc.at(300) != Qfc.at_noneq_3Tvib((300, 300, 300), 300)

    # Below we redefine Evib so that Evib + Erot = E
    define_Evib_as_sum_of_Evibi(Qfc.df)

    # New Test:
    assert np.isclose(Qfc.at(300), Qfc.at_noneq(300, 300), rtol=0.001)
    assert np.isclose(Qfc.at(300), Qfc.at_noneq_3Tvib((300, 300, 300), 300), rtol=0.001)


# @pytest.mark.fast   # (very fast only once the cached database has been generated, else decently fast)
def test_recompute_Q_from_QvibQrot_Dunham_Evib3_Evib12Erot(
    verbose=True, warnings=True, *args, **kwargs
):
    """Calculate vibrational and rotational partition functions:

    - with Dunham expansions. Evib, Erot = (Evib3, Evib1+Evib2+Erot)
    - under nonequilibrium

    Calculate total rovibrational partition function, and compare

    Test if partition function can be recomputed correctly from vibrational
    populations and rotational partition function (note that we are in a coupled
    case so partition function is not simply the product of Qvib, Qrot)
    """

    iso = 1
    Tvib = 1500
    Trot = 300

    S = Molecules["CO2"][iso]["X"]

    Qf = PartFunc_Dunham(
        S,
        use_cached=True,
        group_energy_modes_in_2T_model={"CO2": (["Evib3"], ["Evib1", "Evib2", "Erot"])},
    )

    Q = Qf.at_noneq(Tvib, Trot)
    _, Qvib, dfQrot = Qf.at_noneq(Tvib, Trot, returnQvibQrot=True)
    if verbose:
        printm("Q", Q)
    if verbose:
        printm("Qvib", Qvib)

    # 1) Test Q vs Q recomputed from Qrot, Qvib

    # Recompute Qtot
    df = dfQrot
    Q2 = ((df.gvib * exp(-df.Evib * hc_k / Tvib)) * df.Qrot).sum()
    # Todo: non Boltzmann case

    assert np.isclose(Q, Q2)

    if verbose:
        printm("Tested Q vs recomputed from (Qvib, Qrot) are the same: OK")

    return True


# @pytest.mark.fast   # (very fast only once the cached database has been generated, else decently fast)
def test_recompute_Q_from_QvibQrot_Dunham_Evib123_Erot(
    verbose=True, warnings=True, *args, **kwargs
):
    """Calculate vibrational and rotational partition functions:

    - with Dunham expansions. Evib, Erot = (Evib1+Evib2+Evib3, Erot)
    - under nonequilibrium

    Calculate total rovibrational partition function, and compare

    Test if partition function can be recomputed correctly from vibrational
    populations and rotational partition function (note that we are in a coupled
    case so partition function is not simply the product of Qvib, Qrot)
    """

    iso = 1
    Tvib = 1500
    Trot = 300

    S = Molecules["CO2"][iso]["X"]

    Qf = PartFunc_Dunham(
        S,
        use_cached=True,
        group_energy_modes_in_2T_model={"CO2": (["Evib1", "Evib2", "Evib3"], ["Erot"])},
    )

    Q = Qf.at_noneq(Tvib, Trot)
    _, Qvib, dfQrot = Qf.at_noneq(Tvib, Trot, returnQvibQrot=True)
    if verbose:
        printm("Q", Q)
    if verbose:
        printm("Qvib", Qvib)

    # 1) Test Q vs Q recomputed from Qrot, Qvib

    # Recompute Qtot
    df = dfQrot
    Q2 = ((df.gvib * exp(-df.Evib * hc_k / Tvib)) * df.Qrot).sum()
    # Todo: non Boltzmann case

    assert np.isclose(Q, Q2)

    if verbose:
        printm("Tested Q vs recomputed from (Qvib, Qrot) are the same: OK")

    return True


# @pytest.mark.needs_config_file
# @pytest.mark.fast
# @pytest.mark.needs_db_CDSD_HITEMP
# def test_recompute_Q_from_QvibQrot_CDSD_PCN(verbose=True, warnings=True, *args, **kwargs):
#    '''
#    Calculate vibrational and rotational partition functions:
#
#    - in CDSD with (p,c,N) convention for vibrational levels
#    - under nonequilibrium
#
#    Calculate total rovibrational partition function, and compare
#
#    Test if partition function can be recomputed correctly from vibrational
#    populations and rotational partition function (note that we are in a coupled
#    case so partition function is not simply the product of Qvib, Qrot)
#    '''
#
#    from radis.misc.config import getDatabankEntries
#
#    iso = 1
#
#    try:
#        energies = getDatabankEntries('CDSD-HITEMP-PCN')['levels']
#        levelsfmt = getDatabankEntries('CDSD-HITEMP-PCN')['levelsfmt']
#
#        Tvib = 1500
#        Trot = 300
#
#        Qf = PartFuncCO2_CDSDcalc(energies[iso],
#                                  isotope=iso,
#                                  use_cached=True,
#                                  levelsfmt=levelsfmt)
#        Q = Qf.at_noneq(Tvib, Trot)
#        _, Qvib, dfQrot = Qf.at_noneq(Tvib, Trot, returnQvibQrot=True)
#        if verbose:
#            printm('Q', Q)
#        if verbose:
#            printm('Qvib', Qvib)
#
#        # 1) Test Q vs Q recomputed from Qrot, Qvib
#
#        # Recompute Qtot
#        df = dfQrot
#        Q2 = ((df.gvib*exp(-df.Evib*hc_k/Tvib))*df.Qrot).sum()
#        # Todo: non Boltzmann case
#
#        assert np.isclose(Q, Q2)
#
#        if verbose:
#            printm('Tested Q vs recomputed from (Qvib, Qrot) are the same: OK')
#
#        return True
#
#    except DatabankNotFound as err:
#        assert IgnoreMissingDatabase(err, __file__, warnings)


@pytest.mark.needs_config_file
@pytest.mark.fast
@pytest.mark.needs_db_CDSD_HITEMP
def test_recompute_Q_from_QvibQrot_CDSD_PC(
    verbose=True, warnings=True, *args, **kwargs
):
    """
    Calculate vibrational and rotational partition functions:

    - in CDSD with (p,c) convention for vibrational levels
    - under nonequilibrium

    Recompute total partition function, and compare

    Test if partition function can be recomputed correctly from vibrational
    populations and rotational partition function (note that we are in a coupled
    case so partition function is not simply the product of Qvib, Qrot)
    """

    from radis.misc.config import getDatabankEntries

    iso = 1

    energies = getDatabankEntries("CDSD-HITEMP-PC")["levels"]
    levelsfmt = getDatabankEntries("CDSD-HITEMP-PC")["levelsfmt"]

    Tvib = 1500
    Trot = 300

    Qf = PartFuncCO2_CDSDcalc(
        energies[iso], isotope=iso, use_cached=True, levelsfmt=levelsfmt
    )
    Q = Qf.at_noneq(Tvib, Trot)
    _, Qvib, dfQrot = Qf.at_noneq(Tvib, Trot, returnQvibQrot=True)
    if verbose:
        printm("Q", Q)
    if verbose:
        printm("Qvib", Qvib)

    # 1) Test Q vs Q recomputed from Qrot, Qvib

    # Recompute Qtot
    df = dfQrot
    Q2 = ((df.gvib * exp(-df.Evib * hc_k / Tvib)) * df.Qrot).sum()
    # Todo: non Boltzmann case

    assert np.isclose(Q, Q2)

    if verbose:
        printm("Tested Q vs recomputed from (Qvib, Qrot) are the same: OK")

    return True


# @pytest.mark.fast            # (fast only once the cached database has been generated)


def test_Q_1Tvib_vs_Q_3Tvib(T=1500, verbose=True, warnings=True, *args, **kwargs):
    """Test if partition function calculated in 1-Tvib mode returns the same
    result as partition function calculated in 3-Tvib mode
    """

    b = True

    # input
    M = "CO2"
    I = 1  # isotope
    S = Molecules[M][I]["X"]

    Qf = PartFunc_Dunham(S, use_cached=True)
    df = Qf.df

    # First make sure energies match
    if not (df.Evib == df.Evib1 + df.Evib2 + df.Evib3).all():
        b *= False
        if warnings:
            printm("WARNING in test_Q_1Tvib_vs_Q_3Tvib: Evib != Evib1 + Evib2 + Evib3")

    # Then calculate Q vs Q3T
    Q = Qf.at_noneq(T, T)
    if verbose:
        printm("Q", Q)

    Q3T = Qf.at_noneq_3Tvib((T, T, T), T)
    if verbose:
        printm("Q3T", Q3T)

    assert np.isclose(Q, Q3T)
    if verbose:
        printm("Tested Q in 1-Tvib vs Q in 3-Tvib modes (T={0:.0f}K): OK".format(T))

    return True


@pytest.mark.fast
def test_Morse_Potential_effect_CO(
    T=3000, rtol=1e-4, verbose=True, warnings=True, *args, **kwargs
):
    """Quantify effect of calculating upper levels near dissociation limit
    with Morse Potential

    Returns True if difference is less than rtol
    """

    vmax = 11
    vmax_morse = 48
    jmax = 300
    iso = 1
    S = Molecules["CO"][iso]["X"]

    # TODO: remove vmax, vmax_morse from Dunham method. Only use the one in ElecState
    db = PartFunc_Dunham(S, vmax=vmax, vmax_morse=0, Jmax=jmax, use_cached=False)
    Q_nomorse = db.at(T)

    db = PartFunc_Dunham(
        S, vmax=vmax, vmax_morse=vmax_morse, Jmax=jmax, use_cached=False
    )
    Q_morse = db.at(T)

    if verbose:
        printm("Morse vs no Morse potential (T={0}K)".format(T))
        printm("Q_morse: {0:.3f}".format(Q_morse))
        printm("Q_nomorse: {0:.3f}".format(Q_nomorse))
        printm("Difference: {0:.4f}%".format(abs(Q_nomorse - Q_morse) / Q_morse * 100))

    assert abs(Q_nomorse - Q_morse) / Q_morse < rtol


def test_levels_regeneration(verbose=True, warnings=True, *args, **kwargs):
    """Test that we regenerate levels file automatically if manually changed

    see https://github.com/radis/radis/issues/90

    """

    # from warnings import catch_warnings, filterwarnings
    def run_example():

        from radis.test.utils import (
            define_Evib_as_sum_of_Evibi,
            discard_lines_with_na_levels,
        )

        setup_test_line_databases(
            verbose=True
        )  # add HITEMP-CO2-HAMIL-TEST in ~/radis.json if not there

        sf = SpectrumFactory(
            wavenum_min=2283.7,
            wavenum_max=2285.1,
            wstep=0.001,
            cutoff=1e-30,
            path_length=0.1,
            mole_fraction=400e-6,
            isotope=[1],
            verbose=2,
        )
        sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
        sf.load_databank(
            "HITEMP-CO2-HAMIL-TEST",
            load_energies=True,  # Loading Energy Database
            db_use_cached=True,  # important to test Cache file here
        )

        # Now generate vibrational energies for a 2-T model
        levels = sf.parsum_calc["CO2"][1]["X"].df
        define_Evib_as_sum_of_Evibi(levels)
        discard_lines_with_na_levels(sf)

        # Calculate populations using the non-equilibrium module:
        # This will crash the first time because the Levels Database is just a fragment and does not include all levels.
        try:
            sf.non_eq_spectrum(300, 300)
        except AssertionError:  # expected
            sf.df0.dropna(inplace=True)

        s = sf.non_eq_spectrum(300, 300)
        s.plot()

    # run calculation (SpectrumFactory)
    run_example()
    # get the time when .levels file was last modified
    levels_last_modification = getmtime(
        getTestFile(r"co2_cdsd_hamiltonian_fragment.levels")
    )
    # get the time when .levels.h5 file was last modified
    cache_last_modification = getmtime(
        getTestFile(r"co2_cdsd_hamiltonian_fragment.levels.h5")
    )

    # change the time when .levels file was last modified
    stinfo = os.stat(getTestFile(r"co2_cdsd_hamiltonian_fragment.levels"))
    access_time = stinfo.st_atime
    os.utime(
        getTestFile(r"co2_cdsd_hamiltonian_fragment.levels"),
        (levels_last_modification + 1, access_time + 1),
    )

    # check if the change was successful
    levels_last_modification_again = getmtime(
        getTestFile(r"co2_cdsd_hamiltonian_fragment.levels")
    )

    assert levels_last_modification_again > levels_last_modification

    # run calculations once again to see if the .levels.h5 (cache file) is regenerated
    run_example()

    cache_last_modification_again = getmtime(
        getTestFile(r"co2_cdsd_hamiltonian_fragment.levels.h5")
    )
    assert cache_last_modification_again > cache_last_modification


def test_tabulated_partition_functions(
    verbose=True, plot=True, rtol=1e-2, *args, **kwargs
):
    """Test just-in-time tabulated partition return the same results as
    full summation within 0.5%  (adjust value with ``rtol``)

    We compute up to 10,000 K (beyond the validity range of most spectroscopic
    parameters)

    Run with verbose=True to check the accuracy and calibrate the
    :py:attr:`radis.levels.partfunc.RovibParFuncCalculator.N_bins_scaling` function"""

    from radis.db.molecules import CO2_X_626
    from radis.levels.partfunc import PartFunc_Dunham, PartFuncHAPI

    Z_sum = PartFunc_Dunham(CO2_X_626, mode="full summation")
    Z_tab = PartFunc_Dunham(
        CO2_X_626, mode="tabulation", verbose=3 if verbose else False
    )

    # For reference, start by comparing full-summation partition functions with TIPS tabulated partition functions
    Z_tips = PartFuncHAPI(M=2, I=1)  # CO  # isotope
    for T in [296, 3000, 5000]:
        z1 = Z_sum.at(T)
        z2 = Z_tips.at(T)
        accuracy = max(z1, z2) / min(z1, z2) - 1
        # accuracy_dict[
        if verbose:
            print(
                "full-sum & TIPS tabulated partition function of CO2 at {0}K close within {1:.2%}".format(
                    T, accuracy
                )
            )
        assert (
            accuracy < 4 * rtol
        )  # doesnt make sense to be extremely accurate in our tabulation if the full-summation doesnt match TIPS anyway

    # accuracy_dict = {1:[], 2:[], 3:[]}  # accuracy as a function of number of temperatures (helps how to scale N_bins)

    def compare_partition_functions(Z1, Z2, *T):
        z1 = Z1(*T)
        z2 = Z2(*T)
        # accuracy_dict[
        if verbose:
            accuracy = max(z1, z2) / min(z1, z2) - 1
            print(
                "full-sum & jit-tabulated partition function of CO2 at {0}K close within {1:.2%}".format(
                    T, accuracy
                )
            )
        assert np.isclose(z1, z2, rtol=rtol)

    # Equilibrium (same < 0.2%)
    compare_partition_functions(Z_sum.at, Z_tab.at, 296)
    compare_partition_functions(Z_sum.at, Z_tab.at, 3000)
    compare_partition_functions(Z_sum.at, Z_tab.at, 5000)
    compare_partition_functions(Z_sum.at, Z_tab.at, 10000)

    # Nonequilibrium (same << 0.1%)

    #  ... Compare with Partition function computed from PartFunc_Dunham
    compare_partition_functions(Z_sum.at_noneq, Z_tab.at_noneq, 296, 296)
    compare_partition_functions(Z_sum.at_noneq, Z_tab.at_noneq, 1000, 300)
    compare_partition_functions(Z_sum.at_noneq, Z_tab.at_noneq, 1000, 3000)
    compare_partition_functions(Z_sum.at_noneq, Z_tab.at_noneq, 3000, 3000)
    compare_partition_functions(Z_sum.at_noneq, Z_tab.at_noneq, 5000, 5000)
    compare_partition_functions(Z_sum.at_noneq, Z_tab.at_noneq, 10000, 10000)

    # Nonequilibrium 3 Tvib (same << 0.1%)

    #  ... Compare with Partition function computed from PartFunc_Dunham
    compare_partition_functions(
        Z_sum.at_noneq_3Tvib, Z_tab.at_noneq_3Tvib, (296, 296, 296), 296
    )
    compare_partition_functions(
        Z_sum.at_noneq_3Tvib, Z_tab.at_noneq_3Tvib, (1000, 1000, 2000), 300
    )
    compare_partition_functions(
        Z_sum.at_noneq_3Tvib, Z_tab.at_noneq_3Tvib, (1000, 1000, 3500), 3000
    )
    compare_partition_functions(
        Z_sum.at_noneq_3Tvib, Z_tab.at_noneq_3Tvib, (3000, 3000, 3000), 3000
    )
    compare_partition_functions(
        Z_sum.at_noneq_3Tvib, Z_tab.at_noneq_3Tvib, (5000, 5000, 5000), 5000
    )
    compare_partition_functions(
        Z_sum.at_noneq_3Tvib, Z_tab.at_noneq_3Tvib, (10000, 10000, 10000), 10000
    )

    # ... change Grid :
    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        import matplotlib.pyplot as plt

        plt.ion()

        Tvib = 3000
        Trot_arr = np.linspace(300, 3000, 10)
        plt.figure()
        plt.plot(Trot_arr, [Z_sum.at_noneq(Tvib, T) for T in Trot_arr])
        for N_bins in [10, 100, 1000]:
            Z_tab.N_bins == N_bins
            plt.plot(
                Trot_arr,
                [Z_tab.at_noneq(Tvib, T) for T in Trot_arr],
                "--",
                label=N_bins,
            )
        plt.legend()


def test_parsum_mode_in_factory(verbose=True, plot=True, *args, **kwargs):
    """Test Partition function modes in SpectrumFactory

    using :py:meth:`~radis.spectrum.spectrum.Spectrum.print_perf_profile`

    ::
        # sf.params.parsum_mode = 'full summation'   # default

        full summation profiler :
            spectrum_calculation      0.197s ████████████████
                check_line_databank             0.002s
                check_non_eq_param              0.000s
                reinitialize                    0.006s
                    copy_database                   0.003s
                    memory_usage_warning            0.003s
                    reset_population                0.000s
                calc_noneq_population           0.119s █████████
                    part_function                   0.020s █
                    population                      0.099s ████████
                scaled_non_eq_linestrength      0.002s
                    map_part_func                   0.000s
                    corrected_population_se         0.002s
                calc_emission_integral          0.012s
                calc_lineshift                  0.000s
                calc_hwhm                       0.013s █
                generate_wavenumber_arrays      0.000s
                calc_line_broadening            0.040s ███
                    precompute_LDM_lineshapes       0.002s
                    LDM_Initialized_vectors         0.000s
                    LDM_closest_matching_line       0.000s
                    LDM_Distribute_lines            0.001s
                    LDM_convolve                    0.037s ██
                calc_other_spectral_quan        0.004s
                generate_spectrum_obj           0.000s

    ::
        # sf.params.parsum_mode = 'tabulation'

        tabulation profiler :
            spectrum_calculation      0.104s ████████████████
                check_line_databank             0.002s
                check_non_eq_param              0.001s
                reinitialize                    0.004s
                    copy_database                   0.001s
                    memory_usage_warning            0.003s
                    reset_population                0.000s
                calc_noneq_population           0.026s ███
                    part_function                   0.019s ██
                    population                      0.007s █
                scaled_non_eq_linestrength      0.003s
                    map_part_func                   0.000s
                    corrected_population_se         0.003s
                calc_emission_integral          0.010s █
                calc_lineshift                  0.001s
                calc_hwhm                       0.017s ██
                generate_wavenumber_arrays      0.001s
                calc_line_broadening            0.035s █████
                    precompute_LDM_lineshapes       0.001s
                    LDM_Initialized_vectors         0.000s
                    LDM_closest_matching_line       0.001s
                    LDM_Distribute_lines            0.000s
                    LDM_convolve                    0.032s ████
                    others                          0.001s
                calc_other_spectral_quan        0.003s
                generate_spectrum_obj           0.000s
                others                          0.002s

    So calculation of populations is ~5x faster and the spectrum calculation itself
    is 2x faster.

    """
    from radis import SpectrumFactory

    wmin, wmax = 2284, 2285
    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        molecule="CO2",
        isotope="1",
        truncation=2.5,  # cm-1
        medium="air",
        path_length=10.32,  # cm
        wstep=0.001,
        verbose=0,
    )
    conditions = {
        "mole_fraction": 0.2,
        "pressure": 1,
        "Ttrans": 2000,
        "Trot": 2100,
        "Tvib": 1500,
    }

    # ... initialize tests :
    sf.params.parsum_mode = "full summation"  # default
    sf.load_databank("HITEMP-CO2-TEST")
    sf._add_bands()
    _ = sf.non_eq_spectrum(**conditions)  # initialize energies, etc.
    # ... compare:
    s_HITEMP = sf.non_eq_spectrum(**conditions, name=sf.params.parsum_mode)
    s_HITEMP.print_perf_profile()

    # Now, again with tabulation:
    # ... init :
    sf.params.parsum_mode = "tabulation"
    sf.load_databank("HITEMP-CO2-TEST")
    sf._add_bands()
    _ = sf.non_eq_spectrum(**conditions)  # initialize energies, first tabulation
    # ... compare:
    s_HITEMP2 = sf.non_eq_spectrum(**conditions, name=sf.params.parsum_mode)
    s_HITEMP2.print_perf_profile()

    if plot:
        import matplotlib.pyplot as plt

        plt.ion()
        from radis import plot_diff

        plot_diff(s_HITEMP2, s_HITEMP, method="ratio")

    from radis import get_residual

    assert get_residual(s_HITEMP, s_HITEMP2, "abscoeff") < 6e-5


def _run_testcases(verbose=True, warnings=True, *args, **kwargs):

    # Test 0: delete all cached energies
    test_delete_all_cached_energies(verbose=verbose, warnings=warnings)

    # Test 1: test cache mechanism
    test_cache_file_generation_and_update(verbose=verbose, warnings=warnings)

    # Test 2: compare calculated PartFunc for CO to HAPI
    test_calculatedQ_match_HAPI_CO()
    # Test 2b: compare for many molecules and isotopes
    test_calculatedQ_match_HAPI()

    # Test 3: compare calculated PartFunc to the tabulated one with CDSD
    # Test 4: compare calculated PartFunc to harcoded references
    test_CDSD_calc_vs_tab(verbose=verbose, warnings=warnings)
    test_CDSD_calc_vs_ref(warnings=warnings)
    test_reduced_CDSD_calc_vs_tab(verbose=verbose, warnings=warnings)
    test_reduced_CDSD_calc_noneq(verbose=verbose, warnings=warnings)

    # Test 5a, 5b: recompute Q from QvibQrot
    test_recompute_Q_from_QvibQrot_Dunham_Evib123_Erot(
        verbose=verbose, warnings=warnings
    )
    test_recompute_Q_from_QvibQrot_Dunham_Evib3_Evib12Erot(
        verbose=verbose, warnings=warnings
    )
    test_recompute_Q_from_QvibQrot_CDSD_PC(verbose=verbose, warnings=warnings)
    # test_recompute_Q_from_QvibQrot_CDSD_PCN(verbose=verbose, warnings=warnings)  # ignore in released version

    # Test 6:
    test_Q_1Tvib_vs_Q_3Tvib(verbose=verbose, warnings=warnings)

    # Test 7:
    test_Morse_Potential_effect_CO(verbose=verbose, warnings=warnings)

    # Test 8: Regenerates levels file if it's manually changed
    test_levels_regeneration(verbose=verbose, warnings=True, *args, **kwargs)

    # Test 9 : tabulation
    test_tabulated_partition_functions(verbose=verbose, *args, **kwargs)
    test_parsum_mode_in_factory(verbose=verbose, *args, **kwargs)

    return True


if __name__ == "__main__":
    # printm("Testing parfunc: {0}".format(_run_testcases()))
    printm("Testing partfunc.py:", pytest.main(["test_partfunc.py", "--pdb"]))
