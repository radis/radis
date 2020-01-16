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

from __future__ import absolute_import, unicode_literals, division, print_function
from radis.levels.partfunc import PartFunc_Dunham, PartFuncHAPI
from radis.levels.partfunc_cdsd import PartFuncCO2_CDSDtab, PartFuncCO2_CDSDcalc
from radis.phys.constants import hc_k
from radis.misc.utils import DatabankNotFound
from radis.test.utils import IgnoreMissingDatabase
from radis.misc.printer import printm
from radis.misc.cache_files import DeprecatedFileError
from radis.db.molecules import Molecules
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp
import os
from os.path import basename, exists
import pytest

fig_prefix = basename(__file__) + ": "


# %% Test

# never add @pytest.mark.fast so we don't delete cached files for 'fast' tests
def test_delete_all_cached_energies(verbose=True, warnings=True, *args, **kwargs):
    """ Doesnt really test anything, but cleans all cached energy levels """

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
    """ Test that cache file process works correctly, using CO as an example.

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
    with pytest.raises(DeprecatedFileError):  # DeprecatedFileError is expected
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
    """ Test 1: compare calculated PartFunc to the tabulated one """

    from radis.misc.config import getDatabankEntries

    try:

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

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


@pytest.mark.fast
def test_calculatedQ_match_HAPI_CO(
    vmax=11, jmax=300, plot=False, verbose=True, *args, **kwargs
):
    """ Tested that Q ab_initio (Dunham) match HAPI for CO at different temperatures"""

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

    hapi = PartFuncHAPI(M=5, I=1,)  # CO  # isotope

    us = []
    hap = []
    T = np.linspace(300, 3000)
    for Ti in T:
        us.append(db.at(Ti))
        hap.append(hapi.at(Ti))

    if plot:
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
    """ Tested that Q ab_initio (Dunham) match HAPI for different molecules 
    and isotopes"""

    # molecule, isotope, temperature, absolute tolerance
    for molecule, iso, T, atol in [
        ("CO2", 1, 300, 0.9),
        ("CO2", 1, 1000, 2.0),
        ("CO2", 1, 3000, 2000),
        ("CO2", 2, 300, 1.8),
        ("CO2", 2, 1000, 25),
        ("CO2", 2, 3000, 4962),
        ("CO", 1, 300, 0.16),
        ("CO", 1, 1000, 1.9),
        ("CO", 1, 3000, 27),
        ("CO", 2, 300, 0.31),
        ("CO", 2, 1000, 4.1),
        ("CO", 2, 3000, 56.9),
        ("CO", 3, 300, 0.16),
        ("CO", 3, 1000, 2.1),
        ("CO", 3, 3000, 28.6),
    ]:

        S = Molecules[molecule][iso]["X"]

        # Dont use cached: force recalculating
        db = PartFunc_Dunham(S)

        from radis.io.hitran import get_molecule_identifier

        hapi = PartFuncHAPI(M=get_molecule_identifier(molecule), I=iso)

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
    """ Test partition functions calculated with CDSD energy levels against 
    hardcoded values """

    from radis.misc.config import getDatabankEntries

    iso = 1

    try:
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
            Qf.at_noneq(3000, 3000, overpopulation={"(0,1)": 3}, returnQvibQrot=True)[
                0
            ],
            120053.34252537244,
            rtol=0.001,
        )

        if verbose:
            printm("Tested Q_CDSD values are correct : OK")

        return True

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


# @pytest.mark.fast   # (very fast only once the cached database has been generated, else decently fast)
def test_recompute_Q_from_QvibQrot_Dunham_Evib3_Evib12Erot(
    verbose=True, warnings=True, *args, **kwargs
):
    """     Calculate vibrational and rotational partition functions:
        
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
    """     Calculate vibrational and rotational partition functions:
        
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

    try:
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

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


# @pytest.mark.fast            # (fast only once the cached database has been generated)


def test_Q_1Tvib_vs_Q_3Tvib(T=1500, verbose=True, warnings=True, *args, **kwargs):
    """ Test if partition function calculated in 1-Tvib mode returns the same
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
    """ Quantify effect of calculating upper levels near dissociation limit
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


# from warnings import catch_warnings, filterwarnings


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

    return True


if __name__ == "__main__":
    printm("Testing parfunc: {0}".format(_run_testcases()))

    verbose = True
    warnings = True

#    test_CDSD_calc_vs_tab(verbose=verbose, warnings=warnings)
#    test_recompute_Q_from_QvibQrot_CDSD_PC(verbose=verbose, warnings=warnings)
#    test_recompute_Q_from_QvibQrot_CDSD_PCN(verbose=verbose, warnings=warnings)
