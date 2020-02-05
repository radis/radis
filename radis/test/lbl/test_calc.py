# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 09:36:25 2017

@author: erwan

Examples
--------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)

>>> pytest -m fast
          

-------------------------------------------------------------------------------

  
"""

from __future__ import unicode_literals, print_function, absolute_import, division
from radis.lbl import SpectrumFactory, LevelsList
from radis.lbl.calc import calc_spectrum
from radis.misc.printer import printm
from radis.phys.blackbody import sPlanck
from radis.misc.utils import DatabankNotFound
from radis.test.utils import IgnoreMissingDatabase
import numpy as np
from os.path import basename
import pytest

fig_prefix = basename(__file__) + ": "

# %% Tests
# ----------------


@pytest.mark.fast
def test_sPlanck_conversions(verbose=True, *args, **kwargs):

    if verbose:
        printm("Testing sPlanck conversions: ")

    s_cm = sPlanck(1000, 10000, T=1500, eps=0.3)
    I_cm2cm = s_cm.get("radiance_noslit", Iunit="mW/cm2/sr/cm_1")[1]
    I_cm2nm = s_cm.get("radiance_noslit", Iunit="mW/cm2/sr/nm")[1]

    s_nm = sPlanck(1000, 10000, T=1500, eps=0.3)
    I_nm2nm = s_nm.get("radiance_noslit", Iunit="mW/cm2/sr/nm")[1]
    I_nm2cm = s_nm.get("radiance_noslit", Iunit="mW/cm2/sr/cm_1")[1]

    assert np.allclose(I_cm2cm, I_nm2cm)
    assert np.allclose(I_nm2nm, I_cm2nm)


# @pytest.mark.needs_config_file
# @pytest.mark.needs_db_HITEMP_CO2_DUNHAM
@pytest.mark.needs_connection
def test_calc_spectrum(verbose=True, plot=True, warnings=True, *args, **kwargs):
    """ Basic example, used as a non-regression test

    Notes
    -----
    
    How long it tooks to calculate this Spectrum?

    Performance test on old NeQ package, with the [CDSD-HITEMP-JMIN] databank.
    See the caveats in the E. Pannier "Limits of CO2 NonEquilibrium Models" paper.
    (just used here as a performance monitoring)
    
    - neq 0.9.20: 18.7s

    - neq 0.9.20*: 15.4s   (removed 2nd loop on 1st isotope because of groupby().apply())

    - neq 0.9.20**: 11.7s  (after replacing fill_Evib with map() ) 

    - neq 0.9.21: 9.4s     (improve Qrot / nrot fetching performance) 

    - neq 0.9.22: 8.4s     
    
    Starting from RADIS 1.0.1, the test is run on [HITRAN-2016]_, which 
    is not valid for these temperatures but can be more conveniently 
    downloaded automatically and thus executed everytime with [Travis]_
    
    (we also expect the test to be much faster than above, but that's just 
    because the database is smaller!)
    
    - radis 0.9.20 : 2.49 s    on [HITRAN-2016]
                     4.05 s    on [CDSD-HITEMP-JMIN]
    
    """

    if verbose:
        printm("Testing calc_spectrum match reference")

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        import matplotlib.pyplot as plt

        plt.ion()

    try:
        s = calc_spectrum(
            wavelength_min=4165,
            wavelength_max=4200,
            #                          databank='CDSD-HITEMP-JMIN',
            databank="fetch",  # not appropriate for these temperatures, but convenient for automatic testing
            Tgas=300,
            Tvib=1700,
            Trot=1550,
            path_length=0.1,
            mole_fraction=0.5,
            molecule="CO2",
            isotope="1,2",
            wstep=0.01,
            cutoff=1e-25,
            use_cached=True,
            medium="vacuum",
            verbose=verbose,
            warnings={
                "MissingSelfBroadeningWarning": "ignore",
                "NegativeEnergiesWarning": "ignore",
                "HighTemperatureWarning": "ignore",
            },
        )
        s.apply_slit((2, 2.5), "nm", shape="trapezoidal")

        if plot:
            s.plot(wunit="nm")

        w, I = s.get("radiance", wunit="nm")
        w_ref = w[::100]
        # Compare against hardcoded results (neq 0.9.22, 28/06/18)
        #        I_ref = np.array([0.28694463, 0.29141711, 0.32461613, 0.32909566, 0.21939511, 0.18606445,
        #                          0.19740763, 0.16948599, 0.16780345, 0.15572173, 0.16770853, 0.14966064,
        #                          0.13041356, 0.11751016, 0.10818072, 0.11592531, 0.04666677, 0.00177108,
        #                          0.00069339])
        # Harcoded results changed for RADIS  with the change of
        # database (HITEMP-2010 -> HITRAN-2016) and of Tvib model
        # CDSD with (P,C,Jmin,N) in CDSD polyad -> RADIS built-in constants)
        #        I_ref = np.array([ 0.29148768,  0.29646856,  0.32999337,  0.32249701,  0.2078451 ,
        #                          0.18974631,  0.2019285 ,  0.17346687,  0.17211401,  0.15939359,
        #                          0.17240575,  0.15395179,  0.13374185,  0.11997065,  0.10858693,
        #                          0.11114162,  0.04575873,  0.00163863,  0.00062654])
        # Updated again in RADIS 0.9.20 (19/08/19) to account for the use of DLM (difference
        # not significant)
        #        I_ref = np.array([ 0.29060991,  0.29756722,  0.32972058,  0.3206278 ,  0.20696867,
        #                           0.19218358,  0.20155747,  0.17336405,  0.17218653,  0.1589136 ,
        #                           0.17110649,  0.15403513,  0.13376804,  0.11932659,  0.10882006,
        #                           0.11112725,  0.0458288 ,  0.00247956,  0.00144128])
        # Updated again in RADIS 0.9.20 (02/09/19) with switch to tabulated Q(Tref)
        I_ref = np.array(
            [
                0.29048064,
                0.29743104,
                0.32955513,
                0.32047172,
                0.20688813,
                0.19210952,
                0.20148265,
                0.17330909,
                0.17213373,
                0.15887159,
                0.17106096,
                0.15400039,
                0.13374285,
                0.11930822,
                0.10880631,
                0.11111394,
                0.04582291,
                0.00247955,
                0.00144128,
            ]
        )
        if plot:
            plt.plot(w_ref, I_ref, "or", label="ref")
            plt.legend()
        assert np.allclose(I[::100], I_ref, atol=1e-6)

        return True

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


# @pytest.mark.needs_db_CDSD_HITEMP_PCN
@pytest.mark.needs_connection
def test_calc_spectrum_overpopulations(
    verbose=True, plot=False, warnings=True, *args, **kwargs
):
    """ Non-regression test: 
        
    Example using overpopulation of the 001 asymmetric stretch first level of CO2,
    which is written (p,c,N) = (3,1,4) in [CDSD-4000]_ notation
    
    Notes
    -----
    
    In old Neq package (before RADIS):
    
    the test uses a CDSD-PCN notation for vibrational energy assignation, i.e,
    Evib = minimal energy of a (p,c,N) polyad. See the discussion on the implications
    in the E. Pannier "Limits of CO2 NonEquilibrium models" paper. 
    Better use the assignation scheme suggested in the paper. 
    But it's okay here as a non-regression test.
    
    Starting from RADIS 1.0.1, the test is run on [HITRAN-2016]_, which 
    is not valid for these temperatures but can be more conveniently 
    downloaded automatically and thus executed everytime with [Travis]_
    
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        import matplotlib.pyplot as plt

        plt.ion()

    try:
        s = calc_spectrum(
            wavelength_min=4165,
            wavelength_max=4200,
            #                          databank='CDSD-HITEMP-PCN',
            databank="fetch",  # not appropriate for these temperatures, but convenient for automatic testings
            Tgas=300,
            Tvib=1700,
            Trot=1550,
            #                          overpopulation={'(3,1,4)': 3},  # 00'0'1 in spectroscopic notation
            overpopulation={"(0,0,0,1)": 3},  # 00'0'1 in spectroscopic notation
            path_length=0.1,
            mole_fraction=0.5,
            molecule="CO2",
            isotope="1,2",
            wstep=0.01,
            cutoff=1e-25,
            use_cached=True,
            medium="vacuum",
            verbose=verbose,
            warnings={
                "MissingSelfBroadeningWarning": "ignore",
                "NegativeEnergiesWarning": "ignore",
                "HighTemperatureWarning": "ignore",
            },
        )
        s.apply_slit((2, 2.5), "nm", shape="trapezoidal")

        if plot:
            s.plot()

        w, I = s.get("radiance", wunit="nm")
        w_ref = w[::100]
        # Compare against hardcoded results (neq 0.9.22, 28/06/18)
        #        I_ref = np.array([0.61826008, 0.65598262, 0.79760003, 0.7958013 , 0.5792486 ,
        #                          0.56727691, 0.60361258, 0.51549598, 0.51012651, 0.47133131,
        #                          0.50770568, 0.45093953, 0.39129824, 0.35125324, 0.32238316,
        #                          0.34542781, 0.13908073, 0.00506012, 0.00189535])
        # Harcoded results changed for RADIS v1.0.1  with the change of
        # database (HITEMP-2010 -> HITRAN-2016) and of Tvib model
        # CDSD with (P,C,Jmin,N) in CDSD polyad -> RADIS built-in constants)
        #
        #        I_ref = np.array([ 0.62299838,  0.66229013,  0.81037059,  0.79899315,  0.57215806,
        #                          0.57626389,  0.61424273,  0.52454807,  0.5200812 ,  0.47920924,
        #                          0.51843533,  0.46058817,  0.3983277 ,  0.35582979,  0.32095204,
        #                          0.32821575,  0.13525543,  0.00469489,  0.00174166])
        # Updated again in RADIS 0.9.20 (16/08/19) to account for the use of DLM (difference
        # not significant)
        #        I_ref = np.array([ 0.62134142,  0.66722021,  0.81016539,  0.79387937,  0.56974945,
        #                           0.58280035,  0.6120114 ,  0.52319075,  0.5193041 ,  0.47686282,
        #                           0.51374777,  0.46022548,  0.3979033 ,  0.3534643 ,  0.32129239,
        #                           0.32786479,  0.1351593 ,  0.0068877 ,  0.00387545])
        # Updated again in RADIS 0.9.20 (02/09/19) with switch to tabulated Q(Tref)
        I_ref = np.array(
            [
                0.62109562,
                0.66695661,
                0.80983176,
                0.79356445,
                0.56958189,
                0.58264143,
                0.61185167,
                0.52307454,
                0.51919288,
                0.47677519,
                0.51365307,
                0.46015383,
                0.39785172,
                0.35342697,
                0.32126465,
                0.32783797,
                0.13514737,
                0.00688769,
                0.00387544,
            ]
        )
        if plot:
            plt.plot(w_ref, I_ref, "or", label="ref")
            plt.legend()
            s.plot_populations()

        assert np.allclose(I[::100], I_ref, atol=1e-6)

        if verbose:
            printm("Test overpopulations: OK")

        return True

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


@pytest.mark.needs_config_file
@pytest.mark.needs_db_CDSD_HITEMP_PC
# @pytest.mark.needs_connection
def test_all_calc_methods(
    verbose=True, plot=False, warnings=True, rtol=1e-3, *args, **kwargs
):
    """ Test same spectrum for 3 different calculation variants (equilibrium, 
    non-equilibrium, per band and recombine 
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        import matplotlib.pyplot as plt

        plt.ion()

    try:

        Tgas = 1500

        iso = 1
        sf = SpectrumFactory(
            wavelength_min=4170,
            wavelength_max=4175,
            mole_fraction=1,
            path_length=0.025,
            cutoff=1e-25,
            molecule="CO2",
            isotope=iso,
            db_use_cached=True,
            lvl_use_cached=True,
            verbose=verbose,
        )
        sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
        sf.warnings["NegativeEnergiesWarning"] = "ignore"
        sf.warnings["HighTemperatureWarning"] = "ignore"
        #        sf.fetch_databank()   # uses HITRAN: not really valid at this temperature, but runs on all machines without install
        sf.load_databank("CDSD-HITEMP-PC")

        s_bands = sf.non_eq_bands(Tvib=Tgas, Trot=Tgas)
        lvl = LevelsList(sf.parsum_calc["CO2"][iso]["X"], s_bands, sf.params.levelsfmt)
        s_bd = lvl.non_eq_spectrum(Tvib=Tgas, Trot=Tgas)

        s_nq = sf.non_eq_spectrum(Tvib=Tgas, Trot=Tgas)
        s_eq = sf.eq_spectrum(Tgas=Tgas)

        #
        if plot:
            fig = plt.figure(fig_prefix + "Compare all calc methods")
            s_bd.plot(nfig=fig.number, color="b", lw=5, label="from bands code")
            s_nq.plot(nfig=fig.number, lw=3, label="non eq code")
            s_eq.plot(nfig=fig.number, lw=2, color="r", label="equilibrum code")
            plt.legend()

        assert np.isclose(s_bd.get_power(), s_nq.get_power(), rtol=rtol)
        assert np.isclose(s_bd.get_power(), s_eq.get_power(), rtol=rtol)

        if verbose:
            printm(
                "Eq == non-eq:\t",
                np.isclose(s_eq.get_power(), s_nq.get_power(), rtol=rtol),
            )
            printm(
                "Bands == Non-eq:\t",
                np.isclose(s_bd.get_power(), s_nq.get_power(), rtol=rtol),
            )

        if verbose:
            printm("Test all methods comparison: OK")

        return True

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


@pytest.mark.needs_connection
def test_eq_vs_noneq_isotope(verbose=True, plot=False, warnings=True, *args, **kwargs):
    """ Test same spectrum for 2 different calculation codes (equilibrium, 
    non-equilibrium) in the presence of isotopes 
    
    Notes
    -----
    
    On the old NeQ package the test used [HITEMP-2010]_
    
    Starting from RADIS 1.0.1, the test is run on [HITRAN-2016]_, which 
    is not valid for these temperatures but can be more conveniently 
    downloaded automatically and thus executed everytime with [Travis]_
    
    """

    try:
        Tgas = 1500

        sf = SpectrumFactory(
            wavelength_min=4250,
            wavelength_max=4350,
            mole_fraction=1,
            path_length=1,
            cutoff=1e-25,
            molecule="CO2",
            isotope="1,2",
            db_use_cached=True,
            verbose=verbose,
        )
        sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
        sf.warnings["NegativeEnergiesWarning"] = "ignore"
        sf.warnings["HighTemperatureWarning"] = "ignore"
        sf.fetch_databank()  # uses HITRAN: not really valid at this temperature, but runs on all machines without install
        #        sf.load_databank('HITEMP-CO2-DUNHAM')
        s_nq = sf.non_eq_spectrum(Tvib=Tgas, Trot=Tgas, name="Non-eq")
        s_eq = sf.eq_spectrum(Tgas=Tgas, name="Eq")

        rtol = 5e-3  # 2nd isotope calculated with placeholder energies
        match_eq_vs_non_eq = s_eq.compare_with(
            s_nq, spectra_only="abscoeff", rtol=rtol, plot=plot
        )
        match_eq_vs_non_eq *= s_eq.compare_with(
            s_nq, spectra_only="radiance_noslit", rtol=rtol, plot=plot
        )

        if verbose:
            printm(
                "Tested eq vs non-eq (<{0:.1f}% error) with isotopes: {1}".format(
                    rtol * 100, bool(match_eq_vs_non_eq)
                )
            )

        assert match_eq_vs_non_eq

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


def _run_testcases(plot=False, verbose=True, warnings=True, *args, **kwargs):

    # Test sPlanck and conversion functions
    test_sPlanck_conversions()

    # Test calc_spectrum function
    test_calc_spectrum()

    # Test calc_spectrum with overpopulation
    test_calc_spectrum_overpopulations(
        verbose=verbose, plot=plot, warnings=warnings, *args, **kwargs
    )

    # Compare all calc methods
    test_all_calc_methods(
        verbose=verbose, plot=plot, warnings=warnings, *args, **kwargs
    )

    # Compare same spectrum with two calculation methods
    test_eq_vs_noneq_isotope(
        verbose=verbose, plot=plot, warnings=warnings, *args, **kwargs
    )

    return True


# --------------------------
if __name__ == "__main__":

    printm("Testing calc.py: ", _run_testcases(verbose=True))
