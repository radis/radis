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

from os.path import basename

import numpy as np
import pytest

from radis.lbl import LevelsList, SpectrumFactory
from radis.lbl.calc import calc_spectrum
from radis.misc.printer import printm
from radis.phys.blackbody import sPlanck

fig_prefix = basename(__file__) + ": "

# %% Tests
# ----------------


@pytest.mark.fast
def test_sPlanck_conversions(verbose=True, *args, **kwargs):

    if verbose:
        printm("Testing sPlanck conversions: ")

    s_cm = sPlanck(1000, 10000, T=1500, eps=0.3)
    I_cm2cm = s_cm.get("radiance_noslit", Iunit="mW/cm2/sr/cm-1")[1]
    I_cm2nm = s_cm.get("radiance_noslit", Iunit="mW/cm2/sr/nm")[1]

    s_nm = sPlanck(1000, 10000, T=1500, eps=0.3)
    I_nm2nm = s_nm.get("radiance_noslit", Iunit="mW/cm2/sr/nm")[1]
    I_nm2cm = s_nm.get("radiance_noslit", Iunit="mW/cm2/sr/cm-1")[1]

    assert np.allclose(I_cm2cm, I_nm2cm)
    assert np.allclose(I_nm2nm, I_cm2nm)


# @pytest.mark.needs_config_file
# @pytest.mark.needs_db_HITEMP_CO2_DUNHAM
@pytest.mark.needs_connection
def test_calc_spectrum(verbose=True, plot=True, warnings=True, *args, **kwargs):
    """Basic example, used as a non-regression test.

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

    Starting from RADIS 1.0.1, the test is run on [HITRAN-2020]_, which
    is not valid for these temperatures but can be more conveniently
    downloaded automatically and thus executed everytime with `Travis CI <https://travis-ci.com/radis/radis>`_

    (we also expect the test to be much faster than above, but that's just
    because the database is smaller!)

    - radis 0.9.20 : 2.49 s    on [HITRAN-2020]
                     4.05 s    on [CDSD-HITEMP-JMIN]
    - radis 0.12.2 : 1.1s      on [HITRAN-2020]   (most time spent in the Evib, Erot look-up)

    """
    if verbose:
        printm("Testing calc_spectrum match reference")

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        import matplotlib.pyplot as plt

        plt.ion()

    s = calc_spectrum(
        wavelength_min=4165,
        wavelength_max=4200,
        # wavenum_min=1e7/4200,
        # wavenum_max=1e7/4165,
        #                          databank='CDSD-HITEMP-JMIN',
        databank="hitran",  # not appropriate for these temperatures, but convenient for automatic testing
        # databank="HITRAN-CO2-TEST",  # to use, but only has 1 isotope. TODO add new test file with 2 isotopes
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
        verbose=3,
        optimization="simple",
        broadening_method="fft",
        neighbour_lines=5,  # previously: broadening_max_width [FWHM] = 10
        warnings={
            "MissingSelfBroadeningWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "HighTemperatureWarning": "ignore",
        },
    )
    s.apply_slit((2, 2.5), "nm", shape="trapezoidal")

    if plot:
        s.plot(wunit="nm", Iunit="mW/cm2/sr/nm")

    w, I = s.get("radiance", wunit="nm", Iunit="mW/cm2/sr/nm")
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
    # Updated again in RADIS 0.9.20 (19/08/19) to account for the use of LDM (difference
    # not significant)
    #        I_ref = np.array([ 0.29060991,  0.29756722,  0.32972058,  0.3206278 ,  0.20696867,
    #                           0.19218358,  0.20155747,  0.17336405,  0.17218653,  0.1589136 ,
    #                           0.17110649,  0.15403513,  0.13376804,  0.11932659,  0.10882006,
    #                           0.11112725,  0.0458288 ,  0.00247956,  0.00144128])
    # Updated again in RADIS 0.9.20 (02/09/19) with switch to tabulated Q(Tref)
    # I_ref = np.array([         0.29048064, 0.29743104, 0.32955513, 0.32047172, 0.20688813,
    #                           0.19210952, 0.20148265, 0.17330909, 0.17213373, 0.15887159,
    #                           0.17106096, 0.15400039, 0.13374285, 0.11930822, 0.10880631,
    #                           0.11111394, 0.04582291, 0.00247955, 0.00144128])
    # Updated again on (05/08/20) with implementation of optimized weights:
    # I_ref = np.array([0.29043204, 0.29740738, 0.32954171, 0.32045394, 0.20680637,
    #                           0.19205883, 0.2014279 , 0.17322236, 0.17206767, 0.15879478,
    #                           0.17107564, 0.15400038, 0.13372559, 0.11929585, 0.10881116,
    #                           0.11111882, 0.04581152, 0.00247154, 0.00143631])
    # Updated with RADIS 0.9.26+ (13/12/20) and switch to latest HAPI version
    # therefore TIPS 2017 (instead of TIPS 2011) which changed CO2 partition functions
    # Q(Tref=300) 291.9025 --> 291.0405
    # I_ref = np.array([0.29057002, 0.29755271, 0.32971832, 0.32062057, 0.20689238,
    #       0.19213799, 0.20150788, 0.17328113, 0.17212415, 0.15883971,
    #       0.17112437, 0.15403756, 0.13375255, 0.1193155 , 0.10882585,
    #       0.11113302, 0.04581781, 0.00247155, 0.00143631])
    # Update on 17/06/21 (Radis 0.9.30) with the new HITRAN-2020 lines for CO2:
    # I_ref = np.array([0.29009315, 0.29732781, 0.32959567, 0.32061864, 0.20588816,
    #         0.19179586, 0.20101157, 0.17255937, 0.17163118, 0.15821026,
    #         0.17120068, 0.15403946, 0.13356296, 0.11907705, 0.10890721,
    #         0.11158052, 0.04583138, 0.00232573, 0.00136583])
    # Update on 05/11/21 (Radis 0.10.4) with use of tabulated instead of ab-initio
    # partition functions for Qref  (0.02% constant difference)
    #
    # Update on 14/11/21: after switching from HITRAN fetched by astroquery (partial range)
    # to HITRAN fetched by HAPI (full database) > atol=1e-6 fails, because of
    # numerical errors in how the two databases are stored. Switched to rtol=1e-3
    w_ref = np.array(
        [
            4197.60321744,
            4195.84148905,
            4194.08123884,
            4192.32246493,
            4190.56516548,
            4188.80933862,
            4187.05498252,
            4185.30209532,
            4183.55067518,
            4181.80072025,
            4180.05222871,
            4178.30519870,
            4176.55962842,
            4174.81551601,
            4173.07285967,
            4171.33165756,
            4169.59190786,
            4167.85360877,
            4166.11675846,
        ]
    )

    I_ref = np.array(
        [
            0.29003961,
            0.29727294,
            0.32953484,
            0.32055948,
            0.20585016,
            0.19176046,
            0.20097447,
            0.17252753,
            0.17159951,
            0.15818106,
            0.17116909,
            0.15401103,
            0.13353831,
            0.11905508,
            0.10888712,
            0.11155993,
            0.04582293,
            0.0023253,
            0.00136557,
        ]
    )

    if plot:
        plt.plot(w_ref, I_ref, "or", label="ref")
        plt.legend()

    # [71:-71] because of the shift introduced in 0.9.30 where convolved
    # arrays are not cropped; but np.nan values are added
    assert np.allclose(w[71:-71][::100], w_ref, atol=1e-6)
    assert np.allclose(I[71:-71][::100], I_ref, rtol=1e-3)

    return True


# @pytest.mark.needs_db_CDSD_HITEMP_PCN
@pytest.mark.needs_connection
def test_calc_spectrum_overpopulations(
    verbose=True, plot=False, warnings=True, *args, **kwargs
):
    """Non-regression test.

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

    Starting from RADIS 1.0.1, the test is run on [HITRAN-2020]_, which
    is not valid for these temperatures but can be more conveniently
    downloaded automatically and thus executed everytime with `Travis CI <https://travis-ci.com/radis/radis>`_

    """
    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        import matplotlib.pyplot as plt

        plt.ion()

    s = calc_spectrum(
        wavelength_min=4165,
        wavelength_max=4200,
        #                          databank='CDSD-HITEMP-PCN',
        databank="hitran",  # not appropriate for these temperatures, but convenient for automatic testings
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
        optimization="simple",
        broadening_method="fft",  # For this particular test case
        neighbour_lines=5,  # previously: broadening_max_width [FWHM] = 10
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
    # [71:-71] because of the shift introduced in 0.9.30 where convolved
    # arrays are not cropped; but np.nan values are added
    w_ref = w[71:-71][::100]
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
    # Updated again in RADIS 0.9.20 (16/08/19) to account for the use of LDM (difference
    # not significant)
    #        I_ref = np.array([ 0.62134142,  0.66722021,  0.81016539,  0.79387937,  0.56974945,
    #                           0.58280035,  0.6120114 ,  0.52319075,  0.5193041 ,  0.47686282,
    #                           0.51374777,  0.46022548,  0.3979033 ,  0.3534643 ,  0.32129239,
    #                           0.32786479,  0.1351593 ,  0.0068877 ,  0.00387545])
    # Updated again in RADIS 0.9.20 (02/09/19) with switch to tabulated Q(Tref)
    #       I_ref = np.array([  0.62109562,0.66695661,0.80983176,0.79356445,0.56958189,
    #                           0.58264143,0.61185167,0.52307454,0.51919288,0.47677519,
    #                           0.51365307,0.46015383,0.39785172,0.35342697,0.32126465,
    #                           0.32783797,0.13514737,0.00688769,0.00387544])
    # Updated again on (05/08/20) with implementation of optimized weights:

    #        I_ref = np.array([0.62097252, 0.66685971, 0.80982863, 0.7935332 , 0.56939115,
    #                       0.58255747, 0.61175655, 0.52287059, 0.51905438, 0.47659305,
    #                       0.51375266, 0.46019418, 0.39782806, 0.35340763, 0.32128853,
    #                       0.32785594, 0.13511584, 0.00686547, 0.00386199])

    # Updated with RADIS 0.9.26+ (13/12/20) and switch to latest HAPI version
    # therefore TIPS 2017 (instead of TIPS 2011) which changed CO2 partition functions
    # Q(Tref=300) 291.9025 --> 291.0405

    # I_ref = np.array([0.62123487, 0.667141  , 0.81018478, 0.79386946, 0.56957012,
    #    0.58272738, 0.61192736, 0.52299488, 0.51917337, 0.4766868 ,
    #    0.51385405, 0.46027087, 0.39788325, 0.35344755, 0.32131818,
    #    0.32788457, 0.13512857, 0.00686548, 0.00386199])

    # Update on 17/06/20 (Radis 0.9.30) with the new HITRAN-2020 lines for CO2
    # array([0.62009366, 0.66644973, 0.81033866, 0.79357758, 0.56690276,
    #        0.58248594, 0.61106305, 0.52129611, 0.51819065, 0.47520791,
    #        0.51455743, 0.46064655, 0.39761512, 0.35293948, 0.32169967,
    #        0.3293075 , 0.13525292, 0.00645308, 0.00366828])
    # Update on 05/11/21 (Radis 0.10.4) with use of tabulated instead of ab-initio
    # partition functions for Qref  (0.02% constant difference)
    #
    # Update on 14/11/21: after switching from HITRAN fetched by astroquery (partial range)
    # to HITRAN fetched by HAPI (full database) > atol=1e-6 fails, because of
    # numerical errors in how the two databases are stored. Switched to rtol=1e-3

    I_ref = np.array(
        [
            0.61997922,
            0.66632674,
            0.81018912,
            0.79343113,
            0.56679814,
            0.58237845,
            0.61095028,
            0.52119991,
            0.51809502,
            0.47512021,
            0.51446247,
            0.46056154,
            0.39754174,
            0.35287435,
            0.3216403,
            0.32924673,
            0.13522796,
            0.00645189,
            0.0036676,
        ]
    )
    if plot:
        plt.plot(w_ref, I_ref, "or", label="ref")
        plt.legend()
        s.plot_populations()

    # [71:-71] because of the shift introduced in 0.9.30 where convolved
    # arrays are not cropped; but np.nan values are added
    assert np.allclose(I[71:-71][::100], I_ref, rtol=1e-3)

    if verbose:
        printm("Test overpopulations: OK")

    return True


# @pytest.mark.needs_config_file
# @pytest.mark.needs_db_CDSD_HITEMP_PC
# @pytest.mark.needs_connection
# def test_all_calc_methods(
#    verbose=True, plot=False, warnings=True, rtol=1e-3, *args, **kwargs
# ):
#    """ Test same spectrum for 3 different calculation variants (equilibrium,
#    non-equilibrium, per band and recombine
#    """
#
#    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
#        import matplotlib.pyplot as plt
#
#        plt.ion()
#
#    Tgas = 1500
#
#    iso = 1
#    sf = SpectrumFactory(
#        wavelength_min=4170,
#        wavelength_max=4175,
#        mole_fraction=1,
#        path_length=0.025,
#        cutoff=1e-25,
#        molecule="CO2",
#        isotope=iso,
#        db_use_cached=True,
#        lvl_use_cached=True,
#        verbose=verbose,
#    )
#    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
#    sf.warnings["NegativeEnergiesWarning"] = "ignore"
#    sf.warnings["HighTemperatureWarning"] = "ignore"
#    #        sf.fetch_databank()   # uses HITRAN: not really valid at this temperature, but runs on all machines without install
#    sf.load_databank("CDSD-HITEMP-PC")
#
#    s_bands = sf.non_eq_bands(Tvib=Tgas, Trot=Tgas)
#    lvl = LevelsList(sf.parsum_calc["CO2"][iso]["X"], s_bands, sf.params.levelsfmt)
#    s_bd = lvl.non_eq_spectrum(Tvib=Tgas, Trot=Tgas)
#
#    s_nq = sf.non_eq_spectrum(Tvib=Tgas, Trot=Tgas)
#    s_eq = sf.eq_spectrum(Tgas=Tgas)
#
#    #
#    if plot:
#        fig = plt.figure(fig_prefix + "Compare all calc methods")
#        s_bd.plot(nfig=fig.number, color="b", lw=5, label="from bands code")
#        s_nq.plot(nfig=fig.number, lw=3, label="non eq code")
#        s_eq.plot(nfig=fig.number, lw=2, color="r", label="equilibrum code")
#        plt.legend()
#
#    assert np.isclose(s_bd.get_power(), s_nq.get_power(), rtol=rtol)
#    assert np.isclose(s_bd.get_power(), s_eq.get_power(), rtol=rtol)
#
#    if verbose:
#        printm(
#            "Eq == non-eq:\t",
#            np.isclose(s_eq.get_power(), s_nq.get_power(), rtol=rtol),
#        )
#        printm(
#            "Bands == Non-eq:\t",
#            np.isclose(s_bd.get_power(), s_nq.get_power(), rtol=rtol),
#        )
#
#    if verbose:
#        printm("Test all methods comparison: OK")
#
#    return True


def test_all_calc_methods_CO2pcN(
    verbose=True, plot=False, warnings=True, rtol=1e-3, *args, **kwargs
):
    """Test same spectrum for 3 different calculation variants (equilibrium,
    non-equilibrium, per band and recombine

    Uses CO2 Levels database where the energy partitioning is done as follow:

        2 nonequilibrium modes
        Evib is the minimum of a "p,c,N" group
        Erot = E - Evib

    This corresponds to the levelsfmt = 'cdsd-pcN' in
    :data:`~radis.lbl.loader.KNOWN_LVLFORMAT`
    """

    from radis.misc.config import getDatabankEntries
    from radis.test.utils import (
        define_Evib_as_min_of_polyad,
        discard_lines_with_na_levels,
        setup_test_line_databases,
    )

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        import matplotlib.pyplot as plt

        plt.ion()
    # %%
    Tgas = 500

    iso = 1
    sf = SpectrumFactory(
        wavenum_min=2284,
        wavenum_max=2285,
        truncation=2.5,  # TODO @EP: crashes with 0.15?
        mole_fraction=1,
        path_length=0.025,
        cutoff=1e-25,
        molecule="CO2",
        isotope=iso,
        verbose=verbose,
        export_lines=True,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.warnings["HighTemperatureWarning"] = "ignore"

    # Preparation:

    setup_test_line_databases()

    # Generate a Levels database with p,c,N energy partitioning
    # ... copy "HITEMP-CO2-HAMIL-TEST" info
    database_kwargs = getDatabankEntries("HITEMP-CO2-HAMIL-TEST")
    # ... adapt it to 'cdsd-pcN' mode
    del database_kwargs["info"]
    database_kwargs["levelsfmt"] = "cdsd-pcN"
    # ... load the new database
    sf.load_databank(**database_kwargs, load_energies=True, load_columns="noneq")

    # Now, define Evib:
    Q_calc = sf.parsum_calc["CO2"][1]["X"]
    Q_calc.df = define_Evib_as_min_of_polyad(Q_calc.df, keys=["p", "c", "N"])

    # With this Evib definition, clean the Lines database from where Evib is not defined
    # (because Levels does not exist in the reduced, test Level Database)
    discard_lines_with_na_levels(sf)

    # %%---------------------
    # Ready, let's start the tests:

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

    assert np.isclose(
        s_bd.get_integral("abscoeff"), s_nq.get_integral("abscoeff"), rtol=rtol
    )
    assert np.isclose(
        s_bd.get_integral("abscoeff"), s_eq.get_integral("abscoeff"), rtol=rtol
    )
    assert np.isclose(
        s_nq.get_integral("abscoeff"), s_eq.get_integral("abscoeff"), rtol=rtol
    )

    # TODO @EP: assertion fail in emission. This is due to the slight shift
    # in intensity also observed in the Planck test (test_base.py::test_optically_thick_limit_1iso()).
    #    assert np.isclose(s_bd.get_power(), s_nq.get_power(), rtol=rtol)
    #    assert np.isclose(s_bd.get_power(), s_eq.get_power(), rtol=rtol)
    #    assert np.isclose(s_nq.get_power(), s_eq.get_power(), rtol=rtol)

    return True


@pytest.mark.needs_connection
def test_eq_vs_noneq_isotope(verbose=True, plot=False, warnings=True, *args, **kwargs):
    """Test same spectrum for 2 different calculation codes (equilibrium,
    non-equilibrium) in the presence of isotopes

    Notes
    -----

    On the old NeQ package the test used [HITEMP-2010]_

    Starting from RADIS 1.0.1, the test is run on [HITRAN-2020]_, which
    is not valid for these temperatures but can be more conveniently
    downloaded automatically and thus executed everytime with `Travis CI <https://travis-ci.com/radis/radis>`_

    """

    Tgas = 1500

    sf = SpectrumFactory(
        wavelength_min=4250,
        wavelength_max=4350,
        mole_fraction=1,
        path_length=1,
        cutoff=1e-25,
        molecule="CO2",
        isotope="1,2",
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.warnings["HighTemperatureWarning"] = "ignore"
    sf.fetch_databank(
        "hitran", load_columns="noneq"
    )  # uses HITRAN: not really valid at this temperature, but runs on all machines without install
    s_nq = sf.non_eq_spectrum(Tvib=Tgas, Trot=Tgas, name="Non-eq")
    s_eq = sf.eq_spectrum(Tgas=Tgas, name="Eq")

    rtol = 7e-3  # 2nd isotope calculated with placeholder energies
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


@pytest.mark.needs_connection
def test_calc_spectrum_multiple_molecules(
    verbose=True, plot=True, warnings=True, *args, **kwargs
):
    """Test calculations with multiple molecules

    Note: try to keep the same wavelength ranges for each of the multi-molecule
    tests, so that databases are only downloaded once, and cached!"""

    s_co = calc_spectrum(
        wavelength_min=4165,
        wavelength_max=5000,
        Tgas=1000,
        path_length=0.1,
        mole_fraction=1,
        isotope={"CO": "1,2,3"},
        verbose=verbose,
    )

    s_co2 = calc_spectrum(
        wavelength_min=4165,
        wavelength_max=5000,
        Tgas=1000,
        path_length=0.1,
        mole_fraction=1,
        isotope={"CO2": "1,2"},
        verbose=verbose,
    )

    s_both = calc_spectrum(
        wavelength_min=4165,
        wavelength_max=5000,
        Tgas=1000,
        path_length=0.1,
        mole_fraction=1,
        isotope={"CO2": "1,2", "CO": "1,2,3"},
        verbose=verbose,
    )
    if plot:
        s_both.plot(wunit="nm")

    # Check calculation went fine:
    assert set(s_both.conditions["molecule"]) == set(["CO2", "CO"])

    # Compare
    from radis.los.slabs import MergeSlabs

    assert s_both.compare_with(MergeSlabs(s_co, s_co2), plot=False)

    return True


@pytest.mark.needs_connection
def test_calc_spectrum_multiple_molecules_otherinputs(
    verbose=True, plot=True, warnings=True, *args, **kwargs
):
    """Test calculations with differnet kind of inputs for multiple molecules

    Note: try to keep the same wavelength ranges for each of the multi-molecule
    tests, so that databases are only downloaded once, and cached!"""

    # Give molecule:
    s = calc_spectrum(
        wavelength_min=4165,
        wavelength_max=5000,
        Tgas=1000,
        path_length=0.1,
        molecule=["CO2", "CO"],
        mole_fraction=1,
        isotope={"CO2": "1,2", "CO": "1,2,3"},
        verbose=verbose,
    )
    assert set(s.conditions["molecule"]) == set(["CO2", "CO"])

    # Give isotope only
    s = calc_spectrum(
        wavelength_min=4165,
        wavelength_max=5000,
        Tgas=1000,
        path_length=0.1,
        isotope={"CO2": "1,2", "CO": "1,2,3"},
        verbose=verbose,
    )
    assert set(s.conditions["molecule"]) == set(["CO2", "CO"])

    # Give mole fractions only
    s = calc_spectrum(
        wavelength_min=4165,
        wavelength_max=5000,
        Tgas=1000,
        path_length=0.1,
        mole_fraction={"CO2": 0.2, "CO": 0.8},
        isotope="1,2",
        verbose=verbose,
    )
    assert set(s.conditions["molecule"]) == set(["CO2", "CO"])

    return True


# @pytest.mark.needs_config_file
# @pytest.mark.needs_db_HITEMP_CO2_DUNHAM
@pytest.mark.needs_connection
def test_calc_spectrum_multiple_molecules_inputerror(
    verbose=True, plot=True, warnings=True, *args, **kwargs
):
    """Test calculations with multiple molecules

    Note: try to keep the same wavelength ranges for each of the multi-molecule
    tests, so that databases are only downloaded once, and cached!"""

    # Contradictory:
    with pytest.raises(ValueError):
        calc_spectrum(
            wavelength_min=4165,
            wavelength_max=5000,
            Tgas=1000,
            path_length=0.1,
            molecule=["CO2"],  # contradictory
            mole_fraction=1,
            isotope={"CO2": "1,2", "CO": "1,2,3"},
            verbose=verbose,
        )

    # Partial:
    with pytest.raises(ValueError):
        calc_spectrum(
            wavelength_min=4165,
            wavelength_max=5000,
            Tgas=1000,
            path_length=0.1,
            molecule=["CO2", "CO"],  # contradictory
            mole_fraction=1,
            isotope={"CO2": "1,2"},  # unclear for CO
            verbose=verbose,
        )

    return True


@pytest.mark.needs_connection
def test_calc_spectrum_multiple_molecules_wstep_auto(
    verbose=True, plot=True, warnings=True, *args, **kwargs
):
    """Tests multiple molecules spectrum for wstep = 'auto'
    and checks that minimum wstep value is selected with
    resample = "intersect"""
    from radis import calc_spectrum

    # Merging the CO, CO2 spectrum itself in calc_spectrum
    s = calc_spectrum(
        wavelength_min=4165,
        wavelength_max=5000,  # cm-1
        isotope="1,2,3",
        pressure=1.01325,  # bar
        Tgas=700,  # K
        mole_fraction={"CO": 0.1, "CO2": 0.1},
        path_length=1,  # cm
        wstep="auto",
        databank="hitran",  # or use 'hitemp'
        verbose=verbose,
    )

    # Check calculation went fine:
    assert set(s.conditions["molecule"]) == set(["CO2", "CO"])
    assert s.get_conditions()["wstep"] in ("N/A", 0.013)


def test_check_wavelength_range(verbose=True, warnings=True, *args, **kwargs):
    """Check that input wavelength is correctly taken into account.
    See https://github.com/radis/radis/issues/214
    """
    if verbose:
        printm("Testing calc_spectrum wavelength range")

    wstep = 0.01

    s = calc_spectrum(
        wavelength_min=4348,  # nm
        wavelength_max=5000,
        molecule="CO",
        isotope="1,2,3",
        pressure=1.01325,  # bar
        Tvib=1700,  # K
        Trot=1700,  # K
        databank="HITRAN-CO-TEST",
        wstep=wstep,
    )
    w, I = s.get("radiance_noslit", wunit="nm", Iunit="mW/sr/cm2/nm")

    assert np.isclose(w.min(), 4348, atol=wstep)
    assert np.isclose(w.max(), 5000, atol=wstep)

    return True


def test_non_air_diluent_calc(verbose=True, plot=False, warnings=True, *args, **kwargs):
    from radis import plot_diff

    if verbose:
        printm("Regenerating database with all params, and calculating non_air diluent")

    # Regenerating database of CO and calculating non_air_diluent spectrum
    s1 = calc_spectrum(
        wavenum_min=2245,
        wavenum_max=2250,
        Tgas=700,
        path_length=0.1,
        molecule="CO",
        mole_fraction=0.2,
        isotope=1,
        databank="hitran",
        use_cached="regen",
        diluent={"CO2": 0.2, "H2O": 0.1, "air": 0.5},
        export_lines=True,
    )

    # Calculating spectrum with air as diluent
    s2 = calc_spectrum(
        wavenum_min=2245,
        wavenum_max=2250,
        Tgas=700,
        path_length=0.1,
        molecule="CO",
        mole_fraction=0.2,
        isotope=1,
        databank="hitran",
        export_lines=True,
    )

    assert s1.get_conditions()["diluents"] == {"CO2": 0.2, "H2O": 0.1, "air": 0.5}
    assert s2.get_conditions()["diluents"] == {"air": 0.8}

    hwhm_voigt_s1 = s1.lines.hwhm_voigt
    hwhm_voigt_s2 = s2.lines.hwhm_voigt

    assert (hwhm_voigt_s2 < hwhm_voigt_s1).all()
    assert np.isclose(hwhm_voigt_s2[0], 0.022718389546218788)
    assert np.isclose(hwhm_voigt_s1[0], 0.02530070148135749)

    if plot:
        plot_diff(
            s1,
            s2,
            method=["diff"],
            label1="Diluent CO2, H2O, air",
            label2="Diluent air",
            title="Diluent CO2:0.2, H2O:0.1, air:0.5",
        )

    return True


def test_diluent_invalid(verbose=True, plot=False, *args, **kwargs):

    with pytest.raises(KeyError) as err:
        calc_spectrum(
            wavenum_min=2245,
            wavenum_max=2250,
            Tgas=700,
            path_length=0.1,
            molecule="CO",
            mole_fraction=0.2,
            isotope=1,
            databank="hitran",
            diluent="CO",
        )

    assert (
        "is being called as molecule and diluent, please remove it from diluent."
        in str(err.value)
    )


def _run_testcases(plot=True, verbose=True, warnings=True, *args, **kwargs):

    # Test sPlanck and conversion functions
    test_sPlanck_conversions()

    # Test calc_spectrum function
    test_calc_spectrum()

    # Test calc_spectrum with overpopulation
    test_calc_spectrum_overpopulations(
        verbose=verbose, plot=plot, warnings=warnings, *args, **kwargs
    )

    # Compare all calc methods
    #    test_all_calc_methods_CO2(
    #        verbose=verbose, plot=plot, warnings=warnings, *args, **kwargs
    #    )
    test_all_calc_methods_CO2pcN(
        verbose=verbose, plot=plot, warnings=warnings, *args, **kwargs
    )

    # Compare same spectrum with two calculation methods
    test_eq_vs_noneq_isotope(
        verbose=verbose, plot=plot, warnings=warnings, *args, **kwargs
    )

    # Run test for multiple molecules
    test_calc_spectrum_multiple_molecules()
    test_calc_spectrum_multiple_molecules_otherinputs()
    test_calc_spectrum_multiple_molecules_inputerror()
    test_calc_spectrum_multiple_molecules_wstep_auto()

    test_check_wavelength_range()
    test_non_air_diluent_calc()

    return True


# --------------------------
if __name__ == "__main__":

    printm("Testing calc.py: ", _run_testcases(verbose=True))
