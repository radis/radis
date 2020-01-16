# -*- coding: utf-8 -*-
"""
Created on Mon May  7 17:34:52 2018

@author: erwan
"""

from __future__ import absolute_import, unicode_literals, division, print_function
from radis.lbl import SpectrumFactory
from radis.misc.utils import DatabankNotFound
from radis.misc.printer import printm
from radis.test.utils import IgnoreMissingDatabase, setup_test_line_databases
from radis.phys.convert import cm2nm
import pytest
import numpy as np
import astropy.units as u
from radis.misc.progress_bar import ProgressBar
from radis import get_residual, sPlanck
import radis
import matplotlib.pyplot as plt


@pytest.mark.fast
def test_populations(plot=True, verbose=True, warnings=True, *args, **kwargs):
    """ See wavelength difference in air and vacuum """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    try:
        if verbose:
            printm(">>> _test_media_line_shift")

        setup_test_line_databases()  # add HITRAN-CO-TEST in ~/.radis if not there

        sf = SpectrumFactory(
            wavelength_min=4500,
            wavelength_max=4600,
            wstep=0.001,
            parallel=False,
            bplot=False,
            cutoff=1e-30,
            path_length=0.1,
            mole_fraction=400e-6,
            isotope=[1],
            db_use_cached=True,
            medium="vacuum",
            broadening_max_width=10,
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

        # Compare with factory
        # Plot populations:
        with pytest.raises(ValueError):
            sf.plot_populations()  # no isotope given: error expected
        if plot:
            sf.plot_populations("rovib", isotope=1)
            plt.close()  # no need to keep it open, we just tested the function

        # Test calculated quantities are there
        assert hasattr(sf.df1, "Qref")
        assert hasattr(sf.df1, "Qvib")
        assert hasattr(sf.df1, "Qrotu")
        assert hasattr(sf.df1, "Qrotl")

        # Test hardcoded populations
        assert np.isclose(pops["rovib"]["n"].iloc[0], 0.0091853446840826653)
        assert np.isclose(pops["rovib"]["n"].iloc[1], 0.027052543988733215)
        assert np.isclose(pops["rovib"]["n"].iloc[2], 0.04345502115897712)

        return True

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


# @pytest.mark.needs_connection
def test_optically_thick_limit_1iso(verbose=True, plot=True, *args, **kwargs):
    """ Test that we find Planck in the optically thick limit 
    
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
    DEBUG_MODE = radis.DEBUG_MODE
    radis.DEBUG_MODE = True

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
            broadening_max_width=1,
            export_populations=False,  #'vib',
            export_lines=False,
            isotope=1,
            use_cached=True,
            wstep=wstep,
            pseudo_continuum_threshold=0,
            pressure=P,
            verbose=False,
        )
        #        sf.fetch_databank('astroquery')
        sf.warnings["NegativeEnergiesWarning"] = "ignore"
        sf.load_databank("HITEMP-CO2-TEST")
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

                nfig = "test_opt_thick_limit_1iso {0}".format(s.name)
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
            printm("Tested optically thick limit is Planck (1 isotope): OK")

    finally:
        # Reset DEBUG_MODE
        radis.DEBUG_MODE = DEBUG_MODE


# @pytest.mark.needs_connection
# @pytest.mark.needs_db_HITEMP_CO2_DUNHAM
# @pytest.mark.needs_connection
def test_optically_thick_limit_2iso(verbose=True, plot=True, *args, **kwargs):
    """ Test that we find Planck in the optically thick limit 
    
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
    DEBUG_MODE = radis.DEBUG_MODE
    radis.DEBUG_MODE = True

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
            broadening_max_width=1,
            export_populations=False,  #'vib',
            export_lines=False,
            molecule="CO2",
            isotope=[1, 2],
            db_use_cached=True,
            wstep=wstep,
            pseudo_continuum_threshold=0,
            pressure=P,
            verbose=False,
        )
        #        sf.fetch_databank('astroquery')
        sf.load_databank("HITEMP-CO2-TEST")
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
        radis.DEBUG_MODE = DEBUG_MODE


def _run_testcases(verbose=True, plot=True):

    test_populations(plot=plot, verbose=verbose)
    test_optically_thick_limit_1iso(plot=plot, verbose=verbose)
    test_optically_thick_limit_2iso(plot=plot, verbose=verbose)


if __name__ == "__main__":
    _run_testcases()
