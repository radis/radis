# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 21:37:59 2017

@author: erwan

Reproduce the plasma torch experiment of Packan 2003 in in atmospheric air

Conditions
--------

There is 6 m of room air between the plasma and the detector. The temperature
and humidity in the laboratory were recorded during experiment and calibration.
In both cases the measured ambient temperature was 25°C, and the relative humidity
was 42% (10% uncertainty), which corresponds to a mole fraction of H2O of
1.3+/-0.1 x 1E-2. We used 1.4e-2 because it gave better agreement for absorption
on the fundamental bands of NO. The experiments were done in about 1997, and the
room-air CO2 concentration was assumed to be 330 ppm.

Temperatures and concentrations of all slabs are stored in
`test_compare_torch_CO2_conditions_JTHT2003.dat`



"""

from time import time

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from radis.lbl import SpectrumFactory
from radis.los import MergeSlabs, SerialSlabs
from radis.misc.printer import printm
from radis.phys.convert import nm2cm
from radis.spectrum import experimental_spectrum, get_residual
from radis.test.utils import getValidationCase


@pytest.mark.needs_config_file
@pytest.mark.needs_db_HITEMP_CO2_DUNHAM
@pytest.mark.needs_db_HITEMP_CO_DUNHAM
def test_compare_torch_CO2(
    verbose=True, plot=False, save=False, warnings=True, use_cache=True, *args, **kwargs
):
    """Reproduce the plasma torch experiment of Packan 2003 in in atmospheric air

    Notes
    -----

    Thresholds parameters are reduced for a faster calculation time. For instance:

    - truncation should be 25 rather than 10
    - cutoff should be 1e-27 rather than 1e-25
    - wstep should be 0.01 or even 0.008 rather than 0.1

    Performance:

    - neq==0.9.20: Finished test_compare_torch_CO2 in 758.640324s

    - neq==0.9.21: Finished test_compare_torch_CO2 in 672s

    - neq==0.9.21*: (with ParallelFactory) Finished test_compare_torch_CO2 in 298s

    - neq==0.9.22: (Parallel + continuum) Finished in 65s
      RADIS 0.9.9 == neq 0.9.24

    - RADIS 0.9.26: Finished test_compare_torch_CO2 in 57s (LDM, no continuum)

    - RADIS 0.9.28: ParallelFactory is removed. Finished test_compare_torch_CO2 in 45s

    Reference
    --------

    Packan et al 2003 "Measurement and Modeling of OH, NO, and CO Infrared Radiation
    at 3400 K", JTHT

    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    t0 = time()

    # %% Conditions

    wlmin = 4100
    wlmax = 4900

    wmin = nm2cm(wlmax)
    wmax = nm2cm(wlmin)
    wstep = 0.1

    # Load concentration profile
    conds = pd.read_csv(
        getValidationCase(
            "test_compare_torch_CO2_data/test_compare_torch_CO2_conditions_JTHT2003.dat"
        ),
        comment="#",
        delim_whitespace=True,
    )
    slab_width = 0.05  # cm,  WARNING. Hardcoded (look up the table)

    # %% Calculate slabs

    # CO2
    sf = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=None,
        path_length=None,
        molecule="CO2",
        isotope="1,2",
        wstep=wstep,
        export_lines=False,  # saves some memory
        export_populations=False,  # saves some memory
        cutoff=1e-25,
        truncation=10,
        # pseudo_continuum_threshold=0.01, # use pseudo-continuum, no LDM. Note : 56s on 20/08.
        # optimization=None,
        pseudo_continuum_threshold=0,  # use LDM, no pseudo-continuum. Note : 84s on 20/08
        optimization="min-RMS",
        verbose=False,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    #        # Init database: only if we want to save all spectra being calculated
    # (discarded for the moment)
    #        if use_cache:
    #            sf.init_database('test_compare_torch_CO2_SpectrumDatabase_HITEMP_1e25',
    #                             add_info=['Tgas'])
    sf.init_databank(
        "HITEMP-CO2-DUNHAM",
        db_use_cached=use_cache,
    )  # saves some memory
    #    sf.init_databank('CDSD-4000')

    # CO
    sfco = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=None,
        path_length=None,
        molecule="CO",
        isotope="1,2",
        wstep=wstep,
        export_lines=False,  # saves some memory
        export_populations=False,  # saves some memory
        cutoff=1e-25,
        truncation=10,
        optimization=None,
        verbose=False,
    )
    sfco.warnings["MissingSelfBroadeningWarning"] = "ignore"
    # Init database: only if we want to save all spectra being calculated
    # (discarded for the moment)
    #        if use_cache:
    #            sfco.init_database(
    #                'test_compare_torch_CO_SpectrumDatabase_HITEMP_1e25', add_info=['Tgas'])
    sfco.init_databank("HITEMP-CO-DUNHAM", db_use_cached=use_cache)

    # Calculate all slabs
    # .. CO2 at equilibrium
    slabsco2 = []
    for Tgas, xCO2 in zip(conds.T_K, conds.CO2):
        slabsco2.append(
            sf.eq_spectrum(Tgas=Tgas, mole_fraction=xCO2, path_length=slab_width)
        )
    # .. CO at equilibrium
    slabsco = []
    for Tgas, xCO in zip(conds.T_K, conds.CO):
        slabsco.append(
            sfco.eq_spectrum(Tgas=Tgas, mole_fraction=xCO, path_length=slab_width)
        )

    # Room absorption
    s0 = sf.eq_spectrum(Tgas=300, mole_fraction=330e-6, path_length=600)
    # Warning: This slab includes the CO2 emission from the 300 K air
    # within the spectrometer. But experimentally it is substracted
    # by there the chopper (this is corrected below)
    # ... see RADIS 2019 paper for more details

    # Add radiance for everyone if not calculated (ex: loaded from database)
    for s in slabsco2 + slabsco + [s0]:
        s.update()

    # Merge CO + CO2 for each slab
    slabstot = [MergeSlabs(s, sco) for (s, sco) in zip(slabsco2, slabsco)]

    # %% Line-of-sight

    # Solve RTE along the line of sight
    # --------
    # two semi profile, + room absortion
    line_of_sight = slabstot[1:][::-1] + slabstot + [s0]
    stot = SerialSlabs(*line_of_sight)
    #        stot = SerialSlabs(*slabstot[1:][::-1], *slabstot, s0)  # Python 3 syntax only

    # Also calculate the contribution of pure CO
    # ---------
    s0tr = s0.copy()  # transmittance only
    s0tr.conditions["thermal_equilibrium"] = False
    s0tr._q["radiance_noslit"] *= 0  # hack to remove CO2 emission
    s0tr._q["emisscoeff"] *= 0  # hack to remove CO2 emission
    line_of_sight = slabsco[1:][::-1] + slabsco + [s0tr]
    sco = SerialSlabs(*line_of_sight)
    #        sco = SerialSlabs(*slabsco[1:][::-1], *slabsco, s0tr)  # Python 3 syntax only

    # Generate Slit
    # ------
    disp = 4  # dispersion conversion function (nm/mm)
    s_in = 1  # entrance slit (mm)
    s_out = 2.8  # exit slit (mm)
    top = (s_out - s_in) * disp
    base = (s_out + s_in) * disp
    slit_function = (top, base)  # (nm)

    # Final plot unit
    unit = "µW/cm2/sr"
    norm_by = "max"  # should be 'max' if unit ~ 'µW/cm2/sr',
    # 'area' if unit ~ 'µW/cm2/sr/nm'

    # Convolve with Slit
    # -------
    sco.apply_slit(slit_function, unit="nm", norm_by=norm_by, shape="trapezoidal")
    stot.apply_slit(slit_function, unit="nm", norm_by=norm_by, shape="trapezoidal")

    # Remove emission from within the spectrometer (substracted
    # by the chopper experimentaly)
    # -------
    s0spectro = s0.copy()
    s0spectro.rescale_path_length(75 * 4)  # spectrometer length
    s0spectro.apply_slit(slit_function, unit="nm", norm_by=norm_by, shape="trapezoidal")
    _, I0spectro = s0spectro.get("radiance", Iunit=unit)

    wtot, Itot = stot.get("radiance", wunit="nm", Iunit=unit)
    stot_corr = experimental_spectrum(
        wtot,
        Itot - I0spectro,  # hack to remove
        wunit="nm",  # in air by default
        # emission from within
        # spectrometer
        Iunit=unit,
    )

    # %% Compare with experiment data

    # Plot experimental data
    # ------
    exp = pd.read_csv(
        getValidationCase(
            "test_compare_torch_CO2_data/test_compare_torch_CO2_spectrum_JTHT2003.dat"
        ),
        skiprows=5,
        delim_whitespace=True,
    )
    exp.w /= 10  # Angstrom to nm
    exp = exp[(exp.w > wlmin) & (exp.w < wlmax)]
    sexp = experimental_spectrum(exp.w, exp.I, wunit="nm", Iunit="mW/cm2/sr")

    if plot:

        sexp.plot("radiance", wunit="nm", Iunit=unit, lw=2, label="Packan 2003")

        # Plot calculated
        # ----------

        stot.plot(
            "radiance",
            wunit="nm",
            Iunit=unit,
            nfig="same",
            color="r",
            ls="--",
            label="All: CO$_2$+CO",
        )

        # Plot internal spectrometer emission
        s0spectro.plot(
            "radiance",
            wunit="nm",
            Iunit=unit,
            nfig="same",
            color="b",
            label="CO$_2$ in Spectrometer",
        )

        stot_corr.plot(
            "radiance",
            wunit="nm",
            Iunit=unit,
            nfig="same",
            color="r",
            lw=2,
            label="All-CO$_2$ in Spectro",
            zorder=10,
        )

        sco.plot(
            "radiance",
            wunit="nm",
            Iunit=unit,
            nfig="same",
            color="g",
            ls="-",
            label="CO",
        )

        plt.ylim((0, 9))
        plt.legend(loc="upper right")

        if save:
            plt.savefig("test_compare_torch_CO2_brd20.png")
            plt.savefig("test_compare_torch_CO2_brd20k.pdf")

    assert get_residual(stot_corr, sexp, "radiance") < 0.02

    if verbose:
        print(("Finished test_compare_torch_CO2 in {0:.0f}s".format(time() - t0)))


if __name__ == "__main__":
    printm("test_compare_torch_CO2:", test_compare_torch_CO2(plot=True, use_cache=True))
