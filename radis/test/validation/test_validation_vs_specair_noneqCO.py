# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 21:10:31 2017

@author: erwan

Plot Absorbance (cm-1) of CO with our code, database HITRAN 2012 against Specair,
under non-equilibrium conditions
"""

from os.path import join

import matplotlib.pyplot as plt
import numpy as np
import pytest

from radis import SpectrumFactory, load_spec, plot_diff
from radis.misc.printer import printm
from radis.test.utils import getValidationCase, setup_test_line_databases


@pytest.mark.fast
def test_validation_vs_specair(rtol=1e-2, verbose=True, plot=False, *args, **kwargs):
    """Test RADIS output on CO IR bands against SPECAIR

    Test is only performed on integrals of absorption coefficient

    RADIS doesnt actually match Specair exactly, but this is due to line intensity
    differences (Specair has no rovibrational specific intensities) rather than
    differences in populations calculations, as evidenced by the partition functions
    comparison in the RADIS presentation article.

    """

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    # %% Specair calculation
    # -----------

    specair_300_300 = load_spec(
        getValidationCase(
            join(
                "test_validation_vs_specair_noneqCO_data",
                "specair_CO_IR_Tvib300_Trot300.spec",
            )
        ),
        binary=True,
    )
    specair_300_2000 = load_spec(
        getValidationCase(
            join(
                "test_validation_vs_specair_noneqCO_data",
                "specair_CO_IR_Tvib300_Trot2000.spec",
            )
        ),
        binary=True,
    )
    specair_2000_300 = load_spec(
        getValidationCase(
            join(
                "test_validation_vs_specair_noneqCO_data",
                "specair_CO_IR_Tvib2000_Trot300.spec",
            )
        ),
        binary=True,
    )

    article_version = False  # just for the article

    # %% Compare with RADIS
    # ----------

    wstep = 0.002
    pl = SpectrumFactory(
        wavelength_min=4400,
        wavelength_max=4900,
        mole_fraction=1,
        pressure=0.01,  # bar
        path_length=1,  # we dont care for abscoeff anyway
        cutoff=1e-30,
        wstep=wstep,
        isotope=1,  # '1,2,3',
        medium="vacuum",
        export_lines=True,
    )  # 0.2)
    pl.warnings["MissingSelfBroadeningWarning"] = "ignore"

    if article_version:
        pl.load_databank("HITEMP-CO-DUNHAM")
    else:
        pl.load_databank("HITRAN-CO-TEST")
        # Available on all systems, convenient for fast testing, but you
        # will be missing some small lines for T ~ 2000 K .
        print(
            "Using HITRAN: small lines will be missing at T ~ 2000K. Use HITEMP if you want them"
        )

    s_300_300 = pl.eq_spectrum(300, name="RADIS")

    s_2000_300 = pl.non_eq_spectrum(Tvib=2000, Trot=300, Ttrans=300, name="RADIS")

    s_300_2000 = pl.non_eq_spectrum(Tvib=300, Trot=2000, Ttrans=2000, name="RADIS")

    # %% Test

    # Compare integrals

    b1 = np.isclose(
        specair_300_300.get_integral("abscoeff"),
        s_300_300.get_integral("abscoeff"),
        rtol=rtol,
    )
    b1 *= np.isclose(
        specair_2000_300.get_integral("abscoeff"),
        s_2000_300.get_integral("abscoeff"),
        rtol=rtol,
    )
    b1 *= np.isclose(
        specair_300_2000.get_integral("abscoeff"),
        s_300_2000.get_integral("abscoeff"),
        rtol=rtol,
    )

    # Compare partition functions to hardcoded values
    b2 = np.isclose(s_2000_300.lines.attrs["Q"], 139, atol=1)  # Specair: 139
    b2 *= np.isclose(s_300_2000.lines.attrs["Q"], 727, atol=1)  # Specair: 727

    if verbose:
        printm(
            ">>> comparing RADIS vs SPECAIR on CO: integrals of abscoeff is are close"
            + " to within {0:.1f}%: {1} ({2:.1f}%, {3:.1f}%, {4:.1f}%)".format(
                rtol * 100,
                bool(b1),
                abs(
                    specair_300_300.get_integral("abscoeff")
                    / s_300_300.get_integral("abscoeff")
                    - 1
                )
                * 100,
                abs(
                    specair_2000_300.get_integral("abscoeff")
                    / s_2000_300.get_integral("abscoeff")
                    - 1
                )
                * 100,
                abs(
                    specair_300_2000.get_integral("abscoeff")
                    / s_300_2000.get_integral("abscoeff")
                    - 1
                )
                * 100,
            )
        )
        printm(
            ">>> comparing RADIS vs SPECAIR on CO: partition functions "
            + "are equal to round error: {0}".format(bool(b2))
        )

    if plot:
        plot_diff(
            specair_300_300,
            s_300_300,
            title=r"T$_\mathregular{vib}$ 300 K, T$_\mathregular{rot}$ 300 K",
            diff_window=int(
                0.02 // wstep
            ),  # compensate for small shifts in both codes. we're comparing intensities here.
            lw_multiplier=1,  # 0.75,
            wunit="nm_vac",
            plot_medium=True,
        )
        plt.xlim((4500, 4900))
        if article_version:
            plt.savefig("out/test_validation_vs_specair_noneqCO_Tvib300_Trot300.png")
            plt.savefig("out/test_validation_vs_specair_noneqCO_Tvib300_Trot300.pdf")

        plot_diff(
            specair_2000_300,
            s_2000_300,
            title=r"T$_\mathregular{vib}$ 2000 K, T$_\mathregular{rot}$ 300 K",
            diff_window=int(
                0.02 // wstep
            ),  # compensate for small shifts in both codes. we're comparing intensities here.
            lw_multiplier=1,  # 0.75,
            wunit="nm_vac",
            plot_medium=True,
        )
        plt.xlim((4500, 4900))
        if article_version:
            plt.savefig("out/test_validation_vs_specair_noneqCO_Tvib2000_Trot300.png")
            plt.savefig("out/test_validation_vs_specair_noneqCO_Tvib2000_Trot300.pdf")

        plot_diff(
            specair_300_2000,
            s_300_2000,
            title=r"T$_\mathregular{vib}$ 300 K, T$_\mathregular{rot}$ 2000 K",
            diff_window=int(
                0.02 // wstep
            ),  # compensate for small shifts in both codes. we're comparing intensities here.
            lw_multiplier=1,  # 0.75,
            wunit="nm_vac",
            plot_medium=True,
        )
        plt.xlim((4500, 4900))
        if article_version:
            plt.savefig("out/test_validation_vs_specair_noneqCO_Tvib300_Trot2000.png")
            plt.savefig("out/test_validation_vs_specair_noneqCO_Tvib300_Trot2000.pdf")

    return bool(b1 * b2)


if __name__ == "__main__":
    printm(
        "test_validation_vs_specair_noneqCO:",
        test_validation_vs_specair(3e-2, plot=True, verbose=True),
    )
