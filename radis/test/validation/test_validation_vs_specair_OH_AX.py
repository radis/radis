# -*- coding: utf-8 -*-
"""
Validate RADIS OH(A-X) electronic spectrum against Specair reference data.

Three tests:
- Test 1: RADIS equilibrium vs Specair equilibrium (integrated radiance regression)
- Test 2: RADIS eq_spectrum vs non_eq_spectrum at same T (self-consistency)
- Test 3: RADIS non-eq vs Specair non-eq (integrated radiance regression)
"""

from os.path import join

import numpy as np
import pandas as pd
import pytest

from radis import Spectrum, SpectrumFactory, plot_diff
from radis.misc.printer import printm
from radis.test.utils import getValidationCase


def _integrate_radiance(s, wunit="nm", Iunit="mW/cm2/sr/cm-1"):
    """Get integrated radiance from a Spectrum, filtering NaN and sorting."""
    w, I = s.get("radiance", wunit=wunit, Iunit=Iunit)
    mask = ~np.isnan(I)
    w, I = w[mask], I[mask]
    idx = np.argsort(w)
    return np.trapezoid(I[idx], w[idx])


def _load_Specair_csv(csv_path):
    """Load a Specair reference CSV and return a Spectrum object.

    Handles unit conversion: Case 1 CSV is in W/cm^2/nm/sr,
    Cases 2 & 3 are in mW/cm^2/sr/cm^-1.
    """
    data = pd.read_csv(csv_path)
    header = ",".join(data.columns)
    w = data.iloc[:, 0].values  # nm
    I = data.iloc[:, 1].values

    # Filter NaN and sort
    mask = ~np.isnan(I)
    w, I = w[mask], I[mask]
    idx = np.argsort(w)
    w, I = w[idx], I[idx]

    # Unit conversion for Case 1
    if "W_cm2_nm_sr" in header:
        I = I * w**2 / 1e4  # W/cm^2/nm/sr -> mW/cm^2/sr/cm^-1

    return Spectrum.from_array(
        w, I, "radiance", wunit="nm", unit="mW/cm2/sr/cm-1", name="Specair"
    )


def _setup_OH_factory():
    """Create a SpectrumFactory for the OH A-X band (300-340 nm)."""
    sf = SpectrumFactory(
        wavelength_min=300,
        wavelength_max=340,
        molecule="OH",
        isotope="1",
        pressure=1,
        path_length=1,
        mole_fraction=1e-3,
        wstep=0.001,
        verbose=False,
        self_absorption=True,
    )
    sf.warnings["ElectronicSpectraWarning"] = "ignore"
    sf.fetch_databank("exomol", "MoLLIST-OH", load_energies=True)
    return sf


@pytest.mark.needs_connection
def test_OH_AX_eq_vs_specair(verbose=True, plot=False):
    """Test 1: RADIS equilibrium OH(A-X) vs Specair equilibrium at 400K.

    Regression test against known RADIS/Specair integrated radiance ratio.
    """
    s_specair = _load_Specair_csv(
        getValidationCase(
            join(
                "test_validation_vs_specair_OH_AX_data",
                "radiance_spectrum.csv",
            )
        )
    )

    sf = _setup_OH_factory()
    s_radis = sf.non_eq_spectrum(
        Ttrans=400, Trot=400, Tvib=400, Telec=10000, name="RADIS"
    )
    s_radis.apply_slit(0.1)

    ratio = _integrate_radiance(s_radis) / _integrate_radiance(s_specair)

    if verbose:
        printm(
            f">>> OH(A-X) eq 400K: RADIS/Specair integrated radiance ratio = {ratio:.4f}"
        )

    if plot:
        plot_diff(
            s_radis,
            s_specair,
            "radiance",
            wunit="nm",
            Iunit="mW/cm2/sr/cm-1",
            method="diff",
        )

    assert np.isclose(
        ratio, 1.23, rtol=0.03
    ), f"Regression: eq ratio {ratio:.4f} deviates from expected 1.23"


@pytest.mark.needs_connection
def test_OH_AX_eq_vs_noneq_self_consistency(verbose=True, plot=False):
    """Test 2: RADIS eq_spectrum vs non_eq_spectrum at same temperature.

    At thermal equilibrium (Tvib=Trot=Ttrans=400K, no Telec), both code
    paths should produce identical spectra within numerical error.
    """
    sf = _setup_OH_factory()

    s_eq = sf.eq_spectrum(Tgas=400, name="RADIS-eq")
    s_noneq = sf.non_eq_spectrum(Tvib=400, Trot=400, Ttrans=400, name="RADIS-noneq")

    I_eq = s_eq.get_integral("radiance_noslit")
    I_noneq = s_noneq.get_integral("radiance_noslit")
    ratio = I_eq / I_noneq

    if verbose:
        printm(
            f">>> OH(A-X) self-consistency: eq/noneq radiance_noslit ratio = {ratio:.4f}"
        )

    if plot:
        plot_diff(s_eq, s_noneq, "radiance_noslit")

    assert np.isclose(
        ratio, 1.0, rtol=3e-3
    ), f"eq vs noneq ratio {ratio:.4f} should be ~1.0 at same temperature"


@pytest.mark.needs_connection
def test_OH_AX_noneq_vs_specair(verbose=True, plot=False):
    """Test 3: RADIS non-eq OH(A-X) vs Specair non-eq.

    Regression test against known RADIS/Specair integrated radiance ratios.
    The discrepancy is documented and likely due to database differences.

    Known ratios:
      Case 2 (Trot=500, Tvib=2000, Telec=10000): ~1.22
      Case 3 (Trot=2000, Tvib=500, Telec=15000): ~1.12
    """
    sf = _setup_OH_factory()

    # Case 2
    s_Specair_c2 = _load_Specair_csv(
        getValidationCase(
            join(
                "test_validation_vs_specair_OH_AX_data",
                "radiance_spectrum_Tgas300_Trot500_Tvib2000_Telec10000.csv",
            )
        )
    )
    s_radis_c2 = sf.non_eq_spectrum(
        Ttrans=300, Trot=500, Tvib=2000, Telec=10000, name="RADIS-case2"
    )
    s_radis_c2.apply_slit(0.1)

    ratio2 = _integrate_radiance(s_radis_c2) / _integrate_radiance(s_Specair_c2)

    # Case 3
    s_Specair_c3 = _load_Specair_csv(
        getValidationCase(
            join(
                "test_validation_vs_Specair_OH_AX_data",
                "radiance_spectrum_Tgas300_Trot2000_Tvib500_Telec15000.csv",
            )
        )
    )
    s_radis_c3 = sf.non_eq_spectrum(
        Ttrans=300, Trot=2000, Tvib=500, Telec=15000, name="RADIS-case3"
    )
    s_radis_c3.apply_slit(0.1)

    ratio3 = _integrate_radiance(s_radis_c3) / _integrate_radiance(s_Specair_c3)

    if verbose:
        printm(
            f">>> OH(A-X) non-eq Case 2: RADIS/Specair ratio = {ratio2:.4f} (expected ~1.22)"
        )
        printm(
            f">>> OH(A-X) non-eq Case 3: RADIS/Specair ratio = {ratio3:.4f} (expected ~1.12)"
        )

    if plot:
        plot_diff(
            s_radis_c2,
            s_Specair_c2,
            "radiance",
            wunit="nm",
            Iunit="mW/cm2/sr/cm-1",
            method="diff",
        )
        plot_diff(
            s_radis_c3,
            s_Specair_c3,
            "radiance",
            wunit="nm",
            Iunit="mW/cm2/sr/cm-1",
            method="diff",
        )

    assert np.isclose(
        ratio2, 1.22, rtol=0.03
    ), f"Regression: Case 2 ratio {ratio2:.4f} deviates from expected 1.22"
    assert np.isclose(
        ratio3, 1.12, rtol=0.03
    ), f"Regression: Case 3 ratio {ratio3:.4f} deviates from expected 1.12"


if __name__ == "__main__":
    test_OH_AX_eq_vs_specair(plot=True, verbose=True)
    test_OH_AX_eq_vs_noneq_self_consistency(plot=True, verbose=True)
    test_OH_AX_noneq_vs_specair(plot=True, verbose=True)
