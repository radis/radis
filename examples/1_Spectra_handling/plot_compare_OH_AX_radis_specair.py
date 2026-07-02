# -*- coding: utf-8 -*-
"""
===============================================================
Compare OH(A-X) electronic spectra: RADIS vs Specair
===============================================================

Auto-download the MoLLIST-OH line database from ExoMol and compute
the OH A-X electronic band (300-340 nm) under equilibrium and
non-equilibrium conditions. Compare with Specair reference data.

.. note::
    Electronic spectra computation in RADIS is a work in progress.
    RADIS currently shows ~10-20% discrepancy in integrated radiance
    compared to Specair, likely due to database differences.
    See :py:mod:`radis.test.validation.test_validation_vs_Specair_OH_AX`
    for the full validation test suite.

"""

from os.path import join

import numpy as np
import pandas as pd

from radis import Spectrum, SpectrumFactory, plot_diff
from radis.test.utils import getValidationCase


def _load_Specair_csv(csv_path):
    """Load a Specair reference CSV and return a Spectrum object."""
    data = pd.read_csv(csv_path)
    header = ",".join(data.columns)
    w = data.iloc[:, 0].values
    I = data.iloc[:, 1].values

    mask = ~np.isnan(I)
    w, I = w[mask], I[mask]
    idx = np.argsort(w)
    w, I = w[idx], I[idx]

    if "W_cm2_nm_sr" in header:
        I = I * w**2 / 1e4

    return Spectrum.from_array(
        w, I, "radiance", wunit="nm", unit="mW/cm2/sr/cm-1", name="Specair"
    )


def _integrate_radiance(s, wunit="nm", Iunit="mW/cm2/sr/cm-1"):
    """Get integrated radiance, filtering NaN and sorting."""
    w, I = s.get("radiance", wunit=wunit, Iunit=Iunit)
    mask = ~np.isnan(I)
    w, I = w[mask], I[mask]
    idx = np.argsort(w)
    return np.trapezoid(I[idx], w[idx])


# %% Setup SpectrumFactory
# -------------------------

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


# %% Plot 1: Equilibrium RADIS vs Specair (400 K)
# -------------------------------------------------

s_Specair_eq = _load_Specair_csv(
    getValidationCase(
        join("test_validation_vs_Specair_OH_AX_data", "radiance_spectrum.csv")
    )
)

s_radis_eq = sf.non_eq_spectrum(
    Ttrans=400, Trot=400, Tvib=400, Telec=10000, name="RADIS"
)
s_radis_eq.apply_slit(0.1)

fig, [ax0, ax1] = plot_diff(
    s_radis_eq,
    s_Specair_eq,
    "radiance",
    wunit="nm",
    Iunit="mW/cm2/sr/cm-1",
    method="diff",
)
ratio = _integrate_radiance(s_radis_eq) / _integrate_radiance(s_Specair_eq)
ax0.set_title(f"Equilibrium 400 K: RADIS vs Specair (ratio={ratio:.2f})")


# %% Plot 2: Self-consistency check (eq vs noneq at same T)
# -----------------------------------------------------------

s_eq = sf.eq_spectrum(Tgas=400, name="RADIS-eq")
s_noneq = sf.non_eq_spectrum(Tvib=400, Trot=400, Ttrans=400, name="RADIS-noneq")

fig2, [ax2, ax3] = plot_diff(s_eq, s_noneq, "radiance_noslit")
ratio2 = s_eq.get_integral("radiance_noslit") / s_noneq.get_integral("radiance_noslit")
ax2.set_title(f"Self-consistency: eq vs noneq at 400 K (ratio={ratio2:.4f})")


# %% Plot 3: Non-equilibrium Case 2 (Trot=500K, Tvib=2000K, Telec=10000K)
# -------------------------------------------------------------------------

s_Specair_c2 = _load_Specair_csv(
    getValidationCase(
        join(
            "test_validation_vs_Specair_OH_AX_data",
            "radiance_spectrum_Tgas300_Trot500_Tvib2000_Telec10000.csv",
        )
    )
)

s_radis_c2 = sf.non_eq_spectrum(
    Ttrans=300, Trot=500, Tvib=2000, Telec=10000, name="RADIS"
)
s_radis_c2.apply_slit(0.1)

fig3, [ax4, ax5] = plot_diff(
    s_radis_c2,
    s_Specair_c2,
    "radiance",
    wunit="nm",
    Iunit="mW/cm2/sr/cm-1",
    method="diff",
)
ratio3 = _integrate_radiance(s_radis_c2) / _integrate_radiance(s_Specair_c2)
ax4.set_title(f"Non-eq Trot=500 Tvib=2000 Telec=10000 (ratio={ratio3:.2f})")


# %% Plot 4: Non-equilibrium Case 3 (Trot=2000K, Tvib=500K, Telec=15000K)
# -------------------------------------------------------------------------

s_Specair_c3 = _load_Specair_csv(
    getValidationCase(
        join(
            "test_validation_vs_Specair_OH_AX_data",
            "radiance_spectrum_Tgas300_Trot2000_Tvib500_Telec15000.csv",
        )
    )
)

s_radis_c3 = sf.non_eq_spectrum(
    Ttrans=300, Trot=2000, Tvib=500, Telec=15000, name="RADIS"
)
s_radis_c3.apply_slit(0.1)

fig4, [ax6, ax7] = plot_diff(
    s_radis_c3,
    s_Specair_c3,
    "radiance",
    wunit="nm",
    Iunit="mW/cm2/sr/cm-1",
    method="diff",
)
ratio4 = _integrate_radiance(s_radis_c3) / _integrate_radiance(s_Specair_c3)
ax6.set_title(f"Non-eq Trot=2000 Tvib=500 Telec=15000 (ratio={ratio4:.2f})")
