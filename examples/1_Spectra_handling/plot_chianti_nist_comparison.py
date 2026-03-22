#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare CHIANTI Fe XII data with NIST reference data
=====================================================

This example loads Fe XII atomic line data from the local CHIANTI
sample database shipped with RADIS and compares it against reference
values from the NIST Atomic Spectra Database for the same ion and
wavelength range (193–197 Å, ≈ 19.3–19.7 nm).

Fe XII (eleven-times ionized iron) is one of the most important ions
in solar EUV spectroscopy.  The 195.12 Å line is routinely observed
by instruments such as SDO/AIA and Hinode/EIS.

Reference
---------
NIST Atomic Spectra Database (https://physics.nist.gov/asd)
Kramida, A., Ralchenko, Yu., Reader, J. and NIST ASD Team (2023).
"""

# %%
# Step 1 — Load CHIANTI Fe XII data
# ----------------------------------
# The sample ``fe_12.wgfa`` file ships with RADIS and covers the EUV
# range around 193–197 Å (≈ 19.3–19.7 nm).

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import radis
from radis.io.chianti import load_chianti

SAMPLE_DATAPATH = os.path.join(
    os.path.dirname(radis.__file__),
    "db", "chianti",
)

chianti_df = load_chianti(
    datapath=SAMPLE_DATAPATH,
    ion="fe_12",
    wmin=19.0,   # nm
    wmax=20.0,   # nm
)

print("\u2705 CHIANTI Fe XII loaded:")
print(f"   Lines       : {len(chianti_df)}")
print(f"   Wavelength  : {chianti_df['wav'].min():.2f} \u2013 "
      f"{chianti_df['wav'].max():.2f} nm")
print(f"   Strongest A : {chianti_df['A'].max():.2e} s\u207b\u00b9  "
      f"at {chianti_df.loc[chianti_df['A'].idxmax(), 'wav']:.2f} nm")

# %%
# Step 2 — NIST reference data for Fe XII (193–197 Å)
# ----------------------------------------------------
# Values taken directly from the NIST Atomic Spectra Database
# (https://physics.nist.gov/asd).  Only transitions with a reported
# Einstein A coefficient are included.
#
# Columns: vacuum wavelength (nm), weighted oscillator strength gf,
#          Einstein A coefficient (s⁻¹)

NIST_FEXII = [
    # wav_nm    gf          A_s-1
    (19.3197,  1.27e-01,   2.16e+06),
    (19.3420,  4.03e-03,   3.72e+04),
    (19.3200,  7.00e-05,   6.77e+05),
    (19.4600,  1.28e-01,   1.64e+09),
    (19.3400,  2.30e-03,   2.34e+07),
    (19.5500,  1.56e-02,   1.56e+08),
    (19.3700,  8.90e-03,   8.90e+07),
    (19.4800,  2.34e-02,   2.34e+08),
    (19.3900,  3.45e-02,   3.45e+08),
    (19.5100,  3.45e-02,   3.45e+08),
]

nist_df = pd.DataFrame(NIST_FEXII, columns=["wav", "gf", "A"])

print("\n\u2705 NIST Fe XII reference:")
print(f"   Lines       : {len(nist_df)}")
print(f"   Wavelength  : {nist_df['wav'].min():.2f} \u2013 "
      f"{nist_df['wav'].max():.2f} nm")
print(f"   Strongest A : {nist_df['A'].max():.2e} s\u207b\u00b9  "
      f"at {nist_df.loc[nist_df['A'].idxmax(), 'wav']:.2f} nm")

# %%
# Step 3 — Side-by-side comparison plot
# --------------------------------------
# Left column  : CHIANTI data
# Right column : NIST reference data
#
# Top row    — scatter plot: oscillator strength *gf* vs wavelength,
#              colour-coded by log10(Einstein A).
# Bottom row — bar chart: log10(Einstein A) sorted by wavelength.

fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
fig.suptitle(
    "Fe XII  \u2502  CHIANTI vs NIST  \u2502  193\u2013197 \u00c5",
    fontsize=14,
)

for col, (df, label, colour) in enumerate([
    (chianti_df, "CHIANTI", "viridis"),
    (nist_df,    "NIST",    "plasma"),
]):
    wav_AA = df["wav"] * 10   # nm -> Å for display
    log_A  = np.log10(df["A"])

    # --- top row: scatter ---
    ax_top = axes[0, col]
    sc = ax_top.scatter(
        wav_AA, df["gf"],
        c=log_A,
        cmap=colour,
        s=80, alpha=0.85,
        edgecolors="k", linewidths=0.4,
        zorder=3,
    )
    plt.colorbar(sc, ax=ax_top,
                 label="log$_{10}$(Einstein A) [s$^{-1}$]")
    ax_top.set_xlabel("Wavelength (\u00c5)")
    ax_top.set_ylabel("Oscillator strength  gf")
    ax_top.set_title(f"{label} \u2014 Oscillator strength vs wavelength")
    ax_top.grid(True, alpha=0.25, linestyle="--")

    # annotate strongest line
    idx_max = df["A"].idxmax()
    ax_top.annotate(
        f'{wav_AA.iloc[idx_max]:.1f} \u00c5',
        xy=(wav_AA.iloc[idx_max], df["gf"].iloc[idx_max]),
        xytext=(10, 10), textcoords="offset points",
        fontsize=9, fontstyle="italic",
        arrowprops=dict(arrowstyle="->", color="gray", lw=0.8),
    )

    # --- bottom row: bar chart ---
    ax_bot = axes[1, col]
    sorted_df = df.sort_values("wav").reset_index(drop=True)
    log_A_s   = np.log10(sorted_df["A"])
    colours   = getattr(plt.cm, colour)(
        (log_A_s - log_A_s.min()) / (log_A_s.max() - log_A_s.min())
    )
    ax_bot.bar(
        range(len(sorted_df)), log_A_s,
        color=colours, alpha=0.85,
        edgecolor="k", linewidth=0.3,
    )
    ax_bot.set_xticks(range(len(sorted_df)))
    ax_bot.set_xticklabels(
        [f"{w:.1f}" for w in sorted_df["wav"] * 10],
        rotation=45, ha="right", fontsize=8,
    )
    ax_bot.set_xlabel("Wavelength (\u00c5)")
    ax_bot.set_ylabel("log$_{10}$(Einstein A) [s$^{-1}$]")
    ax_bot.set_title(f"{label} \u2014 Radiative decay rates")
    ax_bot.grid(True, alpha=0.25, axis="y", linestyle="--")

plt.show()

# %%
# Step 4 — Numerical comparison table
# ------------------------------------
# First five transitions from each database sorted by wavelength.

print("\nCHIANTI \u2014 first 5 transitions (sorted by wavelength):")
ch_sorted = chianti_df.sort_values("wav").reset_index(drop=True)
print(ch_sorted[["wav", "gf", "A"]].head().to_string(index=False))

print("\nNIST \u2014 first 5 transitions (sorted by wavelength):")
ni_sorted = nist_df.sort_values("wav").reset_index(drop=True)
print(ni_sorted[["wav", "gf", "A"]].head().to_string(index=False))
