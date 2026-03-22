#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare CHIANTI Fe XII data with NIST/KURUCZ reference data
============================================================

This example demonstrates how to load and visualize atomic line data
from the CHIANTI database for Fe XII using RADIS. The sample database
file is included in the RADIS repository.

Fe XII (eleven-times ionized iron) is one of the most important ions
in solar EUV spectroscopy. The 195.12 Å line is routinely observed
by instruments such as SDO/AIA and Hinode/EIS.
"""

# %%
# Load Fe XII data from the sample CHIANTI database shipped with RADIS.
# The sample ``fe_12.wgfa`` file contains transitions in the EUV range
# around 193-197 Å (≈ 19.3-19.7 nm), which is the primary wavelength
# range used by solar observatories for coronal diagnostics.

import os

import matplotlib.pyplot as plt
import numpy as np

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

print("\u2705 CHIANTI Data Summary:")
print(f"   Total lines: {len(chianti_df)}")
if len(chianti_df) > 0:
    print(f"   Wavelength range: {chianti_df['wav'].min():.2f}"
          f" \u2013 {chianti_df['wav'].max():.2f} nm")
    print(f"   Strongest transition (max A): "
          f"{chianti_df.loc[chianti_df['A'].idxmax(), 'wav']:.2f} nm"
          f"  (A = {chianti_df['A'].max():.2e} s\u207b\u00b9)")

# %%
# Visualize Fe XII atomic line data.
#
# **Top panel** - Each point is one Fe XII transition. The colour
# encodes the radiative decay rate (Einstein A) on a log scale.
# Brighter points correspond to stronger transitions that dominate
# the observed solar EUV spectrum.
#
# **Bottom panel** - Einstein A coefficients span several orders of
# magnitude. The bar chart shows log10(A) for each transition,
# with x-tick labels giving the wavelength in Å.

fig, axes = plt.subplots(2, 1, figsize=(10, 9), constrained_layout=True)

wav_angstrom = chianti_df["wav"] * 10   # nm -> Å for display

# --- Top panel: scatter plot ---
ax1 = axes[0]
scatter1 = ax1.scatter(
    wav_angstrom,
    chianti_df["gf"],
    c=np.log10(chianti_df["A"]),
    cmap="viridis",
    s=80,
    alpha=0.85,
    edgecolors="k",
    linewidths=0.4,
    zorder=3,
)
ax1.set_xlabel("Wavelength (\u00c5)")
ax1.set_ylabel("Oscillator strength  gf")
ax1.set_title("CHIANTI Fe XII \u2014 Oscillator strength vs wavelength")
ax1.grid(True, alpha=0.25, linestyle="--")
cbar1 = plt.colorbar(scatter1, ax=ax1, pad=0.02)
cbar1.set_label("log$_{10}$(Einstein A) [s$^{-1}$]")

# Annotate the strongest line
idx_max = chianti_df["A"].idxmax()
ax1.annotate(
    f'{wav_angstrom.iloc[idx_max]:.1f} \u00c5',
    xy=(wav_angstrom.iloc[idx_max], chianti_df["gf"].iloc[idx_max]),
    xytext=(10, 10), textcoords="offset points",
    fontsize=9, fontstyle="italic",
    arrowprops=dict(arrowstyle="->", color="gray", lw=0.8),
)

# --- Bottom panel: bar chart ---
ax2 = axes[1]
sorted_df = chianti_df.sort_values("wav").reset_index(drop=True)
log_A = np.log10(sorted_df["A"])
colours = plt.cm.viridis(
    (log_A - log_A.min()) / (log_A.max() - log_A.min())
)
ax2.bar(
    range(len(sorted_df)),
    log_A,
    color=colours,
    alpha=0.85,
    edgecolor="k",
    linewidth=0.3,
)
ax2.set_ylabel("log$_{10}$(Einstein A) [s$^{-1}$]")
ax2.set_title("CHIANTI Fe XII \u2014 Radiative decay rates")
ax2.grid(True, alpha=0.25, axis="y", linestyle="--")
ax2.set_xticks(range(len(sorted_df)))
ax2.set_xticklabels(
    [f"{w:.1f}" for w in sorted_df["wav"] * 10],
    rotation=45, ha="right", fontsize=8,
)
ax2.set_xlabel("Wavelength (\u00c5)")

plt.show()

# %%
# Data table - first five transitions sorted by wavelength.
#
# Columns: vacuum wavelength (nm), weighted oscillator strength *gf*,
# and Einstein A coefficient (s-1).

print("\nFirst 5 transitions (sorted by wavelength):")
print(sorted_df[["wav", "gf", "A"]].head().to_string(index=False))

print(f"\nTotal transitions loaded: {len(chianti_df)}")
print(f"Mean log10(A): {np.log10(chianti_df['A']).mean():.2f}")
print(f"Median wavelength: {chianti_df['wav'].median():.2f} nm"
      f" ({chianti_df['wav'].median() * 10:.1f} \u00c5)")
