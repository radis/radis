#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare CHIANTI Fe XII data with NIST/KURUCZ reference data
============================================================

This example demonstrates how to load and visualize atomic line data
from the CHIANTI database for Fe XII using RADIS.
The sample database file is included in the RADIS repository.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from radis.io.chianti import load_chianti

# Path to the sample CHIANTI database included in RADIS
import radis
SAMPLE_DATAPATH = os.path.join(
    os.path.dirname(radis.__file__),
    "db", "chianti"
)

# Load Fe XII data from sample database
chianti_df = load_chianti(
    datapath=SAMPLE_DATAPATH,
    ion='fe_12',
    wmin=0.0,
    wmax=40.0
)

print("\n✅ CHIANTI Data Summary:")
print(f"   Total lines: {len(chianti_df)}")

# Create visualization
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Plot 1: Wavelength vs gf (oscillator strength)
ax1 = axes[0]
scatter1 = ax1.scatter(
    chianti_df['wav']*10,
    chianti_df['int'],
    c=np.log10(chianti_df['A']),
    cmap='viridis',
    s=50,
    alpha=0.8
)
ax1.set_xlabel('Wavelength (Å)', fontsize=12)
ax1.set_ylabel('Oscillator Strength (gf)', fontsize=12)
ax1.set_title('CHIANTI Fe XII — Oscillator Strength vs Wavelength', fontsize=14)
ax1.grid(True, alpha=0.3)
cbar1 = plt.colorbar(scatter1, ax=ax1)
cbar1.set_label('log10(Einstein A) [s-1]', fontsize=10)

# Plot 2: Line index vs A coefficient
ax2 = axes[1]
line_indices = np.arange(len(chianti_df))
ax2.bar(
    line_indices,
    np.log10(chianti_df['A']),
    color='steelblue',
    alpha=0.8,
)
ax2.set_xlabel('Line Index', fontsize=12)
ax2.set_ylabel('log10(Einstein A) [s-1]', fontsize=12)
ax2.set_title('CHIANTI Fe XII — Radiative Decay Rates', fontsize=14)
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.show()

print("\nFirst 5 transitions:")
print(chianti_df[['wav', 'int', 'A']].head())
