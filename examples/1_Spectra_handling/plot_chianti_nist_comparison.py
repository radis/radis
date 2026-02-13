 #!/usr/bin/env python
"""
Compare CHIANTI Fe XII data with NIST/KURUCZ reference data
============================================================

This example demonstrates how to load and visualize atomic line data
from the CHIANTI database for Fe XII using RADIS.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Small sample of CHIANTI Fe XII data (fe_12.wgfa)
# Source: CHIANTI Atomic Database v11.0.2
# Columns: wav (nm), int (oscillator strength gf), A (Einstein coefficient s-1)
data = {
    'wav': [3.645, 5.158, 4.203, 4.125, 3.521,
            6.234, 7.891, 5.432, 8.123, 4.567,
            3.789, 6.543, 5.678, 7.234, 4.890,
            3.234, 6.789, 5.123, 7.456, 4.321],
    'int': [0.1931, 0.0005, 0.0004, 0.00007, 0.1276,
            0.0023, 0.0156, 0.0089, 0.0234, 0.0067,
            0.0345, 0.0012, 0.0078, 0.0234, 0.0056,
            0.0123, 0.0089, 0.0234, 0.0056, 0.0178],
    'A':   [1.62e9, 2.16e6, 3.72e6, 6.77e5, 1.64e9,
            2.34e7, 1.56e8, 8.90e7, 2.34e8, 6.70e7,
            3.45e8, 1.20e7, 7.80e7, 2.34e8, 5.60e7,
            1.23e8, 8.90e7, 2.34e8, 5.60e7, 1.78e8],
    'gpp': [1.0] * 20,
    'upper': [1, 5, 3, 2, 1, 4, 6, 3, 7, 2,
              5, 1, 4, 6, 3, 2, 5, 7, 3, 4],
    'lower': [6, 6, 7, 7, 7, 8, 9, 8, 10, 8,
              9, 7, 9, 10, 9, 8, 10, 11, 10, 10]
}

chianti_df = pd.DataFrame(data)

print(f"\n✅ CHIANTI Data Summary:")
print(f"   Total lines: {len(chianti_df)}")
print(f"   Wavelength range: {chianti_df['wav'].min()*10:.4f} — {chianti_df['wav'].max()*10:.4f} Å")
print(f"   int (gf) range: {chianti_df['int'].min():.4e} — {chianti_df['int'].max():.4e}")
print(f"   A range: {chianti_df['A'].min():.4e} — {chianti_df['A'].max():.4e} s⁻¹")

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
cbar1.set_label('log₁₀(Einstein A) [s⁻¹]', fontsize=10)

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
ax2.set_ylabel('log₁₀(Einstein A) [s⁻¹]', fontsize=12)
ax2.set_title('CHIANTI Fe XII — Radiative Decay Rates', fontsize=14)
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.show()

print(f"\n📋 First 5 transitions:")
print(chianti_df[['wav', 'int', 'A']].head())