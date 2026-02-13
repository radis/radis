#!/usr/bin/env python
"""
Compare CHIANTI Fe XII data with NIST/KURUCZ reference data
"""

import matplotlib.pyplot as plt
import numpy as np
from radis.io.chianti import load_chianti

# Load CHIANTI Fe XII data
print("Loading CHIANTI Fe XII data...")
chianti_df = load_chianti(wmin=0.0, wmax=400.0, ion='fe_12')  # Full range
if len(chianti_df) == 0:
    print("❌ No CHIANTI data found! Check your CHIANTI path.")
    exit(1)

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
    alpha=0.8
)
ax2.set_xlabel('Line Index', fontsize=12)
ax2.set_ylabel('log₁₀(Einstein A) [s⁻¹]', fontsize=12)
ax2.set_title('CHIANTI Fe XII — Radiative Decay Rates', fontsize=14)
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('chianti_nist_comparison.png', dpi=150)
print(f"\n✅ Saved: chianti_nist_comparison.png")
plt.show()

print(f"\n📊 First 5 transitions:")
print(chianti_df[['wav', 'int', 'A']].head())