# -*- coding: utf-8 -*-
"""
=========================================================
Compare NO cross-sections from HITRAN and ExoMol
=========================================================

Auto-download, calculate, and compare NO spectra from the HITRAN and GEISA databases

"""
import astropy.units as u

from radis import calc_spectrum, plot_diff

conditions = {
    "wmin": 1900 / u.cm,
    "wmax": 2100 / u.cm,
    "molecule": "NO",
    "isotope": "1",
    "pressure": 1.01325,  # bar
    "Tgas": 1000,  # K
    "mole_fraction": 0.1,
    "path_length": 1,  # cm
    "verbose": True,
    "neighbour_lines": 20,  # we account for the effect on neighbour_lines by computing ``20cm-1``
    "wstep": 0.0074,
}

s_geisa = calc_spectrum(**conditions, databank="geisa", name="GEISA")
s_hitran = calc_spectrum(**conditions, databank="hitran", name="HITRAN")

fig, [ax0, ax1] = plot_diff(
    s_geisa,
    s_hitran,
    "xsection",
    yscale="log",
)

# Adjust diff plot to be in linear scale
ax1.set_yscale("linear")
ax0.set_ylim(ymin=1e-24)  # to adjust the range of display

# The two spectra are different near 2000 cm-1. You can still check that the
# overal integrated cross-sections are similar (within 2% in March 2025).
# We verify this :
print(
    f"Ratio of integrated area = {s_geisa.get_integral('xsection')/s_hitran.get_integral('xsection'):.2f}"
)
