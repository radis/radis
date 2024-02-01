# -*- coding: utf-8 -*-
"""
=========================================================
Compare CO cross-sections from HITRAN, HITEMP, GEISA, and ExoMol
=========================================================

Auto-download and calculate CO spectrum from the HITRAN, HITEMP, GEISA, and ExoMol databases.
ExoMol references multiple databases for CO. Here we do not
use the ExoMol recommended database (see :py:func:`~radis.io.exomol.get_exomol_database_list`)
but we use the HITEMP database hosted on ExoMol.

Output should be similar, but no! By default these ExoMol does not provide
broadening coefficients for air. However, the Einstein coefficients & linestrengths
should be the same, therefore the integrals under the lines should be similar.

We verify this below :

For further exploration, ExoMol and HITEMP lines can be downloaded and accessed separately using
:py:func:`~radis.io.exomol.fetch_exomol` and :py:func:`~radis.io.hitemp.fetch_hitemp`

"""
import astropy.units as u

from radis import calc_spectrum, plot_diff

conditions = {
    "wmin": 2062 / u.cm,
    "wmax": 2093 / u.cm,
    "molecule": "CO",
    "isotope": "1",
    "pressure": 1.01325,  # bar
    "Tgas": 1000,  # K
    "mole_fraction": 0.1,
    "path_length": 1,  # cm
    "verbose": True,
    "neighbour_lines": 20,  # we account for the effect on neighbour_lines by computing ``20cm-1``
}
#%% Geisa VS HITRAN
s_geisa = calc_spectrum(**conditions, databank="geisa", name="GEISA")
s_hitran = calc_spectrum(
    **conditions,
    databank="hitran",
    name="HITRAN",
)
fig, [ax0, ax1] = plot_diff(s_geisa, s_hitran, "xsection", yscale="log")
# Adjust diff plot to be in linear scale
ax1.set_yscale("linear")
ax0.set_ylim(ymax=ax0.get_ylim()[1] * 10)  # more space for legend

# Note: these two spectra are alike.

#%% ExoMol VS HITEMP
s_exomol = calc_spectrum(
    **conditions, databank="exomol", name="ExoMol's HITEMP (default broadening)"
)
s_hitemp = calc_spectrum(
    **conditions,
    databank="hitemp",
    name="HITEMP (Air broadened)",
)
fig, [ax0, ax1] = plot_diff(s_exomol, s_hitemp, "xsection", yscale="log")
# Adjust diff plot to be in linear scale
ax1.set_yscale("linear")
ax0.set_ylim(ymax=ax0.get_ylim()[1] * 10)  # more space for legend

# Broadening coefficients are different in these databases, so lineshapes
# end up being very different; however the areas under the lines should be the same.
# We verify this :
print(
    f"Ratio of integrated area = {s_exomol.get_integral('xsection')/s_hitemp.get_integral('xsection'):.2f}"
)
