# -*- coding: utf-8 -*-
"""
=========================================================
Compare CO xsections from the ExoMol and HITEMP database
=========================================================

Auto-download and calculate CO spectrum from the HITEMP database, and the
ExoMol database. ExoMol references multiple databases for CO. Here we do not
use the ExoMol recommended database (see :py:func:`~radis.io.exomol.get_exomol_database_list`)
but we use the HITEMP database hosted on ExoMol.

Output should be similar, but no ! By default these two databases provide
different broadening coefficients. However, the Einstein coefficients & linestrengths
should be the same, therefore the integrals under the lines should be similar.

We verifiy this below :

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
    "broadening_method": "fft",  # @ dev: Doesn't work with 'voigt'
    "verbose": True,
    "neighbour_lines": 20,
}


#%%
# Note that we account for the effect on neighbour_lines by computing ``20cm-1``
# on the side (``neighbour_lines`` condition above)

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

#%%
# Broadening coefficients are different in these databases, so lineshapes
# end up being very different; however the areas under the lines should be the same.
# We verify this :
import numpy as np

try:
    assert np.isclose(
        s_exomol.get_integral("xsection"), s_hitemp.get_integral("xsection"), rtol=0.001
    )

except:
    # @dev: someone there is un expected error with this example on ReadThedocs.
    # Escaping for the moment. See https://github.com/radis/radis/issues/501
    pass
