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

For further exploratino, ExoMol and HITEMP lines can be downloaded and accessed separately using
:py:func:`~radis.io.exomol.fetch_exomol` and :py:func:`~radis.io.hitemp.fetch_hitemp`

"""
import astropy.units as u

from radis import calc_spectrum, plot_diff

conditions = {
    "wmin": 2002 / u.cm,
    "wmax": 2300 / u.cm,
    "molecule": "CO",
    "isotope": "2",
    "pressure": 1.01325,  # bar
    "Tgas": 1000,  # K
    "mole_fraction": 0.1,
    "path_length": 1,  # cm
    "broadening_method": "fft",  # @ dev: Doesn't work with 'voigt'
    "verbose": True,
}

s_exomol = calc_spectrum(
    **conditions, databank="exomol", name="ExoMol's HITEMP (H2 broadened)"
)
s_hitemp = calc_spectrum(
    **conditions,
    databank="hitemp",
    name="HITEMP (Air broadened)",
)
plot_diff(s_exomol, s_hitemp, "xsection")

#%% Broadening coefficients are different but areas under the lines should be the same :
import numpy as np

assert np.isclose(
    s_exomol.get_integral("xsection"), s_hitemp.get_integral("xsection"), rtol=0.001
)
