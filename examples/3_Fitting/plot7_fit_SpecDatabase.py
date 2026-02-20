# -*- coding: utf-8 -*-
"""
================================================================================
Example #7: Database fitting using SpecDatabase.fit_spectrum
================================================================================

This example demonstrates the automatic residual plotting feature of
:py:meth:`~radis.tools.database.SpecDatabase.fit_spectrum`.
"""

import os
import shutil

import numpy as np

from radis import SpectrumFactory

# 1. Setup paths and clean previous database
db_name = "radis_gallery_db"
if os.path.exists(db_name):
    shutil.rmtree(db_name)

# 2. Create a database using the inbuilt init_database feature
sf = SpectrumFactory(wavenum_min=2100, wavenum_max=2150, molecule="CO", verbose=0)
sf.load_databank("HITRAN-CO-TEST")
db = sf.init_database(db_name)

print("Precomputing database...")
for T in np.arange(1000, 2000, 100):
    for x in np.arange(0.01, 0.2, 0.02):
        sf.eq_spectrum(Tgas=T, mole_fraction=x)

# 3. Fit an experimental spectrum and plot residuals automatically
s_exp = sf.eq_spectrum(Tgas=1500, mole_fraction=0.1)

# Add noise to mimic an experimental spectrum (2% of max radiance)
s_exp = s_exp.take("radiance_noslit")
noise = np.random.normal(0, 0.02 * s_exp.max().value, len(s_exp))
s_exp = s_exp + noise

# Toggle plot=True/False to see the difference!
# This will find Tgas and mole_fraction as varying conditions and plot them.
db.fit_spectrum(s_exp, var="radiance_noslit", plot=True)
