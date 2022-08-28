# -*- coding: utf-8 -*-
"""
.. _example_explore_database_vaex:

================================
Explore Line Database Parameters
================================

Database will be downloaded automatically and can be edited locally.

The :ref:`Download HITEMP Database example <example_download_hitemp>` showed
how to download the HITEMP database under Pandas format.

RADIS can also be used to explore and visualize larger databases, using the
underlying :py:mod:`vaex` library.

"""

from radis.io.hitemp import fetch_hitemp

df = fetch_hitemp("CO2", output="vaex", load_wavenum_min=2150, load_wavenum_max=2450)
print(f"{len(df)} lines in HITEMP CO2; 2150 - 2450 cm-1")

#%%
# Note the use of `output='vaex'` in :py:func:`~radis.io.hitemp.fetch_hitemp` above.
# The returned DataFrame is a Vaex DataFrame.
# Loading times takes only few tens of milliseconds even for the largest HITEMP or ExoMol
# databases
#
# We can also use Vaex graph functions.
# See Vaex vizualisations : https://vaex.readthedocs.io/en/latest/guides/advanced_plotting.html#
#
# For instance, we plot a :py:func:`~vaex.viz.DataFrameAccessorViz.heatmap` showing Self and Air-broadening coefficients

import matplotlib.pyplot as plt

plt.figure()
df.viz.heatmap("airbrd", "selbrd", limits="99%")

# TODO / idea : compare CO2 broadening with Air broadening for HITRAN database ?
#%%
# Or below we plot the number of lines using Vaex's :py:func:`~vaex.viz.DataFrameAccessorViz.histogram`

plt.figure()
df.viz.histogram("wav", shape=1000)
plt.yscale("log")
