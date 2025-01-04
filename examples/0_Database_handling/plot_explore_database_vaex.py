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

In this example, we plot the dependance of the broadening coefficients to the rotational number and the line density.

"""
import matplotlib.pyplot as plt

from radis.io.hitran import fetch_hitran

#%% Adapt this example if vaex is installed or not
# June 2024: vaex is not compatible with python>=3.11, see https://github.com/radis/radis/pull/656 and the WIP solution here https://github.com/radis/radis/pull/698
try:
    import vaex

    output = "vaex"
    vaex  # fix linting
except ImportError:
    output = "pytables"  # standard pandas format

df = fetch_hitran("CO2", output=output, load_wavenum_min=100, load_wavenum_max=10000)
print(f"{len(df)} lines in HITEMP CO2; 2150 - 2450 cm-1")

#%%
if output == "vaex":
    # Note the use of `output='vaex'` in :py:func:`~radis.io.hitemp.fetch_hitemp` above.
    # The returned DataFrame is a Vaex DataFrame.
    # Loading times takes only few tens of milliseconds even for the largest HITEMP or ExoMol
    # databases
    #
    # We can also use Vaex graph functions.
    # See Vaex vizualisations : https://vaex.readthedocs.io/en/latest/guides/advanced_plotting.html#
    #

    plt.figure()
    df.viz.heatmap("jl", "airbrd", limits="99%")

    # Or below we plot the number of lines using Vaex's :py:meth:`~vaex.viz.DataFrameAccessorViz.histogram
    plt.figure()
    df.viz.histogram("wav", shape=1000)
    plt.yscale("log")

elif output == "pytables":

    df.plot.hexbin(x="jl", y="airbrd", gridsize=30)

    # Or below we plot the number of lines using pandas' histogram
    plt.figure()
    df["wav"].plot.hist(bins=100, logy=True)
