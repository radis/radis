# -*- coding: utf-8 -*-
"""
.. _example_download_hitemp:

============================
Download the HITEMP database
============================

Database will be downloaded automatically and can be edited locally.

To compute a spectrum with the HITEMP database, see the
:ref:`Calculate a HITEMP spectrum example <example_calculate_hitemp>`

By default the database is returned as a Pandas DataFrame.
To explore huge databases (like CO2, CH4 or H2O) that do not fit in RAM,
RADIS allows you to use a Vaex DataFrame instead (out-of-RAM).
See the :ref:`Explore Database with Vaex example <example_explore_database_vaex>`
"""

from radis.io.hitemp import fetch_hitemp

df = fetch_hitemp("OH")
print(df.columns)


#%%
# Returns:
# ::
#    Index(['id', 'iso', 'wav', 'int', 'A', 'airbrd', 'selbrd', 'El', 'Tdpair',
#           'Pshft', 'globu', 'globl', 'locu', 'locl', 'ierr', 'iref', 'lmix', 'gp',
#           'gpp'],
#          dtype='object')
#
# Columns are described in :py:attr:`~radis.io.hitran.columns_2004`
#
# A specific can be retrieved with the same :py:func:`~radis.io.hitemp.fetch_hitemp`
# function. The already downloaded database will be used:
#

df = fetch_hitemp("OH", load_wavenum_min=31500, load_wavenum_max=33000, isotope="1")
df.plot("wav", "int")
