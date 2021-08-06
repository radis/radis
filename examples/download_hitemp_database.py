# -*- coding: utf-8 -*-
"""
============================
Download the HITEMP database
============================

Database will be downloaded automatically and can be edited locally.

To compute a spectrum with the HITEMP database, see the
:ref:`Calculate a HITEMP spectrum example <example_calculate_hitemp>`

"""

from radis.io.hitemp import fetch_hitemp

df = fetch_hitemp("CO")
print(df.columns)


"""
Returns:
::
    Index(['id', 'iso', 'wav', 'int', 'A', 'airbrd', 'selbrd', 'El', 'Tdpair',
       'Pshft', 'ierr', 'iref', 'lmix', 'gp', 'gpp', 'Fu', 'branch', 'jl',
       'syml', 'Fl', 'vu', 'vl'],
      dtype='object')

Columns are described in :py:attr:`~radis.io.hitran.column_2004`
"""
