# -*- coding: utf-8 -*-
"""
================================================
Partition Functions from spectroscopic constants
================================================

RADIS can calculate equilibrium and non-LTE Partition Functions from a given
set of spectroscopic constants using a Dunham expansion.

Calculations use the the :py:class:`~radis.levels.partfunc.PartFunc_Dunham` class

Default spectroscopic constants and spectroscopic models used are given in
:ref:`default spectroscopic constants <label_db_spectroscopic_constants>`.
You can also use your own set of spectroscopic constants.

See Also
--------
:py:class:`~radis.levels.partfunc.PartFuncHAPI`

"""

from radis.db.molecules import Molecules
from radis.levels.partfunc import PartFunc_Dunham

isotope = 1
electronic_state = "X"
S = Molecules["CO"][isotope][electronic_state]

# Equilibrium partition functions :
Qf = PartFunc_Dunham(S)
print(Qf.at(T=3000))  # K

# Nonequilibrium partition functions :
print(Qf.at_noneq(Tvib=2000, Trot=1000))  # K
