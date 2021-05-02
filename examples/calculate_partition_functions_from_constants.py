# -*- coding: utf-8 -*-
"""
================================================
Partition Functions from spectroscopic constants
================================================

By default and for equilibrium calculations, RADIS calculates Partition Functions
using the TIPS program through [HAPI]_. These partition functions can be retrieved
with the :py:class:`~radis.levels.partfunc.PartFunc_Dunham` class::

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
