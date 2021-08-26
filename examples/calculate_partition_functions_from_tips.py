# -*- coding: utf-8 -*-
"""
=============================
Partition Functions from TIPS
=============================

By default and for equilibrium calculations, RADIS calculates Partition Functions
using the TIPS program ([TIPS-2020]_) through [HAPI]_. These partition functions can be retrieved
with the :py:class:`~radis.levels.partfunc.PartFuncTIPS` class::

See Also
--------
:py:class:`~radis.levels.partfunc.PartFunc_Dunham`
"""

from radis.db.classes import get_molecule_identifier
from radis.levels.partfunc import PartFuncHAPI

M = get_molecule_identifier("N2O")
iso = 1

Q = PartFuncHAPI(M, iso)
print(Q.at(T=1500))
