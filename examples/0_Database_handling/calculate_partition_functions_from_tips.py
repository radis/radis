# -*- coding: utf-8 -*-
"""
=============================
Partition Functions from TIPS
=============================

By default and for equilibrium calculations, RADIS calculates Partition Functions
using the TIPS program ([TIPS-2020]_) through [HAPI]_. These partition functions can be retrieved
with the :py:class:`~radis.levels.partfunc.PartFuncTIPS` class. See also :py:class:`~radis.levels.partfunc.PartFunc_Dunham`.
"""

from radis.db.classes import get_molecule_identifier
from radis.levels.partfunc import PartFuncHAPI

M = get_molecule_identifier("N2O")
iso = 1

Q = PartFuncHAPI(M, iso)
print(Q.at(T=1500))

"""
RADIS can also be used to compute Partition Functions from rovibrational energies
calculated with the built-in
:ref:`spectroscopic constants <label_db_spectroscopic_constants>`.

The calculation uses the :py:meth:`~radis.levels.partfunc.RovibParFuncCalculator.at`
method of the :py:class:`~radis.levels.partfunc.PartFunc_Dunham` class,
which reads :py:class:`~radis.db.molecules.Molecules`.
"""

from radis.db.molecules import Molecules
from radis.levels.partfunc import PartFunc_Dunham

iso = 1
electronic_state = "X"
S = Molecules["CO2"][iso][electronic_state]
Qf = PartFunc_Dunham(S)
print(Qf.at(T=3000))  # K

"""
Nonequilibrium partition functions can also be computed with
:py:meth:`~radis.levels.partfunc.RovibParFuncCalculator.at_noneq` ::
"""
print(Qf.at_noneq(Tvib=2000, Trot=1000))  # K

"""
:py:meth:`~radis.levels.partfunc.RovibParFuncCalculator.at_noneq`
can also return the vibrational partition function
and the table of rotational partition functions for each vibrational
state.
"""
