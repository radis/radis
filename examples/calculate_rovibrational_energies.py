# -*- coding: utf-8 -*-
"""
================================
Calculate Rovibrational Energies
================================

RADIS can simply be used to calculate the rovibrational energies of molecules, using the
built-in :ref:`spectroscopic constants <label_db_spectroscopic_constants>` (or your own!).
See the :py:func:`~radis.db.molecules.getMolecule` function,
and the :py:data:`~radis.db.molecules.Molecules` list containing all :py:class:`~radis.db.classes.ElectronicState`
objects.

"""

from radis import getMolecule

# Here we get the energy of the v=6, J=3 level of the 2nd isotope of CO::

CO = getMolecule("CO", 2, "X")
print(CO.Erovib(6, 3))
