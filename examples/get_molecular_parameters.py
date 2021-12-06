# -*- coding: utf-8 -*-
"""
========================
Get Molecular Parameters
========================

Retrieve isotopologue abundances or molecular mass, based on HITRAN data,
using :py:class:`~radis.db.molparam.MolParams`

See also how to use custom (non-terrestrial) abundances in
:ref:`the Custom Abundances example <example_custom_abundances>`

"""

#%% Access the parameters:
from radis.db.molparam import MolParams

molpar = MolParams()

print("CO2 abundances for the first 2 isotopes :")
print(molpar.get("CO2", 1, "abundance"))  # 1 for isotopologue number
print(molpar.get("CO2", 2, "abundance"))  # 1 for isotopologue number


print("H2O molar mass for the first 2 isotopes")
print(molpar.get("H2O", 1, "mol_mass"))  # 1 for isotopologue number
print(molpar.get("CO2", 2, "mol_mass"))  # 1 for isotopologue number


print("Tabulated CO partition function at 296 K for the first 2 isotopes")
print(molpar.get("CO", 1, "Q_296K"))  # 1 for isotopologue number
print(molpar.get("CO", 2, "Q_296K"))  # 1 for isotopologue number


#%% All parameters are stored in a DataFrame:
print("All parameters: ", list(molpar.df.columns))
print(molpar.df)
