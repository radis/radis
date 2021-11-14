# -*- coding: utf-8 -*-
"""
.. _example_custom_abundances:

=====================
Use Custom Abundances
=====================

Custom isotopologue abundances can be defined, to model non-terrestrial
atmospheres.

Abundances are read and set directly in the :py:class:`~radis.lbl.factory.SpectrumFactory`
using the :py:meth:`~radis.lbl.loader.DatabankLoader.get_abundance`
and :py:meth:`~radis.lbl.loader.DatabankLoader.set_abundance` methods.

Below, we compute a CO2 spectrum with custom abundances for the two main
terrestrial isotopes (12C-16O2 ; 13C-16O2)

See Also
--------
:py:class:`~radis.db.molparam.MolParams`

"""

from radis.test.utils import setup_test_line_databases

setup_test_line_databases()  # creates "HITEMP-CO2-TEST" for this example

from radis import SpectrumFactory

sf = SpectrumFactory(
    2284.2,
    2284.6,
    wstep=0.001,  # cm-1
    pressure=20 * 1e-3,  # bar
    mole_fraction=400e-6,
    molecule="CO2",
    isotope="1,2",
    verbose=False,
)
sf.load_databank("HITEMP-CO2-TEST")
#%%To explicitely identify the isotopes we can use the molparam attribute the Factory
print(sf.molparam.get("CO2", 1, "isotope_name"))  # >> (12C)(16O)2
print(sf.molparam.get("CO2", 2, "isotope_name"))  # >> (13C)(16O)2


#%%Print the default abundance of the CO2 isotopes, compute a spectrum
print("Abundance of CO2[1,2]", sf.get_abundance("CO2", [1, 2]))
sf.eq_spectrum(2000).plot("abscoeff")

#%% Set the abundance of CO2(626) to 0.8; and the abundance of CO2(636) to 0.2 (arbitrary):
sf.set_abundance("CO2", [1, 2], [0.8, 0.2])
print("New abundance of CO2[1,2]", sf.get_abundance("CO2", [1, 2]))

sf.eq_spectrum(2000).plot("abscoeff", nfig="same")
