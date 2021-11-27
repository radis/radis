# -*- coding: utf-8 -*-
"""Definition of molecules and list of spectroscopic constants
"""

from .molecules import Molecules, getMolecule


# %% Get list of supported molecules
def _get_supported_molecules_equilibrium():
    """Molecules supported in RADIS equilibrium calculations Basically, all
    [HITRAN-2020]_ species."""
    # Add all HITRAN species
    # Add ExoMol species
    from .classes import EXOMOL_MOLECULES, HITRAN_MOLECULES

    return list(set(HITRAN_MOLECULES).union(set(EXOMOL_MOLECULES)))


def _get_supported_molecules_nonequilibrium():
    """Molecules supported in RADIS non equilibrium calculations without need
    for extra databases.

    Basically, molecules whose spectroscopic constants are built-in
    RADIS database (see radis.db)
    """

    return list(Molecules.keys())


MOLECULES_LIST_EQUILIBRIUM = (
    __supported_molecules_equilibrium__
) = _get_supported_molecules_equilibrium()
""" list: molecules that can be calculated in RADIS at equilibrium.
All [HITRAN-2020]_ and [ExoMol-2020]_ species are available.

Absorption coefficient calculated with RADIS at 300 K, 1 atm are shown for all
[HITRAN-2020]_ molecules in the :ref:`HITRAN spectra page <label_examples_hitran_spectra>` .

- 1 	``'H2O'`` : 	Water 	(`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/0%20-%20H2O%20infrared%20spectrum.png>`__)
- 2 	``'CO2'`` : 	Carbon Dioxide    (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/1%20-%20CO2%20infrared%20spectrum.png>`__)
- 3 	``'O3'`` : 	Ozone  (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/2%20-%20O3%20infrared%20spectrum.png>`__)
- 4 	``'N2O'`` : 	Nitrogen oxide 	  (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/3%20-%20N2O%20infrared%20spectrum.png>`__)
- 5 	``'CO'`` : 	Carbon Monoxide    (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/4%20-%20CO%20infrared%20spectrum.png>`__)
- 6 	``'CH4'`` : 	Methane   (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/5%20-%20CH4%20infrared%20spectrum.png>`__)
- 7 	``'O2'`` : 	Oxygen
- 8 	``'NO'`` : 	Nitric Oxide   (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/7%20-%20NO%20infrared%20spectrum.png>`__)
- 9 	``'SO2'`` : 	Sulfur Dioxide    (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/8%20-%20SO2%20infrared%20spectrum.png>`__)
- 10 	``'NO2'`` : 	Nitrogen Dioxide     (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/9%20-%20NO2%20infrared%20spectrum.png>`__)
- 11 	``'NH3'`` : 	Ammonia  (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/10%20-%20NH3%20infrared%20spectrum.png>`__)
- 12 	``'HNO3'`` : 	Nitric Acid     (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/11%20-%20HNO3%20infrared%20spectrum.png>`__)
- 13 	``'OH'`` : 	Hydroxyl  (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/12%20-%20OH%20infrared%20spectrum.png>`__)
- 14 	``'HF'`` : 	Hydrogen Fluoride     (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/13%20-%20HF%20infrared%20spectrum.png>`__)
- 15 	``'HCl'`` : 	Hydrogen Chloride    (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/14%20-%20HCl%20infrared%20spectrum.png>`__)
- 16 	``'HBr'`` : 	Hydrogen Bromide     (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/15%20-%20HBr%20infrared%20spectrum.png>`__)
- 17 	``'HI'`` : 	Hydrogen Iodide   (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/16%20-%20HI%20infrared%20spectrum.png>`__)
- 18 	``'ClO'`` : 	Chlorine Monoxide    (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/17%20-%20ClO%20infrared%20spectrum.png>`__)
- 19 	``'OCS'`` : 	Carbonyl Sulfide     (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/18%20-%20OCS%20infrared%20spectrum.png>`__)
- 20 	``'H2CO'`` : 	Formaldehyde    (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/19%20-%20H2CO%20infrared%20spectrum.png>`__)
- 21 	``'HOCl'`` : 	Hypochlorous Acid   (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/20%20-%20HOCl%20infrared%20spectrum.png>`__)
- 22 	``'N2'`` : 	Nitrogen
- 23 	``'HCN'`` : 	Hydrogen Cyanide
- 24 	``'CH3Cl'`` : 	Methyl Chloride    (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/23%20-%20CH3Cl%20infrared%20spectrum.png>`__)
- 25 	``'H2O2'`` : 	Hydrogen Peroxide   (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/24%20-%20H2O2%20infrared%20spectrum.png>`__)
- 26 	``'C2H2'`` : 	Acetylene   (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/25%20-%20C2H2%20infrared%20spectrum.png>`__)
- 27 	``'C2H6'`` : 	Ethane  (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/26%20-%20C2H6%20infrared%20spectrum.png>`__)
- 28 	``'PH3'`` : 	Phosphine    (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/27%20-%20PH3%20infrared%20spectrum.png>`__)
- 29 	``'COF2'`` : 	Carbonyl Fluoride   (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/28%20-%20COF2%20infrared%20spectrum.png>`__)
- 30 	``'SF6'`` : 	Sulfur Hexafluoride
- 31 	``'H2S'`` : 	Hydrogen Sulfide     (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/30%20-%20H2S%20infrared%20spectrum.png>`__)
- 32 	``'HCOOH'`` : 	Formic Acid    (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/31%20-%20HCOOH%20infrared%20spectrum.png>`__)
- 33 	``'HO2'`` : 	Hydroperoxyl     (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/32%20-%20HO2%20infrared%20spectrum.png>`__)
- 34 	``'O'`` : 	Oxygen Atom
- 35 	``'ClONO2'`` : 	Chlorine Nitrate
- 36 	``'NO+'`` : 	Nitric Oxide Cation  (`spectrum <https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/35%20-%20NO%2B%20infrared%20spectrum.png>`__)
- 37 	``'HOBr'`` : 	Hypobromous Acid
- 38 	``'C2H4'`` : 	Ethylene
- 39 	``'CH3OH'`` : 	Methanol
- 40 	``'CH3Br'`` : 	Methyl Bromide
- 41 	``'CH3CN'`` : 	Acetonitrile
- 42 	``'CF4'`` : 	CFC-14
- 43 	``'C4H2'`` : 	Diacetylene
- 44 	``'HC3N'`` : 	Cyanoacetylene
- 45 	``'H2'`` : 	Hydrogen
- 46 	``'CS'`` : 	Carbon Monosulfide
- 47 	``'SO3'`` : 	Sulfur trioxide
- 48 	``'C2N2'`` : 	Cyanogen
- 49 	``'COCl2'`` : 	Phosgene

The code to calculate these spectra is also available on the :ref:`example page <label_examples_hitran_spectra>`.

See Also
--------

:py:data:`~radis.db.molecules.Molecules`,
:py:func:`~radis.db.molecules.getMolecule`

"""

MOLECULES_LIST_NONEQUILIBRIUM = (
    __supported_molecules_nonequilibrium__
) = _get_supported_molecules_nonequilibrium()
""" list: molecules that can be calculated in RADIS at nonequilibrium.
Spectroscopic constants to calculate energy levels are needed.

RADIS features some built-in :ref:`spectroscopic constants <label_db_spectroscopic_constants>`
for the following species ([HITRAN-2020]_ nomenclature):

- 2 	``'CO2'`` : 	Carbon Dioxide
- 5 	``'CO'`` : 	Carbon Monoxide


See Also
--------

:py:data:`~radis.db.molecules.Molecules`,
:py:func:`~radis.db.molecules.getMolecule`
"""


__all__ = [
    "MOLECULES_LIST_EQUILIBRIUM",
    "MOLECULES_LIST_NONEQUILIBRIUM",
    "Molecules",
    "getMolecule",
]
