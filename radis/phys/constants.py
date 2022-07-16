# -*- coding: utf-8 -*-
"""Physical constants.

Note : more units available in :
- scipy.constants
- sympy.physics.units


-------------------------------------------------------------------------------
"""

# %% SI units

eV = 1.602176634e-19  # J
"""float: electron volt (J)

https://physics.nist.gov/cgi-bin/cuu/Value?e
"""

h = 6.62607015e-34
"""float: Planck constant (m2.kg.s-1 = J.s)

https://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=planck"""

k_b = 1.380649e-23
"""float: Boltzmann constant (m2.kg.s-2.K-1)

https://physics.nist.gov/cgi-bin/cuu/Value?k|search_for=boltzmann
"""

c = 2.99792458e8
"""float: light velocity (m/s)"""

a_0 = 5.29177e-11
"""float: First Bohr radius (m)"""

m_e = 9.10938291e-31
"""float: Electron mass (kg)"""

Ry = 13.60569253
"""float: Rydberg unit of energy (eV)"""

Av = 6.02214076e23
"""float: Avogadro constant (molecules / mol)

https://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=avogadro"""
Na = Av  # just an alias

eps_0 = 8.854187817e-12
"""float: Vacuum permittivity (Farads / m)"""

hc_k = h * c / k_b * 100  #
"""float:  cm to K conversion (~ 1.44 K/cm), used in Boltman factors"""


# %% HITRAN (CGS) units

k_b_CGS = 1.380648813e-16
"""float: Boltzman constant (erg / K)"""

c_CGS = 2.99792458e10
"""float: Light velocity (cm / s)"""

h_CGS = 6.626196e-27
"""float: Planck constant (erg * s)"""
