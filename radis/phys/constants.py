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

eV2wn = 8065.54  # 1[eV]=8065.54[cm^-1] # from: https://github.com/HajimeKawahara/exojax/blob/78466cef0170ee1a2768b6a6f7b7c911d715c1bd/src/exojax/utils/constants.py#L35

h = 6.62607015e-34
"""float: Planck constant (m2.kg.s-1 = J.s)

https://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=planck"""

e_CGS = 4.80320450e-10
"""[esu]=[dyn^0.5*cm] #elementary charge"""

m_e_CGS = 9.10938356e-28
"""[g] !electron mass"""

k_b = 1.3806490e-23
"""float: Boltzmann constant (m2.kg.s-2.K-1)

https://physics.nist.gov/cgi-bin/cuu/Value?k|search_for=boltzmann
"""

c = 2.99792458e8
"""float: light velocity (m/s)

https://physics.nist.gov/cgi-bin/cuu/Value?c
"""

a_0 = 5.29177e-11
"""float: First Bohr radius (m)"""

m_e = 9.1093837015e-31
"""float: Electron mass (kg)

https://physics.nist.gov/cgi-bin/cuu/Value?me|search_for=electron+mass"""

Ry = 13.605693122994
"""float: Rydberg unit of energy (eV)

https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev|search_for=rydberg"""

Av = 6.02214076e23
"""float: Avogadro constant (molecules / mol)

https://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=avogadro"""
Na = Av  # just an alias

eps_0 = 8.8541878128e-12
"""float: Vacuum permittivity (Farads / m)

https://physics.nist.gov/cgi-bin/cuu/Value?ep0|search_for=vacuum+permittivity"""

hc_k = h * c / k_b * 100  #
"""float:  cm to K conversion (~ 1.44 K/cm), used in Boltzmann factors"""


# %% HITRAN (CGS) units

k_b_CGS = 1.380648813e-16
"""float: Boltzmann constant (erg / K)"""

c_CGS = 2.99792458e10
"""float: Light velocity (cm / s)"""

h_CGS = 6.626196e-27
"""float: Planck constant (erg * s)"""
