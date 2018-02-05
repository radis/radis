# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 19:43:04 2015

@author: Erwan

Physical constants

Note : more units available in :
- scipy.constants
- sympy.physics.units

"""

# %% SI units

# eV in J
eV = 1.60217657e-19  # J

# Planck constant
h = 6.62606957e-34  # m2.kg.s-1 = J.s

# Boltzman constant
k_b = 1.3806488e-23  # m2.kg.s-2.K-1

# Light velocity
c = 2.99792458e8  # m/s

# First Bohr radius
a_0 = 5.29177e-11  # m

# Electron mass
m_e = 9.10938291e-31  # kg

# Rydberg unit of energy
Ry = 13.60569253  # eV

# Avogadro constant
Na = Av = 6.02214129e23  # mol-1

# Vacuum permittivity
eps_0 = 8.854187817e-12  # Farads / m

# cm to K conversion (used in Boltman factors)
hc_k = h*c/k_b*100          #   ~ 1.44 cm.K


# %% HITRAN (CGS) units

# Boltzman constant
k_b_CGS = 1.380648813E-16 # erg/K, CGS

# Light velocity
c_CGS = 2.99792458e10 # cm/s, CGS

# Planck constant
h_CGS = 6.626196e-27 # erg*s, CGS
