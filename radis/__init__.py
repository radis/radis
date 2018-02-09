# -*- coding: utf-8 -*-
"""

RADIS 

A code to simulate infrared spectra of molecules.

---

RADIS is nonequilibrium emission and absorption line-by-line code, for use by 
infrared spectroscopic that want to compare line databases, or experimentalist 
that want to fit their experimental line-of-sight spectra.

Written as a general purpose radiative solver, the code is built around the HITRAN, 
HITEMP and CDSD databases for molecules in their electronic ground state. Energy 
levels are read from tabulated databases or calculated from Dunham developments. 
Boltzmann, Treanor, and state specific vibrational distributions can be 
generated. A modular architecture makes it possible to add new species without 
modifications to the core code. Thus far, CO2, CO are featured for non-equilibrium 
calculations, and all species present in the HITRAN database are featured for 
equilibrium calculations. To fit experimental spectra, RADIS includes a line 
survey tool, an interface with a look-up database to improve fitting convergence 
times, and a multi-slab module with a radiative transfer equation solver to 
reproduce line-of-sight experiments. Validation cases against existing spectral 
codes and experimental results from various plasma sources are presented. 

The code will soon be available on this repository under under GNU General Public 
License v3.0

"""

from __future__ import absolute_import, division, print_function, unicode_literals

# %% Debug mode
DEBUG_MODE = False
# change this at runtime with 
# >>> radis.DEBUG_MODE = True
# Use the printdbg() function in radis.misc, typically with:
# >>> if __debug__: printdbg(...)
# so that printdbg are removed by the Python preprocessor when running in 
# optimize mode:   
# >>> python -O *.py

from .spectrum import *        # Spectrum object
from .lbl import *             # line-by-line module 
from .los import *             # line-of-sight module
from .tools import *           # slit, database, line survey, etc.