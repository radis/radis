# RADIS 

A code to simulate infrared spectra of molecules.

---

- Docs: http://radis.readthedocs.io/

---

RADIS is nonequilibrium emission and absorption line-by-line code, for use by 
infrared spectroscopic that want to compare line databases, or experimentalist 
that want to fit their experimental line-of-sight spectra.

Written as a general purpose radiative solver, the code is built around the HITRAN, 
HITEMP and CDSD databases for molecules in their electronic ground state. Energy 
levels are read from tabulated databases or calculated from Dunham developments. 
Boltzmann, Treanor, and state specific vibrational distributions can be 
generated. A modular architecture makes it possible to add new species without 
modifications to the core code. Thus far, CO<sub>2</sub>, CO are featured for non-equilibrium 
calculations, and all species present in the HITRAN database are featured for 
equilibrium calculations. To fit experimental spectra, RADIS includes a line 
survey tool, an interface with a look-up database to improve fitting convergence 
times, and a multi-slab module with a radiative transfer equation solver to 
reproduce line-of-sight experiments. Validation cases against existing spectral 
codes and experimental results from various plasma sources are included. 

The code will soon be available on this repository under under GNU General Public 
License v3.0
