*****
RADIS
***** 

A code to simulate infrared spectra of molecules.

RADIS is nonequilibrium emission and absorption line-by-line code, for use by 
infrared spectroscopic that want to compare line databases, or experimentalist 
that want to fit their experimental line-of-sight spectra.

- Docs: http://radis.readthedocs.io/
- Source: https://github.com/radis
- Pypi: https://pypi.python.org/pypi/radis

.. image:: https://img.shields.io/pypi/v/radis.svg
    :target: https://pypi.python.org/pypi/radis
    :alt: PyPI

.. image:: https://img.shields.io/travis/radis/radis.svg
    :target: https://travis-ci.org/radis/radis
    :alt: Continuous Integration
    
.. image:: https://codecov.io/gh/radis/radis/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/radis/radis
  
.. image:: https://readthedocs.org/projects/climt/badge/
    :target: https://radis.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
    
Description
-----------

Written as a general purpose radiative solver, the code is built around the HITRAN, 
HITEMP and CDSD databases for molecules in their electronic ground state. Energy 
levels are read from tabulated databases or calculated from Dunham developments. 
Boltzmann, Treanor, and state specific vibrational distributions can be 
generated. A modular architecture makes it possible to add new species without 
modifications to the core code. Thus far, |CO2|, CO are featured for non-equilibrium 
calculations, and all species present in the HITRAN database are featured for 
equilibrium calculations. To fit experimental spectra, RADIS includes a line 
survey tool, an interface with a look-up database to improve fitting convergence 
times, and a multi-slab module with a radiative transfer equation solver to 
reproduce line-of-sight experiments. Validation cases against existing spectral 
codes and experimental results from various plasma sources are included. 

The code will soon be fully available on this repository under 
`GNU LESSER GENERAL PUBLIC LICENSE (v3) <./LICENSE>`_

.. |CO2| replace:: CO\ :sub:`2`\ 
