
.. image:: https://img.shields.io/pypi/v/radis.svg
    :target: https://pypi.python.org/pypi/radis
    :alt: PyPI

.. image:: https://img.shields.io/travis/radis/radis.svg
    :target: https://travis-ci.org/radis/radis
    :alt: Tests
    
.. image:: https://codecov.io/gh/radis/radis/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/radis/radis
    :alt: Coverage
  
.. image:: https://readthedocs.org/projects/climt/badge/
    :target: https://radis.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
  

=====
RADIS
=====

A code to simulate infrared spectra of molecules.

.. warning::
    Deployment is still in progress, not all modules are available yet. The
    documentation will be updated as more modules are included. 

Documentation
-------------

User guide, install procedure and examples are available on the RADIS Website:

    http://radis.readthedocs.io/


Description
-----------
    
RADIS is nonequilibrium emission and absorption line-by-line code, for use 
by infrared spectroscopic that want to compare line databases, or experimentalist 
that want to fit their experimental line-of-sight spectra.

Written as a general purpose radiative solver, the code is built around the HITRAN, 
HITEMP and CDSD databases for molecules in their electronic ground state. Energy 
levels are read from tabulated databases or calculated from Dunham developments. 
Boltzmann, Treanor, and state specific vibrational distributions can be generated. 
A modular architecture makes it possible to add new species without modifications 
to the core code. Thus far, |CO2|, CO are featured for non-equilibrium calculations, 
and all species present in the HITRAN database are featured for equilibrium 
calculations. To fit experimental spectra, RADIS includes a line survey tool, an 
interface with a look-up database to improve fitting convergence times, and a 
multi-slab module with a radiative transfer equation solver to reproduce line-of-sight 
experiments. Validation cases against existing spectral codes and experimental 
results from various plasma sources are included.


License
-------

The code will soon be fully available on this repository under 
`GNU LESSER GENERAL PUBLIC LICENSE (v3) <./LICENSE>`_


Links
-----

- Documentation: http://radis.readthedocs.io/
- Source files: https://github.com/radis
- PyPi project: https://pypi.python.org/pypi/radis
- Test status: https://travis-ci.org/radis/radis
- Test coverage: https://codecov.io/gh/radis/radis


.. |CO2| replace:: CO\ :sub:`2`

