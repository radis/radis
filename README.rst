
*****************************************
`RADIS <https://radis.readthedocs.io/>`__
*****************************************

.. image:: https://img.shields.io/pypi/v/radis.svg
    :target: https://pypi.python.org/pypi/radis
    :alt: PyPI

.. image:: https://img.shields.io/badge/License-LGPL3-blue.svg
    :target: ./License
    :alt: License

.. image:: https://zenodo.org/badge/doi/10.1016/j.jqsrt.2018.09.027.svg
    :target: https://linkinghub.elsevier.com/retrieve/pii/S0022407318305867
    :alt: Article
    


.. image:: https://img.shields.io/travis/radis/radis.svg
    :target: https://travis-ci.org/radis/radis
    :alt: Tests
    
.. image:: https://codecov.io/gh/radis/radis/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/radis/radis
    :alt: Coverage
  
.. image:: https://readthedocs.org/projects/radis/badge/
    :target: https://radis.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb
    :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb

    

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/radis-radiation/community
    :alt: Gitter

    
    
[![Gitter]()](https://gitter.im/PlasmaPy/Lobby)
    
    
Description
-----------
    
RADIS is a nonequilibrium infrared emission and absorption line-by-line code.

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

The code is available for use and modifications on `GitHub <https://github.com/radis/radis>`__
under a `GNU LESSER GENERAL PUBLIC LICENSE (v3) <https://github.com/radis/radis/blob/master/LICENSE>`__,
i.e., that modifications must remain public. 

Documentation
-------------

User guide, install procedure and examples are available on the RADIS Website:

    http://radis.readthedocs.io/



Examples
--------

A `3-temperature fit <http://radis.readthedocs.io/en/latest/#multi-temperature-fit>`_ built on top of RADIS. 

.. figure:: https://raw.githubusercontent.com/radis/radis-examples/master/docs/multi-temperature-fit.gif

Interactive RADIS sessions can be run directly from the Browser, without installing Python. 
For instance, run the Quick Start example online:

.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb
    :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb


More examples can be found in the `documentation <http://radis.readthedocs.io/>`_ or in the 
`RADIS Inteactive examples <https://github.com/radis/radis-examples>`_ project. 


License
-------

The code is available on this repository under 
`GNU LESSER GENERAL PUBLIC LICENSE (v3) <./LICENSE>`_

.. image:: https://img.shields.io/badge/License-LGPL3-blue.svg
    :target: ./License
    :alt: License


Links
-----

RADIS:

- Documentation: http://radis.readthedocs.io/

.. image:: https://readthedocs.org/projects/radis/badge/
    :target: https://radis.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

- Source Code: https://github.com/radis/radis
- Article: https://linkinghub.elsevier.com/retrieve/pii/S0022407318305867

.. image:: https://zenodo.org/badge/doi/10.1016/j.jqsrt.2018.09.027.svg
    :target: https://linkinghub.elsevier.com/retrieve/pii/S0022407318305867
    :alt: Article

And also:

- Test Status: https://travis-ci.org/radis/radis

.. image:: https://img.shields.io/travis/radis/radis.svg
    :target: https://travis-ci.org/radis/radis
    :alt: Tests
    
- Test Coverage: https://codecov.io/gh/radis/radis

.. image:: https://codecov.io/gh/radis/radis/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/radis/radis
    :alt: Coverage
  
- PyPi Repository: https://pypi.org/project/radis/

.. image:: https://img.shields.io/pypi/v/radis.svg
    :target: https://pypi.python.org/pypi/radis
    :alt: PyPI

- Interactive Examples: https://github.com/radis/radis-examples

.. |CO2| replace:: CO\ :sub:`2`


Support
-------

If encountering any problem, first refer to the list of known 
<Issues `https://github.com/radis/radis/issues?utf8=%E2%9C%93&q=is%3Aissue`>__ on GitHub.
We appreciate your feedback and suggestions!

For any question, please join the discussion channel on Gitter:

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/radis-radiation/community
    :alt: Gitter

    
--------


.. image:: https://github.com/radis/radis/blob/master/docs/radis_ico.png
    :target: https://radis.readthedocs.io/
    :scale: 50 %
    :alt: RADIS logo
