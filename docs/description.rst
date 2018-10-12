=====
RADIS
=====

A nonequilibrium infrared emission and absorption line-by-line code.

Description
-----------
    
Written as a general purpose radiative solver, the code is built around the [HITRAN-2016]_, 
[HITEMP-2010]_ and [CDSD-4000]_ databases for molecules in their electronic ground state. Energy 
levels are read from tabulated databases or calculated from Dunham developments. 
Boltzmann, Treanor, and state specific vibrational distributions can be generated. 
A modular architecture makes it possible to add new species without modifications 
to the core code. Thus far, |CO2|, CO are featured for non-equilibrium calculations, 
and all species present in the HITRAN database are featured for equilibrium 
calculations. To fit experimental spectra, RADIS includes a 
:py:class:`~radis.tools.line_survey.LineSurvey` tool, an 
interface with a look-up :py:class:`~radis.tools.database.SpecDatabase` 
to improve fitting convergence times, and a 
multi-slab module with a radiative transfer equation solver to reproduce line-of-sight 
experiments. Validation cases against existing spectral codes and experimental 
results from various plasma sources are included.

The code is available for use and modifications on `GitHub <https://github.com/radis/radis>`__
under a `GNU LESSER GENERAL PUBLIC LICENSE (v3) <https://github.com/radis/radis/blob/master/LICENSE>`__,
i.e., that modifications must remain public. 

.. |CO2| replace:: CO\ :sub:`2`

