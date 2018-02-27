=====
RADIS
=====

A code to simulate infrared spectra of molecules.

RADIS is nonequilibrium emission and absorption line-by-line code, for use 
by infrared spectroscopic that want to compare line databases, or experimentalist 
that want to fit their experimental line-of-sight spectra.

- Docs: http://radis.readthedocs.io/
- Source: https://github.com/radis/radis
- PyPi: https://pypi.python.org/pypi/radis

.. warning::
    Deployment is still in progress, not all modules are available yet. The
    documentation will be updated as more modules are included. 

Description
-----------
    
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

The code will soon be fully available on this repository under 
`GNU LESSER GENERAL PUBLIC LICENSE (v3) <https://github.com/radis/radis/blob/master/LICENSE>`_

*Example of input produced by the* :class:`~radis.tools.line_survey.LineSurvey` *tool.*

.. raw:: html

    <iframe id="igraph" src="https://plot.ly/~erwanp/6/" width="700" height="450" seamless="seamless" scrolling="no"></iframe>
	
.. |CO2| replace:: CO\ :sub:`2`
.. |H2O| replace:: H\ :sub:`2`\ O