.. RADIS documentation master file, created by
   sphinx-quickstart on Tue Feb  6 03:32:15 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====
RADIS
=====

.. |logo_png| image:: radis_ico.png

A code to simulate infrared spectra of molecules.

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

===============
Getting Started
===============

.. toctree::
   :maxdepth: 2

   install


==================
User documentation
==================


The line-by-line (LBL) module
-----------------------------

This is the core of RADIS: it calculates the spectral densities for a homogeneous
slab of gas, and returns a Spectrum object. 

.. toctree::
   :maxdepth: 2
   
   lbl/index

   
The line-of-sight (LOS) module
------------------------------

This module takes several Spectrum objects as an input and combines then along the 
line-of-sight (SerialSlabs) or at the same spatial position (MergeSlabs), to 
reproduce line-of-sight experiments 

.. toctree::
   :maxdepth: 2
   
   los/index

   
The Spectrum class
------------------

This module contains the :class:`~radis.spectrum.spectrum.Spectrum` object itself, with several methods that can be 
applied after the Spectrum was calculated: rescale, apply instrumental slit function, 
store or retrieve from a Spectrum database. 

.. toctree::
   :maxdepth: 2
   
   spectrum/spectrum
   spectrum/howto

   
Tools
------------------

Different tools to work with Spectrum objects

.. toctree::
   :maxdepth: 2
   
   tools/line_survey
   tools/database



==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`




.. |CO2| replace:: CO\ :sub:`2`
.. |H2O| replace:: H\ :sub:`2`\ O
