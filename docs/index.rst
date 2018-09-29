.. RADIS documentation master file, created by
   sphinx-quickstart on Tue Feb  6 03:32:15 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |logo_png| image:: radis_ico.png

.. include:: description.rst

===============
Getting Started
===============

.. toctree::
   :maxdepth: 2

   install


==================
User documentation
==================


Line-by-line (LBL) module
-------------------------

This is the core of RADIS: it calculates the spectral densities for a homogeneous
slab of gas, and returns a :class:`~radis.spectrum.spectrum.Spectrum` object. Calculations
are performed within the :class:`~radis.lbl.factory.SpectrumFactory` class. 

.. toctree::
   :maxdepth: 2
   
   lbl/index

   
Line-of-sight (LOS) module
--------------------------

This module takes several :class:`~radis.spectrum.spectrum.Spectrum` objects 
as input and combines then along the line-of-sight (:func:`~radis.los.slabs.SerialSlabs`) 
or at the same spatial position (:func:`~radis.los.slabs.MergeSlabs`), to reproduce 
line-of-sight experiments. The module allows combination of Spectra such as::

    s_line_of_sight = (s_plasma_CO2 // s_plasma_CO) > (s_room_absorption) 

.. toctree::
   :maxdepth: 2
   
   los/index

   
The Spectrum class
------------------

This module contains the :class:`~radis.spectrum.spectrum.Spectrum` object itself, with several methods that can be 
applied after the Spectrum was calculated: rescale, apply instrumental slit function, 
store or retrieve from a Spectrum database, plot or compare with another Spectrum object. 

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


Interfaces
----------

RADIS includes parsers and interfaces to read and return data in different formats: 

.. toctree::
   :maxdepth: 2
   
   io/parsers
   io/databases
   io/thermo


===============
Developer Guide
===============

   
.. toctree::
   :maxdepth: 2
   
   dev/test
   
   

=====================
Access Module Methods
=====================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



