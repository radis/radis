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

Examples
--------


Multi Temperature Fit
~~~~~~~~~~~~~~~~~~~~~

A 3 temperature fitting example reproducing the validation case of Klarenaar 2017 [1]_, who calculated a transmittance
spectrum from the initial data of Dang 1973 [2]_, with a 1 rotational temperature + 
3 vibrational temperature (Treanor distributions) model. 

.. [1] Klarenaar et al 2017, "Time evolution of vibrational temperatures in a CO2 glow 
       discharge measured with infrared absorption spectroscopy" doi/10.1088/1361-6595/aa902e

.. [2] Dang et al 1982, "Detailed vibrational population distributions in a CO2 laser 
        discharge as measured with a tunable diode laser" doi/10.1007/BF00694640

|CO2| Energies are calculated from Dunham developments in an uncoupled harmonic oscillator - rigid rotor model. 
The example is based on one of `RADIS validation cases <https://github.com/radis/radis/tree/master/radis/test/validation>`_.
It makes use of the RADIS `Spectrum <http://radis.readthedocs.io/en/latest/#the-spectrum-class>`_
class and the associated compare and load functions

.. only:: html

   .. figure:: https://raw.githubusercontent.com/radis/radis-examples/master/docs/multi-temperature-fit.gif

      3 temperature fit in RADIS


More examples can be found in the `RADIS examples <https://github.com/radis/radis-examples>`_ project. 


Line Survey
~~~~~~~~~~~


Example of input produced by the :class:`~radis.tools.line_survey.LineSurvey` tool:

.. raw:: html

    <iframe id="igraph" src="https://plot.ly/~erwanp/6/" width="700" height="450" seamless="seamless" scrolling="no"></iframe>
	
.. |CO2| replace:: CO\ :sub:`2`
.. |H2O| replace:: H\ :sub:`2`\ O



Interactive Examples
~~~~~~~~~~~~~~~~~~~~

RADIS in-the-browser sessions can be run from the `RADIS examples <https://github.com/radis/radis-examples>`_ project. 
