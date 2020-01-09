.. _label_examples:

========
Examples
========

Many examples scripts are available on the `radis-examples project <https://github.com/radis/radis-examples>`__. 



Line Survey
===========


Example of input produced by the :class:`~radis.tools.line_survey.LineSurvey` tool::

    from radis import SpectrumFactory
    sf = SpectrumFactory(
                         wavenum_min=2380,
                         wavenum_max=2400,
                         mole_fraction=400e-6,
                         path_length=100,  # cm
                         isotope=[1],
                         ) 
    sf.load_databank('HITRAN-CO2-TEST')
    s = sf.eq_spectrum(Tgas=1500)
    s.apply_slit(0.5)
    s.line_survey(overlay='radiance_noslit', barwidth=0.01)


.. raw:: html

    <iframe id="igraph" src="https://plot.ly/~erwanp/6/" width="650" height="420" seamless="seamless" scrolling="no"></iframe>
	
.. |CO2| replace:: CO\ :sub:`2`
.. |H2O| replace:: H\ :sub:`2`\ O



RADIS in-the-browser
====================

RADIS in-the-browser sessions can be run from the `RADIS examples <https://github.com/radis/radis-examples>`_ project.
No installation needed, you don't even need Python on your computer. 

For example, run the Quick Start :ref:`first example <label_first_example>` by clicking on the link below:

.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb
    :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb

Or start a bare RADIS online session:
    
.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb
    :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb

The full list can be found on the `RADIS Interactive Examples <https://github.com/radis/radis-examples#interactive-examples>`_ project. 


Get rovibrational energies
==========================

RADIS can simply be used to calculate the rovibrational energies of molecules, using the 
built-in :ref:`spectroscopic constants <label_db_spectroscopic_constants>`. 
See the :py:func:`~radis.db.molecules.getMolecule` function,  
and the :py:data:`~radis.db.molecules.Molecules` list containing all :py:class:`~radis.db.classes.ElectronicState` 
objects. 

Here we get the energy of the asymmetric mode of CO2::

    from radis import getMolecule
    CO2 = getMolecule('CO2', 1, 'X')
    print(CO2.Erovib(0, 0, 0, 1, 0))
    >>> 2324.2199999

Here we get the energy of the v=6, J=3 level of the 2nd isotope of CO::

    CO = getMolecule('CO', 2, 'X')
    print(CO.Erovib(6, 3))
    >>> 12218.8130906978


.. _label_examples_multitemperature_fit:

Multi Temperature Fit
=====================

A `3 temperature fitting example <https://github.com/radis/radis-examples/tree/master/multi-temperature-fit>`__ .
reproducing the validation case of Klarenaar 2017 [1]_, who calculated a transmittance
spectrum from the initial data of Dang 1973 [2]_, with a 1 rotational temperature + 
3 vibrational temperature (Treanor distributions) model. 

|CO2| Energies are calculated from Dunham developments in an uncoupled harmonic oscillator - rigid rotor model. 
The example is based on one of 
`RADIS validation cases <https://github.com/radis/radis/blob/master/radis/test/validation/test_CO2_3Tvib_vs_klarenaar.py>`_.
It makes use of the RADIS `Spectrum <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`_
class and the associated compare and load functions

.. only:: html

   .. figure:: https://raw.githubusercontent.com/radis/radis-examples/master/docs/multi-temperature-fit.gif
      :alt: CO2 multi temperature fitting
      :target: https://raw.githubusercontent.com/radis/radis-examples/master/docs/multi-temperature-fit.gif

The fitting script can be found in the 
`radis-examples project <https://github.com/radis/radis-examples/tree/master/multi-temperature-fit>`__ .

.. [1] Klarenaar et al 2017, "Time evolution of vibrational temperatures in a CO2 glow 
       discharge measured with infrared absorption spectroscopy" doi/10.1088/1361-6595/aa902e

.. [2] Dang et al 1982, "Detailed vibrational population distributions in a CO2 laser 
        discharge as measured with a tunable diode laser" doi/10.1007/BF00694640


 .. _label_examples_hitran_spectra:
 
HITRAN spectra
==============

The absorption coefficient of all HITRAN species (see :py:data:`~radis.io.MOLECULES_LIST_EQUILIBRIUM`)
is calculated in `plot_all_hitran_spectra.py <https://github.com/radis/radis-examples/blob/master/hitran_spectra/plot_all_hitran_spectra.py>`__ 
at 300 K, 1 atm for the first isotope.

For instance:

- 1 	``'H2O'`` : 	Water absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='H2O',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/0%20-%20H2O%20infrared%20spectrum.png
   :width: 600
   :alt: Water H2O infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/0%20-%20H2O%20infrared%20spectrum.png

- 2 	``'CO2'`` : 	Carbon Dioxide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='CO2',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')


.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/1%20-%20CO2%20infrared%20spectrum.png
   :width: 600
   :alt: Carbon Dioxide CO2 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/1%20-%20CO2%20infrared%20spectrum.png

- 3 	``'O3'`` : 	Ozone absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='O3',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')


.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/2%20-%20O3%20infrared%20spectrum.png
   :width: 600
   :alt: Ozone O3 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/2%20-%20O3%20infrared%20spectrum.png



- 4 	``'N2O'`` : 	Nitrogen oxide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='N2O',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')


.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/3%20-%20N2O%20infrared%20spectrum.png
   :width: 600
   :alt: Nitrogen oxide N2O infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/3%20-%20N2O%20infrared%20spectrum.png



- 5 	``'CO'`` : 	Carbon Monoxide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='CO',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')


.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/4%20-%20CO%20infrared%20spectrum.png
   :width: 600
   :alt: Carbon Monoxide CO infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/4%20-%20CO%20infrared%20spectrum.png


- 6 	``'CH4'`` : 	Methane absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='CH4',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/5%20-%20CH4%20infrared%20spectrum.png
   :width: 600
   :alt: Methane CH4 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/5%20-%20CH4%20infrared%20spectrum.png


- 7 	``'O2'`` : 	Oxygen absorption coefficient at 300 K : no lines for ``isotope='1'`` (symmetric!)
- 8 	``'NO'`` : 	Nitric Oxide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='NO',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/7%20-%20NO%20infrared%20spectrum.png
   :width: 600
   :alt: Nitric Oxide NO infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/7%20-%20NO%20infrared%20spectrum.png


- 9 	``'SO2'`` : 	Sulfur Dioxide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='SO2',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/8%20-%20SO2%20infrared%20spectrum.png
   :width: 600
   :alt: Sulfur Dioxide SO2 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/8%20-%20SO2%20infrared%20spectrum.png


- 10 	``'NO2'`` : 	Nitrogen Dioxide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='NO2',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/9%20-%20NO2%20infrared%20spectrum.png
   :width: 600
   :alt: Nitrogen Dioxide NO2 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/9%20-%20NO2%20infrared%20spectrum.png


- 11 	``'NH3'`` : 	Ammonia absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='NH3',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/10%20-%20NH3%20infrared%20spectrum.png
   :width: 600
   :alt: Ammonia NH3 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/10%20-%20NH3%20infrared%20spectrum.png


- 12 	``'HNO3'`` : 	Nitric Acid absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='HNO3',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/11%20-%20HNO3%20infrared%20spectrum.png
   :width: 600
   :alt: Nitric Acid HNO3 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/11%20-%20HNO3%20infrared%20spectrum.png


- 13 	``'OH'`` : 	Hydroxyl absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='OH',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/12%20-%20OH%20infrared%20spectrum.png
   :width: 600
   :alt: Hydroxyl OH infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/12%20-%20OH%20infrared%20spectrum.png


- 14 	``'HF'`` : 	Hydrogen Fluoride absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='HF',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/13%20-%20HF%20infrared%20spectrum.png
   :width: 600
   :alt: Hydrogen Fluoride HF infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/13%20-%20HF%20infrared%20spectrum.png


- 15 	``'HCl'`` : 	Hydrogen Chloride absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='HCl',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/14%20-%20HCl%20infrared%20spectrum.png
   :width: 600
   :alt: Hydrogen Chloride HCl infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/14%20-%20HCl%20infrared%20spectrum.png


- 16 	``'HBr'`` : 	Hydrogen Bromide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='HBr',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/15%20-%20HBr%20infrared%20spectrum.png
   :width: 600
   :alt: Hydrogen Bromide HBr infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/15%20-%20HBr%20infrared%20spectrum.png


- 17 	``'HI'`` : 	Hydrogen Iodide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='HI',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/16%20-%20HI%20infrared%20spectrum.png
   :width: 600
   :alt: Hydrogen Iodide HI infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/16%20-%20HI%20infrared%20spectrum.png


- 18 	``'ClO'`` : 	Chlorine Monoxide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='ClO',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/17%20-%20ClO%20infrared%20spectrum.png
   :width: 600
   :alt: Chlorine Monoxide ClO infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/17%20-%20ClO%20infrared%20spectrum.png


- 19 	``'OCS'`` : 	Carbonyl Sulfide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='OCS',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/18%20-%20OCS%20infrared%20spectrum.png
   :width: 600
   :alt: Carbonyl Sulfide OCS infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/18%20-%20OCS%20infrared%20spectrum.png


- 20 	``'H2CO'`` : 	Formaldehyde absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='H2CO',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/19%20-%20H2CO%20infrared%20spectrum.png
   :width: 600
   :alt: Formaldehyde H2CO infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/19%20-%20H2CO%20infrared%20spectrum.png


- 21 	``'HOCl'`` : 	Hypochlorous Acid absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='HOCl',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/20%20-%20HOCl%20infrared%20spectrum.png
   :width: 600
   :alt: Hypochlorous Acid HOCl infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/20%20-%20HOCl%20infrared%20spectrum.png


- 22 	``'N2'`` : 	Nitrogen absorption coefficient at 300 K : no lines for ``isotope='1'`` (symmetric!)
- 23 	``'HCN'`` : 	Hydrogen Cyanide absorption coefficient at 300 K : not calculated
- 24 	``'CH3Cl'`` : 	Methyl Chloride absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='CH3Cl',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/23%20-%20CH3Cl%20infrared%20spectrum.png
   :width: 600
   :alt: Methyl Chloride CH3Cl infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/23%20-%20CH3Cl%20infrared%20spectrum.png


- 25 	``'H2O2'`` : 	Hydrogen Peroxide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='H2O2',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/24%20-%20H2O2%20infrared%20spectrum.png
   :width: 600
   :alt: Hydrogen Peroxide H2O2 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/24%20-%20H2O2%20infrared%20spectrum.png


- 26 	``'C2H2'`` : 	Acetylene absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='C2H2',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/25%20-%20C2H2%20infrared%20spectrum.png
   :width: 600
   :alt: Acetylene C2H2 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/25%20-%20C2H2%20infrared%20spectrum.png


- 27 	``'C2H6'`` : 	Ethane absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='C2H6',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/26%20-%20C2H6%20infrared%20spectrum.png
   :width: 600
   :alt: Ethane C2H6 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/26%20-%20C2H6%20infrared%20spectrum.png


- 28 	``'PH3'`` : 	Phosphine absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='PH3',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/27%20-%20PH3%20infrared%20spectrum.png
   :width: 600
   :alt: Phosphine PH3 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/27%20-%20PH3%20infrared%20spectrum.png


- 29 	``'COF2'`` : 	Carbonyl Fluoride absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='COF2',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/28%20-%20COF2%20infrared%20spectrum.png
   :width: 600
   :alt: Carbonyl Fluoride COF2 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/28%20-%20COF2%20infrared%20spectrum.png


- 30 	``'SF6'`` : 	Sulfur Hexafluoride absorption coefficient at 300 K : not calculated
- 31 	``'H2S'`` : 	Hydrogen Sulfide absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='H2S',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/30%20-%20H2S%20infrared%20spectrum.png
   :width: 600
   :alt: Hydrogen Sulfide H2S infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/30%20-%20H2S%20infrared%20spectrum.png


- 32 	``'HCOOH'`` : 	Formic Acid absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='HCOOH',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/31%20-%20HCOOH%20infrared%20spectrum.png
   :width: 600
   :alt: Formic Acid HCOOH infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/31%20-%20HCOOH%20infrared%20spectrum.png


- 33 	``'HO2'`` : 	Hydroperoxyl absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='HO2',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/32%20-%20HO2%20infrared%20spectrum.png
   :width: 600
   :alt: Hydroperoxyl HO2 infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/32%20-%20HO2%20infrared%20spectrum.png


- 34 	``'O'`` : 	Oxygen Atom absorption coefficient at 300 K : not calculated
- 35 	``'ClONO2'`` : 	Chlorine Nitrate absorption coefficient at 300 K : not calculated
- 36 	``'NO+'`` : 	Nitric Oxide Cation absorption coefficient at 300 K ::

    s = calc_spectrum(wavelength_min=1000, 
                      wavelength_max=20000,
                      Tgas=300,
                      pressure=1,
                      molecule='NO+',
                      lineshape_optimization=None,
                      cutoff=1e-23,
                      isotope='1')
    s.plot('abscoeff', wunit='nm')
 

.. image:: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/35%20-%20NO%2B%20infrared%20spectrum.png
   :width: 600
   :alt: Nitric Oxide Cation NO+ infrared absorption coefficient
   :target: https://raw.githubusercontent.com/radis/radis-examples/master/hitran_spectra/out/35%20-%20NO%2B%20infrared%20spectrum.png


- 37 	``'HOBr'`` : 	Hypobromous Acid absorption coefficient at 300 K : not calculated
- 38 	``'C2H4'`` : 	Ethylene absorption coefficient at 300 K : not calculated
- 39 	``'CH3OH'`` : 	Methanol absorption coefficient at 300 K : not calculated
- 40 	``'CH3Br'`` : 	Methyl Bromide absorption coefficient at 300 K : not calculated
- 41 	``'CH3CN'`` : 	Acetonitrile absorption coefficient at 300 K : not calculated
- 42 	``'CF4'`` : 	CFC-14 absorption coefficient at 300 K : not calculated
- 43 	``'C4H2'`` : 	Diacetylene absorption coefficient at 300 K : not calculated
- 44 	``'HC3N'`` : 	Cyanoacetylene absorption coefficient at 300 K : not calculated
- 45 	``'H2'`` : 	Hydrogen absorption coefficient at 300 K : not calculated
- 46 	``'CS'`` : 	Carbon Monosulfide absorption coefficient at 300 K : not calculated
- 47 	``'SO3'`` : 	Sulfur trioxide absorption coefficient at 300 K : not calculated
- 48 	``'C2N2'`` : 	Cyanogen absorption coefficient at 300 K : not calculated
- 49 	``'COCl2'`` : 	Phosgene absorption coefficient at 300 K : not calculated





CH4 Full Spectrum Benchmark
===========================

Here we reproduce the full spectrum (0.001 - 11500 cm-1) of Methane for a ``broadening_max_width`` 
corresponding to about 50 HWHMs, as in the Benchmark case of [HAPI]_, Table 7, Methane_III,
also featured in the [RADIS-2018]_ article ::

    from radis import SpectrumFactory
    
    benchmark_line_brd_ratio = 50    # “WavenumberWingHW”/HWHMs
    dnu = 0.01         # step in HAPI Benchmark article
    molecule = 'CH4'
    wavenum_min = 0.001
    wavenum_max = 11505
    pressure_bar = 1.01315
    T = 296
    isotopes = [1, 2, 3, 4]
    
    sf = SpectrumFactory(wavenum_min=wavenum_min,
                         wavenum_max=wavenum_max,
                         isotope=isotopes,  #'all',
                         verbose=2,
                         wstep=dnu,     # depends on HAPI benchmark. 
                         cutoff=1e-23,  
                         broadening_max_width=5.73,  # Corresponds to WavenumberWingHW/HWHM=50 in HAPI
                         molecule=molecule,
                         )
    sf.fetch_databank('astroquery', load_energies=False)
    
    s = sf.eq_spectrum(Tgas=T, pressure=pressure_bar)
    s.plot()


The comparison in terms of performance with HAPI can be found in the ``radis/test/benchmark/radis_vs_hapi_CH4_full_spectrum.py`` 
case::

    cd radis
    python radis/test/benchmark/radis_vs_hapi_CH4_full_spectrum.py 

Using the different :ref:`Performance <label_lbl_performance>` optimisations available in RADIS, 
the calculation is typically 100 times faster in RADIS::

    >>> Calculated with HAPI in 157.41s
    >>> Calculated with RADIS in 1.65s

