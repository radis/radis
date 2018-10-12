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

Run the example in `radis/test/validation/test_CO2_3Tvib_vs_klarenaar.py`

Line Survey
~~~~~~~~~~~


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



Interactive Examples
~~~~~~~~~~~~~~~~~~~~

RADIS in-the-browser sessions can be run from the `RADIS examples <https://github.com/radis/radis-examples>`_ project. 
