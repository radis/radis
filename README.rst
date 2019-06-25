
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

*****************************************
`RADIS <https://radis.readthedocs.io/>`__
*****************************************

A nonequilibrium infrared emission and absorption line-by-line code &
a post-processing library to compare experimental and calculated spectra.

User guide, install procedure and examples are available on the `RADIS Website <http://radis.readthedocs.io/>`__:

.. image:: https://readthedocs.org/projects/radis/badge/
    :target: https://radis.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


===============
Getting Started
===============

Install
-------

Assuming you have Python installed with the `Anaconda <https://www.anaconda.com/download/>`_ distribution just use::

    pip install radis 
    
**That's it!** You can now run your first example below.
If you encounter any issue, or to upgrade the package later, please refer to the 
`detailed installation procedure <https://radis.readthedocs.io/en/latest/install.html#label-install>`__ . 

Quick Start
-----------


Calculate a CO equilibrium spectrum from the HITRAN database, using the
`calc_spectrum <https://radis.readthedocs.io/en/latest/source/radis.lbl.calc.html#radis.lbl.calc.calc_spectrum>`__ function. Output is a 
`Spectrum object <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`__: ::

    from radis import calc_spectrum
    s = calc_spectrum(1900, 2300,         # cm-1
                      molecule='CO',
                      isotope='1,2,3',
                      pressure=1.01325,   # bar
                      Tgas=700,           # K
                      mole_fraction=0.1, 
                      path_length=1,      # cm
                      )
    s.apply_slit(0.5, 'nm')       # simulate an experimental slit
    s.plot('radiance')

.. figure:: https://radis.readthedocs.io/en/latest/_images/co_spectrum_700K.png
    :scale: 60 %

Calculate a CO *nonequilibrium* spectrum from the HITRAN database
(on your first call, this will calculate and cache the CO(X) rovibrational
energies): ::

    s2 = calc_spectrum(1900, 2300,         # cm-1
                      molecule='CO',
                      isotope='1,2,3',
                      pressure=1.01325,   # bar
                      Tvib=700,           # K
                      Trot=300,           # K
                      mole_fraction=0.1, 
                      path_length=1,      # cm
                      )
    s2.apply_slit(0.5, 'nm')
    s2.plot('radiance', nfig='same')    # compare with previous
    
The Quick Start examples automatically download the line databases from `HITRAN-2016 <https://radis.readthedocs.io/en/latest/bibliography.html#hitran-2016>`__, which is valid for temperatures below 700 K. 
For *high temperature* cases, you may need to use other line databases such as 
`HITEMP-2010 <https://radis.readthedocs.io/en/latest/bibliography.html#hitemp-2010>`__ (typically T < 2000 K) or `CDSD-4000 <https://radis.readthedocs.io/en/latest/bibliography.html#cdsd-4000>`__ (T < 5000 K). These databases must be described in a ``~/.radis`` 
`Configuration file <https://radis.readthedocs.io/en/latest/lbl/index.html#label-lbl-config-file>`__. 

More complex `examples <https://radis.readthedocs.io/en/latest/examples.html#label-examples>`__ will require to use the `SpectrumFactory <https://radis.readthedocs.io/en/latest/source/radis.lbl.factory.html#radis.lbl.factory.SpectrumFactory>`__
class, which is the core of RADIS line-by-line calculations. 
`calc_spectrum <https://radis.readthedocs.io/en/latest/source/radis.lbl.calc.html#radis.lbl.calc.calc_spectrum>`__ is a wrapper to `SpectrumFactory <https://radis.readthedocs.io/en/latest/source/radis.lbl.factory.html#radis.lbl.factory.SpectrumFactory>`__
for the simple cases. 

Experimental spectra can be loaded using the `experimental_spectrum <https://radis.readthedocs.io/en/latest/source/radis.spectrum.models.html#radis.spectrum.models.experimental_spectrum>`__ function 
and compared with the `plot_diff <https://radis.readthedocs.io/en/latest/source/radis.spectrum.compare.html#radis.spectrum.compare.plot_diff>`__ function. For instance::

    from numpy import loadtxt
    from radis import experimental_spectrum, plot_diff
    w, I = loadtxt('my_file.txt').T    # assuming 2 columns 
    sexp = experimental_spectrum(w, I, Iunit='mW/cm2/sr/nm')
    plot_diff(sexp, s)    # comparing with previously spectrum 's' calculated previously 

Typical output of `plot_diff <https://radis.readthedocs.io/en/latest/source/radis.spectrum.compare.html#radis.spectrum.compare.plot_diff>`__:

.. image:: https://radis.readthedocs.io/en/latest/_images/cdsd4000_vs_hitemp_3409K.svg
    :scale: 60 %
    :target: https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#compare-two-spectra
    :alt: https://radis.readthedocs.io/en/latest/_images/cdsd4000_vs_hitemp_3409K.svg

Refer to the `Examples <https://radis.readthedocs.io/en/latest/examples.html#label-examples>`__ section for more examples, and to the  
`Spectrum page <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`__ for more post-processing functions. 

In the browser (no installation needed!)
----------------------------------------

Alternatively, you can also run RADIS directly in the browser with the  
`RADIS Interactive Examples <https://github.com/radis/radis-examples#interactive-examples>`_ project. 
For instance, run the Quick Start example on the link below:

.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb
    :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb

Or start a bare RADIS online session:
    
.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb
    :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb


---------------------------------------------------------------------

===============
Developer Guide
===============

Architecture
------------

RADIS internals are described in the `Developer Guide <https://radis.readthedocs.io/en/latest/developer.html>`__ :

.. image:: https://radis.readthedocs.io/en/latest/_images/RADIS_flow_chart.svg
     :target:   https://radis.readthedocs.io/en/latest/dev/architecture.html#label-dev-architecture
     :alt: https://radis.readthedocs.io/en/latest/_images/RADIS_flow_chart.svg


License
-------

The code is available on this repository under 
`GNU LESSER GENERAL PUBLIC LICENSE (v3) <./LICENSE>`_

.. image:: https://img.shields.io/badge/License-LGPL3-blue.svg
    :target: ./License
    :alt: License



Support
-------

If encountering any problem, first refer to the list of known 
`Issues <https://github.com/radis/radis/issues?utf8=%E2%9C%93&q=is%3Aissue>`__ on GitHub.
We appreciate your feedback and suggestions!

For any question, please join the discussion channel on Gitter:

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/radis-radiation/community
    :alt: Gitter


---------------------------------------------------------------------

==========
References
==========

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





Other Spectroscopic tools
-------------------------

Similar packages or softwares you could be interested in (please reference your own if not there!) : 

- `specutil <https://github.com/astropy/specutils>`__: a Python package for spectral analysis in astronomy 
- `pyspeckit <https://github.com/pyspeckit/pyspeckit>`__: a python spectroscopic toolkit 
- `rampy <https://github.com/charlesll/rampy>`__: a Python package for spectral data processing (IR, Raman, XAS...) 
- `scikit-spectra <https://github.com/hugadams/scikit-spectra>`__: Python pandas-based toolkit for explorative spectroscopy, in particular UVVis spectroscopic data. 
- `WrightTools <https://joss.theoj.org/papers/a82637112ac3e03df961d4494bc927d4>`__: a Python package for multidimensional spectroscopy 
- `spectools <https://pyhdust.readthedocs.io/en/latest/spectools.html#module-pyhdust.spectools>`__: Python tools of the BeACoN group
- `SpectroscoPyx <https://github.com/PlasmaPy/SpectroscoPyx>`__: a Python package for spectroscopy
- `Spectragryph <https://www.effemm2.de/spectragryph/index.html>`__: software for FTIR / organic spectroscopy 

And in general the list of `GitHub spectroscopy related packages <https://github.com/topics/spectroscopy>`__




--------

.. |CO2| replace:: CO\ :sub:`2`


.. image:: https://github.com/radis/radis/blob/master/docs/radis_ico.png
    :target: https://radis.readthedocs.io/
    :scale: 50 %
    :alt: RADIS logo
