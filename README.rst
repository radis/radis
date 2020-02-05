 | Users: |badge_pypi|  |badge_pypistats| |badge_article| |badge_docs| |badge_binder| |badge_gitter| 
 | Devs : |badge_contributors| |badge_travis| |badge_coverage| |badge_license| 
    
*****************************************
`RADIS <https://radis.readthedocs.io/>`__
*****************************************

RADIS is a fast line-by-line code for high resolution infrared molecular spectra (emission / absorption, 
equilibrium / nonequilibrium). 

Includes post-processing tools to compare experimental spectra and spectra calculated 
with RADIS or other spectral codes.

User guide, install procedure and examples are available on the `RADIS Website <http://radis.readthedocs.io/>`__:

|badge_docs|


===============
Getting Started
===============

Install
-------

Assuming you have Python installed with the `Anaconda <https://www.anaconda.com/download/>`_ distribution just use::

    pip install radis 
    
**That's it!** You can now run your first example below.
If you encounter any issue, or to upgrade the package later, please refer to the 
`detailed installation procedure <https://radis.readthedocs.io/en/latest/dev/install.html#label-install>`__ . 

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
    
The Quick Start examples automatically download the line databases from `HITRAN-2016 <https://radis.readthedocs.io/en/latest/references/references.html#hitran-2016>`__, which is valid for temperatures below 700 K. 
For *high temperature* cases, you may need to use other line databases such as 
`HITEMP-2010 <https://radis.readthedocs.io/en/latest/references/references.html#hitemp-2010>`__ (typically T < 2000 K) or `CDSD-4000 <https://radis.readthedocs.io/en/latest/references/references.html#cdsd-4000>`__ (T < 5000 K). These databases must be described in a ``~/.radis`` 
`Configuration file <https://radis.readthedocs.io/en/latest/lbl/lbl.html#label-lbl-config-file>`__. 

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
    plot_diff(sexp, s)    # comparing with spectrum 's' calculated previously 

Typical output of `plot_diff <https://radis.readthedocs.io/en/latest/source/radis.spectrum.compare.html#radis.spectrum.compare.plot_diff>`__:

.. image:: docs/spectrum/cdsd4000_vs_hitemp_3409K.svg
    :scale: 60 %
    :target: https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#compare-two-spectra
    :alt: https://radis.readthedocs.io/en/latest/_images/cdsd4000_vs_hitemp_3409K.svg

Refer to the `Examples <https://radis.readthedocs.io/en/latest/examples/examples.html>`__ section for more examples, and to the  
`Spectrum page <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html>`__ for more post-processing functions. 

In the browser (no installation needed!)
----------------------------------------

Alternatively, you can also run RADIS directly in the browser with the  
`RADIS Interactive Examples <https://github.com/radis/radis-examples#interactive-examples>`_ project. 
For instance, run the Quick Start example on the link below:

.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb
    :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb

Or start a bare RADIS online session:

|badge_binder|


Cite
----

Article is available at |badge_article| 

For reproducibility, do not forget to cite the line database used, and the spectroscopic constants 
if running nonquilibrium  calculations. See `How to cite? <https://radis.readthedocs.io/en/latest/references/references.html#cite>`__

---------------------------------------------------------------------

===============
Developer Guide
===============

Architecture
------------

RADIS internals are described in the `Developer Guide <https://radis.readthedocs.io/en/latest/developer.html>`__

.. image:: https://radis.readthedocs.io/en/latest/_images/RADIS_flow_chart.svg
     :target:   https://radis.readthedocs.io/en/latest/dev/architecture.html#label-dev-architecture
     :alt: https://radis.readthedocs.io/en/latest/_images/RADIS_flow_chart.svg


License
-------

The code is available on this repository under 
`GNU LESSER GENERAL PUBLIC LICENSE (v3) <./LICENSE>`_   |badge_license|


Support
-------

If encountering any problem, first refer to the list of known 
`Issues <https://github.com/radis/radis/issues?utf8=%E2%9C%93&q=is%3Aissue>`__ on GitHub.
We appreciate your feedback and suggestions!

For any question, please join the discussion channel on Gitter: |badge_gitter|
or the `Q&A Forum <https://groups.google.com/forum/#!forum/radis-radiation>`__ 



---------------------------------------------------------------------

==========
References
==========

Links
-----

- Documentation: |badge_docs|

- Help: |badge_gitter| `Q&A forum <https://groups.google.com/forum/#!forum/radis-radiation>`__

- Article: |badge_article|

- Source Code: |badge_stars| |badge_contributors| |badge_license|

- Test Status: |badge_travis| |badge_coverage|
 
- PyPi Repository: |badge_pypi|  |badge_pypistats|

- Interactive Examples: `radis_examples <https://github.com/radis/radis-examples>`__ |badge_examples| |badge_binder|




Other Spectroscopic tools
-------------------------

See `awesome-spectra <https://github.com/erwanp/awesome-spectra>`__   |badge_awesome_spectra|

--------

.. image:: https://github.com/radis/radis/blob/master/docs/radis_ico.png
    :target: https://radis.readthedocs.io/
    :scale: 50 %
    :alt: RADIS logo
    
    
    
    

.. |CO2| replace:: CO\ :sub:`2`

.. |badge_docs| image:: https://readthedocs.org/projects/radis/badge/
                :target: https://radis.readthedocs.io/en/latest/?badge=latest
                :alt: Documentation Status

.. |badge_article| image:: https://zenodo.org/badge/doi/10.1016/j.jqsrt.2018.09.027.svg
                   :target: https://linkinghub.elsevier.com/retrieve/pii/S0022407318305867
                   :alt: Article

.. |badge_stars| image:: https://img.shields.io/github/stars/radis/radis.svg?style=social&label=Star
                :target: https://github.com/radis/radis/stargazers
                :alt: GitHub
   
.. |badge_contributors| image:: https://img.shields.io/github/contributors/radis/radis.svg
                        :target: https://github.com/radis/radis/stargazers
                        :alt: Contributors

.. |badge_license| image:: https://img.shields.io/badge/License-LGPL3-blue.svg
                   :target: ./License.md
                   :alt: License

.. |badge_travis| image:: https://img.shields.io/travis/radis/radis.svg
                  :target: https://travis-ci.com/radis/radis
                  :alt: Tests

.. |badge_coverage| image:: https://codecov.io/gh/radis/radis/branch/master/graph/badge.svg
                    :target: https://codecov.io/gh/radis/radis
                    :alt: Coverage

.. |badge_pypi| image:: https://img.shields.io/pypi/v/radis.svg
                :target: https://pypi.python.org/pypi/radis
                :alt: PyPI

.. |badge_pypistats| image:: https://img.shields.io/pypi/dw/radis.svg
                     :target: https://pypistats.org/packages/radis
                     :alt: Downloads

.. |badge_examples| image:: https://img.shields.io/github/stars/radis/radis-examples.svg?style=social&label=Star
                :target: https://github.com/radis/radis-examples
                :alt: Examples

.. |badge_awesome_spectra| image:: https://img.shields.io/github/stars/erwanp/awesome-spectra.svg?style=social&label=Star
                           :target: https://github.com/erwanp/awesome-spectra
                           :alt: Examples

.. |badge_binder| image:: https://mybinder.org/badge.svg 
                  :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb
                  :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb

.. |badge_gitter| image:: https://badges.gitter.im/Join%20Chat.svg
                  :target: https://gitter.im/radis-radiation/community
                  :alt: Gitter
    
