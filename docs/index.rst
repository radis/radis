.. RADIS documentation master file, created by
   sphinx-quickstart on Tue Feb  6 03:32:15 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |logo_png| image:: radis_ico.png

*****
RADIS
*****

RADIS is a fast :ref:`line-by-line code <label_line_by_line>` for high resolution infrared molecular spectra (emission / absorption,
equilibrium / nonequilibrium).

It also includes :ref:`post-processing tools <label_spectrum>` to compare experimental spectra and spectra calculated
with RADIS, or with other spectral codes.

===============
Getting Started
===============

Install
=======

Assuming you have Python installed with the `Anaconda <https://www.anaconda.com/download/>`_ distribution you can use 'pip'::

    pip install radis

or 'conda' via the conda-forge channel::

    conda install radis -c conda-forge

**That's it!** You can now run your first example below.

Had a problem or want something different?

- refer to the :ref:`detailed installation procedure <label_install>`
- try the RADIS-app, our GUI interface :ref:`🌱 RADIS Online <label_radis_online>` !
- try RADIS-lab, our Jupyter notebooks :ref:`🔬 RADIS-lab <label_radis_lab>` !

.. _label_first_example:

Quick Start
===========

Calculate a CO equilibrium spectrum from the [HITRAN-2020]_ database, using the
:py:func:`~radis.lbl.calc.calc_spectrum` function. Lines are downloaded automatically::

    from radis import calc_spectrum
    s = calc_spectrum(1900, 2300,         # cm-1
                      molecule='CO',
                      isotope='1,2,3',
                      pressure=1.01325,   # bar
                      Tgas=700,           # K
                      mole_fraction=0.1,
                      path_length=1,      # cm
                      databank='hitran',  # or 'hitemp', 'geisa', 'exomol'
                      )
    s.apply_slit(0.5, 'nm')       # simulate an experimental slit
    s.plot('radiance')

.. figure:: examples/co_spectrum_700K.png
    :scale: 60 %

=======
Content
=======

:doc:`auto_examples/index`
    Explore the capabilities of RADIS with our examples.
:doc:`features/features`
    More details about RADIS features, what RADIS can and can not provide.
:doc:`lbl/lbl`
    Very detailed section on the Line-by-line (LBL) module and how spectra are generated in RADIS (databanks, line profiles, etc.)
:doc:`spectrum/spectrum`
    Manual around the :class:`~radis.spectrum.spectrum.Spectrum` object, and how to use it to post-process spectra.
:doc:`los/los`
    The Line-of-sight (LOS) module for users dealing with line-of-sight experiments. The module allows combination of Spectra such as::

      s_line_of_sight = (s_plasma_CO2 // s_plasma_CO) > (s_room_absorption)
:doc:`hitran-spectra/hitran-spectra`
    The HITRAN spectra database, plotted for every molecule at ambient conditions.
:doc:`online/online`
    Detailed documentation on the two RADIS online interfaces.
:doc:`dev/developer`
    For developers: how to contribute to RADIS.
:doc:`references/references`
    References, our published research work, license, and recordings of our conferences.
:doc:`modules`
    Detailed documentation of all RADIS modules and classes.

---------------------------------------------------------------------

Cite
====

RADIS is built on the shoulders of many state-of-the-art packages and databases. If using RADIS
to compute spectra, make sure you cite all of them, for proper reproducibility and acknowledgement of
the work ! See :ref:`How to cite? <label_cite>`

---------------------------------------------------------------------

* :ref:`modindex`


.. toctree::
   :maxdepth: 2
   :hidden:

   auto_examples/index
   features/features
   lbl/lbl
   spectrum/spectrum
   los/los
   hitran-spectra/hitran-spectra
   online/online
   dev/developer
   references/references
   modules

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: COUCOU
   user_guide/test1/index
   user_guide/test2/index




---------------------------------------------------------------------

`Q&A Forum <https://groups.google.com/forum/#!forum/radis-radiation>`__

|badge_pypi|  |badge_pypistats| |badge_article| |badge_docs| |badge_license| |badge_contributors| |badge_travis| |badge_coverage| |badge_binder| |badge_gitter| |badge_slack|

|badge_stars|


.. |badge_docs| image:: https://readthedocs.org/projects/radis/badge/
                :target: https://radis.readthedocs.io/en/latest/?badge=latest
                :alt: Documentation Status

.. |badge_article| image:: https://zenodo.org/badge/doi/10.1016/j.jqsrt.2018.09.027.svg
                   :target: https://linkinghub.elsevier.com/retrieve/pii/S0022407318305867
                   :alt: Article

.. |badge_stars| image:: https://img.shields.io/github/stars/radis/radis.svg?style=social&label=GitHub
                :target: https://github.com/radis/radis/stargazers
                :alt: GitHub

.. |badge_contributors| image:: https://img.shields.io/github/contributors/radis/radis.svg
                        :target: https://github.com/radis/radis/graphs/contributors
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

.. |badge_slack| image:: https://img.shields.io/badge/slack-join-green.svg?logo=slack
                  :target: https://radis.github.io/slack-invite/
                  :alt: Slack

