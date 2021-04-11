==========
References
==========

.. include:: constants.rst

.. include:: bibliography.rst

Licence
=======

The code is available for use and modifications on `GitHub <https://github.com/radis/radis>`__
under a `GNU LESSER GENERAL PUBLIC LICENSE (v3) <https://github.com/radis/radis/blob/master/LICENSE>`__,
i.e. modifications must remain public and under LGPLv3.


.. _label_db_cite:


Cite
====

RADIS is built on the shoulders of many state-of-the-art packages and databases. If using RADIS
for your work, **cite all of them that made it possible** :

Line-by-line algorithm :

- cite the line-by-line code as [RADIS-2018]_ |badge_article1|
- if using the default optimization, cite the new spectral synthesis algorithm [Spectral-Synthesis-Algorithm]_ |badge_article2|
- for reproducibility, mention the RADIS version number. Get your version with :py:func:`~radis.get_version`
  (latest version available is |badge_pypi|) ::

    import radis
    radis.get_version()

Database and database retrieval algorithms :

- cite the Line Databases used (for instance, [HITRAN-2016]_, [HITEMP-2010]_ or [CDSD-4000]_ ).
- if downloading [HITRAN-2016]_ directly with ``fetch('hitran')``, cite [HAPI]_ which is the underlying
  interface that makes it work !
- if running nonequilibrium calculations, mention the reference of the spectroscopic constants used
  to calculate the energies (for instance, the
  :ref:`RADIS built-in constants <label_db_spectroscopic_constants>`)


Research Work
=============

Research papers using RADIS and the associated algorithms :

- Papers citing |badge_article1| : https://scholar.google.fr/scholar?cites=5826973099621508256

- Papers citing |badge_article2| : https://scholar.google.fr/scholar?cites=17363432006874800849


Useful Links
============

RADIS:

- Documentation: |badge_docs|

- Help: |badge_slack| |badge_gitter| `Q&A forum <https://groups.google.com/forum/#!forum/radis-radiation>`__

- RADIS Articles: |badge_article1| |badge_article2|

- Source Code: |badge_stars| |badge_contributors| |badge_license|

- Test Status: |badge_travis| |badge_coverage|

- PyPi Repository: |badge_pypi|  |badge_pypistats|

- Interactive Examples: |badge_examples| |badge_binder|

- Try online : :ref:`🌱 RADIS Lab <label_radis_online>`

.. include:: similar_tools.rst


.. |badge_docs| image:: https://readthedocs.org/projects/radis/badge/
                :target: https://radis.readthedocs.io/en/latest/?badge=latest
                :alt: Documentation Status

.. |badge_article1| image:: https://zenodo.org/badge/doi/10.1016/j.jqsrt.2018.09.027.svg
                   :target: https://linkinghub.elsevier.com/retrieve/pii/S0022407318305867
                   :alt: Main Article

.. |badge_article2| image:: https://zenodo.org/badge/doi/10.1016/j.jqsrt.2020.107476.svg
                   :target: https://linkinghub.elsevier.com/retrieve/pii/S0022407320310049
                   :alt: Spectral Synthesis Algorithm

.. |badge_stars| image:: https://img.shields.io/github/stars/radis/radis.svg?style=social&label=GitHub
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

.. |badge_examples| image:: https://img.shields.io/github/stars/radis/radis-examples.svg?style=social&label=Examples
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
