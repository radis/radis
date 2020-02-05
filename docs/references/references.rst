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

Cite
====

If using RADIS for your work, cite: 

- the line-by-line code as [RADIS-2018]_ |badge_article|  
- the :ref:`version number <https://radis.readthedocs.io/en/latest/source/radis.html#radis.get_version>`__ 
  (latest version available is |badge_pypi|) ::
  
    import radis
    radis.get_version()

- the Line Databases used (for instance, [HITRAN-2016]_, [HITEMP-2010]_ or [CDSD-4000]_ ).
- if running nonequilibrium calculations, mention the reference of the spectroscopic constants used 
  to calculate the energies (for instance, the 
  :ref:`RADIS built-in constants <label_db_spectroscopic_constants>`)


Useful Links
------------

Links
-----

RADIS:

- Documentation: |badge_docs|

- Help: |badge_gitter| `Q&A forum <https://groups.google.com/forum/#!forum/radis-radiation>`__

- Article: |badge_article|

- Source Code: |badge_stars| |badge_contributors| |badge_license|

- Test Status: |badge_travis| |badge_coverage|
 
- PyPi Repository: |badge_pypi|  |badge_pypistats|

- Interactive Examples: `radis_examples <https://github.com/radis/radis-examples>`__ |badge_examples| |badge_binder|


Other Spectroscopic Tools:

.. include:: similar_tools.rst 



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
    