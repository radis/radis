========
Features
========

.. include:: description.rst


Features
========


RADIS is both an infrared :ref:`line-by-line code <label_line_by_line>`
and a :ref:`post-processing library <label_spectrum>`.
It includes:

- Absorption and emission spectra of all [HITRAN-2020]_ and [ExoMol-2020]_ species under equilibrium calculations (:py:data:`~radis.db.MOLECULES_LIST_EQUILIBRIUM`)
- Absorption and emission spectra of CO2 and CO for non-LTE calculations (see :py:data:`~radis.db.MOLECULES_LIST_NONEQUILIBRIUM` )
- Different Line Databases: support of [HITRAN-2020]_, [HITEMP-2010]_ and [CDSD-4000]_ line databases (see :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`)
- Calculation of :ref:`Rovibrational Energies of molecules <label_examples_rovibrational_energies>`.
- Calculation of equilibrium and nonequilibrium :ref:`Partition Functions <label_examples_partition_functions>`.
- Spatially heterogeneous spectra (see :ref:`see line-of-sight <label_los_index>`)
- Post-processing tools to load and :ref:`compare with experimental spectra <label_spectrum_howto_compare>` (see :ref:`the Spectrum object <label_spectrum>`)
- A :ref:`Line Survey <label_spectrum_linesurvey>` tool to identify which lines correspond to a spectral feature.

RADIS does *not* include, so far:

- Line-mixing effects and speed-dependant lineshapes. [HAPI]_ is a Python alternative that does it.
- Collisional-induced absorption (CIA) or emission.
- Electronic states other than electronic ground states
- Hamiltonian calculations (a private module for CO2 is available `on request <mailto:erwan.pannier@gmail.com>`__)
- Raman spectra (contribute in `#43 <https://github.com/radis/radis/issues/43>`__)

RADIS also features:

- :ref:`High Performances <label_lbl_performance>`: spectra are calculated up to several orders of magnitude faster than equivalent line-by-line codes.
- In-the-browser calculations (no install needed) : see :ref:`🌱 RADIS Online <label_radis_online>`.
- Automatic download of the latest HITRAN and HITEMP databases with :py:func:`~radis.lbl.calc.calc_spectrum`
- Automatic testing and continuous integration tools for a reliable :ref:`Open-source Development <label_developer_guide>`.

Remarks and request for features can be done on `GitHub <https://github.com/radis/radis/issues>`__ ,
on the `Q&A forum <https://groups.google.com/forum/#!forum/radis-radiation>`__ or on the Gitter community chat:

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/radis-radiation/community
    :alt: Gitter

.. image:: https://img.shields.io/badge/slack-join-green.svg?logo=slack
    :target: https://radis.github.io/slack-invite/
    :alt: Slack


Use Cases
---------

Use RADIS to:

- Quickly compare different line databases:
  Various line database formats are supported by RADIS, and can be easily switched
  to be compared. See the list of supported line databases formats:
  :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`
  and refer to the :ref:`Configuration file <label_lbl_config_file>` on how to use them.

  See the comparison of two CO2 spectra calculated with [HITEMP-2010]_ and [CDSD-4000]_
  below:

  .. image:: https://radis.readthedocs.io/en/latest/_images/cdsd4000_vs_hitemp_3409K.svg
      :alt: https://radis.readthedocs.io/en/latest/_images/cdsd4000_vs_hitemp_3409K.svg
      :scale: 100 %


- Use the RADIS post-processing methods with the calculation results of another spectral code. For instance,
  `pySpecair <https://spectralfit.gitlab.io/specair/>`__, the Python interface to `SPECAIR <http://www.specair-radiation.net/>`__,
  uses the RADIS :py:class:`~radis.spectrum.spectrum.Spectrum` object for post-processing
  (see :ref:`How to generate a Spectrum? <label_howto_generate_spectrum>`)


Refer to the :ref:`Examples <label_examples>` section for more examples, or to
the `RADIS Interactive Examples <https://github.com/radis/radis-examples#interactive-examples>`_ project.




See the :ref:`Architecture <label_dev_architecture>` section for an overview of the RADIS calculation
steps.


Line Databases
==============

List of supported line databases formats: :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT` :

- [HITRAN-2016]_
- [HITRAN-2020]_
- [HITEMP-2010]_
- [CDSD-4000]_
- [ExoMol-2020]_

For download and configuration of line databases, see the :ref:`Line Databases section <label_line_databases>`



Interfaces
==========

RADIS includes parsers and interfaces to read and return data in different formats:

.. include:: thermo.rst


New features
============

RADIS is open-source, so everyone can `contribute <https://github.com/radis/radis/contribute>`__
to the code development. Read the :ref:`Developer Guide <label_developer_guide>` to get started.

You can also suggest or vote for new features below:

.. image:: https://feathub.com/radis/radis?format=svg
   :target: https://feathub.com/radis/radis
