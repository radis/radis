========
General description
========

.. include:: description.rst


Features
========


RADIS is both an infrared :ref:`line-by-line code <label_line_by_line>`
and a :ref:`post-processing library <label_spectrum>`.
It includes:

- Absorption and emission spectra of all [HITRAN-2020]_ and [ExoMol-2020]_ species under equilibrium calculations (:py:data:`~radis.db.MOLECULES_LIST_EQUILIBRIUM`)
- Absorption and emission spectra of CO2 and CO for non-LTE calculations (see :py:data:`~radis.db.MOLECULES_LIST_NONEQUILIBRIUM` )
- A common API to Different Line Databases: support of [HITRAN-2020]_, [HITEMP-2010]_, [CDSD-4000]_, [ExoMol-2020]_, [GEISA-2020]_ line databases (see :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`). The API is also compatible with the :py:mod:`exojax` code.
- Calculation of rovibrational energies of molecules, see :ref:`database handling <label_example_database_handling>`.
- Calculation of equilibrium and nonequilibrium partition functions, see :ref:`database handling <label_example_database_handling>`.
- Spatially heterogeneous spectra, see :ref:`see line-of-sight <label_los_index>`
- Post-processing tools to load and :ref:`compare with experimental spectra <label_spectrum_howto_compare>`
- A :ref:`Line Survey <label_spectrum_linesurvey>` tool to identify which lines correspond to a spectral feature.
- A :ref:`Matlab implementation <label_matlab_access>`)

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

New features
============

RADIS is open-source, so everyone can contribute
to the code development or suggest new features in our `GitHub page <https://github.com/radis/radis>`__.
Read the :ref:`Developer Guide <label_developer_guide>` to get started.

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/radis-radiation/community
    :alt: Gitter

.. image:: https://img.shields.io/badge/slack-join-green.svg?logo=slack
    :target: https://radis.github.io/slack-invite/
    :alt: Slack
