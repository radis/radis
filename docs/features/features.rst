========
Features
========

.. include:: description.rst


Features
========


RADIS is both an infrared :ref:`line-by-line code <label_line_by_line>` 
and a :ref:`post-processing library <label_spectrum>`. 
It includes: 

- all [HITRAN-2016]_ species for equilibrium calculations (see :py:data:`~radis.io.MOLECULES_LIST_EQUILIBRIUM`)
- CO2 and CO for nonequilibrium calculations (see :py:data:`~radis.io.MOLECULES_LIST_NONEQUILIBRIUM` )
- different line databases: support of [HITRAN-2016]_, [HITEMP-2010]_ and [CDSD-4000]_ line databases (see :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`) 
- fast calculations: up to several orders of magnitude faster than equivalent line-by-line codes (see :ref:`Performance <label_lbl_performance>`)
- spatially heterogeneous spectra (see :ref:`see line-of-sight <label_los_index>`)
- post processing tools to load and :ref:`compare with experimental spectra <label_spectrum_howto_compare>` (see :ref:`the Spectrum object <label_spectrum>`)
- a line survey tool to identify which lines correspond to a spectral feature (see :ref:`Line Survey <label_spectrum_linesurvey>`).
- automatic testing and continuous integration tools for a reliable open-source development (see the :ref:`Developer Guide <label_developer_guide>`)
- in-the-browser calculations (no install needed: see the `RADIS Interactive Examples <https://github.com/radis/radis-examples#interactive-examples>`_ project)

RADIS does not include, so far: 

- line-mixing effects and speed-dependant lineshapes. [HAPI]_ is a Python alternative that does it. 
- collisional-induced absorption (CIA) or emission. 
- electronic states other than electronic ground states
- Hamiltonian calculations (a private module for CO2 is available `on request <mailto:erwan.pannier@gmail.com>`__)
- Raman spectra (contribute in `#43 <https://github.com/radis/radis/issues/43>`__)

Remarks and request for features can be done on `GitHub <https://github.com/radis/radis/issues>`__ ,
on the `Q&A forum <https://groups.google.com/forum/#!forum/radis-radiation>`__ or on the Gitter community chat: 

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/radis-radiation/community
    :alt: Gitter


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

  .. image:: spectrum/cdsd4000_vs_hitemp_3409K.*
      :alt: https://radis.readthedocs.io/en/latest/_images/cdsd4000_vs_hitemp_3409K.svg


- Use the RADIS post-processing methods with the calculation results of another spectral code. For instance, 
  `pySpecair <https://spectralfit.gitlab.io/specair/>`__, the Python interface to `SPECAIR <http://www.specair-radiation.net/>`__,
  uses the RADIS :py:class:`~radis.spectrum.spectrum.Spectrum` object for post-processing 
  (see :ref:`How to generate a Spectrum? <label_howto_generate_spectrum>`)
  
  
Refer to the :ref:`Examples <label_examples>` section for more examples, or to 
the `RADIS Interactive Examples <https://github.com/radis/radis-examples#interactive-examples>`_ project. 




See the :ref:`Architecture <label_dev_architecture>` section for an overview of the RADIS calculation 
steps. 


Interfaces
==========

RADIS includes parsers and interfaces to read and return data in different formats: 

.. include:: databases.rst

.. include:: thermo.rst
