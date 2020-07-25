.. _label_line_databases:

Line databases
==============

List of supported line databases formats: :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`

HITRAN
------

RADIS can automatically fetch HITRAN lines using the `Astroquery <https://astroquery.readthedocs.io>`_ 
module. This is done by specifying ``databank=='fetch'`` in :py:func:`~radis.lbl.calc.calc_spectrum`
or by using the :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank` method in 
:class:`~radis.lbl.factory.SpectrumFactory`. 
Refer to :py:func:`~radis.io.query.fetch_astroquery` for more information.

You can also download the HITRAN databases files locally: 

- HITRAN can be downloaded from https://hitran.org/lbl/. Expect
  ~80 Mb for CO2, or 50 Mb for H2O. Cite with [HITRAN-2016]_. 

.. note::

    RADIS has parsers to read line databases in Pandas dataframes. 
    This can be useful if you want to edit the database. 
    see :func:`~radis.io.hitran.hit2df`

    There are also functions to get HITRAN molecule ids, and vice-versa:
    :func:`~radis.io.hitran.get_molecule`, :func:`~radis.io.hitran.get_molecule_identifier`


HITEMP
------

RADIS can read files from the HITEMP 2010 database, however files have to be 
downloaded manually.

- HITEMP-2010 files can be downloaded from https://hitran.org/hitemp/. Expect
  ~3 Gb for CO2 or ~10 Gb for H2O. Cite with [HITEMP-2010]_ 

The ``~/.radis`` is then used to properly handle the line databases 
on the User environment. See the :ref:`Configuration file <label_lbl_config_file>` section, as well as 
the :py:mod:`radis.misc.config` module and the :py:func:`~radis.misc.config.getDatabankList` 
function for more information. 

CDSD-4000
---------

RADIS can read files from the CDSD-4000 database, however files have to be 
downloaded manually.

- CDSD-4000 files can be downloaded from ftp://ftp.iao.ru/pub/. Expect ~50 Gb for all CO2. 
  Cite with [CDSD-4000]_. 

The ``~/.radis`` is then used to properly handle the line databases 
on the User environment. See the :ref:`Configuration file <label_lbl_config_file>` section, as well as 
the :py:mod:`radis.misc.config` module and the :py:func:`~radis.misc.config.getDatabankList` 
function for more information. 

.. note::

    See :func:`~radis.io.cdsd.cdsd2df` for the conversion to a Pandas DataFrame. 
