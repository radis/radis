
Line databases
--------------

List of supported line databases formats: :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`

HITRAN
''''''

RADIS can automatically fetch HITRAN lines using the `Astroquery <https://astroquery.readthedocs.io>`_ 
module. 

This is done by specifying ``databank=='fetch'`` in :py:func:`~radis.lbl.calc.calc_spectrum`
or by using the :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank` method in 
:class:`~radis.lbl.factory.SpectrumFactory`. 

Refer to :py:func:`~radis.io.query.fetch_astroquery` for more information on 
the wrapper to Astroquery. 

RADIS has parsers to read line databases in Pandas dataframes: 
see :func:`~radis.io.hitran.hit2df`

There are also functions to get HITRAN molecule ids, and vice-versa:
:func:`~radis.io.hitran.get_molecule`, :func:`~radis.io.hitran.get_molecule_identifier`


HITEMP / CDSD-4000
''''''''''''''''''

RADIS can read files from the HITEMP-2010 and CDSD-4000 databases 
(see `parsers <https://radis.readthedocs.io/en/latest/io/parsers.html>`__), 
however the user needs to download them manually.

The ``~/.radis`` is then used to properly handle the line databases 
on the User environment. 

See the :ref:`Configuration file <label_lbl_config_file>` section, as well as 
the :py:mod:`radis.misc.config` module and the :py:func:`~radis.misc.config.getDatabankList` 
function for more information. 

See :func:`~radis.io.cdsd.cdsd2df` for the conversion to a Pandas DataFrame. 