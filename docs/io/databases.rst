
**************
Line databases
**************

HITRAN
------

RADIS can automatically fetch HITRAN lines using the `Astroquery <https://astroquery.readthedocs.io>`_ 
module. 

This is done by specifying ``databank=='fetch'`` in :py:func:`~radis.lbl.calc.calc_spectrum`
or by using the :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank` method in 
:class:`~radis.lbl.factory.SpectrumFactory`. 

Refer to :py:func:`~radis.io.query.fetch_astroquery` for more information on 
the wrapper to Astroquery. 


HITEMP / CDSD
-------------

RADIS can read files from the HITEMP-2010 and CDSD-4000 databases 
(see `parsers <https://radis.readthedocs.io/en/latest/io/parsers.html>`__), 
however the user needs to download them manually.

The ``~/.radis`` is then used to properly handle the line databases 
on the User environment. 

See the :py:mod:`radis.misc.config` module and the 
:py:func:`~radis.misc.config.getDatabankList` function for more information. 
