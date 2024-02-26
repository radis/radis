.. _label_line_databases:

Line databases
==============

List of supported line databases formats: :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`.
The databases are stored in a hidden folder ``~/.radisdb`` on the user's computer.
The informations regarding each database are stored in the configuration file ``~/radis.json``, see the :ref:`Configuration file <label_lbl_config_file>` section.

HITRAN / HITEMP
------

HITRAN is designed for atmospheric applications. It is the most widely used line database.
The HITRAN group also supports the HITEMP database, which is designed for high-temperature applications.
The files of the HITEMP database are formatted in the same way as HITRAN files but are much larger.

By default, RADIS automatically fetch HITRAN/HITEMP lines using the `Astroquery <https://astroquery.readthedocs.io>`_
module. This is done by specifying ``databank=='hitran'`` or ``databank=='hitemp'``.
- See the following example for :py:func:`~radis.lbl.calc.calc_spectrum` ::

    from radis import calc_spectrum
    calc_spectrum(
        wavenum_min=2500 / u.cm,
        wavenum_max=4500 / u.cm,
        molecule="OH",
        Tgas=600,
        databank="hitemp",  # test by fetching directly
        verbose=False,
    )

- If you just want to parse the HITEMP files, use :py:func:`~radis.io.hitemp.fetch_hitemp`. You can also refer to :py:func:`~radis.io.query.fetch_astroquery` for more information. ::

    from radis.io.hitemp import fetch_hitemp
    fetch_hitemp("NO")



You can also download the HITRAN/HITEMP databases files locally, but you will have to properly refer to your database in ``~/radis.json``.

.. note::

    RADIS has parsers to read line databases in Pandas dataframes.
    This can be useful if you want to edit the database.
    see :func:`~radis.io.hitran.hit2df`

    There are also functions to get HITRAN molecule ids, and vice-versa:
    :func:`~radis.db.classes.get_molecule`, :func:`~radis.db.classes.get_molecule_identifier`

GEISA
---------

The GEISA database is designed for atmospheric applications. Using GEISA in RADIS is similar to the HITRAN database, simply specify ``databank=='geisa'``.

CDSD-4000
---------

RADIS can read files from the CDSD-4000 database, however files have to be
downloaded manually.

- CDSD-4000 files can be downloaded from ftp://ftp.iao.ru/pub/. Expect ~50 Gb for all CO2.
  Cite with [CDSD-4000]_.
- Tabulated partition functions are available in the ``partition_functions.txt`` file on the
  [CDSD-4000]_ FTP : ftp://ftp.iao.ru/pub/CDSD-4000/  . They can be loaded and interpolated
  with :py:class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDtab`. This can be done automatically
  providing ``parfuncfmt: cdsd`` and ``parfunc = PATH/TO/cdsd_partition_functions.txt`` is given
  in the ``~/radis.json`` configuration file (see the :ref:`Configuration file <label_lbl_config_file>`).

The ``~/radis.json`` is  used to properly handle the line databases
on the User environment. See the :ref:`Configuration file <label_lbl_config_file>` section, as well as
the :py:mod:`radis.misc.config` module and the :py:func:`~radis.misc.config.getDatabankList`
function for more information.

.. note::

    See :func:`~radis.io.cdsd.cdsd2df` for the conversion to a Pandas DataFrame.
