.. _label_line_databases:

Line databases
==============

List of supported line databases formats: :py:data:`~radis.lbl.loader.KNOWN_DBFORMAT`

HITRAN
------

RADIS can automatically fetch HITRAN lines using the `Astroquery <https://astroquery.readthedocs.io>`_
module. This is done by specifying ``databank=='hitran'`` in :py:func:`~radis.lbl.calc.calc_spectrum`
or by using the :py:meth:`~radis.lbl.loader.DatabankLoader.fetch_databank` method in
:class:`~radis.lbl.factory.SpectrumFactory`.
Refer to :py:func:`~radis.io.query.fetch_astroquery` for more information.

You can also download the HITRAN databases files locally:

- HITRAN can be downloaded from https://hitran.org/lbl/. Expect
  ~80 Mb for CO2, or 50 Mb for H2O. Cite with [HITRAN-2020]_.

.. note::

    RADIS has parsers to read line databases in Pandas dataframes.
    This can be useful if you want to edit the database.
    see :func:`~radis.io.hitran.hit2df`

    There are also functions to get HITRAN molecule ids, and vice-versa:
    :func:`~radis.db.classes.get_molecule`, :func:`~radis.db.classes.get_molecule_identifier`


HITEMP
------

RADIS can read files from the HITEMP database.

- HITEMP-2010 files can be downloaded from https://hitran.org/hitemp/. Expect
  ~3 Gb for CO2 or ~10 Gb for H2O. Cite with [HITEMP-2010]_

The ``~/radis.json`` is then used to properly handle the line databases
on the User environment. See the :ref:`Configuration file <label_lbl_config_file>` section, as well as
the :py:mod:`radis.misc.config` module and the :py:func:`~radis.misc.config.getDatabankList`
function for more information.

📣 starting from radis==0.9.28 you can also download HITEMP directly. Example ::

    from radis import calc_spectrum
    calc_spectrum(
        wavenum_min=2500 / u.cm,
        wavenum_max=4500 / u.cm,
        molecule="OH",
        Tgas=600,
        databank="hitemp",  # test by fetching directly
        verbose=False,
    )

- Some HITEMP line databases are pre-configured in the :ref:`RADIS-lab <label_radis_lab>` online environment.
  No install needed !

- If you just want to parse the HITEMP files, use :py:func:`~radis.io.hitemp.fetch_hitemp` ::

    from radis.io.hitemp import fetch_hitemp
    fetch_hitemp("NO")


CDSD-4000
---------

RADIS can read files from the CDSD-4000 database, however files have to be
downloaded manually.

- CDSD-4000 files can be downloaded from ftp://ftp.iao.ru/pub/. Expect ~50 Gb for all CO2.
  Cite with [CDSD-4000]_.
- Tabulated partition functions are availabe in the ``partition_functions.txt`` file on the
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
