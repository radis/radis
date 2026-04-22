.. _label_dev_test:

Test Status
-----------

Test routines are automatically executed by `Travis CI <https://app.travis-ci.com/github/radis/radis/branches>`_
whenever changes are made to the `GitHub <https://github.com/radis/radis>`_ repository.

It is recommended to run these tests locally to identify potential errors before pushing changes.
To run the tests locally, assuming you have cloned the source code
(see the :ref:`Install section <label_install>`), execute the following command in
the ``radis`` directory::

    cd radis
    pytest

The entire test procedure takes approximately 5–10 minutes. You can use pytest filtering keys
to :ref:`run specific tests only <label_dev_select_test>`.

.. _label_dev_select_test:

Selecting Tests
---------------

You can select or exclude specific tests using the ``-m`` option. For example, to run only the fast tests, use::

    pytest -m fast

Commonly used tags include:

- **fast**: Tests that run in under 1 second.
- **needs_connection**: Requires an internet connection.
- **needs_config_file**: Requires the ``~/radis.json`` configuration file.

Additional tags for user-defined [HITEMP-2010]_ databases, as described in the :ref:`Configuration file <label_lbl_config_file>` section:

- **needs_db_HITEMP_CO2_DUNHAM**: Requires the HITEMP-CO2-DUNHAM database in ``~/radis.json``.
- **needs_db_HITEMP_CO_DUNHAM**: Requires the HITEMP-CO-DUNHAM database in ``~/radis.json``.
- **needs_db_CDSD_HITEMP_PC**: Requires the CDSD-HITEMP-PC database in ``~/radis.json``.

The default test routine executed on GitHub actions (see the ``radis/.travis.yml`` file) are::

    pytest -m "fast and not needs_db_CDSD_HITEMP_PCN and not needs_db_CDSD_HITEMP and not needs_db_CDSD_HITEMP_PC"
    pytest -m "not fast and not needs_cuda and not download_large_databases and not needs_db_CDSD_HITEMP and not needs_db_CDSD_HITEMP_PCN and not needs_db_CDSD_HITEMP_PC and not needs_db_HITEMP_CO2_DUNHAM and not needs_db_HITEMP_CO_DUNHAM"


Writing New Tests
-----------------

Any function starting with ``test`` will be automatically detected by pytest. Use `assert` statements within the function to validate object properties. For example::

    def test_positive(a):
        assert a > 0

Tip: Ensure that no figures are opened by default in your test routines, as this may cause the process to hang when executed by pytest. To enforce non-blocking behavior, add the following lines within your test function::

    if plot:
        import matplotlib.pyplot as plt
        plt.ion()  # Prevent Matplotlib from blocking execution during pytest

Refer to previous pull requests for examples of test cases written for new features (e.g., `#495 <https://github.com/radis/radis/pull/495>`_ or `#3697 <https://github.com/statsmodels/statsmodels/issues/3697>`_).

.. _label_dev_test_files:

Test Files
----------

To ensure errors are reproducible, use the test files provided in the
`radis/test/files <https://github.com/radis/radis/tree/develop/radis/test/files>`__
and `radis/test/validation <https://github.com/radis/radis/tree/develop/radis/test/validation>`__ directories.
These files include examples of line databases, spectra, and energy levels.
You can retrieve the paths to these test files using the :py:func:`~radis.test.utils.getTestFile` and
:py:func:`~radis.test.utils.getValidationCase` functions, respectively.

Load a line database file::

    from radis.test.utils import getTestFile
    from radis.api.hitranapi import hit2df
    df = hit2df(getTestFile("hitran_CO_fragment.par"))

    print(df)  # Replace with your test code

    >>> Out:
           id  iso       wav           int             A  ...  gpp  branch  jl  vu  vl
    0   5    1  3.705026  2.354000e-44  2.868000e-10  ...  1.0       1   0   4   4
    ...

Load a Spectrum object::

    from radis.test.utils import getTestFile
    from radis import load_spec
    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.5.spec"))

    print(s)  # Replace with your test code

    >>> Out:
        Spectrum Name:  CO_Tgas1500K_mole_fraction0.5.spec
    Spectral Arrays
    ----------------------------------------
       abscoeff 	(37,870 points)
    ...

Debugging
---------

Use the :py:func:`~radis.misc.debug.printdbg` function in ``radis.misc`` and
the :py:data:`~radis.DEBUG_MODE` global variable for debugging.

Code Coverage
-------------

Code coverage ensures that every line in RADIS is properly tested. View the current code coverage status (click the badge for details):

.. image:: https://codecov.io/gh/radis/radis/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/radis/radis
  :alt: Code Coverage

To view the test coverage report locally, use ``codecov``, which is integrated with pytest through the ``--cov=./`` command::

    pip install codecov pytest-cov
    cd radis/test
    pytest --cov=./

Performance Benchmarks
----------------------

RADIS performance is benchmarked against previous versions in a dedicated project: `radis-benchmark <https://github.com/radis/radis-benchmark>`__.

Results are available at: 🔗 https://radis.github.io/radis-benchmark/

.. image:: http://img.shields.io/badge/benchmarked%20by-asv-blue.svg?style=flat
  :target: https://github.com/radis/radis-benchmark
  :alt: Benchmarks
