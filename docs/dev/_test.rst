.. _label_dev_test:

Test
====


Test status
-----------

Test routines are performed automatically by `Travis CI <https://travis-ci.com/radis/radis>`_
whenever a change is made on the `GitHub <https://github.com/radis/radis>`_ repository.
See the current test status below:

.. image:: https://img.shields.io/travis/radis/radis.svg
    :target: https://travis-ci.com/radis/radis
    :alt: https://travis-ci.com/radis/radis

It is a good practice to perform these tests locally to detect potential
errors before pushing.
To run the tests locally, assuming you cloned the source code
(see the :ref:`Install section <label_install>` ), run the following command in
the ``radis`` directory::

    cd radis
    pytest

The whole test procedure takes 5 - 10 min. You can use pytest filtering keys
to :ref:`run specific tests only <label_dev_select_test>`


Code coverage
-------------

Code coverage makes sure every line in RADIS is properly tested. See
the current code coverage status (click the badge for more details):

.. image:: https://codecov.io/gh/radis/radis/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/radis/radis
  :alt: https://codecov.io/gh/radis/radis



If you want to see the test coverage report locally use ``codecov`` that
is interfaced with pytest through the ``--cov=./`` command::

    pip install codecov pytest-cov
    cd radis/test
    pytest --cov=./

Performance benchmarks
----------------------

RADIS performance is tested against past versions on a dedicated project : `radis-benchmark <https://github.com/radis/radis-benchmark>`__.

Results can be found on : 🔗 https://radis.github.io/radis-benchmark/

.. image:: http://img.shields.io/badge/benchmarked%20by-asv-blue.svg?style=flat
            :target: https://github.com/radis/radis-benchmark
            :alt: Benchmarks


.. _label_dev_select_test:

Select tests
------------

You can select or ignore some tests with the ``-m`` option, for instance
run only the fast tests with::

    pytest -m fast --cov=./

The list of tags is:

- fast : each test should run under 1s
- needs_connection : requires an internet connection
- needs_python3 : a Python-3 only test
- needs_config_file: needs ``~/radis.json`` to be defined.

Plus some tags for user-defined [HITEMP-2010]_ databases, as given in the :ref:`Configuration file <label_lbl_config_file>`
section

- needs_db_HITEMP_CO2_DUNHAM : requires HITEMP-CO2-DUNHAM database in ``~/radis.json``
- needs_db_HITEMP_CO_DUNHAM : requires HITEMP-CO-DUNHAM database in ``~/radis.json``
- needs_db_CDSD_HITEMP_PC : requires CDSD-HITEMP-PC database in ``~/radis.json``

The default test routine run on `Travis CI <https://travis-ci.com/radis/radis>`__
is (see the ``radis/.gitlab-ci.yml`` file)::

    pytest -m "not needs_config_file" --cov=./;

Which ignores all test that rely on the [HITEMP-2010]_ databases, as they cannot (yet) be downloaded
automatically on the test machine.

Write new tests
---------------

Any function starting with ``test`` will be picked by pytest. Use assert
statements within the function to test properties of your objects. Ex::

    def test_positive(a):
        assert a>0

Tips: make sure you don't open any figure by default in your test routines,
else it will be stuck when called by pytest. Or, force a non-blocking behaviour
adding the following lines within your test function::

    if plot:
        import matplotlib.pyplot as plt
        plt.ion()   # dont get stuck with Matplotlib if executing through pytest

You can refer to Pull Request: https://github.com/radis/radis/pull/495 to see
how test cases are written when a new feature is added.

See: https://github.com/statsmodels/statsmodels/issues/3697



.. _label_dev_test_files:

Test files
----------

To make the errors as reproducible as possible, try to use the test files provided in the
Test files are provided in the `radis\test\files <https://github.com/radis/radis/tree/develop/radis/test/files>`__
and `radis\test\validation <https://github.com/radis/radis/tree/develop/radis/test/validation>`__ folders.
They contain examples of line databases, spectra, or energy levels.
The path to these test files can be retrieved using the :py:func:`~radis.test.utils.getTestFile` and
:py:func:`~radis.test.utils.getValidationCase` functions, respectively.

Load a line database file ::

    from radis.test.utils import getTestFile
    from radis.api.hitranapi import hit2df
    df = hit2df(getTestFile("hitran_CO_fragment.par"))

    print(df)  # replace with your test code

    >>> Out:

           id  iso       wav           int             A  ...  gpp  branch  jl  vu  vl
    0   5    1  3.705026  2.354000e-44  2.868000e-10  ...  1.0       1   0   4   4
    1   5    1  3.740024  1.110000e-38  5.999000e-09  ...  1.0       1   0   3   3
    2   5    1  3.775024  9.233000e-34  1.947000e-08  ...  1.0       1   0   2   2
    3   5    1  3.810028  5.706000e-29  4.130000e-08  ...  1.0       1   0   1   1
    4   5    1  3.845033  3.300000e-24  7.207000e-08  ...  1.0       1   0   0   0
    5   5    1  7.409906  1.815000e-43  2.726000e-09  ...  3.0       1   1   4   4
    6   5    1  7.479900  8.621000e-38  5.746000e-08  ...  3.0       1   1   3   3
    7   5    1  7.549901  7.177000e-33  1.867000e-07  ...  3.0       1   1   2   2
    8   5    1  7.619908  4.436000e-28  3.961000e-07  ...  3.0       1   1   1   1

    [9 rows x 16 columns]

Load a Spectrum object ::

    from radis.test.utils import getTestFile
    from radis import load_spec
    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.5.spec"))

    print(s)    # replace with your test code

    >>> Out:

        Spectrum Name:  CO_Tgas1500K_mole_fraction0.5.spec
    Spectral Arrays
    ----------------------------------------
       abscoeff 	(37,870 points)
    Populations Stored
    ----------------------------------------
       CO 		 [1]
    Physical Conditions
    ----------------------------------------
       Tgas                 1500 K
       Trot                 1500 K
       Tvib                 1500 K
       isotope              1
       mole_fraction        0.5
       molecule             CO
       path_length          0.01 cm
       pressure_mbar        1013.25 mbar
       rot_distribution     boltzmann
       self_absorption      True
       state                X
       thermal_equilibrium  True
       vib_distribution     boltzmann
       wavelength_max       4801.3089 nm
       wavelength_min       4401.1999 nm
       wavenum_max          2272.1077 cm-1
       wavenum_min          2082.7654 cm-1
    Computation Parameters
    ----------------------------------------
       Tref                 296 K
       broadening_max_width  10 cm-1
       cutoff               1e-25 cm-1/(#.cm-2)
       db_assumed_sorted    True
       dbformat             hitran
       dbpath               d:/github/radis/radis/test/files/hitran_co_3iso_2000_2300cm.par
       levelsfmt            neq
       parfuncfmt           hapi
       pseudo_continuum_threshold  0
       wavenum_max_calc     2277.1104 cm-1
       wavenum_min_calc     2077.7654 cm-1
       waveunit             cm-1
       wstep                0.005 cm-1
    Information
    ----------------------------------------
       calculation_time     0.14 s
       chunksize            10000000.0
       db_use_cached        True
    ----------------------------------------


Report errors
-------------

If you encounter any error, open an `Issue on GitHub <https://github.com/radis/radis/issues>`__

To simplify the debugging process, provide a code snippet that reproduces
the error. If you need a line database, spectrum, or energy level, try to use one
of the :ref:`test files <label_dev_test_files>`.

Debugging
---------

See the :py:func:`~radis.misc.debug.printdbg` function in ``radis.misc``, and
the :py:data:`~radis.DEBUG_MODE` global variable.
