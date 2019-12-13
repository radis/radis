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
- needs_config_file: needs ``~/.radis`` to be defined.

Plus some tags for user-defined [HITEMP-2010]_ databases, as given in the :ref:`Configuration file <label_lbl_config_file>`
section

- needs_db_HITEMP_CO2_DUNHAM : requires HITEMP-CO2-DUNHAM database in ``~/.radis``
- needs_db_HITEMP_CO_DUNHAM : requires HITEMP-CO-DUNHAM database in ``~/.radis`` 
- needs_db_CDSD_HITEMP_PC : requires CDSD-HITEMP-PC database in ``~/.radis``

The default test routine run on `Travis CI <https://travis-ci.com/radis/radis>`__
is (see the ``radis/.gitlab-ci.yml` file)::

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
        
See: https://github.com/statsmodels/statsmodels/issues/3697




Report errors
-------------

If you encounter any error, open an `Issue on GitHub <https://github.com/radis/radis/issues>`__



Debugging
---------

See the :py:func:`~radis.misc.debug.printdbg` function in ``radis.misc``, and
the :py:data:`~radis.DEBUG_MODE` global variable. 
    