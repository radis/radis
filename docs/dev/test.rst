====
Test
====

Run tests
---------

Test routines are performed automatically by Travis CI whenever a push is 
made to the ``develop`` or ``master`` GitHub branch. See current test status:

.. image:: https://img.shields.io/travis/radis/radis.svg
    :target: https://travis-ci.org/radis/radis
    :alt: Continuous Integration
    

It is a good practice to perform these tests locally to detect potential 
errors before pushing. To run tests in RADIS use ``pytest``::

    cd radis/test
    pytest
    
By convention in RADIS, test functions that run in less than 5s should have
``fast`` in their name. You can use pytest filtering keys to run 
only these. The whole test procedure should last about one or two minutes::

    cd radis/test 
    pytest -k fast 
    

Code coverage 
-------------

Code coverage makes sure every line in RADIS is properly tested. See 
the current code coverage status:
    
.. image:: https://codecov.io/gh/radis/radis/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/radis/radis
  
 
If you want to see the test coverage report locally use ``codecov`` that 
is interfaced with pytest through the ``--cov=./`` command::

    pip install codecov pytest-cov
    cd radis/test
    pytest --cov=./


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
    