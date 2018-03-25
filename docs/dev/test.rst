====
Test
====

Test routines are performed automatically whenever a push is made to the 
develop or master GitHub branch. It is a good practice to perform these 
tests locally to detect potential errors before pushing. 

To run tests in RADIS use ``pytest``::

    cd radis/test
    pytest
    
By convention, test functions that run in less than 5s should have
``fast`` in their name. You can use pytest filtering keys to run 
only these::

    cd radis/test 
    pytest -k fast 
    
If you want to see the test coverage use ``codecov`` that is interfaced
with pytest with the ``--cov=./`` command::

    pip install codecov pytest-cov
    cd radis/test
    pytest --cov=./