# -*- coding: utf-8 -*-
"""
Test that line survey works
"""

#from radis.misc.utils import DatabankNotFound
from radis.test.utils import getTestFile
from radis.tools.database import load_spec, SpecDatabase
import os
from os.path import exists, dirname, basename

def test_database_functions__fast(verbose=True, plot=True, warnings=True, *args, **kwargs):
    ''' Test SpecDatabase functions '''
    
    import pytest
    
    db = SpecDatabase(dirname(getTestFile('.')))

    # Database visualisation methods    
    if verbose:
        print('{0} items in test database: {1}'.format(len(db), db.see(['Tvib', 'Trot'])))
    if plot:
        db.plot('Tvib', 'Trot')

    # Database get methods
    db.get_closest(Tgas=1300, path_length=1)
    s = db.get_unique(Tgas=1500, path_length=0.01, mole_fraction=0.5)  # there should be one only
                                                      # ... note that this is just so we can test 
                                                      # ... get_unique(). If we were to add new
                                                      # ... test cases with matching conditions
                                                      # ... let's add more criteria to keep it unique
    match = db.get(**s.conditions)
    assert len(match) == 1 
    
    # Database add method
    s2 = s.copy()
    s2.conditions['Tgas'] = 0  # make it unique (for testing)
    l = db.add(s2, if_exists_then='increment')
    assert exists(l)
    try:
        assert s2 in db
        # .. ensures that you cant add it twice
        with pytest.raises(ValueError):
            db.add(s2)
    finally:
        os.remove(l)
    db.update(force_reload=True)            # update database 
    assert s2 not in db
        

def _run_testcases(plot=True, verbose=True, *args, **kwargs):

    # Show media line_shift
    test_database_functions__fast(plot=plot, verbose=verbose, *args, **kwargs)

    return True

if __name__ == '__main__':

    print('Testing line survey functions:', _run_testcases(plot=True))