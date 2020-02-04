# -*- coding: utf-8 -*-
"""
Test that :class:`~radis.tools.database.SpecDatabase` works

-------------------------------------------------------------------------------

"""

from __future__ import print_function, absolute_import, division, unicode_literals

# from radis.misc.utils import DatabankNotFound
from radis.test.utils import getTestFile
from radis.tools.database import SpecDatabase, SpecList
import os
from os.path import exists, dirname
import pytest


@pytest.mark.fast
def test_speclist(*args, **kwargs):

    from radis.test.utils import getTestFile

    # Get a list of test files
    spec_list = SpecDatabase(getTestFile(".")).get()

    # Now generate a SpecList
    dbl = SpecList(*spec_list)
    # check methods work:
    dbl.get_closest(Trot=1550)


@pytest.mark.fast
def test_database_functions(
    verbose=True, plot=True, close_plots=True, warnings=True, *args, **kwargs
):
    """ Test SpecDatabase functions """

    if plot:
        import matplotlib.pyplot as plt

        plt.ion()
    if close_plots:
        plt.close("all")

    import pytest

    db = SpecDatabase(dirname(getTestFile(".")))

    # Database visualisation methods
    if verbose:
        print(
            "{0} items in test database: {1}".format(len(db), db.see(["Tvib", "Trot"]))
        )
    if plot:
        db.plot_cond("Tvib", "Trot")

    # Database get methods
    db.get_closest(Tgas=1300, path_length=1)
    s = db.get_unique(
        Tgas=1500, path_length=0.01, mole_fraction=0.5
    )  # there should be one only
    # ... note that this is just so we can test
    # ... get_unique(). If we were to add new
    # ... test cases with matching conditions
    # ... let's add more criteria to keep it unique
    match = db.get(**s.conditions)
    #    assert len(match) == 1
    # TODO: not working in Python 2.7 yet

    # Database add method
    s2 = s.copy()
    s2.conditions["Tgas"] = 0  # make it unique (for testing)

    matching_spectrum_in_db = s2 in db

    l = db.add(
        s2,
        add_info=["Tvib", "Trot"],
        discard=[],
        compress=True,
        if_exists_then="increment",
    )
    assert exists(l)
    try:
        assert s2 in db
        # .. ensures that you cant add it twice
        with pytest.raises(ValueError):
            db.add(s2, add_info=["Tvib", "Trot"])

    # Now remove the Spectrum, update the database and ensures it's not there anymore
    finally:
        os.remove(l)
    db.update(force_reload=True)  # update database
    assert not exists(l)

    # make sure we deleted it properly
    if not matching_spectrum_in_db:
        assert s2 not in db


def test_plot_spec(plot=True, close_plots=True, verbose=True, *args, **kwargs):

    from radis.tools.database import plot_spec
    from radis.test.utils import getTestFile

    if plot:
        import matplotlib.pyplot as plt

        plt.ion()
    if close_plots:
        plt.close("all")

    if plot:
        plot_spec(getTestFile("N2C_specair_380nm.spec"))


def _run_testcases(plot=True, close_plots=False, verbose=True, *args, **kwargs):

    test_speclist(*args, **kwargs)
    test_database_functions(
        plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs
    )
    test_plot_spec(plot=plot, close_plots=close_plots, verbose=verbose, *args, **kwargs)

    return True


if __name__ == "__main__":

    print("Testing database functions:", _run_testcases(plot=True))
