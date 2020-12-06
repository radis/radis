# -*- coding: utf-8 -*-
"""
Test that :class:`~radis.tools.database.SpecDatabase` works

-------------------------------------------------------------------------------

"""


import os
from os.path import dirname, exists

import pytest

# from radis.misc.utils import DatabankNotFound
from radis.test.utils import getTestFile
from radis.tools.database import SpecDatabase, SpecList


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
    assert len(match) == 1

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

    from radis.test.utils import getTestFile
    from radis.tools.database import plot_spec

    if plot:
        import matplotlib.pyplot as plt

        plt.ion()
    if close_plots:
        plt.close("all")

    if plot:
        plot_spec(getTestFile("N2C_specair_380nm.spec"))


def test_save_compressed2(verbose=True, *args, **kwargs):
    "Check if saving a spectrum with compress = 2 does not change something."
    import shutil
    from os.path import join

    from radis import SpecDatabase, calc_spectrum
    from radis.test.utils import setup_test_line_databases

    shutil.rmtree(join(dirname(getTestFile(".")), "newDb"), ignore_errors=True)

    # get the spectrum
    setup_test_line_databases()
    s = calc_spectrum(
        2000,
        2300,  # cm-1
        molecule="CO",
        isotope="1,2,3",
        pressure=1.01325,  # bar
        Tgas=700,  # K
        mole_fraction=0.1,
        path_length=1,  # cm
        verbose=False,
        wstep=0.01,
        medium="vacuum",
        databank="HITRAN-CO-TEST",
    )
    try:

        # load in one databse
        db = SpecDatabase(join(dirname(getTestFile(".")), "newDb"))
        db.add(s, compress=2, if_exists_then="error")

        # simulate an experimentalist who come later and load the spectrum
        db2 = SpecDatabase(join(dirname(getTestFile(".")), "newDb"))
        s_bis = db2.get_unique(Tgas=700)

    finally:
        # we want to make sure this folder is deleted
        shutil.rmtree(join(dirname(getTestFile(".")), "newDb"))

    # we check the loaded spectrum contains less information than the calculated one
    assert not s == s_bis
    assert s_bis.get_vars() == ["abscoeff"]  # only this spectral quantity was stored
    assert s_bis.lines is None
    assert s_bis.conditions is not None  # we kept the metadata
    # now we check if it works
    s_bis.update()
    for var in s.get_vars():
        assert s.compare_with(s_bis, spectra_only=var, plot=False, verbose=verbose)


if __name__ == "__main__":
    import pytest

    # -s is to plot
    pytest.main(["test_database.py", "-s"])
    # test_save_compressed2()
