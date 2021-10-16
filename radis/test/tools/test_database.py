# -*- coding: utf-8 -*-
"""
Test that :class:`~radis.tools.database.SpecDatabase` works

-------------------------------------------------------------------------------

"""


import os
from os.path import dirname, exists

import pytest

from radis.test.utils import getTestFile
from radis.tools.database import SpecDatabase, SpecList


@pytest.mark.fast
def test_speclist(*args, **kwargs):
    """Test :py:class:`~radis.tools.database.SpecList` functions"""

    from radis.test.utils import getTestFile

    # Get a list of test files
    spec_list = SpecDatabase(getTestFile("."), lazy_loading=False).get()

    # Now generate a SpecList
    dbl = SpecList(*spec_list)
    # check methods work:
    dbl.get_closest(Trot=1550)


@pytest.mark.fast
def test_database_functions(
    verbose=True, plot=True, close_plots=True, warnings=True, *args, **kwargs
):
    """Test :py:class:`~radis.tools.database.SpecDatabase` functions"""

    if plot:
        import matplotlib.pyplot as plt

        plt.ion()
    if close_plots:
        plt.close("all")

    db = SpecDatabase(dirname(getTestFile(".")), lazy_loading=False)

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
            db.add(s2, add_info=["Tvib", "Trot"], if_exists_then="error")

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


def test_lazy_loading(verbose=True, *args, **kwargs):
    """Test lazy loading"""

    import numpy as np

    from radis import calc_spectrum

    DBNAME = "test_CO_OH_database"
    db = SpecDatabase(DBNAME, lazy_loading=False)  # create the database

    #%% Simulate a fake database

    # Compute spectra of CO and OH and add them to the database
    T_db = np.linspace(300, 1200, 5)  # Temperature array
    for T in T_db:
        s1 = (
            calc_spectrum(
                wavenum_min=3000,
                wavenum_max=4500,
                molecule="CO",
                Tgas=T,
                databank="hitemp",  # test by fetching directly
                verbose=False,
            )
            .apply_slit(2, "nm")
            .take("radiance")
            # TODO : use https://github.com/radis/radis/issues/135 once implemented
        )
        db.add(s1, store_name=f"CO_{T}K.spec", if_exists_then="ignore")
        s2 = (
            calc_spectrum(
                wavenum_min=3000,
                wavenum_max=4500,
                molecule="OH",
                Tgas=T,
                databank="hitemp",  # test by fetching directly
                verbose=False,
            )
            .apply_slit(2, "nm")
            .take("radiance")
            # TODO : use https://github.com/radis/radis/issues/135 once implemented
        )
        db.add(s2, store_name=f"OH_{T}K.spec", if_exists_then="ignore")

    #%% Access the spectra in the database via the get functions.
    # Check what species are in the database at a given temperature (300)
    s300 = db.get(Tgas=300)  # simple use of get
    if verbose:
        print([s.conditions["molecule"] for s in s300])
    # Check what species are in the database at a given species
    sCO = db.get(molecule="CO")
    if verbose:
        print([s.conditions["Tgas"] for s in sCO])
    # Advanced use of the get function
    sdb = db.get('Tgas>=600 & Tgas<900 & molecule=="CO"')
    if verbose:
        print("Number of spectra found for the given conditions: ", len(sdb))

    #%% Do the same in lazy-loading mode, and compare :

    db2 = SpecDatabase(DBNAME, lazy_loading=True)  # create the database

    #%% Access the spectra in the database via the get functions.
    # Check what species are in the database at a given temperature (300)
    s300_2 = db2.get(Tgas=300)  # simple use of get
    print([s.conditions["molecule"] for s in s300_2])

    assert s300 == s300_2

    # Check what species are in the database at a given species
    sCO_2 = db2.get(molecule="CO")
    if verbose:
        print([s.conditions["Tgas"] for s in sCO_2])

    assert sCO == sCO_2

    # Advanced use of the get function
    sdb_2 = db2.get('Tgas>=600 & Tgas<900 & molecule=="CO"')
    if verbose:
        print("Number of spectra found for the given conditions: ", len(sdb))

    assert sdb == sdb_2


if __name__ == "__main__":
    import pytest

    # -s is to plot
    pytest.main(["test_database.py", "-s"])
