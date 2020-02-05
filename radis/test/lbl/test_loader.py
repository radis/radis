# -*- coding: utf-8 -*-
"""
Created on Mon May  7 17:34:52 2018

@author: erwan
"""

from __future__ import absolute_import, unicode_literals, division, print_function
from radis.lbl import SpectrumFactory
from radis.misc.utils import DatabankNotFound
from radis.misc.printer import printm
from radis.test.utils import IgnoreMissingDatabase, setup_test_line_databases
import matplotlib.pyplot as plt
from os.path import exists
from shutil import rmtree


def test_retrieve_from_database(
    plot=True, verbose=True, warnings=True, *args, **kwargs
):
    """ Test autoretrieve from a database:
        
        first generate an empty :py:class:`~radis.tools.database.SpecDatabase`
        associated to a :py:class`~radis.lbl.factory.SpectrumFactory`, 
        then calculate a first spectrum, then calculate it again and make sure
        it is retrieved from the database
        """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    temp_database_name = "temp_spec_database"

    try:
        if verbose:
            printm(">>> test_retrieve_from_database")

        assert not exists(temp_database_name)

        setup_test_line_databases()  # add HITEMP-CO2-TEST in ~/.radis if not there

        sf = SpectrumFactory(
            2284.2,
            2284.6,
            wstep=0.001,  # cm-1
            pressure=20 * 1e-3,  # bar
            cutoff=0,
            path_length=0.1,
            mole_fraction=400e-6,
            molecule="CO2",
            isotope=[1],
            db_use_cached=True,
            medium="vacuum",
            broadening_max_width=10,
            export_populations="rovib",
            verbose=verbose,
        )
        sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
        sf.init_databank(
            "HITEMP-CO2-TEST"
        )  # unlike load_databank, will automatically be built when needed
        db = sf.init_database(temp_database_name, autoretrieve=True)

        assert len(db) == 0

        # Calculate a first spectrum
        s1 = sf.non_eq_spectrum(2000, 1000)

        assert len(db) == 1
        # Note that the Spectrum in database is not s1 because populations
        # are not stored by default

        # hack: now change the `autoretrieve` status to `force` to make sure
        # the spectrum is retrieved. Else, an error will be raised
        sf.autoretrievedatabase = "force"

        # Calculate spectrum under the same conditions
        s2 = sf.non_eq_spectrum(2000, 1000)

        # Note that s1 == s2 won't work because populations are not stored
        # by default in the database
        assert s1.compare_with(s2, spectra_only=True)

        return True

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)

    finally:
        rmtree(temp_database_name)


def _run_testcases(verbose=True, plot=True):

    test_retrieve_from_database(plot=plot, verbose=verbose)


if __name__ == "__main__":
    _run_testcases()
