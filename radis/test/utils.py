# -*- coding: utf-8 -*-
""" Tools to test RADIS library

Summary
-------

Tools to test RADIS library

Examples
--------

Run all tests::
    
    cd radis/test
    pytest
    
Run only "fast" tests (tests that have a "fast" label, and should be 
a few seconds only)::
    
    cd radis/test
    pytest -m fast

-------------------------------------------------------------------------------

"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os
from radis.misc.config import (
    getDatabankList,
    getDatabankEntries,
    addDatabankEntries,
    diffDatabankEntries,
)
from radis.misc.utils import FileNotFoundError
from radis.misc.printer import printr
from os.path import join, dirname

TEST_FOLDER_PATH = join(dirname(dirname(__file__)), "test")
IGNORE_MISSING_DATABASES = False  # TODO: move in ~/.radis
"""bool: how to deal with missing Line databases during tests. If ``True``, 
print a warning. Else, raise an error. See :func:`~radis.test.utils.IgnoreMissingDatabase`
"""


def getTestFile(file):
    """ Return the full path of a test file. Used by test functions not to
    worry about the project architecture"""

    return join(TEST_FOLDER_PATH, "files", file)


def getValidationCase(file):
    """ Return the full path of a validation case file. Used by test functions not to
    worry about the project architecture"""

    return join(TEST_FOLDER_PATH, "validation", file)


try:  # Python 3.6 only
    getTestFile.__annotations__["file"] = os.listdir(join(TEST_FOLDER_PATH, "files"))
    getValidationCase.__annotations__["file"] = os.listdir(
        join(TEST_FOLDER_PATH, "validation")
    )
except:
    pass

# %% Comparison functions


def testEqual(a, b, info=""):
    if a != b:
        print("Mismatch", info, ":", a, "!=", b)
    return a == b


# %% Test Databases


TEST_DATABASES = {
    "HITRAN-CO2-TEST": {
        "info": "HITRAN 2016 database, CO2, 1 main isotope (CO2-626), bandhead: "
        + "2380-2398 cm-1 (4165-4200 nm)",
        "path": [getTestFile(r"hitran_co2_626_bandhead_4165_4200nm.par")],
        "format": "hitran",
        "parfuncfmt": "hapi",
        "levelsfmt": "radis",
    },
    "HITRAN-CO-TEST": {
        "info": "HITRAN 2016 database, CO, 3 main isotopes (CO-26, 36, 28), "
        + "2000-2300 cm-1",
        "path": [getTestFile(r"hitran_co_3iso_2000_2300cm.par")],
        "format": "hitran",
        "parfuncfmt": "hapi",
        "levelsfmt": "radis",
    },
    "HITEMP-CO2-TEST": {
        "info": "HITEMP-2010, CO2, 3 main isotope (CO2-626, 636, 628), "
        + "2283.7-2285.1 cm-1",
        "path": [getTestFile(r"cdsd_hitemp_09_fragment.txt")],
        "format": "cdsd-hitemp",  # CDSD-HITEMP version (same lines as HITEMP-2010).
        "parfuncfmt": "hapi",
        "levelsfmt": "radis",
    },
}
"""dict: test databases added in the :ref:`Configuration file <label_lbl_config_file>`
by :py:func:`~radis.test.utils.setup_test_line_databases`
"""

# %% Utils to test spec module


def setup_test_line_databases(verbose=True):
    """ Build :py:data:`~radis.test.utils.TEST_DATABASES` and add them in ~/.radis. 
    Generate the file if it  doesnt exist

    In particular:

    - HITRAN-CO2-TEST: CO2, HITRAN 2016, 4165-4200 nm 
    - HITRAN-CO-TEST: CO, HITRAN 2016, 2000-2300 cm-1
    - HITEMP-CO2-TEST: CO2, HITEMP-2010, 2283.7-2285.1 cm-1, 3 isotopes

    These test databases are used to run the different test routines. They can
    obviously be used by Users to run simulations, but we suggest Users to download
    their own line databases files and add them to ~/.radis so they have more control
    on it
    
    See Also
    --------
    
    :ref:`Configuration file <label_lbl_config_file>`

    """
    # TODO: generate large band databases for the main species (let's say CO2,
    # H2O and CH4) and main isotopes by fetching the HITRAN 2016 database.

    # Get list of databases
    try:
        dbnames = getDatabankList()
    except FileNotFoundError:
        dbnames = []

    # %% Add test databases

    def add_to_parser(config, name, dic):
        for k, v in dic.items():
            config[name][k] = v
        if verbose:
            print("Adding '{0}' database in ~/.radis".format(name))

    for dbname, dbentries in TEST_DATABASES.items():

        if dbname in dbnames:  # Check entries are correct
            #            for k
            diff = diffDatabankEntries(
                getDatabankEntries(dbname), dbentries, verbose=False
            )
            if diff is not None:
                raise ValueError(
                    "{0}".format(diff)
                    + "\nIn ~/.radis\n----------\n{0}".format(
                        getDatabankEntries(dbname)
                    )
                    + "\n\nExpected\n---------\n{0}\n\n".format(dbentries)
                    + "Test Database {0} doesnt match expected ".format(dbname)
                    + "entries for key `{0}`. See comparison above. ".format(diff)
                    + "To regenerate test databases just delete the {0} ".format(dbname)
                    + "entry in your ~/.radis"
                )

        else:  # add them (create ~/.radis file if doesnt exist yet)
            addDatabankEntries(dbname, dbentries)

    return


# %% Deal with missing databases
def _failsafe_if_no_db(testcase, *args, **kwargs):
    """finally not implemented?"""
    from radis.misc.utils import DatabankNotFound

    try:
        testcase(*args, **kwargs)
    except DatabankNotFound:
        import sys

        print((sys.exc_info()))
        print(
            (
                "Testing {0}: Database not defined. \n".format(testcase.__name__)
                + "Ignoring the test"
            )
        )
        return True


def IgnoreMissingDatabase(err, file="", warnings=True):
    """ A function to deal with MissingDatabases errors. If :data:`~radis.test.utils.IGNORE_MISSING_DATABASES`
    is ``True``, just print a warning. Else, raise the error
    
    Parameters
    ----------
    
    err: an Error
    
    file: str
        where the error occured. Use ``file=__file__`` on function call
    """
    # TODO: make IGNORE_MISSING_DATABASES a ~/.radis parameter
    if IGNORE_MISSING_DATABASES:
        if warnings:
            import sys

            print(sys.exc_info())
            printr(
                "In {0}: Database not defined: {1}".format(file, err.filename)
                + "\n Ignoring the test"
            )
        return True
    else:
        raise err


if __name__ == "__main__":
    #    run_tests()
    setup_test_line_databases()
