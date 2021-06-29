# -*- coding: utf-8 -*-
"""Tools to test RADIS library.

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


import os
from os.path import dirname, exists, join

from radis.db.utils import getFile
from radis.misc.config import (
    addDatabankEntries,
    diffDatabankEntries,
    getDatabankEntries,
    getDatabankList,
)

TEST_FOLDER_PATH = join(dirname(dirname(__file__)), "test")


def getTestFile(file, force=False):
    """Return the full path of a test file, if it exists. Used by test
    functions not to worry about the project architecture. Using :ref:`test
    files <label_dev_test_files>` is recommended when writing tests.

    Parameters
    ----------

    file: str
        filename. See :ref:`the list of available test files <label_dev_test_files>`

    Returns
    -------

    path: str
        absolute path of ``file`` on the local machine. Raises an error
        if test file not present, unless you use ``force=True``

    Examples
    --------

    ::

        from radis.test.utils import getTestFile
        from radis import load_spec
        load_spec(getTestFile('CO_Tgas1500K_mole_fraction0.01.spec'))

    See Also
    --------

    :py:func:`~radis.test.utils.getValidationCase`
    """

    path = join(TEST_FOLDER_PATH, "files", file)

    if not exists(path) and not force:
        raise FileNotFoundError(
            "Test file `{0}` does not exist. Choose one of: \n- {1} or use force=True".format(
                file, "\n- ".join(os.listdir(join(TEST_FOLDER_PATH, "files")))
            )
        )

    return path


def getValidationCase(file, force=False):
    """Return the full path of a validation case file, if it exists. Used by
    test functions not to worry about the project architecture. Using
    :ref:`validation test files <label_dev_test_files>` is recommended when
    writing validation cases.

    Parameters
    ----------

    file: str
        filename. See :ref:`the list of available validation files <label_dev_test_files>`

    Returns
    -------

    path: str
        absolute path of ``file`` on the local machine. Raises an error
        if validation file not present, unless you use ``force=True``

    Examples
    --------

    Load the reference case from the [Klarenaar2017]_ paper ::

        from radis.test.utils import getValidationCase
        from radis import Spectrum

        s_exp = Spectrum.from_txt(
            getValidationCase(
                join(
                    "test_CO2_3Tvib_vs_klarenaar_data", "klarenaar_2017_digitized_data.csv",
                )
            ),
            "transmittance_noslit",
            wunit="cm-1",
            unit="",
            delimiter=",",
            name="Klarenaar 2017",
        )

    See Also
    --------

    :py:func:`~radis.test.utils.getTestFile`
    """

    path = join(TEST_FOLDER_PATH, "validation", file)

    if not exists(path) and not force:
        raise FileNotFoundError(
            "Validation case `{0}` does not exist. Choose one of: \n- {1} or use force=True".format(
                file, "\n- ".join(os.listdir(join(TEST_FOLDER_PATH, "validation")))
            )
        )

    return path


# Python 3.6+ only
getTestFile.__annotations__["file"] = os.listdir(join(TEST_FOLDER_PATH, "files"))
getValidationCase.__annotations__["file"] = os.listdir(
    join(TEST_FOLDER_PATH, "validation")
)


# %% Convenience function


def test_spectrum(**kwargs):
    """Generate the :ref:`first example spectrum <label_first_example>` with ::

        import radis
        s = radis.test_spectrum()
        s.plot()

    Other Parameters
    ----------------
    kwargs: sent to :py:func:`~radis.lbl.calc.calc_spectrum`


    """
    from radis import calc_spectrum

    conditions = {
        "wavenum_min": 1900,
        "wavenum_max": 2300,
        "molecule": "CO",
        "isotope": "1,2,3",
        "pressure": 1.01325,  # bar
        "Tgas": 700,  # K
        "mole_fraction": 0.1,
        "path_length": 1,  # cm
        "databank": "hitran",
    }

    conditions.update(kwargs)

    s = calc_spectrum(**conditions)
    return s


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
    "HITEMP-CO2-HAMIL-TEST": {
        "info": "HITEMP-2010, CO2, 3 main isotope (CO2-626, 636, 628), "
        + "2283.7-2285.1 cm-1, energies calculated from Tashkun effective hamiltonian",
        "path": [getTestFile(r"cdsd_hitemp_09_fragment.txt")],
        "format": "cdsd-hitemp",  # CDSD-HITEMP version (same lines as HITEMP-2010).
        "parfunc": getFile("CO2", "partition_functions.txt"),
        "parfuncfmt": "cdsd",
        "levels": {1: getTestFile(r"co2_cdsd_hamiltonian_fragment.levels")},
        "levelsfmt": "cdsd-hamil",
    },
}
"""dict: test databases added in the :ref:`Configuration file <label_lbl_config_file>`
by :py:func:`~radis.test.utils.setup_test_line_databases`
"""

# %% Utils to test spec module


def setup_test_line_databases(verbose=True):
    """Build :py:data:`~radis.test.utils.TEST_DATABASES` and add them in
    ~/radis.json. Generate the file if it  doesnt exist.

    In particular:

    - HITRAN-CO2-TEST: CO2, HITRAN 2016, 4165-4200 nm
    - HITRAN-CO-TEST: CO, HITRAN 2016, 2000-2300 cm-1
    - HITEMP-CO2-TEST: CO2, HITEMP-2010, 2283.7-2285.1 cm-1, 3 isotopes
    - HITEMP-CO2-HAMIL-TEST: same as previous, with (some) energy levels computed
      from Tashkun effective Hamiltonian.


    These test databases are used to run the different test routines. They can
    obviously be used by Users to run simulations, but we suggest Users to download
    their own line databases files and add them to ~/radis.json so they have more control
    on it

    Examples
    --------

    Initialize the Line databases::

        from radis import setup_test_line_databases
        setup_test_line_databases()

    Plot a CO2 spectrum at high temperature::

        from radis import calc_spectrum
        calc_spectrum(2284,
                      2285,
                      Tgas=2000,
                      pressure=1,
                      molecule='CO2',
                      isotope=1
                      databank='HITEMP-CO2-TEST').plot()

    Note that 'HITEMP-CO2-TEST' is defined on 2283.7-2285.1 cm-1 only, as
    can be shown by reading the Database information:

        from radis.misc.config import printDatabankEntries
        printDatabankEntries('HITEMP-CO2-TEST')

        >>> 'HITEMP-CO2-TEST':
        >>> {'info': 'HITEMP-2010, CO2, 3 main isotope (CO2-626, 636, 628), 2283.7-2285.1 cm-1',
        >>> 'path': ['/USER/PATH/TO\\radis\\radis\\test\\files\\cdsd_hitemp_09_fragment.txt'],
        >>> 'format': 'cdsd-hitemp'
        >>> 'parfuncfmt': 'hapi'
        >>> 'levelsfmt': 'radis'


    See Also
    --------

    :ref:`Configuration file <label_lbl_config_file>`,
    :py:func:`~radis.misc.config.getDatabankList`,
    :py:func:`~radis.misc.config.printDatabankEntries`
    """
    # TODO: generate large band databases for the main species (let's say CO2,
    # H2O and CH4) and main isotopes by fetching the HITRAN 2016 database.

    # Get list of databases
    try:
        dbnames = getDatabankList()
    except FileNotFoundError:
        dbnames = []

    # %% Add test databases
    for dbname, dbentries in TEST_DATABASES.items():

        if dbname in dbnames:  # Check entries are correct
            #            for k
            diff = diffDatabankEntries(
                getDatabankEntries(dbname), dbentries, verbose=False
            )
            if diff is not None:
                raise ValueError(
                    "{0}".format(diff)
                    + "\nIn ~/radis.json\n----------\n{0}".format(
                        getDatabankEntries(dbname)
                    )
                    + "\n\nExpected\n---------\n{0}\n\n".format(dbentries)
                    + "Test Database {0} doesnt match expected ".format(dbname)
                    + "entries for key `{0}`. See comparison above. ".format(diff)
                    + "To regenerate test databases just delete the {0} ".format(dbname)
                    + "entry in your ~/radis.json"
                )

        else:  # add them (create ~/radis.json file if doesnt exist yet)
            addDatabankEntries(dbname, dbentries)

    return


# %% Edit existing Line databases


def define_Evib_as_sum_of_Evibi(levels):
    """Note that this is arbitrary for a polyatomic molecule. Lookup Pannier,
    Dubuet and Laux 2020 for more.

    We also update Erot to maintain the sum Evib+Erot = E :

    ::

        Evib = Evib1 + Evib2 + Evib3
        Erot = E - Evib    # to be consistent with equilibrium
    """

    levels["Evib"] = levels.Evib1 + levels.Evib2 + levels.Evib3
    levels["Erot"] = levels.E - levels.Evib

    return levels


def define_Evib_as_min_of_polyad(levels, keys):
    """Here we define the vibrational energy as the minimum energy in a polyad.
    Here, the polyad is defined for each combination of ``keys`` Typically,
    ``keys=['p', 'c', 'N']`` or keys=['p', 'c'].

    Rotational energy is the rest::

        Evib = min(E(p,c,j,n) for a given set of (p,c))
        Erot = E - Evib

    .. warning::
        See Pannier, Dubuet & Laux 2020 for a quantitative comparison
        of the different possible methods to define vibrational energy.


    Parameters
    ----------

    sf: SpectrumFactory object
    """

    def fill_EvibErot(grp):
        Evib0 = grp.E.min()
        grp["Evib"] = Evib0
        grp["Erot"] = grp.E - Evib0
        return grp

    levels = levels.groupby(keys).apply(fill_EvibErot)
    levels.reset_index()

    return levels


def discard_lines_with_na_levels(sf):
    """In the test Levels databases, not all levels are given (to save space).
    Consequently, in the Line databases, some lines have N/A levels and cannot
    be calculated at nonequilibrium. This function cleans the line databases
    from such lines by first running a dummy calculation and removing the lines
    where levels were N/A.

    .. warning::
        results from such a calculation are physically wrong. Only use
        to test the functions!

    Parameters
    ----------

    sf: SpectrumFactory
    """

    # Calculate populations using the non-equilibrium module:
    # This will crash the first time because the Levels Database is just a fragment and does not include all levels.
    try:
        sf.non_eq_spectrum(300, 300)
    except AssertionError:  # expected
        sf.df0.dropna(inplace=True)


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


if __name__ == "__main__":
    #    run_tests()
    setup_test_line_databases()
