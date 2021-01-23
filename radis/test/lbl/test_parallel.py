# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 09:38:23 2017

@author: erwan

Use
------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that a 'fast' label)

>>> pytest -m fast

"""

from time import time
from warnings import warn

from radis.lbl.parallel import ParallelFactory, SpectrumFactory
from radis.misc.printer import printm
from radis.test.utils import setup_test_line_databases

# %% Main


def test_parallel(plot=False, verbose=True, warnings=True, *args, **kwargs):
    """Core of the parallel test routine

    Test than :class:`~radis.lbl.parallel.ParallelFactory` works.

    Note that it is not always faster than the sequential mode : with the new
    LDM algorithm (>0.9.20), the line broadening step is only weakly dependant
    on the number of lines and therefore the parallelizatino overhead may
    not be necessary.
    """

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        import matplotlib.pyplot as plt

        plt.ion()

    setup_test_line_databases()

    # Performance test
    Tvib = 1000  # , 1200, 1400]
    Trot = [900, 1200, 1250, 1300, 1400, 1600]
    path_length = 0.1
    mole_fraction = 0.5

    cutoff = 0
    broadening_max_width = 50
    wstep = 0.01

    # .. parallel
    pf = ParallelFactory(
        wavelength_min=4165,
        wavelength_max=4200,
        cutoff=cutoff,
        db_use_cached=True,
        lvl_use_cached=True,
        isotope="1",
        broadening_max_width=broadening_max_width,
        wstep=wstep,
        verbose=verbose,
    )
    pf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    pf.warnings["NegativeEnergiesWarning"] = "ignore"
    pf.warnings["EmptyDatabaseWarning"] = "ignore"
    pf.load_databank("HITRAN-CO2-TEST")

    t0 = time()
    s_parallel = pf.non_eq_spectrum(
        Tvib, Trot, path_length=path_length, mole_fraction=mole_fraction
    )
    t0 = time() - t0

    # ... non parallel
    sf = SpectrumFactory(
        wavelength_min=4165,
        wavelength_max=4200,
        cutoff=cutoff,
        db_use_cached=True,
        lvl_use_cached=True,
        isotope="1,2",
        broadening_max_width=broadening_max_width,
        wstep=wstep,
        verbose=verbose,
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.warnings["EmptyDatabaseWarning"] = "ignore"
    sf.load_databank("HITRAN-CO2-TEST")
    t1 = time()
    s_sequential = []
    for Troti in Trot:
        s_sequential.append(
            sf.non_eq_spectrum(
                Tvib, Troti, path_length=path_length, mole_fraction=mole_fraction
            )
        )
    t1 = time() - t1

    if verbose:
        printm(("\nParallelFactory performance test ({0} cases)".format(len(Trot))))
        printm("-------------------")
        printm(("Parallel version: {0:.1f}s".format(t0)))
        printm(("Normal version: {0:.1f}s".format(t1)))
        printm("-------------------")

    # Assert that spectra are the same
    for i in range(len(Trot)):
        assert s_parallel[i].compare_with(
            s_sequential[i], spectra_only=True, plot=False
        )
    if verbose:
        printm("Ensured that spectra are the same: OK")

    #        # parallel version should be at least faster!
    #        # note: ...  may not be the case if the calculated ranges are too small
    #        # ......... as is the case in the test database. Better use // for spectra
    #        # ......... with calculation times > 5s
    #        assert t0 > t1
    if t0 > t1:
        warn(
            "Parallel version ({0:.1f}s) is slower than the Serial version ({1:.1f}s) !".format(
                t0, t1
            )
        )

    return True  # no plot implemented


def _run_testcases(plot=False, verbose=True, *args, **kwargs):

    test_parallel(plot=plot, verbose=verbose, *args, **kwargs)


if __name__ == "__main__":

    printm(("Tested parallel.py:", _run_testcases(plot=True)))
