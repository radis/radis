# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 10:59:49 2017

@author: erwan

Examples
--------

Run all tests::

    pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that a 'fast' label)::

    pytest -m fast


-------------------------------------------------------------------------------


"""

import os
from os.path import basename, exists

import matplotlib.pyplot as plt

from radis.lbl import SpectrumFactory
from radis.misc.printer import printm
from radis.test.utils import setup_test_line_databases
from radis.tools.database import load_spec

fig_prefix = basename(__file__) + ": "

# %% Test


def test_load_spectrum(plot=False, verbose=True, warnings=True, *args, **kwargs):
    """Test load / save

    Fast version: dont save lines / populations, compare spectra only
    """

    setup_test_line_databases()

    temp_file_name = "_test_database_co2_tempfile.spec"
    assert not exists(temp_file_name)

    try:
        sf = SpectrumFactory(
            wavelength_min=4190,
            wavelength_max=4200,
            mole_fraction=400e-6,
            path_length=0.1,  # cm
            isotope=[1],
            cutoff=1e-20,
            verbose=verbose,
        )
        sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
        sf.load_databank("HITRAN-CO2-TEST")
        s1 = sf.eq_spectrum(Tgas=300)
        s1.apply_slit(2, verbose=False)
        s1.update()
        s2 = load_spec(s1.store(temp_file_name, compress=True))
        s2.update()
        if plot:
            fig = plt.figure(fig_prefix + "Calculated vs stored+retrieved")
            s1.plot("absorbance", nfig=fig.number, lw=3, label="calculated")
            s2.plot(
                "absorbance",
                nfig=fig.number,
                color="r",
                label="stored (compressed) and retrieved",
            )
            plt.legend()

        assert s1.compare_with(s2, spectra_only=True, plot=False)

    finally:  # cleaning
        if exists(temp_file_name):
            os.remove(temp_file_name)

    return True


def test_load_lines_pops(plot=False, verbose=True, warnings=True, *args, **kwargs):
    """Test load / save

    Full version: save lines, populations (hundreds of MB for CO2, much less
    for CO), compare everything
    """

    temp_file_name = "_test_database_co_tempfile.spec"
    assert not exists(temp_file_name)

    try:
        sf = SpectrumFactory(
            wavelength_min=4500,
            wavelength_max=4800,
            mole_fraction=400e-6,
            path_length=0.1,  # cm
            isotope=[1, 2],
            cutoff=1e-20,
            verbose=verbose,
        )
        sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
        sf.load_databank("HITRAN-CO-TEST")
        s1 = sf.non_eq_spectrum(Tvib=300, Trot=300)
        s1.apply_slit(2, verbose=False)
        s1.update()
        s2 = load_spec(
            s1.store(
                temp_file_name,
                discard=[],  # we're actually saving lines and
                # populations here! expect hundreds
                # of megabytes for big molecules
                # :line CO2 ... here CO is fine
                # (by default lines and populations
                # are discarded)
                compress=True  # only removes some spectral
                # quantities, cant change population
                # or lines size
            ),
            binary=True,  # uncompress
        )
        s2.update()

        assert s1.compare_with(s2, spectra_only=False, plot=False)

        # Test with json-tricks directly
        # Does not work yet
    #        from json_tricks import dumps, loads
    #        s2b = loads(dumps(s1))
    #        assert s1.compare_with(s2b, spectra_only=False, plot=False)

    finally:  # cleaning
        if exists(temp_file_name):
            os.remove(temp_file_name)

    return True


def _run_testcases(plot=False, verbose=True, warnings=True, *args, **kwargs):

    test_load_spectrum(plot=plot, verbose=verbose, warnings=warnings, *args, **kwargs)

    test_load_lines_pops(plot=plot, verbose=verbose, warnings=warnings, *args, **kwargs)

    return True


if __name__ == "__main__":

    printm("Testing database.py: ", _run_testcases(plot=True))
