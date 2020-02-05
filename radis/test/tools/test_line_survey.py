# -*- coding: utf-8 -*-
"""
Test that line survey works

-------------------------------------------------------------------------------

"""

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.misc.utils import DatabankNotFound
from radis.test.utils import (
    getTestFile,
    IgnoreMissingDatabase,
    setup_test_line_databases,
)
from radis.misc.printer import printm
from radis import SpectrumFactory
from radis.tools.database import load_spec
import os
from os.path import exists
import pytest


@pytest.mark.fast
def test_line_survey(verbose=True, plot=False, warnings=True, *args, **kwargs):
    """ Test line survey """

    _temp_file = "radis_test_line_survey.html"
    if exists(_temp_file):
        os.remove(_temp_file)

    s = load_spec(getTestFile("CO_Tgas1500K_mole_fraction0.01.spec"), binary=True)
    s.line_survey(overlay="abscoeff", display=plot, filename=_temp_file)

    assert exists(_temp_file)

    if verbose:
        print("test_line_survey: html file was correctly generated")

    if not plot:
        # clean after use
        os.remove(_temp_file)

    return True


@pytest.mark.fast
def test_line_survey_CO2(verbose=True, plot=True, warnings=True, *args, **kwargs):

    setup_test_line_databases()

    try:
        pl = SpectrumFactory(
            wavenum_min=2380,
            wavenum_max=2400,
            #                         wavelength_min=4170,
            #                         wavelength_max=4200,
            mole_fraction=400e-6,
            path_length=100,  # cm
            parallel=False,
            cutoff=1e-30,
            isotope=[1],
            save_memory=True,
            db_use_cached=True,
        )  # 0.2)
        pl.warnings["MissingSelfBroadeningWarning"] = "ignore"
        pl.load_databank("HITRAN-CO2-TEST")
        s = pl.eq_spectrum(Tgas=1500)
        s.apply_slit(0.5)
        if plot:
            s.line_survey(overlay="transmittance", barwidth=0.01)

        if verbose:
            printm("no boolean defined for test_line_survey")
        return True  # test not defined (just testing methods work)

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


def _run_testcases(plot=True, verbose=True, *args, **kwargs):

    # Show media line_shift
    test_line_survey(plot=plot, verbose=verbose, *args, **kwargs)
    test_line_survey_CO2(plot=plot, verbose=verbose, *args, **kwargs)

    return True


if __name__ == "__main__":

    print("Testing line survey functions:", _run_testcases(plot=True))
