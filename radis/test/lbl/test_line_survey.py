# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 10:14:40 2017

@author: erwan

Examples
--------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that have a 'fast' label)

>>> pytest -m fast

"""

from radis.lbl import SpectrumFactory
from radis.misc.printer import printm
from radis.misc.utils import DatabankNotFound
from radis.test.utils import IgnoreMissingDatabase
import pytest


@pytest.mark.fast
def test_line_survey(verbose=True, plot=False, warnings=True, *args, **kwargs):

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
            db_use_cached=True)  # 0.2)
        pl.warnings['MissingSelfBroadeningWarning'] = 'ignore'
        pl.load_databank('HITRAN-CO2-TEST')
        s = pl.eq_spectrum(Tgas=1500)
        s.apply_slit(0.5)
        if plot:
            s.line_survey(overlay='transmittance', barwidth=0.01)

        if verbose:
            printm('no boolean defined for test_line_survey')
        return True  # test not defined (just testing methods work)

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


def _run_testcases(plot=False, verbose=True, *args, **kwargs):

    # Show media line_shift
    test_line_survey(plot=plot, verbose=verbose, *args, **kwargs)

    return True


if __name__ == '__main__':

    printm('Testing line survey functions:', _run_testcases(plot=True))
