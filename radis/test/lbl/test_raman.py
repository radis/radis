# -*- coding: utf-8 -*-
"""

Examples
--------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests (i.e: tests that a 'fast' label)

>>> pytest -m fast

"""

from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
from numpy import allclose, linspace
import matplotlib.pyplot as plt
from os.path import basename
from neq.spec.raman import Raman_Stokes_N2_1vib, Raman_rotational_N2
from neq.misc.printer import printm


def test_validation_Raman1vib_Lo2012(plot=True, verbose=True, debug=False, warnings=True,
                                     close_plots=True, *args, **kwargs):
    ''' 
    Validation on a spectrum computed by Lo2012, p234. The slit function is approxamative 
    explaining the small differences. Still good agreement
    '''

    if plot:  # Make sure matplotlib is interactive so that test are not stuck
        plt.ion()
        if close_plots:
            plt.close('all')

    s = Raman_Stokes_N2_1vib(532, 618, T_vib=2942, nb_r=500, nb_v=40)
    s.apply_slit(0.65, shape='gaussian')
    if plot:
        s.plot(normalize=1, nfig='Raman_Stokes_N2_1vib(532, T_gas=618,T_vib=2942)')

    if verbose:
        printm('no test defined for test_validation_Raman1vib_Lo2012')


def test_validation_Ramanrot_Limbach(plot=False, verbose=True, debug=False, warnings=True,
                                     close_plots=True, *args, **kwargs):
    ''' 
    Test 3 spectra taken from Limbach Thesis, p62. Good agreement for the width, 
    some perturbations may not be taken into account
    '''

    if plot:  # Make sure matplotlib is interactive so that test are not stuck
        plt.ion()
        if close_plots:
            plt.close('all')

    R294 = Raman_rotational_N2(532, 294)
    R2000 = Raman_rotational_N2(532, 2000)
    R4000 = Raman_rotational_N2(532, 4000)

    if plot:
        R294.plot(nfig='Raman_rotational_N2(532,294)')
        R2000.plot(nfig='Raman_rotational_N2(532,2000)')
        R4000.plot(nfig='Raman_rotational_N2(532,4000)')

    if verbose:
        printm('no test defined for test_validation_Ramanrot_Limbach')


def _run_all_tests(plot=True, verbose=True, *args, **kwargs):

    test_validation_Raman1vib_Lo2012(
        plot=plot, verbose=verbose, close_plots=False)
    test_validation_Ramanrot_Limbach(
        plot=plot, verbose=verbose, close_plots=False)


if __name__ == '__main__':
    printm('Test spectrum: ', _run_all_tests(debug=False, plot=True))
