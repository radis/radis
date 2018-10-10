# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 10:18:48 2017

@author: erwan

Examples
--------

Run all tests:

>>> pytest       (in command line, in project folder)

Run only fast tests

>>> pytest -m fast

"""

from radis.lbl import SpectrumFactory
from radis.spectrum.spectrum import calculated_spectrum, transmittance_spectrum
from radis.tools.database import load_spec
from radis.los.slabs import SerialSlabs
from radis.tools.slit import (gaussian_slit, triangular_slit, trapezoidal_slit,
                              import_experimental_slit, convolve_with_slit)
from radis.phys.units import is_homogeneous
from radis.misc.utils import DatabankNotFound
from radis.misc.printer import printm
from radis.test.utils import setup_test_line_databases
from radis.test.utils import IgnoreMissingDatabase
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, linspace, abs, trapz
from os.path import basename
import pytest

fig_prefix = basename(__file__)+': '

# %% --------------------------------------------------------------------------
#  Test cases
# -----------------------------------------------------------------------------


@pytest.mark.needs_db_CDSD_HITEMP
def test_constant_source(verbose=True, plot=True, warnings=True, *args, **kwargs):
    ''' Test transmittance = radiance pour une source égale à 1  '''

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    try:

        wmin = 4100
        wmax = 4600
        path_length = 500

        pl = SpectrumFactory(
            wavelength_min=wmin,
            wavelength_max=wmax,
            db_use_cached=True,
            mole_fraction=400e-6,
            path_length=path_length,
            cutoff=1e-20,
            wstep=0.1, 
            isotope=[1, 2])
        pl.warnings.update({'MissingSelfBroadeningWarning': 'ignore',
                            'NegativeEnergiesWarning': 'ignore',
                            'VoigtBroadeningWarning': 'ignore'})
        pl.load_databank('CDSD-HITEMP')
        s0 = pl.eq_spectrum(Tgas=300)
        w0 = s0.get('radiance_noslit')[0]

        s_1 = calculated_spectrum(w0, np.ones_like(w0), Iunit='mW/cm2/sr/nm')
        #s_1 = theoretical_spectrum(nm2cm(w0), np.ones_like(w0), Iunit='mW/cm2/sr/nm', wunit='cm-1')

        s = SerialSlabs(s_1, s0)
        s.apply_slit(1.5, plot_slit=plot)        # in nm (s ~ s_1 in nm)
        s0.apply_slit(1.5, plot_slit=plot)       # in nm (s0 in cm-1)

        if plot:
            fig = plt.figure(
                fig_prefix+'Transmittance vs radiance for constant source')
            s.plot('radiance', nfig=fig.number,
                   label='radiance constant source', force=True)
            s0.plot('transmittance', nfig=fig.number,
                    color='r', label='transmittance', force=True)
            plt.legend()
            plt.title('transmittance = radiance for constant source')
            # Todo: automatize test output

        if verbose:
            printm('\n>>> _test_constant_source : NOT DEFINED\n')
        return True  # nothing defined yet

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


@pytest.mark.fast
def test_resampling(rtol=1e-2, verbose=True, plot=True, warnings=True, *args, **kwargs):
    ''' Test what happens when a spectrum in nm or cm-1, is convolved
    with a slit function in nm. In particular, slit function is generated
    in the spectrum unit, and spectrum is resampled if not evenly spaced'''

    if verbose:
        printm('Test auto resampling')

    if plot:  # Make sure matplotlib is interactive so that test are not stuck in pytest
        plt.ion()

    try:

        setup_test_line_databases()  # add HITRAN-CO-TEST in neq.rc if not there

        plCO = SpectrumFactory(
            wavenum_min=2230,
            wavenum_max=2260,
            mole_fraction=0.02,
            path_length=100,  # cm
            broadening_max_width=20,  # cm^-1
            wstep=0.02,
            isotope=[1, 2, 3],
            verbose=verbose,
        )
        plCO.warnings['MissingSelfBroadeningWarning'] = 'ignore'
        plCO.load_databank('HITRAN-CO-TEST')
        sCO = plCO.eq_spectrum(Tgas=300)

        w_nm, T_nm = sCO.get('transmittance_noslit', wunit='nm', medium='air')
        w_nm, I_nm = sCO.get('radiance_noslit', wunit='nm',
                             Iunit='mW/cm2/sr/nm')
        sCO_nm = transmittance_spectrum(w_nm, T_nm, wunit='nm',
                                        conditions={'medium': 'air'})  # a new spectrum stored in nm
        # sCO_nm = theoretical_spectrum(w_nm, I_nm, wunit='nm', Iunit='mW/cm2/sr/nm') #  a new spectrum stored in nm

        if plot:
            fig = plt.figure(fig_prefix+'auto-resampling')
            sCO.plot('transmittance_noslit', wunit='cm-1', nfig=fig.number,
                     marker='o', color='k', lw=3, ms=10, label='(stored in cm-1)')
            plt.title('No slit function')
            sCO_nm.plot('transmittance_noslit', wunit='cm-1', nfig=fig.number,
                        marker='o', color='r', label='(stored in nm)')
#            plt.xlim((2246.58, 2247.52))
#            plt.ylim((0.87, 1.01))
            plt.legend()

        slit_function = 0.8
        slit_unit = 'cm-1'
        sCO.apply_slit(slit_function, unit=slit_unit)
        sCO_nm.apply_slit(slit_function, unit=slit_unit)

        if plot:
            fig = plt.figure(fig_prefix+'auto-resampling (after convolution)')
            sCO.plot('transmittance', wunit='cm-1', nfig=fig.number,
                     marker='o', color='k', lw=3, ms=10, label='(stored in cm-1)')
            plt.title('Slit function: {0} {1}'.format(
                slit_function, slit_unit))
            sCO_nm.plot('transmittance', wunit='cm-1', nfig=fig.number,
                        marker='o', color='r', label='(stored in nm)')

#            plt.xlim((2246.58, 2247.52))
#            plt.ylim((0.87, 1.01))
            plt.legend()

        w_conv, T_conv = sCO.get('transmittance', wunit='cm-1')
        w_nm_conv, T_nm_conv = sCO_nm.get('transmittance', wunit='cm-1')

        error = abs((trapz(1-T_conv, w_conv)-trapz(1-T_nm_conv,
                                                   w_nm_conv))/trapz(1-T_nm_conv, w_nm_conv))

        if verbose:
            printm('\n>>> _test_resampling\n')
        if verbose:
            printm('Error between 2 spectra ({0:.2f}%) < {1:.2f}%: {2}'.format(error*100,
                                                                               rtol*100, bool(error < rtol)))
        assert bool(error < rtol)

    except DatabankNotFound as err:
        assert IgnoreMissingDatabase(err, __file__, warnings)


def _run_testcases(plot=False, verbose=True, *args, **kwargs):

    #    test_against_specair_convolution(plot=plot, verbose=verbose,     # Moved in RADIS
    #                                           *args, **kwargs)
    #    test_normalisation_mode(plot=plot, verbose=verbose, *args, **kwargs)         # Moved in RADIS
    #    test_slit_energy_conservation(plot=plot, verbose=verbose, *args, **kwargs)   # Moved in RADIS
    #    test_slit_function_effect(plot=plot, verbose=verbose, *args, **kwargs)       # Moved in RADIS
    test_constant_source(plot=plot, verbose=verbose, *args, **kwargs)
#    test_all_slits(plot=plot, verbose=verbose, *args, **kwargs)       # Moved in RADIS
    test_resampling(plot=plot, verbose=verbose, *args, **kwargs)

    return True


if __name__ == '__main__':
    printm('Testing slit.py: ', _run_testcases(plot=True))
