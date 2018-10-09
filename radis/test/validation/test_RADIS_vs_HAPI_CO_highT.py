# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 21:10:31 2017

@author: erwan

Validation case
------------

Plot absorption coefficient (cm-1) of CO at high temperature (2000 K) with 
RADIS, and compare with calculations from HAPI using the HITRAN database 
    
        Physical Conditions
    ----------------------------------------
       Tgas                 1500 K
       Trot                 1500 K
       Tvib                 1500 K
       isotope              1
       mole_fraction        1
       molecule             CO
       path_length          0.1 cm
       pressure_mbar        20000.0 mbar
       wavelength_max       4901.9607843137255 nm
       wavelength_min       4401.408450704225 nm
       wavenum_max          2272 cm-1
       wavenum_min          2040 cm-1
    Computation Parameters
    ----------------------------------------
       Tref                 296 K
       broadening_max_width  30 cm-1
       cutoff               1e-27 cm-1/(#.cm-2)
       db_assumed_sorted    True
       db_use_cached        False
       dbformat             hitran
       dbpath               #### USER DEPENDANT 
       fillmissinglevelswithzero  False
       levelsfmt            neq
       medium               air
       parfuncfmt           hapi
       rot_distribution     boltzmann
       self_absorption      True
       vib_distribution     boltzmann
       wavenum_max_calc     2287.0 cm-1
       wavenum_min_calc     2025.0 cm-1
       waveunit             cm-1
       wstep                0.005 cm-1
    ----------------------------------------
    

"""

from __future__ import absolute_import
from radis import SpectrumFactory, Spectrum
from radis.io.hapi import (db_begin, fetch, tableList, absorptionCoefficient_Voigt,
                                transmittanceSpectrum)
from neq.test.utils import setup_test_line_databases, printm
from os.path import join, dirname
import pytest


@pytest.mark.fast
def test_validation_case(rtol=1e-2, verbose=True, plot=False, debug=False,
                         *args, **kwargs):
    '''
    Plot absorption coefficient (cm-1) of CO at high temperature (2000 K) with 
    RADIS, and compare with calculations from HAPI using the HITRAN database 
    '''

    setup_test_line_databases()  # add HITRAN-CO-TEST in neq.rc if not there

    # Conditions
    T = 1000
    p = 10
    L = 0.1
    dnu = 0.005
    wmin = 2040  # cm-1
    wmax = 2272  # cm-1
#    broadening_max_width = 6  # cm-1
    broadening_max_width = 30  # cm-1

    # %% HITRAN calculation
    # -----------

    # Generate HAPI database locally
    db_begin(join(dirname(__file__), __file__.replace('.py', '_HAPIdata')))
    if not 'CO' in tableList():   # only if data not downloaded already
        fetch('CO', 5, 1, wmin-broadening_max_width /
              2, wmax+broadening_max_width/2)
        # HAPI doesnt correct for side effects

    # Calculate with HAPI
    nu, coef = absorptionCoefficient_Voigt(SourceTables='CO',
                                           Environment={'T': T,  # K
                                                        'p': p/1.01325,  # atm
                                                        },
                                           WavenumberStep=dnu,
                                           HITRAN_units=False,
                                            GammaL='gamma_self')
    nu, trans = transmittanceSpectrum(nu, coef,
                                      Environment={'l': L,  # cm
                                                   })
    s_hapi = Spectrum.from_array(nu, trans, 'transmittance_noslit', 'cm-1', 'I/I0',
                                 conditions={'Tgas': T},
                                 name='HAPI')

    # %% Calculate with RADIS
    # ----------
    pl = SpectrumFactory(
        wavenum_min=wmin,
        wavenum_max=wmax,
        mole_fraction=1,
        path_length=L,
        wstep=dnu,
        pressure=p,
        broadening_max_width=broadening_max_width,
        isotope=[1])  # 0.2)
    pl.warnings['MissingSelfBroadeningWarning'] = 'ignore'
    pl.load_databank('HITRAN-CO-TEST')
#    s = pl.non_eq_spectrum(Tvib=T, Trot=T, Ttrans=T)
    s = pl.eq_spectrum(Tgas=T)  # , Ttrans=300)
    s.name = 'RADIS'

    # %% Compare
    # also shrinks HAPI range to the valid one
    s_hapi.resample(s.get_wavenumber(), unit='cm-1', energy_threshold=0.1)
    if plot:
        from radis import plot_diff
        plot_diff(s, s_hapi, var='transmittance_noslit',
                  title='{0} bar, {1} K, {2} cm'.format(p, T, L))

    # %%
#    s.plot('abscoeff', lw=3, label='RADIS', medium='vacuum')
    #s_hapi.plot('abscoeff', nfig='same', label='HAPI', medium='vacuum')

    # Compare integrals
    diff = abs(s.get_integral('transmittance_noslit', medium='vacuum') /
               s_hapi.get_integral('transmittance_noslit', medium='vacuum')-1)
    b = diff < rtol

    if verbose:
        printm('Integral difference ({0:.2f}%) < {1:.2f}%: {2}'.format(diff*100,
                                                                       rtol*100, b))

    # swamp console namespace with local variables (only for debug!)
    if debug:
        globals().update(locals())

    return b


if __name__ == '__main__':
    printm('test RADIS_vs_HAPI_CO_highT: ', test_validation_case(plot=True,
                                                                 debug=False))
