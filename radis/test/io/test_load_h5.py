# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 20:36:38 2021

@author: erwan
"""

from astropy import units as u

from radis import calc_spectrum
from radis.test.utils import getTestFile


def test_load_h5(*args, **kwargs):

    calc_spectrum(
        wavenum_min=2245 / u.cm,
        wavenum_max=2255 / u.cm,
        molecule="CO2",
        Tgas=600,
        databank=getTestFile("cdsd_hitemp_09.h5"),
        verbose=3,
    )

    return
