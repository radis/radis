# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 21:58:47 2017

@author: erwan

Formula to calculate most spectral quantities

Are stored here all equations that need to be used by both the factory (to 
generate the spectrum) or the spectrum class (to recompute some quantities
in post processing)

"""

from __future__ import absolute_import
from radis.phys.convert import cm2nm
from radis.phys.blackbody import planck

# ------- Get radiance from emissivity

def calc_radiance(wavenumber, emissivity, Tgas, unit='mW/sr/cm2/nm'):
    ''' Derive radiance (mW/cm2/sr/nm) from the emissivity


    Returns
    -------

    radiance: mW/sr/cm2/nm
    '''

    radiance = emissivity * planck(cm2nm(wavenumber), Tgas, unit=unit)

    return radiance

