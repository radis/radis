# -*- coding: utf-8 -*-
"""

=========================================
Calculate a High Temperature H20 Spectrum
=========================================

This example uses the calc_spectrum function to calculate the H20 spectra at 3000K using the 
HITEMP databank. Wavenumbers are set between 2000 and 4000.


"""

#%%

from radis import calc_spectrum
from astropy import units as u
s = calc_spectrum(wavenum_min = 2000 / u.cm, 
                  wavenum_max = 4000 / u.cm,
                  molecule = 'H2O',
                  isotope = '1,2,3',
                  pressure = 1.01325 * u.bar,
                  mole_fraction = 1,
                  path_length = 1 * u.cm,
                  Tgas=3000, 
                  databank='hitemp',
                  verbose=3,        # adding some extra info for the first computation
                  )
#Without verbose=False this will show all the input parameters. 
#With verbose=2,3,etc... we get increasing number of details about the calculation. 

#%% 
#We then plot the results:
#
import matplotlib.pyplot as plt 

s.plot('transmittance_noslit', wunit='nm')


