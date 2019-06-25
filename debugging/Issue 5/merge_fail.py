# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 23:45:17 2018

@author: erwan

See Issue #22
"""

from radis import load_spec, MergeSlabs
from radis.misc.arrays import count_nans

s1 = load_spec('cruden_2611K_1.spec')
s2 = load_spec('cruden_2611K_2.spec')
assert count_nans(s1._q['emissivity_noslit']) == 0
assert count_nans(s2._q['emissivity_noslit']) == 0

#import radis
#radis.DEBUG_MODE = True

s = MergeSlabs(s1, s2, resample_wavespace='full', out_of_bounds='transparent')

assert count_nans(s._q['emissivity_noslit']) == 0



#%% 

s1 = load_spec('cruden_2611K_1.spec')
s2 = load_spec('cruden_2611K_2.spec')

s1.conditions['thermal_equilibrium'] = False
assert count_nans(s1._q['emissivity_noslit']) == 0
assert count_nans(s2._q['emissivity_noslit']) == 0

s = MergeSlabs(s1, s2, resample_wavespace='full', out_of_bounds='transparent')

assert 'emissivity_noslit' not in s.get_quantities()
