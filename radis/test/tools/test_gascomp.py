# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 12:55:28 2021

@author: erwan

Test the :py:func:`~radis.tools.gascomp.get_eq_mole_fraction` function.
"""

### As of version 0.16, cantera is not required to run RADIS.
### Therefore, the test is commented out.

# def test_get_eq_mole_fraction(*args, **kwargs):
#     from radis.misc.basics import all_in
#     from radis.tools.gascomp import get_eq_mole_fraction

#     gas = get_eq_mole_fraction("CO2:1", 3000, 1e5)

#     assert all_in(["C", "CO", "CO2", "O", "O2"], gas.keys())

# if __name__ == "__main__":
#     test_get_eq_mole_fraction()
