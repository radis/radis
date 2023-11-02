"""
.. _example_run_from_matlab:

Example of using RADIS from Matlab
----------------------------------

All methods of RADIS can be accessed using `py.radis` followed by the method to be accessed

In the lines below, we emulate Matlab syntax to generate the example on this
online documentation, so that lines be copied directly into Matlab.

See the real Matlab file : https://github.com/radis/radis/blob/develop/examples/calc_matlab.m
See https://github.com/radis/radis/pull/547 for more details and screenshots
of Radis running in Matlab directly
"""

#%% Section not needed in Matlab

from radis.test.utils import EmulateMatlab  # this line is not needed in Matlab

py = EmulateMatlab()  # this line is not needed in Matlab

#%% Matlab example
# All lines below can be copied in Matlab directly :

s = py.radis.calc_spectrum(
    1900,
    2300,
    molecule="CO",
    isotope="1,2,3",
    pressure=1.01325,
    Tgas=700,
    mole_fraction=0.1,
    path_length=1,
    databank="hitran",
)
s.apply_slit(0.5, "nm")
