# -*- coding: utf-8 -*-
"""
.. _example_potential_lowering_pfs:

========================
Atomic spectrum #1: Atomic partition functions
========================

The architecture allows multiple sources of partition functions each of which can be used with any atomic databank. The source is specified by the `pfsource` parameter of :py:class:`~radis.lbl.factory.SpectrumFactory`. There are currently 3 sources implemented:
    - 'barklem': an interpolator that uses a table of partition functions from [Barklem-&-Collet-2016]
    - 'kurucz' an interpolator that uses a table of partition functions provided with some Kurucz linelists
    - 'nist' a calculator that calculates the partition function from a table of energy levels provided by NIST
The partition functions from all 3 sources depend on temperature, but those from 'kurucz' also depend on potential lowering, which can be set with the `potential_lowering` parameter of :py:class:`~radis.lbl.factory.SpectrumFactory`.

The default `pfsource` is 'nist'. It can be changed on the fly using the :meth:`~radis.lbl.loader.DatabankLoader.set_atomic_partition_functions` method, and the potential lowering for 'kurucz' can be modified on the fly by changing the `sf.input.potential_lowering` attribute of the :py:class:`~radis.lbl.factory.SpectrumFactory` instance ``sf``. The changes are reflected the next time a partition function is calculated.

Allowable values for `potential_lowering` are usually (in cm-1/Zeff**2): -500, -1000, -2000, -4000, -8000, -16000, -32000.
"""
#%%
import traceback

from radis import calc_spectrum

T_high = 50000  # K
s, sf = calc_spectrum(
    205,
    263,
    wunit="nm",
    species="Y_I",
    Tgas=T_high,
    databank="kurucz",
    pfsource="nist",
    return_factory=True,
    save_memory=False,  # to be able to recalculate spectra with the same SpectrumFactory
)

s.plot("radiance_noslit", wunit="cm-1")

#%%
# We can change the `pfsource` to 'kurucz':
#
sf.set_atomic_partition_functions("kurucz")

#%%
# but if we run anything that attempts to calculate a partition function, we get an error:
#

try:
    sf.eq_spectrum(T_high)
except Exception:
    print(traceback.format_exc())

#%%
# until we specify the potential lowering:
#

sf.input.potential_lowering = (
    -500
)  # any of [-500, -1000, -2000, -4000, -8000, -16000, -32000]
s2 = sf.eq_spectrum(T_high)
s2.plot("radiance_noslit", wunit="cm-1")

#%%
# In this case, 50 000 K is beyond the maximal temperature in 'barklem'.
#
sf.set_atomic_partition_functions("barklem")

try:
    sf.eq_spectrum(T_high)
except ValueError:
    print(traceback.format_exc())
#%%
# In this case, 99 K is within the temperature range of 'barklem':
#
T_low = 99  # K
s2 = sf.eq_spectrum(Tgas=T_low)
s2.plot("radiance_noslit", wunit="cm-1")

#%%
# However, 99 K is *not* within the temperature range of 'kurucz':
#
sf.set_atomic_partition_functions(
    "kurucz"
)  # potential lowering is still -500 from above

try:
    sf.eq_spectrum(Tgas=T_low)
except ValueError:
    print(traceback.format_exc())

#%%
# Specifying a value for the potential lowering that isn't present in the tables just returns a KeyError when attempting to calculate:
#

sf.input.potential_lowering = -1100
try:
    sf.eq_spectrum(4000)
except KeyError:
    print(traceback.format_exc())

#%%
# Partition functions could also be available for a species from one source but not another, e.g. available from 'barklem':
#

s, sf = calc_spectrum(
    39000,
    40000,
    species="Os_II",
    Tgas=4000,
    databank="kurucz",
    pfsource="barklem",
    return_factory=True,
    save_memory=False,
)

s.plot("radiance_noslit", wunit="cm-1")

#%%
# but not from 'kurucz', resulting in a warning when you try to set it as the `pfsource`:
#

sf.set_atomic_partition_functions(
    "kurucz"
)  # setting `potential_lowering` is irrelevant as the tables aren't even available

#%%
# and an error if we attempt to calculate:
#

try:
    sf.eq_spectrum(4000)
except Exception:
    print(traceback.format_exc())

#%%
# An example of the opposite, where 'kurucz' does include a species:
#


def lbfunc1(
    **kwargs,
):  # an arbitrary broadening formula as NIST databank requires `lbfunc`
    return 0.1 * (296 / kwargs["Tgas"]) ** 0.8, None


s, sf = calc_spectrum(
    40000,
    50000,
    species="Y_V",
    Tgas=4000,
    databank="nist",
    lbfunc=lbfunc1,
    cutoff=0,
    pfsource="kurucz",
    potential_lowering=-1000,
    return_factory=True,
    save_memory=False,
)

s.plot("radiance_noslit", wunit="cm-1")

#%%
# but 'barklem' doesn't:
#
sf.set_atomic_partition_functions("barklem")

try:
    sf.eq_spectrum(4000)
except ValueError:
    print(traceback.format_exc())

#%%
# References
# ----------
# .. [Barklem-&-Collet-2016] `"Partition functions and equilibrium constants for diatomic molecules and atoms of astrophysical interest" <https://ui.adsabs.harvard.edu/abs/2016A%2526A...588A..96B>`_
