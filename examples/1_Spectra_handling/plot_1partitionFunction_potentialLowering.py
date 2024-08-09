# -*- coding: utf-8 -*-
"""
.. _example_potential_lowering_pfs:

========================
Atomic spectrum #1: Partition functions and potential lowering
========================

Two sources can be used for the partition function: Kurucz and [Barklem-\&-Collet-2016]_.
These sources have different temperature ranges.
The Kurucz partition function includes a dependence on potential lowering. The
latter can be set with the `potential_lowering` parameter of :py:class:`~radis.lbl.factory.SpectrumFactory`.

By default, the [Barklem-\&-Collet-2016]_ partition functions are employed. The
Kurucz partition functions can be used instead by setting the potential lowering
to a non-zero value. The `sf.input.potential_lowering` attribute of the
:py:class:`~radis.lbl.factory.SpectrumFactory` instance ``sf`` can be changed
on the fly to modify the potential lowering value. The change will be reflected
the next time the partition function interpolator's `._at` method is used,
without re-initialization.

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
    potential_lowering=-500,  # any of [-500, -1000, -2000, -4000, -8000, -16000, -32000]
    return_factory=True,
    save_memory=False,  # to be able to recalculate spectra with the same SpectrumFactory
)

s.plot("radiance_noslit", wunit="cm-1")

#%%
# Setting potential_lowering = None automatically sets [Barklem-\&-Collet-2016]_ as the reference partition function.
# In this case, 50 000 K is beyond the maximal temperature in [Barklem-\&-Collet-2016]_.
#
sf.input.potential_lowering = None

try:
    sf.eq_spectrum(T_high)
except ValueError:
    print(traceback.format_exc())
#%%
# In this case, 99 K is within the temperature range of [Barklem-\&-Collet-2016]_:
#
T_low = 99  # K
s2 = sf.eq_spectrum(Tgas=T_low)
s2.plot("radiance_noslit", wunit="cm-1")

#%%
# However, 99 K is *not* within the temperature range of Kurucz:
#
sf.input.potential_lowering = -500  # to use the Kurucz partition function

try:
    sf.eq_spectrum(Tgas=T_low)
except ValueError:
    print(traceback.format_exc())

#%%
# Specifying a value for the potential lowering that isn't present in the tables just returns a KeyError:
#

sf.input.potential_lowering = -1100
try:
    sf.eq_spectrum(4000)
except KeyError:
    print(traceback.format_exc())

#%%
# Specifying the potential_lowering where the tables aren't available has no effect, and the partition functions from [Barklem-\&-Collet-2016]_ get used instead with just the temperature:
#

s = calc_spectrum(
    39000,
    40000,
    species="Os_II",
    Tgas=4000,
    databank="kurucz",
    potential_lowering=-1100,
)

s.plot("radiance_noslit", wunit="cm-1")

#%%
# Occasionally, the opposite is true, and the partition functions from [Barklem-\&-Collet-2016]_ don't include a species:
#

try:
    calc_spectrum(100, 10000, species="Y_V", Tgas=4000, databank="kurucz", cutoff=0)
except Exception:
    print(traceback.format_exc())

#%%
# but it does have dedicated partition function tables with the Kurucz linelists:
#

s = calc_spectrum(
    100,
    10000,
    species="Y_V",
    Tgas=4000,
    databank="kurucz",
    cutoff=0,
    potential_lowering=-1000,
)

s.plot("radiance_noslit", wunit="cm-1")

#%%
# References
# ----------
# .. [Barklem-\&-Collet-2016] `"Partition functions and equilibrium constants for diatomic molecules and atoms of astrophysical interest" <https://ui.adsabs.harvard.edu/abs/2016A%2526A...588A..96B>`_
