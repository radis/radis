# -*- coding: utf-8 -*-
"""
.. _example_potential_lowering_pfs:

========================
Kurucz partition functions and potential lowering
========================

For atomic species in the Kurucz databank, some species come with dedicated partition function tables provided, which depend on both temperature and potential lowering. The `potential_lowering` parameter of :py:class:`~radis.lbl.factory.SpectrumFactory` provides a way for the user to set this value to make use of these tables in the interpolator.

When such tables are unavailable or `potential_lowering` is `None` (the default), the partition functions are instead interpolated from [Barklem-&-Collet-2016]_ Table 8.

The value can thereafter be changed on the fly by changing the `sf.input.potential_lowering` attribute of the :py:class:`~radis.lbl.factory.SpectrumFactory` instance ``sf``, the result of which is reflected the next time partition function interpolator's `._at` method is used, without any need to re-initialise it.

Allowable values are usually (in cm-1/Zeff**2): -500, -1000, -2000, -4000, -8000, -16000, -32000. This is based on what's been encountered in the partition function tables of the species checked so far, but documentation of the Kurucz linelists is unclear on whether this is the case for all species.

The temperature ranges of the partition functions from [Barklem-&-Collet-2016]_ may differ from those of the dedicated Kurucz tables, so where you can calculate spectra at a particular temperature with the Kurucz tables:
"""
#%%

from radis import calc_spectrum, SpectrumFactory
import traceback

s, sf = calc_spectrum(
    205,
    263,
    wunit='nm',
    species="Y_I",
    Tgas=50000,
    databank="kurucz",
    potential_lowering=-500, # any of [-500, -1000, -2000, -4000, -8000, -16000, -32000]
    return_factory=True,
    save_memory=False, #to be able to recalculate spectra with the same SpectrumFactory
)

s.plot("radiance_noslit", wunit="cm-1")

#%%
# you encounter an error at the same temperature with [Barklem-&-Collet-2016]_:
#
sf.input.potential_lowering = None

try:
    sf.eq_spectrum(50000)
except ValueError:
    print(traceback.format_exc())

#%%
# and where you can calculate spectra at a particular temperature with [Barklem-&-Collet-2016]_:
#

s = sf.eq_spectrum(Tgas=99)
s.plot("radiance_noslit", wunit="cm-1")

#%%
# you encounter an error at the same temperature with the Kurucz tables:
#
sf.input.potential_lowering = -500

try:
    sf.eq_spectrum(Tgas=99)
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
# Specifying the potential_lowering where the tables aren't available has no effect, and the partition functions from [Barklem-&-Collet-2016]_ get used instead with just the temperature:
#

s = calc_spectrum(
    39000,
    40000,
    species="Os_II",
    Tgas=4000,
    databank="kurucz",
    potential_lowering=-1100
)

s.plot("radiance_noslit", wunit="cm-1")

#%%
# Occasionally, the opposite is true, and the partition functions from [Barklem-&-Collet-2016]_ don't include a species:
#

try:
    calc_spectrum(
        100,
        10000,
        species="Y_V",
        Tgas=4000,
        databank="kurucz",
        cutoff=0
    )
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
    potential_lowering=-1000
)

s.plot("radiance_noslit", wunit="cm-1")

#%%
# References
# ----------
# .. [Barklem-&-Collet-2016] `"Partition functions and equilibrium constants for diatomic molecules and atoms of astrophysical interest" <https://ui.adsabs.harvard.edu/abs/2016A%2526A...588A..96B>`_