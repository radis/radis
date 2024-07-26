# -*- coding: utf-8 -*-
"""
.. _example_custom_lorentzian_broadening:

========================
Using a custom broadening function for atomic spectra
========================

The `lbfunc` parameter of :py:class:`~radis.lbl.factory.SpectrumFactory` allows us to specify a custom function to use for the Lorentzian broadening of a spectrum. This is especially useful for handling pressure broadening of atomic lines for which multiple theories exist.

By default, RADIS calculates the total Lorentzian broadening of atomic lines as the sum of the Van der Waals broadening, Stark broadening, and radiation broadening as returned by :py:func:`~radis.lbl.broadening.gamma_vald3` from the broadening parameters for each of these types provided in the databank. The line shift is also calculated as 2/3 times the Van der Waals broadening, and `lbfunc` allows you to specify an alternative line shift as the 2nd return argument if you wish.

`lbfunc` can be changed on the fly by changing the `sf.params.lbfunc` attribute of the :py:class:`~radis.lbl.factory.SpectrumFactory` instance, the result of which is reflected next time the library calculates its Lorentzian broadening.

To make use of it, familiarise yourself with the column names that RADIS assigns to the relevant quantities you intend to use in `lbfunc`.

We take an example of neutral atomic oxygen (O_I) and calculate a spectrum using the broadening formula (and default values) of SpectraPlot (https://spectraplot.com/), and assume there is no line shift.
"""
from radis import SpectrumFactory

def lbfunc1(**kwargs):
    return 0.1*(296/kwargs['Tgas'])**0.8, None

mole_fraction = 0.01

sf = SpectrumFactory(
    12850,
    12870,
    species="O_I",
    pressure=1.01325, # = 1 atm
    diluent={'H':1-mole_fraction-1e-3, 'e-': 1e-3}, #so it all adds up to 1
    mole_fraction=mole_fraction,
    path_length=15,
    lbfunc=lbfunc1
)
sf.fetch_databank('kurucz', parfuncfmt="kurucz")
s1 = sf.eq_spectrum(4000)

#%%
# Now compare the result with that of the default handling of broadening in RADIS by removing the `lbfunc` parameter, recalculating the spectrum without it, and plotting the diff between the results:
#
sf.params.lbfunc = None
s_default = sf.eq_spectrum(4000)

from radis import plot_diff

plot_diff(s1, s_default, label1='s1', label2='s_default')

#%%
# Here's another example adding self broadening, calculated using eq(16) of [Minesi et al 2020]_, to the 3 broadening types already handled by RADIS, and keeping the default handling of the line shift:
#

from radis.phys.constants import k_b_CGS
from radis.lbl.broadening import gamma_vald3

def lbfunc2(df, pressure_atm, mole_fraction, Tgas, diluent, isneutral, **kwargs): # assign variable names to the quantities we use and put the rest in kwargs which can be ignored
    beta = 1e-4 # example
    n_emitter = pressure_atm*1.01325*1e6*mole_fraction / (k_b_CGS*Tgas)
    gamma_self = ((1e7/df['wav'])**2)*beta*n_emitter/2.7e19
    print(gamma_self)
    # copied from default code in RADIS:
    gammma_rad, gamma_stark, gamma_vdw = gamma_vald3(Tgas, pressure_atm*1.01325, 
    df['wav'], df['El'], df['ionE'], df['gamRad'], df['gamSta'], df['gamvdW'], diluent, isneutral)
    print(gammma_rad, gamma_stark, gamma_vdw)
    shift = (1.0/3.0)*2*gamma_vdw #Konjević et al. 2012 §4.1.3.2, neglect stark shift by default
    wl = gammma_rad + gamma_stark + gamma_vdw + gamma_self
    return wl, shift

sf.params.lbfunc = lbfunc2
s2 = sf.eq_spectrum(4000)

plot_diff(s2, s_default, label1='s2', label2='s_default')

#%%
# We can even modify broadening parameters of individual lines:
#

def lbfunc3(df, Tgas, pressure_atm, diluent, isneutral, **kwargs):
    # only for Pandas dataframes:
    df.loc[df['orig_wavelen'] == 777.5388, 'gamvdW'] = -7 # should also result in a different shift
    df.loc[df['orig_wavelen'] == 777.1944, 'gamSta'] = -5
    # copied from default code in RADIS:
    gammma_rad, gamma_stark, gamma_vdw = gamma_vald3(Tgas, pressure_atm*1.01325, df['wav'], df['El'], df['ionE'], df['gamRad'], df['gamSta'], df['gamvdW'], diluent, isneutral)
    shift = (1.0/3.0)*2*gamma_vdw #Konjević et al. 2012 §4.1.3.2
    wl = gammma_rad + gamma_stark + gamma_vdw
    return wl, shift

sf.params.lbfunc = lbfunc3
s3 = sf.eq_spectrum(4000)

plot_diff(s3, s_default, label1='s3', label2='s_default')

#%%
# References
# ----------
# [Minesi et al 2020] `"Fully ionized nanosecond discharges in air: the thermal spark" <https://ui.adsabs.harvard.edu/abs/2020PSST...29h5003M>`_