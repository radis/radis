# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 08:12:00 2025

@author: Nicolas Minesi
"""

from radis import SpectrumFactory, plot_diff

conditions = {
    # 'wavenum_min': 500,
    # 'wavenum_max': 10000,
    "wavenum_min": 2007,
    "wavenum_max": 2009.5,
    "path_length": 0.1,
    "molecule": "CO",
    "mole_fraction": 0.2,
    "isotope": 1,  # or "1,2" or "all"
    "pressure": 1e-5,
    "wstep": 0.002,
    # also measure interpolation time
    "save_memory": False,
    "truncation": 50,  # cm-1; Default value is 50 cm-1
}
Tgas = 700

#%% ############## No optimization used ##############
# %% Using a convolution of a Gaussian and Lorentzian, no optimization
sf = SpectrumFactory(
    **conditions,
    broadening_method="convolve",  # Voigt = numpy.convolve("Gaussian", "Lorentzian")
    optimization=None,
)
sf.fetch_databank("hitran")
s_conv = sf.eq_spectrum(Tgas=Tgas)


# %% Using a polynomial approximation, no optimization
sf = SpectrumFactory(
    **conditions,
    broadening_method="voigt",  # Voigt = polynomial approximation derived from Whithing
    optimization=None,
)
sf.fetch_databank("hitran")
s_poly = sf.eq_spectrum(Tgas=Tgas)

# Compare the spectra
s_conv.name = f"Convolution : {s_conv.c['calculation_time']:.1f}s"
s_poly.name = f"Polynomial : {s_poly.c['calculation_time']:.1f}s"
plot_diff(s_conv, s_poly, "absorbance", yscale="log")

#%% ############## Using LDM optimization (linear/simple) ##############
# # %% Using FFT to compute Voigt
# sf = SpectrumFactory(
#     **conditions,
#     broadening_method = "fft", #Voigt = polynomial approximation derived from Whithing
#     optimization = "simple",
#     )
# sf.fetch_databank("hitran")
# s_fft = sf.eq_spectrum(Tgas=Tgas)
