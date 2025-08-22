# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 08:12:00 2025

@author: Nicolas Minesi
"""

from radis import SpectrumFactory, plot_diff

conditions = {
    # 'wavenum_min': 500,
    # 'wavenum_max': 10000,
    # "wavenum_min": 2007,
    # "wavenum_max": 2009.5,
    # "molecule": "CO",
    "wavenum_min": 2000,
    "wavenum_max": 2400,
    "molecule": "CO2",
    "path_length": 0.1,
    "mole_fraction": 0.2,
    "isotope": 1,  # or "1,2" or "all"
    "pressure": 1e-5,
    "wstep": 0.002,
    # also measure interpolation time
    "save_memory": False,
    "truncation": 0.1,  # cm-1; Default value is 50 cm-1
}
Tgas = 700
databank = "hitemp"
database = "2010"  # "default"
#%% ############## No optimization used ##############
conditions["optimization"] = None
# # Using a convolution of a Gaussian and Lorentzian, no optimization
conditions[
    "broadening_method"
] = "convolve"  # Voigt = numpy.convolve("Gaussian", "Lorentzian")
sf = SpectrumFactory(**conditions)
sf.fetch_databank(databank)
s_conv = sf.eq_spectrum(Tgas=Tgas)

# Using a polynomial approximation, no optimization
conditions[
    "broadening_method"
] = "voigt"  # Voigt = polynomial approximation derived from Whithing
sf = SpectrumFactory(**conditions)
sf.fetch_databank(databank, database=database)
s_poly = sf.eq_spectrum(Tgas=Tgas)
s_poly.name = f"Poly. : {s_poly.c['calculation_time']:.1f}s"

# # Compare the spectra
s_conv.name = f"Convolution : {s_conv.c['calculation_time']:.1f}s"
plot_diff(s_conv, s_poly, "absorbance", yscale="log")

#%% ############## Using LDM optimization (linear/simple) ##############
conditions["optimization"] = "simple"

conditions[
    "broadening_method"
] = "convolve"  # Voigt = numpy.convolve("Gaussian", "Lorentzian")
sf = SpectrumFactory(**conditions)
sf.fetch_databank(databank)
s_conv_optim = sf.eq_spectrum(Tgas=Tgas)
s_conv_optim.name = f"Convolution : {s_conv_optim.c['calculation_time']:.1f}s"


conditions[
    "broadening_method"
] = "voigt"  # Voigt = polynomial approximation derived from Whithing
sf = SpectrumFactory(**conditions)
sf.fetch_databank(databank)
s_poly_optim = sf.eq_spectrum(Tgas=Tgas)
s_poly_optim.name = f"Polynomial : {s_poly_optim.c['calculation_time']:.1f}s"


conditions["broadening_method"] = "fft"  # FFT(Voigt) = FFT(Gaussian) * FFT(Lorentzian)
conditions["truncation"] = None
sf = SpectrumFactory(**conditions)
sf.fetch_databank(databank, database=database)
s_fft_optim = sf.eq_spectrum(Tgas=Tgas)
s_fft_optim.name = f"FFT - linear optim.: {s_fft_optim.c['calculation_time']:.1f}s"


# Compare the spectra
plot_diff(s_conv_optim, s_poly_optim, "absorbance", yscale="log")
plot_diff(s_conv_optim, s_fft_optim, "absorbance", yscale="log")
plot_diff(s_poly, s_fft_optim, "absorbance", yscale="linear")
plot_diff(s_poly_optim, s_fft_optim, "absorbance", yscale="linear")
