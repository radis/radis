import matplotlib.pyplot as plt

from radis.io.kurucz import fetch_kurucz
from radis.lbl.factory import SpectrumFactory

hdf5_path = fetch_kurucz("Fe_I")[0]  # replace 'Fe_I' with atom and ionization_state


# Define the parameters for SpectrumFactory
wavenum_min = 4000  # adjust as needed
wavenum_max = 5000  # adjust as needed
species = "Fe_I"

# Create a spectrum factory
sf = SpectrumFactory(wavenum_min, wavenum_max, species="Fe_I", wstep="auto")
sf.load_databank(path=hdf5_path, format="kurucz", parfuncfmt="kurucz")

# Calculate the spectrum
spectrum = sf.eq_spectrum(Tgas=5000)

# Plot the spectrum
spectrum.plot("radiance_noslit", wunit="nm")

plt.show()
