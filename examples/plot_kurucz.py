from radis.lbl.factory import SpectrumFactory
import matplotlib.pyplot as plt
from radis.io.kurucz import fetch_kurucz
from radis.api.kuruczapi import AdBKurucz

import pandas as pd
import os
import h5py
import numpy as np




hdf5_path = fetch_kurucz('Fe',00)[0] # replace 'Fe', 00 with your atom and ionization_state


# Define the parameters for SpectrumFactory
wavenum_min = 4000 # adjust as needed
wavenum_max = 5000 # adjust as needed
atom = 'Fe'  # adjust as needed
ionization_state = '00'  # adjust as needed

# Create a spectrum factory
sf = SpectrumFactory(wavenum_min, wavenum_max, atom=atom, ionization_state=ionization_state, wstep="auto")
sf.load_databank(path=hdf5_path, format='kurucz')

# Calculate the spectrum
spectrum = sf.eq_spectrum(Tgas=5000)

# Plot the spectrum
spectrum.plot('radiance_noslit', wunit='nm')

plt.show()
