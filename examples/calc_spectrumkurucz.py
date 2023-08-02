from radis.lbl.factory import SpectrumFactory
import matplotlib.pyplot as plt
from radis.api.kuruczapi import AdBKurucz

import pandas as pd
import os
import periodictable
import h5py
import numpy as np

#Move to api.kuruczapi.py ?
def get_kurucz_data(atom, ionization_state):
    kurucz = AdBKurucz(atom,ionization_state)
    atomic_number = getattr(periodictable, atom).number
    ionization_state_str = str(ionization_state).zfill(2)
    kurucz_file = f"gf{atomic_number}{ionization_state_str}.all"
    hdf5_file = f"gf{atomic_number}{ionization_state_str}.hdf5"
    kurucz.url = kurucz.get_url(atomic_number, ionization_state)
    kurucz.hdf5_file = hdf5_file  # Set kurucz.hdf5_file to hdf5_file

    # If hdf5 file exists, read data from it
    if os.path.exists(hdf5_file):
        print("HDF5 file already exists, reading data from it.")
        df = kurucz.read_hdf5(hdf5_file)
    else:
        kuruczf = kurucz.download_file()
        df = kurucz.read_kurucz(kuruczf)
        #print(df)
        print(kurucz.hdf5_file)
        kurucz.store_hdf5(df, kurucz.hdf5_file)
        df = kurucz.read_hdf5(hdf5_file)

    return hdf5_file



hdf5_path = get_kurucz_data('Fe',00) # replace 'Fe', 00 with your atom and ionization_state


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
