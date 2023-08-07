""""

Summary
-----------------------

Kurucz database parser

-----------------------

Defines :py:func:`~radis.io.fetch_kurucz` based on :py:class:`~AdBKurucz`

"""

from radis.api.kuruczapi import AdBKurucz
import mendeleev
import os




def fetch_kurucz(atom, ionization_state):
    kurucz = AdBKurucz(atom,ionization_state)
    atomic_number = getattr(mendeleev, atom).atomic_number
    ionization_state_str = str(ionization_state).zfill(2)
    kurucz_file = f"gf{atomic_number}{ionization_state_str}.all"
    hdf5_file = f"gf{atomic_number}{ionization_state_str}.hdf5"
    kurucz.url = kurucz.get_url(atomic_number, ionization_state)
    kurucz.hdf5_file = hdf5_file  # Set kurucz.hdf5_file to hdf5_file

    # If hdf5 file exists, read data from it
    if os.path.exists(hdf5_file):
        print("HDF5 file already exists, reading data from it.")
        df = kurucz.read_hdf5(hdf5_file)
        kurucz.add_airbrd(df)
    else:
        kuruczf = kurucz.download_file()
        df = kurucz.read_kurucz(kuruczf)
        #print(df)
        print(kurucz.hdf5_file)
        kurucz.store_hdf5(df, kurucz.hdf5_file)
        df = kurucz.read_hdf5(hdf5_file)
        kurucz.add_airbrd(df)

    return hdf5_file,df
