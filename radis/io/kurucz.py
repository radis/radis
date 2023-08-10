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




def fetch_kurucz(species):
    kurucz = AdBKurucz(species)
    atomic_number =kurucz.get_atomic_number(species)
    ionization_state_str = kurucz.get_ionization_state(species)
    kurucz_file = f"gf{atomic_number}{ionization_state_str}.all"
    hdf5_file = f"gf{atomic_number}{ionization_state_str}.hdf5"
    kurucz.url = kurucz.get_url(atomic_number, ionization_state_str)
    kurucz.hdf5_file = hdf5_file  # Set kurucz.hdf5_file to hdf5_file

    # If hdf5 file exists, read data from it
    if os.path.exists(hdf5_file):
        print("HDF5 file already exists, reading data from it.")
        df = kurucz.read_hdf5(hdf5_file)
        kurucz.add_airbrd(df)
    else:
        kuruczf = kurucz.download_file() 
        df=kurucz.read_kurucz(kuruczf)
        #print(df)
        kurucz.store_hdf5(df, kurucz.hdf5_file)
        df = kurucz.read_hdf5(hdf5_file)
        kurucz.add_airbrd(df)

    return hdf5_file,df
