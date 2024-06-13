""""

Summary
-----------------------

Kurucz database parser

-----------------------

Defines :py:func:`~radis.io.fetch_kurucz` based on :py:class:`~AdBKurucz`

"""

import os

from radis.api.hdf5 import DataFileManager
from radis.api.kuruczapi import AdBKurucz,get_atomic_number,get_ionization_state


def fetch_kurucz(species, isotope, engine):
    kurucz = AdBKurucz(species)
    data_manager = DataFileManager(engine=engine)
    atomic_number = f"{get_atomic_number(species):02}"
    ionization_state_str = f"{get_ionization_state(species):02}"
    hdf5_file = f"gf{atomic_number}{ionization_state_str}.hdf5"
    kurucz.url = kurucz.get_url(atomic_number, ionization_state_str)
    kurucz.hdf5_file = hdf5_file  # Set kurucz.hdf5_file to hdf5_file

    # If hdf5 file exists, read data from it
    if os.path.exists(hdf5_file):
        print("HDF5 file already exists, reading data from it.")
    else:
        kuruczf = kurucz.download_file()
        df = kurucz.read_kurucz(kuruczf)
        data_manager.write(
            file=hdf5_file,
            df=df,
            append=False,
            format="table",
            data_columns=df.columns,
        )
    df = data_manager.load(fname=hdf5_file, within=[("iso", isotope)] if isotope is not None else [])
    kurucz.add_airbrd(df)

    return hdf5_file, df
