"""
================================================
Spectrum Database
================================================
RADIS has :py:class:`~radis.tools.database.SpecDatabase` feature used to store and retrieve calculated Spectrums. A path can be specified for SpecDatabase all Spectrums are stored as .spec files which can be loaded
from the SpecDatabase object itself.  A csv file is generated which contains all input and conditional parameters of Spectrum.

RADIS also has :py:meth:`~radis.lbl.loader.DatabankLoader.init_database` feature which initializes the SpecDatabase for the SpectrumFactory and every Spectrum
generated from it will be stored in the SpecDatabase automatically.

You can use :py:meth:`~radis.tools.database.SpecList.plot_cond` to make a 2D plot using the conditions of the Spectrums in the SpecDatabase and use a z_label to plot a heat map based on it.


"""

import numpy as np

from radis import SpectrumFactory
from radis.tools import SpecDatabase

sf = SpectrumFactory(
    wavenum_min=2900,
    wavenum_max=3200,
    molecule="OH",
    broadening_max_width=10,  # cm-1
    medium="vacuum",
    verbose=0,  # more for more details
    pressure=10,
    wstep="auto",
)
sf.fetch_databank("hitemp")

# Generating 3 Spectrums
s1 = sf.eq_spectrum(Tgas=300, path_length=1)
s2 = sf.eq_spectrum(Tgas=450, path_length=1)

# Creating SpecDatabase
my_folder = r"/home/pipebomb/Desktop/SpecDatabase_Test/"
db = SpecDatabase(my_folder, lazy_loading=False)

# Method 1: Creating .spec file
db.add(s1)
db.add(s2)

# Method 2: Using init_database()
sf.init_database(my_folder)

# Generates Spectrum and adds to SpecDatabase automatically
for i in np.arange(600, 3100, 150):
    sf.eq_spectrum(Tgas=i, path_length=1)


# Loading SpecDatabase
db_new = SpecDatabase(my_folder)

# Comparing data conditions of different spectrum from csv generated file
db_new.plot_cond("Tgas", "wstep", "calculation_time")
