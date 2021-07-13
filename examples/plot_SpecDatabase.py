"""
================================================
Spectrum Database
================================================
RADIS has SpecDatabase feature used to store and retrieve calculated Spectrums. A path can be specified for SpecDatabase all Spectrums are stored as .spec files which can be loaded
from the SpecDatabase object itself.  A csv file is generated which contains all input and conditional parameters of Spectrum.

:py:class:`~radis.tools.database.SpecDatabase`
"""

from radis import SpectrumFactory
from radis.tools import SpecDatabase

sf = SpectrumFactory(
    wavenum_min=2384,
    wavenum_max=2400,
    molecule="CO2",
    isotope="1,2",
    broadening_max_width=10,  # cm-1
    medium="vacuum",
    verbose=0,  # more for more details
    wstep="auto",
)
sf.fetch_databank()

# Generating 3 Spectrums
s1 = sf.eq_spectrum(name="Spectum_CO2_400", Tgas=400, path_length=1)
s2 = sf.eq_spectrum(name="Spectum_CO2_450", Tgas=450, path_length=1)
s3 = sf.eq_spectrum(name="Spectum_CO2_500", Tgas=500, path_length=1)

# Creating SpecDatabase
my_folder = r"/home/pipebomb/Desktop/SpecDatabase_Test/"
db = SpecDatabase(my_folder, lazy_loading=False)

# Method 1: Creating .spec file
db.add(s1)
db.add(s2)

# Method 2: Creating .spec file manually
# Note: Doesn't get added to SpecDatabase csv file
s3.store(my_folder + s3.name)

# Loading SpecDatabase
"""
Note: If number of spec files in SpecDatabase isn't equal to number of Spectrums
in the csv files, use  lazy_loading=False.
"""
db_new = SpecDatabase(my_folder, lazy_loading=False)

# Loading all spec files in a list
list_Spectrum = []
for s in db_new:
    list_Spectrum.append(s)

# Generating Plot
for spec in list_Spectrum:
    spec.plot("radiance_noslit", nfig="same")
