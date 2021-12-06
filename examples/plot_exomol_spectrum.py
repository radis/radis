# -*- coding: utf-8 -*-
"""
================================
Calculate a spectrum from ExoMol
================================

Auto-download and compute a SiO spectrum from the ExoMol database ([ExoMol-2020]_)

ExoMol lines can be downloaded and accessed separately using
:py:func:`~radis.io.exomol.fetch_exomol`

"""

from radis import calc_spectrum

s = calc_spectrum(
    1080,
    1320,  # cm-1
    molecule="SiO",
    isotope="1",
    pressure=1.01325,  # bar
    Tgas=1000,  # K
    mole_fraction=0.1,
    path_length=1,  # cm
    broadening_method="fft",  # @ dev: Doesn't work with 'voigt'
    databank=("exomol", "EBJT"),  # Simply use 'exomol' for the recommended database
)
s.apply_slit(1, "cm-1")  # simulate an experimental slit
s.plot("radiance")


#%% See line data:
from radis.io.exomol import fetch_exomol

df = fetch_exomol("SiO", database="EBJT", isotope="1", load_wavenum_max=5000)
print(df)


#%% See the list of recommended databases for the 1st isotope of SiO :
from radis.io.exomol import get_exomol_database_list, get_exomol_full_isotope_name

databases, recommended = get_exomol_database_list(
    "SiO", get_exomol_full_isotope_name("SiO", 1)
)
print("Databases for SiO: ", databases)
print("Database recommended by ExoMol: ", recommended)
