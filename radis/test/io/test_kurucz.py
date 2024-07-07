import pandas as pd
import numpy as np
import pytest

from radis.api.kuruczapi import load_pf_Barklem2016, read_kurucz
from radis.io.kurucz import fetch_kurucz
from radis.test.utils import getTestFile

def test_read_kurucz():
    df = read_kurucz(getTestFile('gf4000.all')) # from: http://kurucz.harvard.edu/atoms/4000/gf4000.all
    a = df.loc[5, ['orig_wavelen', 'loggf', 'species', 'El', 'jlower', 'labellower', 'Eu', 'ju', 'labelupper', 'gamRad', 'gamSta', 'gamvdW', 'iso']]#.to_numpy()
    b = np.array([504403.3610,-7.259,' 40.00', 35456.250, 2.0, 'd3p w3P', 35476.070,3.0, '5s4F6s e5F', 8.41,-5.27,-7.53, 0], dtype='object')
    assert (a == b).all()
    a = df.loc[7, ['orig_wavelen', 'loggf', 'species', 'El', 'jlower', 'labellower', 'Eu', 'ju', 'labelupper', 'gamRad', 'gamSta', 'gamvdW', 'iso']]#.to_numpy()
    b = np.array([328425.5787 ,-6.296,' 40.00',   21943.740,  2.0, 'd4 a5D',       21974.180,  1.0, 'd2sp z3S',    7.22, -5.99, -7.79, 0], dtype='object')
    assert (a == b).all()

# def test_fetch_kurucz():
#     species = "Na_I"
#     hdf5_path, df = fetch_kurucz(species)

#     required_columns = [
#         "A",
#         "wav",
#         "El",
#         "eupper",
#         "gu",
#         "jlower",
#         "ju",
#         # "id",
#         "iso",
#         "gamRad",
#         "gamSta",
#         "gamvdW",
#         # "Tdpair",
#         "wlnmair",
#         "loggf",
#         "species",
#         "labellower",
#         "labelupper",
#         # "ref",
#         # "NLTElower",
#         # "NLTEupper",
#         # "isonum",
#         # "hyperfrac",
#         # "isonumdi",
#         # "isofrac",
#         # "hypershiftlower",
#         # "hypershiftupper",
#         # "hyperFlower",
#         # "hypernotelower",
#         # "hyperFupper",
#         # "hypternoteupper",
#         # "strenclass",
#         # "auto",
#         # "landeglower",
#         # "landegupper",
#         # "isoshiftmA",
#         # "airbrd",
#     ]

#     assert set(required_columns).issubset(
#         df.columns
#     ), "Not all required columns are present in the DataFrame"

@pytest.mark.needs_connection
def test_calc_kurucz_spectrum(verbose=3, plot=True, *args, **kwargs):
    from radis import calc_spectrum

    s = calc_spectrum(
        100,
        10000,  # cm-1
        species="Y++++", # should be converted to Y_V by radis.db.classes.to_conventional_name
        pressure=1.01325,  # bar
        Tgas=4000,  # K
        diluent={'H':1-0.01-1e-3, 'e-': 1e-3},        
        mole_fraction=0.01,
        path_length=15,  # cm
        databank='kurucz',
        verbose=verbose,
        potential_lowering=-500,
        cutoff=0
    )
    if plot:
        s.plot("radiance_noslit")

def test_barklem_pf():
    # Testing load_pf_Barklem2016 function
    pfTdat, pfdat = load_pf_Barklem2016()
    assert isinstance(pfTdat, pd.Series)
    assert isinstance(pfdat, pd.DataFrame)
    assert len(pfTdat) == 42
    assert len(pfdat) == 284

if __name__ == "__main__":
    test_read_kurucz()
    test_calc_kurucz_spectrum()
    test_barklem_pf()