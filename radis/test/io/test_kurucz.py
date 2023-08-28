from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from radis.levels.partfunc import PartFuncKurucz
from radis.api.kuruczapi import AdBKurucz
from radis.io.kurucz import fetch_kurucz

def test_adbkurucz_get_url():
    kurucz = AdBKurucz("Fe_I")

    # Testing get_url function
    url = kurucz.get_url(26, 00)
    assert url == "http://kurucz.harvard.edu/linelists/gfall/gf2600.all"

def test_fetch_kurucz():
    species="Na_I"
    hdf5_path,df = fetch_kurucz(species)


    required_columns = ['A', 'nu_lines', 'wav', 'El', 'eupper', 'gu', 'jlower', 'ju', 'id',
       'iso', 'gamRad', 'gamSta', 'gamvdW', 'Tdpair', 'wlnmair', 'loggf',
       'species', 'labellower', 'labelupper', 'ref', 'NLTElower', 'NLTEupper',
       'isonum', 'hyperfrac', 'isonumdi', 'isofrac', 'hypershiftlower',
       'hypershiftupper', 'hyperFlower', 'hypernotelower', 'hyperFupper',
       'hypternoteupper', 'strenclass', 'auto', 'landeglower', 'landegupper',
       'isoshiftmA', 'airbrd']

    assert set(required_columns).issubset(df.columns), "Not all required columns are present in the DataFrame"

@pytest.mark.needs_connection
def test_adbkurucz_external_data_functions():
    kurucz = AdBKurucz("Na_I")

    # Testing load_pf_Barklem2016 function
    pfTdat, pfdat = kurucz.load_pf_Barklem2016()
    assert isinstance(pfTdat, pd.Series)
    assert isinstance(pfdat, pd.DataFrame)
    assert len(pfTdat) == 42
    assert len(pfdat) == 284
