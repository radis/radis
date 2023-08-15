from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from radis.levels.partfunc import PartFuncKurucz
from radis.api.kuruczapi import AdBKurucz

def test_adbkurucz_get_url():
    kurucz = AdBKurucz("Fe_I")

    # Testing get_url function
    url = kurucz.get_url(26, 00)
    assert url == "http://kurucz.harvard.edu/linelists/gfall/gf2600.all"

@pytest.mark.needs_connection
def test_adbkurucz_external_data_functions():
    kurucz = AdBKurucz("Na_I")

    # Testing load_pf_Barklem2016 function
    pfTdat, pfdat = kurucz.load_pf_Barklem2016()
    assert isinstance(pfTdat, pd.Series)
    assert isinstance(pfdat, pd.DataFrame)
    assert len(pfTdat) == 42
    assert len(pfdat) == 284
