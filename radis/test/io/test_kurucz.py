from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from radis.api.kuruczapi import AdBKurucz


def test_adbkurucz_functions():
    kurucz = AdBKurucz()

    # Testing get_url function
    url = kurucz.get_url(26, 00)
    assert url == "http://kurucz.harvard.edu/linelists/gfall/gf2600.all"

#Unecessary,the method to use is radis air2vacuum
def test_air_to_vac():
    kurucz = AdBKurucz()
    wlair = 5000
    wlvac = kurucz.air_to_vac(wlair)
    assert round(wlvac, 5) == round(5001.39485, 5)


def test_partfn():
    kurucz = AdBKurucz()
    with patch.object(kurucz, "pfdat") as mock_pfdat, patch.object(
        kurucz, "pfTdat"
    ) as mock_pfTdat:
        # Assuming partition function values for corresponding temperatures
        mock_pfdat.loc.__getitem__.return_value = pd.Series(
            [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        )
        # Temperature range from 2000K to 3000K
        mock_pfTdat.values.flatten.return_value = np.array(
            [2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0]
        )
        key = "H_I"
        T = 2500.0
        Q = kurucz.partfn(key, T)
        # As T = 2500 is exactly in the middle of [2400, 2600],
        # the interpolated partition function value should be in the middle of [3.0, 3.5]
        expected_Q = 3.25
        assert Q == expected_Q


@pytest.mark.needs_connection
def test_adbkurucz_external_data_functions():
    kurucz = AdBKurucz()

    # Testing load_pf_Barklem2016 function
    pfTdat, pfdat = kurucz.load_pf_Barklem2016()
    assert isinstance(pfTdat, pd.Series)
    assert isinstance(pfdat, pd.DataFrame)
    assert len(pfTdat) == 42
    assert len(pfdat) == 284

    # Testing load_atomicdata function
    atomic_data = kurucz.load_atomicdata()
    assert isinstance(atomic_data, pd.DataFrame)
    assert "ielem" in atomic_data.columns
    assert "ionizationE1" in atomic_data.columns
