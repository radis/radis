from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from radis.levels.partfunc import PartFuncKurucz
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


def test_partfunckurucz():
    atom = "Ca"
    ionization_state = "00"
    kurucz_data = PartFuncKurucz(atom, ionization_state)

    temperatures = np.array([1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1.5e-01, 2e-01, 3e-01, 5e-01, 7e-01, 1.0, 1.3, 1.7, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 70.0, 100.0, 130.0, 170.0, 200.0, 250.0, 300.0, 500.0, 700.0, 1e03, 1.5e03, 2e03, 3e03, 4e03, 5e03, 6e03, 7e03, 8e03, 9e03, 1e04])
    pf_values = np.array([1.0]*33 + [1.00016, 1.00701, 1.04991, 1.17173, 1.41809, 1.85725, 2.60365, 3.81954, 5.69578])
    
    for T, expected_pf in zip(temperatures, pf_values):
        assert np.isclose(kurucz_data._at(T), expected_pf, atol=1e-4), f"Failed at T={T}"
    
    T_test = 3000  # Temperature at which to test interpolation

    # Compute expected partition function value using np.interp directly on the original data
    expected_pf = np.interp(T_test, temperatures , pf_values)

    assert np.isclose(kurucz_data._at(T_test), expected_pf, atol=1e-4), f"Interpolation failed at T={T_test}"

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
