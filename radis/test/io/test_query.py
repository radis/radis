# -*- coding: utf-8 -*-
"""
Test query functions
"""

import pytest
from radis.io.query import fetch_astroquery

@pytest.mark.needs_connection   # ignored by pytest with argument -m "not needs_connection" 
def test_fetch_astroquery(verbose=True, *args, **kwargs):
    fetch_astroquery('CO2', 1, 2200, 2400, verbose=verbose)
    
    return


