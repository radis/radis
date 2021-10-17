# -*- coding: utf-8 -*-
"""
Test query functions

-------------------------------------------------------------------------------

"""

from os.path import exists, join

import pytest

from radis.misc.warning import DeprecatedFileWarning


# ignored by pytest with argument -m "not needs_connection"
@pytest.mark.needs_connection
def test_fetch_astroquery(verbose=True, *args, **kwargs):
    """Test astroquery"""
    from radis.io.query import fetch_astroquery

    df = fetch_astroquery("CO2", 1, 2200, 2400, verbose=verbose, cache=False)

    assert df.iloc[0].id == 2
    assert df.iloc[0].iso == 1

    return


# ignored by pytest with argument -m "not needs_connection"
@pytest.mark.needs_connection
def test_fetch_astroquery_empty(verbose=True, *args, **kwargs):
    """Test astroquery: get a spectral range where there are no lines"""
    from radis.io.query import fetch_astroquery

    df = fetch_astroquery(
        2, 1, 25000, 50000, verbose=verbose, cache=False
    )  # 200-400 nm

    assert len(df) == 0

    return


# ignored by pytest with argument -m "not needs_connection"
@pytest.mark.needs_connection
def test_fetch_astroquery_cache(verbose=True, *args, **kwargs):
    """Test astroquery cache file generation

    - Cache file is generated
    - Cache file is loaded.
    - Different metadata raises an error
    """
    from astroquery.hitran import Hitran

    from radis.io.query import CACHE_FILE_NAME, fetch_astroquery

    df = fetch_astroquery(
        "CO2",
        1,
        2200,
        2400,
        verbose=verbose,
        cache="regen",  # Delete existing cache file if needed
    )

    # assert cache file was created
    assert exists(
        join(
            Hitran.cache_location,
            CACHE_FILE_NAME.format(
                **{"molecule": "CO2", "isotope": 1, "wmin": 2200, "wmax": 2400}
            ),
        )
    )

    # Try to load with different metadata: expect a DeprecatedFileError
    with pytest.raises(DeprecatedFileWarning):
        df2 = fetch_astroquery(
            "CO2",
            1,
            2200,
            2400,
            verbose=verbose,
            cache="force",  # force load of cache file
            expected_metadata={"not_in_metadata": True},
        )
    # Try to load with correct metadata. Expect to work and return the dataframe
    df2 = fetch_astroquery(
        "CO2",
        1,
        2200,
        2400,
        verbose=verbose,
        cache="force",  # force load of cache file
    )
    assert (df == df2).all().all()

    return


def _run_testcases(verbose=True, *args, **kwargs):

    test_fetch_astroquery(verbose=verbose, *args, **kwargs)
    test_fetch_astroquery_empty(verbose=verbose, *args, **kwargs)
    test_fetch_astroquery_cache(verbose=verbose, *args, **kwargs)

    return True


if __name__ == "__main__":
    print("test_query.py: ", _run_testcases(verbose=True))
