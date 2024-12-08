import pytest
import tempfile
from os.path import join, exists
import os
import pandas as pd
from unittest.mock import patch, MagicMock
from radis.io.query import fetch_astroquery, _fix_astroquery_file_format, CACHE_FILE_NAME
from astroquery.hitran import Hitran


@pytest.mark.needs_connection
def test_fetch_astroquery_success():
    """Test successful data fetching with fetch_astroquery."""
    df = fetch_astroquery("CO2", 1, 2200, 2400, verbose=False, cache=False)
    assert not df.empty
    assert "wav" in df.columns
    assert "int" in df.columns


@pytest.mark.needs_connection
def test_fetch_astroquery_empty_range():
    """Test fetch_astroquery when no lines are found in the specified range."""
    df = fetch_astroquery("CO2", 1, 50000, 60000, verbose=False, cache=False)
    assert df.empty


@pytest.mark.needs_connection
def test_fetch_astroquery_invalid_molecule():
    """Test fetch_astroquery with an invalid molecule name."""
    with pytest.raises(NotImplementedError, match="not supported"):
        fetch_astroquery("INVALID", 1, 2200, 2400, verbose=False, cache=False)



@pytest.mark.needs_connection
def test_fetch_astroquery_cache():
    """Test fetch_astroquery with caching enabled."""

    with tempfile.TemporaryDirectory() as tmp_cache_dir:
        # Patch the environment variable for the astroquery cache
        with patch.dict(os.environ, {"ASTROQUERY_CACHE_DIR": tmp_cache_dir}):
            # Force astroquery to reinitialize the cache location
            Hitran.cache_location = tmp_cache_dir

            # Generate the cache file
            fetch_astroquery("CO2", 1, 2200, 2400, verbose=False, cache="regen")

            # Check if the cache file was created
            cache_file = CACHE_FILE_NAME.format(molecule="CO2", isotope=1, wmin=2200, wmax=2400)
            assert exists(join(tmp_cache_dir, cache_file))

            # Load with correct metadata and verify it returns the same DataFrame
            df1 = fetch_astroquery("CO2", 1, 2200, 2400, verbose=False, cache="force")
            df2 = fetch_astroquery("CO2", 1, 2200, 2400, verbose=False, cache="force")

            assert (df1 == df2).all().all()


@pytest.mark.needs_connection
def test_fix_astroquery_file_format():
    """Test _fix_astroquery_file_format for removing empty lines."""
    # Create a temporary file with empty lines
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_file:
        tmp_file.write("Line 1\n\nLine 2\n\n\nLine 3\n")
        tmp_file_path = tmp_file.name

    try:
        # Call the function to fix the file format
        _fix_astroquery_file_format(tmp_file_path)

        # Read back the file and check if blank lines are removed
        with open(tmp_file_path, "r") as file:
            lines = file.readlines()

        assert lines == ["Line 1\n", "Line 2\n", "Line 3\n"]

    finally:
        # Clean up the temporary file
        if os.path.exists(tmp_file_path):
            os.remove(tmp_file_path)


def test_fix_astroquery_file_format_file_not_found():
    """Test _fix_astroquery_file_format with a non-existent file."""
    with pytest.raises(FileNotFoundError):
        _fix_astroquery_file_format("non_existent_file.txt")


if __name__ == "__main__":
    # Run pytest programmatically when the script is executed directly
    pytest.main(["-v", __file__])
