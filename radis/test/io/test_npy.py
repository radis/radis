from unittest.mock import patch  # , MagicMock

import numpy as np
import pandas as pd
import pytest

from radis.io.npy import npy2df


@pytest.fixture
def mock_data():
    """Fixture to provide sample mock data for the arrays and file paths."""
    return {
        "mock_iso": np.array([1, 2, 3]),
        "mock_v0": np.array([1000.0, 2000.0, 3000.0]),
        "mock_da": np.array([0.1, 0.2, 0.3]),
        "mock_log_2gs": np.array([0.01, 0.02, 0.03]),
        "mock_S0": np.array([10.0, 20.0, 30.0]),
        "mock_El": np.array([0.5, 1.5, 2.5]),
        "mock_log_2vMm": np.array([0.001, 0.002, 0.003]),
        "mock_na": np.array([100, 200, 300]),
        "keywords": {
            "wav": "path/to/v0.npy",
            "int": "path/to/int.npy",
            "Pshft": "path/to/da.npy",
            "log_2gs": "path/to/log_2gs.npy",
            "Tdpair": "path/to/na.npy",
            "El": "path/to/El.npy",
        },
    }


@patch("numpy.load")
def test_npy2df_success(mock_load, mock_data):
    """Test successful conversion of numpy arrays to a DataFrame."""
    mock_load.side_effect = [
        mock_data["mock_iso"],
        mock_data["mock_v0"],
        mock_data["mock_da"],
        mock_data["mock_log_2gs"],
        mock_data["mock_S0"],
        mock_data["mock_El"],
        mock_data["mock_log_2vMm"],
        mock_data["mock_na"],
    ]

    # Call the function
    df = npy2df(mock_data["keywords"], verbose=False)

    # Expected DataFrame
    expected_df = pd.DataFrame(
        {
            "iso": mock_data["mock_iso"],
            "wav": mock_data["mock_v0"],
            "Pshft": mock_data["mock_da"],
            "log_2gs": mock_data["mock_log_2gs"],
            "Tdpair": mock_data["mock_na"],
            "log_2vMm": mock_data["mock_log_2vMm"],
            "int": mock_data["mock_S0"],
            "El": mock_data["mock_El"],
        }
    )

    # Assert the DataFrame is as expected
    pd.testing.assert_frame_equal(df, expected_df)


@patch("numpy.load")
def test_npy2df_file_not_found(mock_load, mock_data):
    """Test that FileNotFoundError is raised when a file is missing."""
    # Configure the mock to raise FileNotFoundError
    mock_load.side_effect = FileNotFoundError("File not found")

    with pytest.raises(FileNotFoundError):
        npy2df(mock_data["keywords"], verbose=False)


@patch("numpy.load")
def test_npy2df_verbose_output(mock_load, mock_data):
    """Test verbose output when verbose=True."""
    # Mock the return values
    mock_load.side_effect = [
        mock_data["mock_iso"],
        mock_data["mock_v0"],
        mock_data["mock_da"],
        mock_data["mock_log_2gs"],
        mock_data["mock_S0"],
        mock_data["mock_El"],
        mock_data["mock_log_2vMm"],
        mock_data["mock_na"],
    ]

    with patch("builtins.print") as mock_print:
        npy2df(mock_data["keywords"], verbose=2)  # Set verbose to 2
        assert mock_print.called
        mock_print.assert_any_call("Loading iso...", end=" ")
        mock_print.assert_any_call("Done!")
