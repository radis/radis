import unittest
from unittest.mock import patch, MagicMock
import numpy as np
import pandas as pd
from radis.io.npy import npy2df  


class TestNpy2Df(unittest.TestCase):
    def setUp(self):
        # Sample mock data for the arrays
        self.mock_iso = np.array([1, 2, 3])
        self.mock_v0 = np.array([1000.0, 2000.0, 3000.0])
        self.mock_da = np.array([0.1, 0.2, 0.3])
        self.mock_log_2gs = np.array([0.01, 0.02, 0.03])
        self.mock_S0 = np.array([10.0, 20.0, 30.0])
        self.mock_El = np.array([0.5, 1.5, 2.5])
        self.mock_log_2vMm = np.array([0.001, 0.002, 0.003])
        self.mock_na = np.array([100, 200, 300])

        # Dictionary of file paths
        self.keywords = {
            "wav": "path/to/v0.npy",
            "int": "path/to/int.npy",
            "Pshft": "path/to/da.npy",
            "log_2gs": "path/to/log_2gs.npy",
            "Tdpair": "path/to/na.npy",
            "El": "path/to/El.npy",
        }

    @patch("numpy.load")
    def test_npy2df_success(self, mock_load):
        """Test successful conversion of numpy arrays to a DataFrame."""

        # Mock the return values of np.load
        mock_load.side_effect = [
            self.mock_iso,
            self.mock_v0,
            self.mock_da,
            self.mock_log_2gs,
            self.mock_S0,
            self.mock_El,
            self.mock_log_2vMm,
            self.mock_na,
        ]

        # Call the function
        df = npy2df(self.keywords, verbose=False)

        # Expected DataFrame
        expected_df = pd.DataFrame(
            {
                "iso": self.mock_iso,
                "wav": self.mock_v0,
                "Pshft": self.mock_da,
                "log_2gs": self.mock_log_2gs,
                "Tdpair": self.mock_na,
                "log_2vMm": self.mock_log_2vMm,
                "int": self.mock_S0,
                "El": self.mock_El,
            }
        )

        # Assert the DataFrame is as expected
        pd.testing.assert_frame_equal(df, expected_df)

    @patch("numpy.load")
    def test_npy2df_file_not_found(self, mock_load):
        """Test that FileNotFoundError is raised when a file is missing."""
        # Configure the mock to raise FileNotFoundError
        mock_load.side_effect = FileNotFoundError("File not found")

        with self.assertRaises(FileNotFoundError):
            npy2df(self.keywords, verbose=False)

    @patch("numpy.load")
    def test_npy2df_verbose_output(self, mock_load):
        """Test verbose output when verbose=True."""
        # Mock the return values
        mock_load.side_effect = [
            self.mock_iso,
            self.mock_v0,
            self.mock_da,
            self.mock_log_2gs,
            self.mock_S0,
            self.mock_El,
            self.mock_log_2vMm,
            self.mock_na,
        ]

        with patch("builtins.print") as mock_print:
            npy2df(self.keywords, verbose=2)  # Set verbose to 2
            self.assertTrue(mock_print.called)
            mock_print.assert_any_call("Loading iso...", end=" ")
            mock_print.assert_any_call("Done!")

if __name__ == "__main__":
    unittest.main()
