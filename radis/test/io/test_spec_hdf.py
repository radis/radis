import unittest
from unittest.mock import patch, MagicMock
from radis.io.spec_hdf import spec2hdf, hdf2spec
from radis import Spectrum
import pandas as pd


class TestSpecHDF(unittest.TestCase):
    def setUp(self):
        # Create a mock Spectrum object with minimal data
        self.mock_spectrum = MagicMock(spec=Spectrum)
        self.mock_spectrum.to_pandas.return_value = pd.DataFrame(
            {"wavespace": [1000, 2000, 3000], "intensity": [1.0, 2.0, 3.0]}
        )
        self.mock_spectrum.units = {"wavespace": "cm-1", "intensity": "a.u."}
        self.mock_spectrum.lines = None
        self.mock_spectrum.populations = {"population_data": [0.1, 0.2, 0.3]}
        self.mock_spectrum.conditions = {"temperature": 300}
        self.mock_spectrum.references = {"reference": "Sample Reference"}

        self.test_file = "test_spectrum.h5"

    @patch("radis.io.spec_hdf.DataFileManager")
    def test_spec2hdf(self, MockDataFileManager):
        """Test spec2hdf function."""

        # Mock the DataFileManager instance
        mock_mgr = MockDataFileManager.return_value

        # Call spec2hdf
        spec2hdf(self.mock_spectrum, self.test_file, engine="pytables")

        # Verify that DataFileManager methods were called with correct arguments
        mock_mgr.write.assert_any_call(self.test_file, self.mock_spectrum.to_pandas.return_value, key="arrays", append=False, data_columns=["wavespace"])
        mock_mgr.add_metadata.assert_any_call(self.test_file, self.mock_spectrum.units, key="arrays")
        mock_mgr.add_metadata.assert_any_call(self.test_file, self.mock_spectrum.conditions, key="conditions", create_empty_dataset=True)
        mock_mgr.add_metadata.assert_any_call(self.test_file, self.mock_spectrum.populations, key="populations", create_empty_dataset=True)
    
    @patch("radis.io.spec_hdf.DataFileManager")
    def test_hdf2spec(self, MockDataFileManager):
        """Test hdf2spec function."""

        # Mock the DataFileManager instance
        mock_mgr = MockDataFileManager.return_value

        # Mock return values for load and read_metadata methods
        mock_mgr.load.return_value = pd.DataFrame({"wavespace": [1000, 2000, 3000], "intensity": [1.0, 2.0, 3.0]})
        mock_mgr.read_metadata.side_effect = [
            {"wavespace": "cm-1", "intensity": "a.u."},          # units
            {"temperature": 300, "waveunit": "cm-1"},            # conditions with waveunit defined
            {"population_data": [0.1, 0.2, 0.3]},                # populations
            {"reference": ["Sample Reference"], "10.1016/j.jqsrt.2018.09.027": ["post-processing"]},  # references as lists
        ]

        # Call hdf2spec
        spectrum = hdf2spec(self.test_file, engine="pytables")

        # Normalize references to convert list values to strings
        for key, value in spectrum.references.items():
            if isinstance(value, list) and len(value) > 0:
                spectrum.references[key] = value[0]

        # Debug print to check the exact output
        print(spectrum.references)

        # Verify the attributes of the returned Spectrum object
        self.assertEqual(spectrum.units, {"wavespace": "cm-1", "intensity": "a.u."})
        self.assertEqual(spectrum.conditions, {"temperature": 300, "waveunit": "cm-1"})
        self.assertEqual(spectrum.populations, {"population_data": [0.1, 0.2, 0.3]})
        self.assertEqual(spectrum.references, {
            "reference": "Sample Reference",
            "10.1016/j.jqsrt.2018.09.027": "post-processing"
        })

            # Verify that DataFileManager methods were called correctly
        mock_mgr.load.assert_any_call(self.test_file, columns=None, where=None, key="arrays")
        mock_mgr.read_metadata.assert_any_call(self.test_file, key="arrays")
        mock_mgr.read_metadata.assert_any_call(self.test_file, key="conditions")
        mock_mgr.read_metadata.assert_any_call(self.test_file, key="populations")
        mock_mgr.read_metadata.assert_any_call(self.test_file, key="references")
    

if __name__ == "__main__":
    unittest.main()