from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from radis import Spectrum
from radis.io.spec_hdf import hdf2spec, spec2hdf


@pytest.fixture
def mock_spectrum():
    """Fixture to create a mock Spectrum object with minimal data."""
    mock = MagicMock(spec=Spectrum)
    mock.to_pandas.return_value = pd.DataFrame(
        {"wavespace": [1000, 2000, 3000], "intensity": [1.0, 2.0, 3.0]}
    )
    mock.units = {"wavespace": "cm-1", "intensity": "a.u."}
    mock.lines = None
    mock.populations = {"population_data": [0.1, 0.2, 0.3]}
    mock.conditions = {"temperature": 300}
    mock.references = {"reference": "Sample Reference"}
    return mock


@pytest.fixture
def test_file():
    """Fixture for the test file name."""
    return "test_spectrum.h5"


@patch("radis.io.spec_hdf.DataFileManager")
def test_spec2hdf(MockDataFileManager, mock_spectrum, test_file):
    """Test spec2hdf function."""
    # Mock the DataFileManager instance
    mock_mgr = MockDataFileManager.return_value

    # Call spec2hdf
    spec2hdf(mock_spectrum, test_file, engine="pytables")

    # Verify that DataFileManager methods were called with correct arguments
    mock_mgr.write.assert_any_call(
        test_file,
        mock_spectrum.to_pandas.return_value,
        key="arrays",
        append=False,
        data_columns=["wavespace"],
    )
    mock_mgr.add_metadata.assert_any_call(test_file, mock_spectrum.units, key="arrays")
    mock_mgr.add_metadata.assert_any_call(
        test_file, mock_spectrum.conditions, key="conditions", create_empty_dataset=True
    )
    mock_mgr.add_metadata.assert_any_call(
        test_file,
        mock_spectrum.populations,
        key="populations",
        create_empty_dataset=True,
    )


@patch("radis.io.spec_hdf.DataFileManager")
def test_hdf2spec(MockDataFileManager, test_file):
    """Test hdf2spec function."""
    # Mock the DataFileManager instance
    mock_mgr = MockDataFileManager.return_value

    # Mock return values for load and read_metadata methods
    mock_mgr.load.return_value = pd.DataFrame(
        {"wavespace": [1000, 2000, 3000], "intensity": [1.0, 2.0, 3.0]}
    )
    mock_mgr.read_metadata.side_effect = [
        {"wavespace": "cm-1", "intensity": "a.u."},  # units
        {"temperature": 300, "waveunit": "cm-1"},  # conditions with waveunit defined
        {"population_data": [0.1, 0.2, 0.3]},  # populations
        {
            "reference": ["Sample Reference"],
            "10.1016/j.jqsrt.2018.09.027": ["post-processing"],
        },  # references as lists
    ]

    # Call hdf2spec
    spectrum = hdf2spec(test_file, engine="pytables")

    # Normalize references to convert list values to strings
    for key, value in spectrum.references.items():
        if isinstance(value, list) and len(value) > 0:
            spectrum.references[key] = value[0]

    # Verify the attributes of the returned Spectrum object
    assert spectrum.units == {"wavespace": "cm-1", "intensity": "a.u."}
    assert spectrum.conditions == {"temperature": 300, "waveunit": "cm-1"}
    assert spectrum.populations == {"population_data": [0.1, 0.2, 0.3]}
    assert spectrum.references == {
        "reference": "Sample Reference",
        "10.1016/j.jqsrt.2018.09.027": "post-processing",
    }

    # Verify that DataFileManager methods were called correctly
    mock_mgr.load.assert_any_call(test_file, columns=None, where=None, key="arrays")
    mock_mgr.read_metadata.assert_any_call(test_file, key="arrays")
    mock_mgr.read_metadata.assert_any_call(test_file, key="conditions")
    mock_mgr.read_metadata.assert_any_call(test_file, key="populations")
    mock_mgr.read_metadata.assert_any_call(test_file, key="references")


# Optional main block for running the test file directly
if __name__ == "__main__":
    pytest.main([__file__])
