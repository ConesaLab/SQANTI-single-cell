import pytest
import pandas as pd
import sys
from io import StringIO
from unittest.mock import patch, mock_open
from src.sqanti_sc import fill_design_table


# Mock Arguments Class
class MockArgs:
    def __init__(self, inDESIGN, input_dir):
        self.inDESIGN = inDESIGN
        self.input_dir = input_dir


@pytest.fixture
def mock_csv_data():
    """Mock valid CSV data with required columns."""
    return "sampleID,file_acc\nsample_1,file_acc_1\nsample_2,file_acc_2"


@pytest.fixture
def mock_args(tmp_path):
    """Mock arguments object."""
    return MockArgs(inDESIGN=tmp_path / "design.csv", input_dir="/mock/input")


def test_fill_design_table_correct_format(mock_csv_data, mock_args):
    """Test function with correctly formatted CSV design file."""
    with patch("builtins.open", mock_open(read_data=mock_csv_data)):
        with patch("pandas.read_csv", return_value=pd.read_csv(StringIO(mock_csv_data))):
            with patch("pandas.DataFrame.to_csv") as mock_to_csv:
                df = fill_design_table(mock_args)

                # Explicitly check expected columns
                expected_columns = {"sampleID", "file_acc", "classification_file", "junction_file"}
                assert set(df.columns) == expected_columns, f"Expected columns {expected_columns}, but got {set(df.columns)}"

                # Check if new columns are created correctly
                assert df["classification_file"].iloc[0] == "/mock/input/file_acc_1/sample_1_classification.txt"
                assert df["junction_file"].iloc[1] == "/mock/input/file_acc_2/sample_2_junctions.txt"

                # Ensure DataFrame is written back to CSV
                mock_to_csv.assert_called_once_with(mock_args.inDESIGN, sep=',', index=False)


def test_fill_design_table_incorrect_format(mock_args):
    """Test function with incorrectly formatted CSV design file (missing required columns)."""
    bad_csv_data = "sampleID\nsample_1\nsample_2"  # Missing 'file_acc' column

    with patch("builtins.open", mock_open(read_data=bad_csv_data)):
        with patch("pandas.read_csv", return_value=pd.read_csv(StringIO(bad_csv_data))):
            with pytest.raises(ValueError, match="ERROR: Missing required columns in design table"):
                fill_design_table(mock_args)
