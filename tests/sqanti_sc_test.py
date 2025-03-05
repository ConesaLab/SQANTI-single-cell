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
    return "sampleID,file_acc\nwtc11_PBcDNA,ENCFF003QZT\nwtc11_ONTdRNA,ENCFF104BNW"

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
                
                # Check if new columns are created correctly
                assert "classification_file" in df.columns
                assert "junction_file" in df.columns
                assert df["classification_file"].iloc[0] == "/mock/input/ENCFF003QZT/wtc11_PBcDNA_classification.txt"
                assert df["junction_file"].iloc[1] == "/mock/input/ENCFF104BNW/wtc11_ONTdRNA_junctions.txt"
                
                # Ensure DataFrame is written back to CSV
                mock_to_csv.assert_called_once_with(mock_args.inDESIGN, sep=',', index=False)

def test_fill_design_table_incorrect_format(mock_args):
    """Test function with incorrectly formatted CSV design file."""
    bad_csv_data = "sampleID\nS1\nS2"  # Only one column
    
    with patch("builtins.open", mock_open(read_data=bad_csv_data)):
        with patch("pandas.read_csv", return_value=pd.read_csv(StringIO(bad_csv_data))):
            with patch("sys.exit") as mock_exit:
                with patch("sys.stderr", new_callable=StringIO) as mock_stderr:
                    fill_design_table(mock_args)
                    
                    # Ensure error message is printed
                    assert "ERROR: is incorrectly formatted" in mock_stderr.getvalue()
                    
                    # Ensure the function exits
                    mock_exit.assert_called_once_with(-1)