import sys
import pytest
import pandas as pd
import os
import subprocess
from unittest.mock import patch, MagicMock, mock_open
import pysam
import glob
import hashlib
from io import StringIO

# Add SQANTI3 to Python path
sqanti3_src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../SQANTI3"))
if sqanti3_src_path not in sys.path:
    sys.path.insert(0, sqanti3_src_path)
from src.commands import run_command
from src.module_logging import qc_logger, update_logger

sqanti_sc_src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../src"))
if sqanti_sc_src_path not in sys.path:
    sys.path.insert(0, sqanti_sc_src_path)
from sqanti_sc import *

# Fixture for creating a mock arguments object
@pytest.fixture
def mock_args(tmpdir):
    class MockArgs:
        def __init__(self):
            # Set default values for all arguments
            self.test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")
            self.inDESIGN = str(tmpdir.join("design.csv"))
            self.input_dir = str(tmpdir)
            self.refGTF = os.path.join(self.test_data_dir, "reference_transcriptome.gtf")
            self.refFasta = os.path.join(self.test_data_dir, "reference_genome.fasta")
            self.mode = "reads"
            self.out_dir = str(tmpdir.join("output_dir"))
            self.report = "pdf"
            self.samtools_cpus = 2
            self.verbose = True
            self.factor = None
            self.SKIPHASH = False
            self.ignore_cell_summary = False
            self.min_ref_len = 0
            self.force_id_ignore = False
            self.genename = False
            self.short_reads = None
            self.SR_bam = None
            self.novel_gene_prefix = None
            self.aligner_choice = "minimap2"
            self.gmap_index = None
            self.sites = "ATAC,GCAG,GTAG"
            self.skipORF = False
            self.orf_input = None
            self.CAGE_peak = None
            self.polyA_motif_list = None
            self.polyA_peak = None
            self.phyloP_bed = None
            self.saturation = False
            self.isoform_hits = False
            self.ratio_TSS_metric = "max"
            self.mapping_cpus = 2
            self.chunks = 2
            self.is_fusion = False
            self.expression = None
            self.coverage = None
            self.window = 20
            self.fl_count = None
            self.isoAnnotLite = False
            self.gff3 = None
            self.version = "0.1.0"
            self.log_level = "INFO"

            # Create dummy design file
            design_content = "sampleID,file_acc\nsample1,file1\nsample2,file2"
            tmpdir.join("design.csv").write(design_content)

            os.makedirs(self.out_dir, exist_ok=True)

    return MockArgs()

def test_fill_design_table(mock_args):
    df = fill_design_table(mock_args)

    assert isinstance(df, pd.DataFrame)
    assert 'classification_file' in df.columns
    assert 'junction_file' in df.columns
    assert df['classification_file'][0] == f"{mock_args.input_dir}/file1/sample1_classification.txt"
    assert df['junction_file'][0] == f"{mock_args.input_dir}/file1/sample1_junctions.txt"

@patch('sqanti_sc.run_command')
@patch('sqanti_sc.os.path.isfile')
@patch('sqanti_sc.glob.glob')
def test_get_files_runSQANTI3_gtf(mock_glob, mock_isfile, mock_run_command, mock_args, capsys):
    # Mock glob and isfile to simulate finding files
    mock_glob.return_value = [os.path.join(mock_args.test_data_dir, "isoforms.gtf")]
    mock_isfile.return_value = True

    # Mock design file DataFrame
    design_df = pd.DataFrame({
        'sampleID': ['sample1'],
        'file_acc': ['file1']
    })

    # Mock pandas read_csv to return the mocked DataFrame
    with patch('sqanti_sc.pd.read_csv', return_value=design_df):
        get_files_runSQANTI3(mock_args, design_df)

    # Capture printed output and check for expected messages
    captured = capsys.readouterr()
    assert "[INFO] Running SQANTI-sc qc for sample" in captured.out

    # Assert that run_command was called
    mock_run_command.assert_called()

@patch('sqanti_sc.subprocess.check_call')
@patch('sqanti_sc.subprocess.call')
@patch('sqanti_sc.os.path.exists')
@patch('sqanti_sc.os.remove')
@patch('sqanti_sc.pd.read_csv')
@patch('sqanti_sc.pd.DataFrame.to_csv')
@patch('builtins.open', new_callable=mock_open, read_data="chr1\tHAVANA\ttranscript\t1\t1000\t.\t+\t.\ttranscript_id \"ENST00000456328.2\";\n")
def test_make_UJC_hash(mock_file, mock_to_csv, mock_read_csv, mock_remove, mock_exists, mock_call, mock_check_call, mock_args):
    # Mock file existence and subprocess calls
    mock_exists.return_value = False
    mock_check_call.return_value = 0
    mock_call.return_value = 0

    # Mock DataFrames for classification and UJC files
    classfile_df = pd.DataFrame({
        "isoform": ["ENST00000456328.2", "ENST00000456329.2"],
        "chr": ["chr1", "chr2"],
        "strand": ["+", "-"],
        "associated_transcript": ["transcript1", "transcript2"]
    })

    ujc_df = pd.DataFrame({
        "isoform": ["ENST00000456328.2", "ENST00000456329.2"],
        "jxn_string": ["chr1_+_100_200", "chr2_-_300_400"]
    })

    # Mock read_csv to return these DataFrames in sequence
    mock_read_csv.side_effect = [classfile_df, ujc_df]

    # Mock input design DataFrame
    df = pd.DataFrame({
        "sampleID": ["sample1"],
        "file_acc": ["file1"]
    })

    # Call the function under test
    make_UJC_hash(mock_args, df)

    # Verify that commands were executed
    assert mock_check_call.call_count == 2
    assert "gtftools" in mock_check_call.call_args_list[0][0][0]
    assert "bedtools" in mock_check_call.call_args_list[1][0][0]

    # Verify file operations with corrected paths
    expected_path = f"{mock_args.out_dir}/file1/sample1"
    mock_remove.assert_any_call(f"{expected_path}tmp_introns.bed")
    mock_remove.assert_any_call(f"{expected_path}tmp_UJC.txt")
    mock_remove.assert_any_call(f"{expected_path}_temp.txt")

    # Verify that to_csv was called once with correct arguments
    mock_to_csv.assert_called_once()
    args, kwargs = mock_to_csv.call_args
    assert kwargs['sep'] == '\t'
    assert kwargs['index'] is False

    # Capture the DataFrame passed to to_csv by mocking its creation earlier in the function
    merged_df = pd.merge(classfile_df, ujc_df, on="isoform", how="left")
    
    # Simulate jxn_string creation for validation
    def create_jxn_string(row):
        if pd.isna(row["jxn_string"]):
            return f"{row['chr']}_{row['strand']}_monoexon_{row['associated_transcript']}"
        return row["jxn_string"]

    merged_df["jxn_string"] = merged_df.apply(create_jxn_string, axis=1)
    merged_df["jxnHash"] = merged_df["jxn_string"].apply(lambda x: hashlib.sha256(x.encode("utf-8")).hexdigest())

    # Verify jxn_string and jxnHash creation
    expected_jxn_strings = [
        "chr1_+_100_200",
        "chr2_-_300_400"
    ]
    
    expected_hashes = [
        hashlib.sha256("chr1_+_100_200".encode("utf-8")).hexdigest(),
        hashlib.sha256("chr2_-_300_400".encode("utf-8")).hexdigest()
    ]

    for i, row in merged_df.iterrows():
        assert row["jxn_string"] == expected_jxn_strings[i]
        assert row["jxnHash"] == expected_hashes[i]

@patch('sqanti_sc.pysam.AlignmentFile')
@patch('sqanti_sc.os.path.isfile')
@patch('sqanti_sc.glob.glob')
@patch('sqanti_sc.pd.read_csv')
def test_add_cell_data_bam(mock_read_csv, mock_glob, mock_isfile, mock_alignment_file, mock_args):
    # Test add_cell_data function with BAM file input
    mock_isfile.return_value = True
    mock_glob.return_value = [os.path.join(mock_args.test_data_dir, "unmapped_reads.bam")]

    mock_bam = MagicMock()
    mock_read = MagicMock()
    mock_read.has_tag.side_effect = lambda tag: tag in ["XM", "CB"]
    mock_read.get_tag.side_effect = lambda tag: {"XM": "UMI1", "CB": "cell1"}[tag]
    mock_read.query_name = "iso1"
    mock_bam.return_value.__enter__.return_value = [mock_read]
    mock_alignment_file.return_value = mock_bam

    classification_df = pd.DataFrame({
        "isoform": ["iso1", "iso2"],
        "chr": ["chr1", "chr2"],
        "strand": ["+", "-"]
    })
    mock_read_csv.return_value = classification_df

    df = pd.DataFrame({
        "sampleID": ["sample1"],
        "file_acc": ["file1"]
    })

    add_cell_data(mock_args, df)

    mock_alignment_file.assert_called_with(os.path.join(mock_args.test_data_dir, "unmapped_reads.bam"), "rb", check_sq=False)

@patch('sqanti_sc.pd.read_csv')
@patch('sqanti_sc.os.path.isfile')
@patch('sqanti_sc.glob.glob')
@patch('sqanti_sc.os.remove')
@patch('sqanti_sc.pd.DataFrame.to_csv')
def test_add_cell_data_metadata(mock_to_csv, mock_remove, mock_glob, mock_isfile, mock_read_csv, mock_args, tmpdir, capsys):
    # Define file paths
    metadata_file = tmpdir.join("file1_cell_association.txt")
    tmp_classification_file = tmpdir.join("output_dir", "file1", "sample1_classification_tmp.txt")  # Adjusted path
    output_dir = tmpdir.join("output_dir", "file1")  # Adjusted path

    # Create correct metadata content
    metadata_content = "id\tcell_barcode\tumi\niso1\tcell1\tumi1\niso2\tcell2\tumi2"
    metadata_file.write(metadata_content)

    # Create the expected directory structure
    output_dir.ensure(dir=True)

    # Mock file existence checks
    mock_isfile.side_effect = lambda path: str(path) in [str(metadata_file), str(tmp_classification_file)]

    # Mock BAM files as empty
    mock_glob.return_value = []

    # Mock classification DataFrame *with* 'RTS_stage' column
    classification_df = pd.DataFrame({
        "isoform": ["iso1", "iso2"],
        "chr": ["chr1", "chr2"],
        "strand": ["+", "-"],
        "RTS_stage": [True, False],
        "predicted_NMD": [False, True],
        "within_CAGE_peak": [True, False],
        "polyA_motif_found": [False, False]
    })

    # Mock metadata DataFrame
    metadata_df = pd.DataFrame({
        "id": ["iso1", "iso2"],
        "cell_barcode": ["cell1", "cell2"],
        "umi": ["umi1", "umi2"]
    })

    # Configure mock_read_csv to return different DataFrames based on input path
    def read_csv_side_effect(filepath, *args, **kwargs):
        if filepath == str(metadata_file):
            return metadata_df
        elif filepath == str(tmp_classification_file):
            return classification_df
        raise ValueError(f"Unexpected file path: {filepath}")

    mock_read_csv.side_effect = read_csv_side_effect

    # Create input dataframe
    df = pd.DataFrame({
        "sampleID": ["sample1"],
        "file_acc": ["file1"]
    })

    # Set mock_args.input_dir to tmpdir
    mock_args.input_dir = str(tmpdir)

    # Execute function
    add_cell_data(mock_args, df)

    # Verify output
    captured = capsys.readouterr()
    assert "UMI and cell barcode information successfully added" in captured.out

    # Verify file removal
    mock_remove.assert_called_once_with(str(tmp_classification_file))

    # Verify to_csv was called
    mock_to_csv.assert_called_once()

@patch('sqanti_sc.subprocess.run')
@patch('sqanti_sc.os.path.isfile')
@patch('sqanti_sc.utilitiesPath', 'utilities')  # Mock the global utilitiesPath
def test_generate_report(mock_isfile, mock_run, mock_args, capsys):
    mock_isfile.return_value = True

    # Add necessary attributes to mock_args
    mock_args.ignore_cell_summary = False
    mock_args.report = "pdf"
    mock_args.input_dir = "input"
    mock_args.out_dir = "output"

    df = pd.DataFrame({
        "sampleID": ["sample1"],
        "file_acc": ["file1"]
    })

    generate_report(mock_args, df)

    # Capture printed output
    captured = capsys.readouterr()

    # Check stdout for the report generation message
    assert "Generating SQANTI3 report for file1" in captured.out  # Match the exact message
    assert "SQANTI3 report successfully generated for file1" in captured.out

    # Verify subprocess.run was called with the correct command
    expected_cmd = (
        "Rscript utilities/SQANTI-sc_reads.R "
        "output/file1/sample1_classification.txt "
        "pdf output/file1/sample1"
    )
    mock_run.assert_called_once_with(expected_cmd, shell=True, check=True)

# Additional tests for edge cases and error handling

def test_fill_design_table_missing_columns(mock_args):
    # Create a design file with missing columns
    design_content = "sampleID\nsample1\nsample2"
    mock_args.inDESIGN = mock_args.input_dir + "/design_missing_columns.csv"
    with open(mock_args.inDESIGN, 'w') as f:
        f.write(design_content)

    with pytest.raises(ValueError, match="Missing required columns in design table"):
        fill_design_table(mock_args)

@patch('sqanti_sc.os.path.isfile')
def test_get_files_runSQANTI3_missing_files(mock_isfile, mock_args, capsys):
    mock_isfile.return_value = False
    
    df = pd.DataFrame({
        "sampleID": ["sample1"],
        "file_acc": ["file1"]
    })

    with pytest.raises(SystemExit) as e:
        get_files_runSQANTI3(mock_args, df)

    assert e.type == SystemExit
    assert e.value.code == -1

    captured = capsys.readouterr()
    assert "[ERROR] Genome or annotation file is missing." in captured.err

@patch('sqanti_sc.pysam.AlignmentFile')
@patch('sqanti_sc.os.path.isfile')
@patch('sqanti_sc.glob.glob')
@patch('sqanti_sc.pd.read_csv')
def test_add_cell_data_missing_tags(mock_read_csv, mock_glob, mock_isfile, mock_alignment_file, mock_args, capsys):
    mock_isfile.return_value = True
    mock_glob.return_value = [os.path.join(mock_args.test_data_dir, "unmapped_reads.bam")]

    mock_bam = MagicMock()
    mock_read = MagicMock()
    mock_read.has_tag.side_effect = lambda tag: False  # Simulate missing tags
    mock_read.query_name = "iso1"
    mock_bam.return_value.__enter__.return_value = [mock_read]
    mock_alignment_file.return_value = mock_bam

    classification_df = pd.DataFrame({
        "isoform": ["iso1"],
        "chr": ["chr1"],
        "strand": ["+"]
    })
    mock_read_csv.return_value = classification_df

    df = pd.DataFrame({
        "sampleID": ["sample1"],
        "file_acc": ["file1"]
    })

    add_cell_data(mock_args, df)

    captured = capsys.readouterr()
    assert "No valid reads found in the BAM or cell association file" in captured.out

@patch('sqanti_sc.pd.read_csv')
@patch('sqanti_sc.os.path.isfile')
@patch('sqanti_sc.glob.glob')
def test_add_cell_data_metadata_missing_columns(mock_glob, mock_isfile, mock_read_csv, mock_args, tmpdir, capsys):
    # Test add_cell_data function with metadata file input but missing mandatory column
    metadata_content = "id\tumi\niso1\tumi1\niso2\tumi2" # No cell_barcode column
    metadata_file = tmpdir.join("file1_cell_association.txt")  # Changed to match the expected filename
    metadata_file.write(metadata_content)

    mock_glob.return_value = []  # No BAM files
    mock_isfile.side_effect = lambda x: x == str(metadata_file)

    classification_df = pd.DataFrame({
        "isoform": ["iso1", "iso2"],
        "chr": ["chr1", "chr2"],
        "strand": ["+", "-"]  # Fixed to match the length of other columns
    })
    mock_read_csv.side_effect = [pd.read_csv(StringIO(metadata_content), sep='\t'), classification_df]

    df = pd.DataFrame({
        "sampleID": ["sample1"],
        "file_acc": ["file1"]
    })

    add_cell_data(mock_args, df)

    captured = capsys.readouterr()
    assert "Cell association file" in captured.err
    assert "does not contain the required columns" in captured.err
