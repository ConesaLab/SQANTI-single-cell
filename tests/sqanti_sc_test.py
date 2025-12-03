import os
import sys
from io import StringIO
from unittest.mock import MagicMock, mock_open, patch
import subprocess

import pandas as pd
import pytest

# Add package paths for imports
sqanti3_src_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../SQANTI3")
)
if sqanti3_src_path not in sys.path:
    sys.path.insert(0, sqanti3_src_path)

sqanti_sc_src_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../src")
)
if sqanti_sc_src_path not in sys.path:
    sys.path.insert(0, sqanti_sc_src_path)

# Now import from the modularized implementation
from qc_io import fill_design_table
from sqanti3_qc_runner import run_sqanti3_qc
from classification_enrichment import annotate_with_ujc_hash, annotate_with_cell_metadata
from qc_reports import generate_report, generate_multisample_report
from qc_pipeline import main as pipeline_main
from cell_metrics import calculate_metrics_per_cell


@pytest.fixture
def mock_args(tmpdir):
    """Fixture for creating a mock arguments object."""
    class MockArgs:
        def __init__(self):
            # Set default values for all arguments
            self.test_data_dir = os.path.join(
                os.path.dirname(__file__), "test_data"
            )
            self.inDESIGN = str(tmpdir.join("design.csv"))
            self.input_dir = str(tmpdir)
            self.refGTF = os.path.join(
                self.test_data_dir, "reference_transcriptome.gtf"
            )
            self.refFasta = os.path.join(
                self.test_data_dir, "reference_genome.fasta"
            )
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
            self.version = "0.1.1"
            self.log_level = "INFO"

            # Create dummy design file and output directory
            design_content = "sampleID,file_acc\nsample1,file1"
            tmpdir.join("design.csv").write(design_content)
            os.makedirs(self.out_dir, exist_ok=True)

    return MockArgs()


@patch('qc_io.glob.glob')
@patch('qc_io.os.path.isfile')
def test_fill_design_table(mock_isfile, mock_glob, mock_args, tmpdir):
    """Test discovery of cell association files (BAM and text)."""
    bam_file_path = os.path.join(mock_args.input_dir, "file1.bam")
    assoc_file_path = os.path.join(
        mock_args.input_dir, "file2_cell_association.txt"
    )

    # For file1, glob finds a bam. For file2, glob finds nothing,
    # isfile finds assoc txt.
    mock_glob.side_effect = [[bam_file_path], []]
    mock_isfile.side_effect = lambda path: path == assoc_file_path

    # Update design to have two samples
    design_content = "sampleID,file_acc\nsample1,file1\nsample2,file2"
    with open(mock_args.inDESIGN, 'w') as f:
        f.write(design_content)

    df = fill_design_table(mock_args)

    assert isinstance(df, pd.DataFrame)
    assert 'cell_association' in df.columns
    assert df['cell_association'][0] == os.path.abspath(bam_file_path)
    assert df['cell_association'][1] == os.path.abspath(assoc_file_path)


def test_fill_design_table_preserve_existing(mock_args, tmpdir):
    """If cell_association_file is already set, it should not be overwritten."""
    # Prepare a design with pre-filled association path
    predefined = os.path.join(mock_args.input_dir, "pre_filled.tsv")
    design_content = (
        "sampleID,file_acc,cell_association\n"
        f"sample1,file1,{predefined}\n"
    )
    with open(mock_args.inDESIGN, 'w') as f:
        f.write(design_content)

    df = fill_design_table(mock_args)
    assert df.loc[0, 'cell_association'] == predefined


@patch('sqanti3_qc_runner.run_command')
@patch('sqanti3_qc_runner.os.path.isfile')
@patch('sqanti3_qc_runner.glob.glob')
def test_run_sqanti3_qc_gtf(
        mock_glob, mock_isfile, mock_run_command, mock_args, capsys):
    """Test running the QC step with a GTF file."""
    # Mock glob and isfile to simulate finding files
    gtf_path = os.path.join(mock_args.test_data_dir, "isoforms.gtf")
    mock_glob.return_value = [gtf_path]
    mock_isfile.return_value = True

    # Mock design file DataFrame
    design_df = pd.DataFrame({
        'sampleID': ['sample1'],
        'file_acc': ['file1']
    })

    # Call function directly with provided DataFrame
    run_sqanti3_qc(mock_args, design_df)

    # Capture printed output and check for expected messages
    captured = capsys.readouterr()
    assert "[INFO] Running SQANTI-sc qc for" in captured.out

    # Assert that run_command was called
    mock_run_command.assert_called()


def test_run_sqanti3_qc_fusion_requires_orf(mock_args):
    """In fusion mode without --orf_input and not skipping ORF, raise error."""
    design_df = pd.DataFrame({
        'sampleID': ['sample1'],
        'file_acc': ['file1']
    })
    # Create minimal environment so the function reaches fusion check
    mock_args.is_fusion = True
    mock_args.skipORF = False
    mock_args.orf_input = None
    mock_args.input_dir = os.getcwd()
    # Ensure reference files exist
    mock_args.refFasta = __file__
    mock_args.refGTF = __file__

    # Place a fake gtf file to match glob
    fake_gtf = os.path.join(mock_args.input_dir, 'file1.gtf')
    try:
        with open(fake_gtf, 'w') as f:
            f.write("gtf")
        with patch('sqanti3_qc_runner.glob.glob', return_value=[fake_gtf]):
            with patch('sqanti3_qc_runner.run_command'):
                with pytest.raises(ValueError, match='--orf_input must be provided'):
                    run_sqanti3_qc(mock_args, design_df)
    finally:
        try:
            os.remove(fake_gtf)
        except Exception:
            pass

@patch('classification_enrichment.subprocess.check_call')
@patch('classification_enrichment.os.path.exists')
@patch('classification_enrichment.os.remove')
@patch('classification_enrichment.os.rename')
@patch('classification_enrichment.pd.read_csv')
@patch('classification_enrichment.pd.DataFrame.to_csv')
@patch('builtins.open', new_callable=mock_open)
def test_make_UJC_hash(
        mock_open, mock_to_csv, mock_read_csv, mock_rename, mock_remove,
        mock_exists, mock_check_call, mock_args):
    """Test the UJC hashing function."""
    # Mock file existence and subprocess calls
    mock_exists.return_value = False
    mock_check_call.return_value = 0

    # Mock DataFrames
    class_df_cols = {
        "isoform": ["PB.1.1", "PB.2.1"], "chr": ["chr1", "chr1"],
        "strand": ["+", "+"], "structural_category": ["FSM", "FSM"],
        "associated_gene": ["geneA", "geneB"],
        "associated_transcript": ["txA", "txB"], "ref_length": [100, 200],
        "length": [100, 200]
    }
    class_df = pd.DataFrame(class_df_cols)

    class_df_for_ujc = class_df[
        ["isoform", "chr", "strand", "associated_transcript"]
    ]

    ujc_df = pd.DataFrame({
        "isoform": ["PB.1.1", "PB.2.1"],
        "jxn_string": ["chr1_+_100_200", pd.NA]
    })

    # Side effect for pd.read_csv:
    # 1. Read initial classification file (for UJC calculation)
    # 2. Read UJC temp file
    # 3. Read full, original classification file (for merging)
    mock_read_csv.side_effect = [class_df_for_ujc, ujc_df, class_df]

    # Mock input design DataFrame
    df = pd.DataFrame({"sampleID": ["sample1"], "file_acc": ["file1"]})

    # Call the function
    annotate_with_ujc_hash(mock_args, df)

    # Verify that gffread, gtftools and bedtools were called
    assert mock_check_call.call_count == 3
    assert "gffread" in mock_check_call.call_args_list[0][0][0]
    assert "gtftools" in mock_check_call.call_args_list[1][0][0]
    assert "bedtools" in mock_check_call.call_args_list[2][0][0]

    # Verify that the final merged file is written correctly
    mock_to_csv.assert_called_once()
    args_to_csv, kwargs_to_csv = mock_to_csv.call_args
    assert args_to_csv[0].endswith("_classification_tmp.txt")


@patch('classification_enrichment.pd.read_csv')
@patch('classification_enrichment.pd.DataFrame.to_csv', autospec=True)
@patch('classification_enrichment.os.path.isfile')
@patch('classification_enrichment.os.remove')
@patch('classification_enrichment.os.rename')
@patch('classification_enrichment.pysam.AlignmentFile')
def test_add_cell_data_reads_mode_bam(
        mock_alignment_file, mock_rename, mock_remove, mock_isfile,
        mock_to_csv, mock_read_csv, mock_args):
    """Test adding cell data in 'reads' mode with a BAM file."""
    # Test read mode with a BAM file
    mock_args.mode = "reads"
    bam_path = os.path.join(mock_args.input_dir, "file1.bam")
    design_df = pd.DataFrame({
        "sampleID": ["sample1"], "file_acc": ["file1"],
        "cell_association": [bam_path]
    })

    # Mocks for file system and pysam
    mock_isfile.return_value = True
    mock_read = MagicMock()
    mock_read.has_tag.side_effect = lambda tag: tag in ["XM", "CB"]
    mock_read.get_tag.side_effect = \
        lambda tag: {"XM": "UMI1", "CB": "CELL1"}[tag]
    mock_read.query_name = "iso1"
    mock_alignment_file.return_value.__enter__.return_value = [mock_read]

    # Mock for classification and junction files
    class_df = pd.DataFrame({"isoform": ["iso1", "iso2_dup1"]})
    junc_df = pd.DataFrame({"isoform": ["iso1", "iso2"]})
    mock_read_csv.side_effect = [class_df, junc_df]

    # Use a side effect to capture the dataframes passed to to_csv
    saved_dfs = []

    def to_csv_side_effect(self, path, *args, **kwargs):
        saved_dfs.append(self.copy())
    mock_to_csv.side_effect = to_csv_side_effect

    # Run function
    annotate_with_cell_metadata(mock_args, design_df)

    # Assertions
    mock_alignment_file.assert_called_with(bam_path, "rb", check_sq=False)
    assert mock_to_csv.call_count == 2  # Called for class file and junc file

    # Check that the classification df has CB and UMI
    final_class_df = saved_dfs[0]
    assert "CB" in final_class_df.columns
    assert "UMI" in final_class_df.columns
    assert final_class_df.loc[
        final_class_df.isoform == 'iso1', 'CB'
    ].iloc[0] == 'CELL1'

    # Check that the junctions df has CB
    final_junc_df = saved_dfs[1]
    assert "CB" in final_junc_df.columns


@patch('classification_enrichment.pd.read_csv')
@patch('classification_enrichment.pd.DataFrame.to_csv', autospec=True)
@patch('classification_enrichment.os.path.isfile')
@patch('classification_enrichment.os.remove')
@patch('classification_enrichment.os.rename')
def test_add_cell_data_isoform_mode_tsv(
        mock_rename, mock_remove, mock_isfile, mock_to_csv,
        mock_read_csv, mock_args):
    """Test adding cell data in 'isoforms' mode with a TSV file."""
    # Test isoform mode with a TSV file
    mock_args.mode = "isoforms"
    assoc_path = os.path.join(mock_args.input_dir, "file1_assoc.tsv")
    design_df = pd.DataFrame({
        "sampleID": ["sample1"], "file_acc": ["file1"],
        "cell_association": [assoc_path]
    })

    # Mocks for file system and dataframes
    mock_isfile.return_value = True
    assoc_df = pd.DataFrame({
        "pbid": ["iso1", "iso2"],
        "cell_barcodes": ["CELL1,CELL2", "CELL3"]
    })
    class_df = pd.DataFrame({
        "isoform": ["iso1", "iso2"],
        "RTS_stage": [True, False], "predicted_NMD": [False, True],
        "within_CAGE_peak": [True, True], "polyA_motif_found": [False, True]
    })
    mock_read_csv.side_effect = [assoc_df, class_df]

    # Use a side effect to capture the dataframe passed to to_csv
    saved_dfs = []

    def to_csv_side_effect(self, path, *args, **kwargs):
        saved_dfs.append(self.copy())
    mock_to_csv.side_effect = to_csv_side_effect

    # Run function
    annotate_with_cell_metadata(mock_args, design_df)

    # Assertions
    assert mock_to_csv.call_count == 1  # Not called for junctions in this mode

    final_class_df = saved_dfs[0]
    assert "CB" in final_class_df.columns
    assert "UMI" not in final_class_df.columns  # No UMI in isoform mode
    assert final_class_df.loc[
        final_class_df.isoform == 'iso1', 'CB'
    ].iloc[0] == 'CELL1,CELL2'
    assert final_class_df.loc[
        final_class_df.isoform == 'iso1', 'RTS_stage'
    ].iloc[0] == 'TRUE'


@patch('qc_reports.subprocess.run')
@patch('qc_reports.os.path.isfile')
def test_generate_multisample_report_skips_when_insufficient(mock_isfile, mock_run, mock_args, capsys):
    df = pd.DataFrame({"sampleID": ["s1"], "file_acc": ["f1"]})
    generate_multisample_report(mock_args, df)
    captured = capsys.readouterr()
    assert "fewer than 2".lower() in captured.out.lower()
    assert not mock_run.called


@patch('qc_reports.subprocess.run')
@patch('qc_reports.os.path.isfile')
def test_generate_multisample_report_runs_with_two(mock_isfile, mock_run, mock_args):
    df = pd.DataFrame({
        "sampleID": ["s1", "s2"],
        "file_acc": ["f1", "f2"]
    })
    # Fake two cell summary files present
    def _isfile_side_effect(path):
        return path.endswith('_SQANTI_cell_summary.txt.gz')
    mock_isfile.side_effect = _isfile_side_effect

    generate_multisample_report(mock_args, df)
    assert mock_run.called
    # Quick assert command includes both files
    cmd = mock_run.call_args[0][0]
    assert '--files' in cmd


@patch('qc_pipeline.build_parser')
@patch('qc_pipeline.run_sqanti3_qc')
@patch('qc_pipeline.annotate_with_ujc_hash')
@patch('qc_pipeline.annotate_with_cell_metadata')
@patch('qc_pipeline.calculate_metrics_per_cell')
@patch('qc_pipeline.generate_report')
def test_pipeline_main_smoke(
        mock_generate_report,
        mock_calc,
        mock_add,
        mock_hash,
        mock_get,
        mock_build_parser,
        tmpdir):
    class DummyArgs:
        def __init__(self, write_per_cell_outputs=False):
            self.inDESIGN = str(tmpdir.join('design.csv'))
            self.input_dir = str(tmpdir)
            self.out_dir = str(tmpdir)
            self.refGTF = __file__
            self.refFasta = __file__
            self.mode = 'reads'
            self.report = 'pdf'
            self.SKIPHASH = True
            self.write_per_cell_outputs = write_per_cell_outputs
            self.log_level = 'INFO'
            self.run_clustering = False
    # Write a minimal design
    tmpdir.join('design.csv').write('sampleID,file_acc,cell_association\ns1,f1,assoc.txt\n')

    class DummyParser:
        def parse_args(self_inner):
            return DummyArgs(write_per_cell_outputs=False)

    mock_build_parser.return_value = DummyParser()
    pipeline_main()  # should run through without errors
    assert mock_get.called
    assert mock_add.called
    assert mock_calc.called
    assert mock_generate_report.called


def test_calculate_metrics_divzero_safety(tmpdir, mock_args):
    """Minimal file setup to ensure function writes summary even with zero denoms."""
    # Arrange: create minimal classification file with no CB → empty valid set
    out_dir = str(tmpdir.join('output_dir'))
    os.makedirs(out_dir, exist_ok=True)
    mock_args.out_dir = out_dir
    file_acc = 'f1'
    sampleID = 's1'
    os.makedirs(os.path.join(out_dir, file_acc), exist_ok=True)
    prefix = os.path.join(out_dir, file_acc, sampleID)
    # classification with columns but CB empty
    cls = pd.DataFrame({
        'isoform': ['i1'], 'CB': [''], 'structural_category': ['full-splice_match'],
        'associated_gene': [''], 'associated_transcript': [''],
        'exons': [2], 'length': [100], 'ref_length': [200], 'chrom': ['chr1'],
    })
    cls.to_csv(f"{prefix}_classification.txt", sep='\t', index=False)
    # junc not required, but create empty
    pd.DataFrame({'isoform': []}).to_csv(f"{prefix}_junctions.txt", sep='\t', index=False)

    df = pd.DataFrame({'sampleID': [sampleID], 'file_acc': [file_acc]})
    # Act
    calculate_metrics_per_cell(mock_args, df)
    # Assert
    out_summary = f"{prefix}_SQANTI_cell_summary.txt.gz"
    assert os.path.isfile(out_summary)


@patch('qc_reports.subprocess.run')
@patch('qc_reports.os.path.isfile')
@patch('qc_reports.utilitiesPath', 'utilities')
def test_generate_report(mock_isfile, mock_run, mock_args, capsys):
    """Test report generation command."""
    mock_isfile.return_value = True
    mock_args.mode = "reads"
    df = pd.DataFrame({"sampleID": ["sample1"], "file_acc": ["file1"]})

    generate_report(mock_args, df)

    captured = capsys.readouterr()
    assert "SQANTI3 report generated for file1" in captured.out

    # Verify subprocess.run was called with the correct new command
    class_file = f"{mock_args.out_dir}/file1/sample1_classification.txt"
    junc_file = f"{mock_args.out_dir}/file1/sample1_junctions.txt"
    prefix = f"{mock_args.out_dir}/file1/sample1"
    # The script appends --cell_summary when the file exists
    cell_summary = f"{prefix}_SQANTI_cell_summary.txt.gz"
    expected_cmd = (
        f"Rscript utilities/SQANTI-sc_reads.R "
        f"{class_file} {junc_file} {mock_args.report} {prefix} "
        f"{mock_args.mode} --cell_summary {cell_summary}"
    ).strip()

    # Extract the actual command, handle potential whitespace differences
    actual_cmd = " ".join(mock_run.call_args[0][0].split())

    assert actual_cmd == expected_cmd


@patch('qc_reports.subprocess.run')
@patch('qc_reports.os.path.isfile')
@patch('qc_reports.utilitiesPath', 'utilities')
def test_generate_report_without_cell_summary(mock_isfile, mock_run, mock_args, capsys):
    """When the cell summary file is absent, the command should not include --cell_summary."""
    mock_args.mode = "reads"
    df = pd.DataFrame({"sampleID": ["sample1"], "file_acc": ["file1"]})

    class_file = f"{mock_args.out_dir}/file1/sample1_classification.txt"
    junc_file = f"{mock_args.out_dir}/file1/sample1_junctions.txt"
    prefix = f"{mock_args.out_dir}/file1/sample1"
    cell_summary = f"{prefix}_SQANTI_cell_summary.txt.gz"

    def _isfile_side_effect(path):
        # class and junction exist, cell summary does not
        return path in {class_file, junc_file}

    mock_isfile.side_effect = _isfile_side_effect

    generate_report(mock_args, df)

    actual_cmd = " ".join(mock_run.call_args[0][0].split())
    assert "--cell_summary" not in actual_cmd
    expected_cmd = (
        f"Rscript utilities/SQANTI-sc_reads.R "
        f"{class_file} {junc_file} {mock_args.report} {prefix} "
        f"{mock_args.mode}"
    ).strip()
    assert actual_cmd == expected_cmd


@patch('qc_reports.subprocess.run')
@patch('qc_reports.os.path.isfile')
@patch('qc_reports.utilitiesPath', 'utilities')
def test_generate_report_subprocess_error(mock_isfile, mock_run, mock_args, capsys):
    """If the report subprocess fails, an error message should be printed."""
    mock_args.mode = "reads"
    df = pd.DataFrame({"sampleID": ["sample1"], "file_acc": ["file1"]})

    # Everything exists, but subprocess.run raises
    mock_isfile.return_value = True
    mock_run.side_effect = subprocess.CalledProcessError(1, "Rscript")

    generate_report(mock_args, df)

    captured = capsys.readouterr()
    assert "Error generating report" in captured.out


def test_fill_design_table_missing_columns(mock_args):
    """Test error handling for design file with missing columns."""
    # Create a design file with missing columns
    design_content = "sampleID\nsample1\nsample2"
    mock_args.inDESIGN = mock_args.input_dir + "/design_missing_columns.csv"
    with open(mock_args.inDESIGN, 'w') as f:
        f.write(design_content)

    with pytest.raises(ValueError, match="Missing required columns"):
        fill_design_table(mock_args)


@patch('sqanti3_qc_runner.os.path.isfile')
def test_run_sqanti3_qc_missing_files(mock_isfile, mock_args, capsys):
    """Test error handling for missing reference files."""
    mock_isfile.return_value = False

    df = pd.DataFrame({
        "sampleID": ["sample1"],
        "file_acc": ["file1"]
    })

    with pytest.raises(SystemExit) as e:
        run_sqanti3_qc(mock_args, df)

    assert e.type == SystemExit
    assert e.value.code == -1

    captured = capsys.readouterr()
    assert "[ERROR] Genome or annotation file is missing." in captured.err


@patch('classification_enrichment.pysam.AlignmentFile')
@patch('classification_enrichment.os.path.isfile')
@patch('classification_enrichment.pd.read_csv')
@patch('classification_enrichment.pd.DataFrame.to_csv')
def test_add_cell_data_missing_tags(
        mock_to_csv, mock_read_csv, mock_isfile, mock_alignment_file,
        mock_args, capsys):
    """Test handling of BAM files missing CB/XM tags."""
    mock_args.mode = "reads"
    bam_path = os.path.join(mock_args.input_dir, "file1.bam")
    design_df = pd.DataFrame({
        "sampleID": ["sample1"], "file_acc": ["file1"],
        "cell_association": [bam_path]
    })

    mock_isfile.return_value = True
    mock_read = MagicMock()
    mock_read.has_tag.return_value = False  # Simulate missing tags
    mock_read.query_name = "iso1"
    mock_alignment_file.return_value.__enter__.return_value = [mock_read]

    annotate_with_cell_metadata(mock_args, design_df)

    captured = capsys.readouterr()
    assert "No cell data for file1. Skipping." in captured.out


@patch('classification_enrichment.pd.read_csv')
@patch('classification_enrichment.os.path.isfile')
def test_annotate_with_cell_metadata_missing_columns(
        mock_isfile, mock_read_csv, mock_args, capsys):
    """Test error handling for metadata file with missing columns."""
    mock_args.mode = 'reads'
    assoc_path = os.path.join(mock_args.input_dir, "file1_assoc.tsv")
    design_df = pd.DataFrame({
        "sampleID": ["sample1"], "file_acc": ["file1"],
        "cell_association": [assoc_path]
    })

    metadata_content = "id\tumi\niso1\tumi1"  # No cell_barcode column
    mock_isfile.return_value = True
    mock_read_csv.return_value = pd.read_csv(
        StringIO(metadata_content), sep='\t'
    )

    annotate_with_cell_metadata(mock_args, design_df)

    captured = capsys.readouterr()
    assert "missing required columns" in captured.err
