import os
import sys
from io import StringIO
from unittest.mock import MagicMock, mock_open, patch
import subprocess

import sys
import platform
# --- MOCK FOR WINDOWS TESTING ---
# SQANTI3 and some pandas versions check for POSIX platforms.
if sys.platform == 'win32':
    sys.platform = 'linux'
    import collections
    os.uname = lambda: collections.namedtuple('Uname', ['sysname', 'nodename', 'release', 'version', 'machine'])('Linux', 'localhost', '5.10.0', '1', 'x86_64')
    # Mock pysam so tests can run without the C extension compilation on win32
    sys.modules['pysam'] = MagicMock()
# --------------------------------

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
from classification_enrichment import annotate_with_ujc_hash, annotate_with_cell_metadata, annotate_with_sample_support
from qc_reports import generate_report, generate_multisample_report
from qc_pipeline import main as pipeline_main
from cell_metrics import calculate_metrics_per_cell
from sc_clustering import prepare_anndata, run_clustering_analysis


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
            self.genename = False
            self.short_reads = None
            self.SR_bam = None
            self.novel_gene_prefix = None
            self.aligner_choice = "minimap2"
            self.gmap_index = None
            self.sites = "ATAC,GCAG,GTAG"
            self.include_ORF = False
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
            self.min_cov = 1
            self.ratio_TSS_threshold = 2.0
            self.ref_cov_min_pct = 45.0

            # Clustering-related attributes
            self.run_clustering = False
            self.normalization = 'log1p'
            self.n_neighbors = 15
            self.n_pc = 30
            self.resolution = 0.5
            self.n_top_genes = 2000
            self.clustering_method = 'leiden'
            self.n_clusters = 10

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
    """In fusion mode without --orf_input and including ORF, raise error."""
    design_df = pd.DataFrame({
        'sampleID': ['sample1'],
        'file_acc': ['file1']
    })
    # Create minimal environment so the function reaches fusion check
    mock_args.is_fusion = True
    mock_args.include_ORF = True
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


def test_annotate_with_sample_support(mock_args, tmpdir):
    """recount3 sidecar sample-support folds into classification correctly.

    Exercises the real (unmocked) coordinate join: sidecar start/end ==
    junctions.txt genomic_start_coord/genomic_end_coord (zero offset), the
    per-isoform 'weakest junction' reduction, junctions absent from the sidecar
    counting as 0-sample support, and monoexon isoforms staying NA.
    """
    mock_args.mode = "isoforms"

    # Coverage file + its sidecar sit next to each other; only the sidecar is read.
    coverage = tmpdir.join("blood_recount3.SJ.out.tab")
    coverage.write("dummy\n")
    tmpdir.join("blood_recount3.SJ.out.tab.sample_support.tsv").write(
        "chrom\tstart\tend\tstrand\tn_samples\tpct_samples\n"
        "chr1\t100\t200\t+\t40\t80.0\n"
        "chr1\t300\t400\t+\t10\t20.0\n"
        "chr1\t500\t600\t+\t30\t60.0\n"
    )

    # SQANTI3 outputs live under out_dir/<file_acc>/<sampleID>_*.txt
    sample_dir = tmpdir.join("output_dir", "f1").ensure(dir=True)
    sample_dir.join("s1_junctions.txt").write(
        "isoform\tchrom\tstrand\tgenomic_start_coord\tgenomic_end_coord\tsample_with_cov\n"
        "iso1\tchr1\t+\t100\t200\t1\n"
        "iso1\tchr1\t+\t300\t400\t1\n"   # iso1 weakest: 10 samples
        "iso2\tchr1\t+\t100\t200\t1\n"
        "iso2\tchr1\t+\t500\t600\t1\n"   # iso2 weakest: 30 samples
        "iso3\tchr1\t+\t700\t800\t1\n"   # absent from sidecar -> 0 samples
    )
    class_file = sample_dir.join("s1_classification.txt")
    class_file.write(
        "isoform\tstructural_category\tmin_sample_cov\tmin_cov\n"
        "iso1\tFSM\t1\t5\n"
        "iso2\tNIC\t1\t3\n"
        "iso3\tNNC\t1\t2\n"
        "iso4\tFSM\tNA\tNA\n"          # monoexon: no junctions -> NA
    )

    design_df = pd.DataFrame({
        "sampleID": ["s1"], "file_acc": ["f1"], "coverage": [str(coverage)],
    })

    annotate_with_sample_support(mock_args, design_df)

    # --- junctions.txt: sample_with_cov overwritten, sample_pct added --------
    jout = pd.read_csv(str(sample_dir.join("s1_junctions.txt")), sep="\t", dtype=str,
                       keep_default_na=False)
    assert "sample_with_cov" in jout.columns   # overwritten, not a new column
    assert "sample_pct" in jout.columns        # genuinely new
    assert "sample_cov" not in jout.columns    # no redundant alias
    row_100 = jout[(jout["genomic_start_coord"] == "100") & (jout["isoform"] == "iso1")]
    assert row_100["sample_with_cov"].iloc[0] == "40"
    assert row_100["sample_pct"].iloc[0] == "80.0000"
    row_abs = jout[jout["genomic_start_coord"] == "700"]
    assert row_abs["sample_with_cov"].iloc[0] == "0"

    # --- classification.txt gets per-isoform min columns --------------------
    # keep_default_na=False so the literal 'NA' strings the pipeline writes stay strings
    out = pd.read_csv(str(class_file), sep="\t", dtype=str,
                      keep_default_na=False).set_index("isoform")
    assert "min_sample_pct" in out.columns
    assert out.loc["iso1", "min_sample_cov"] == "10"     # min(40, 10)
    assert out.loc["iso1", "min_sample_pct"] == "20.0000"
    assert out.loc["iso2", "min_sample_cov"] == "30"     # min(40, 30)
    assert out.loc["iso2", "min_sample_pct"] == "60.0000"
    assert out.loc["iso3", "min_sample_cov"] == "0"      # junction not in sidecar
    assert out.loc["iso3", "min_sample_pct"] == "0.0000"
    assert out.loc["iso4", "min_sample_cov"] == "NA"     # monoexon
    assert out.loc["iso4", "min_sample_pct"] == "NA"


def test_annotate_with_sample_support_noop_without_sidecar(mock_args, tmpdir):
    """No sidecar next to the coverage file -> classification left untouched."""
    mock_args.mode = "isoforms"
    coverage = tmpdir.join("star.SJ.out.tab")
    coverage.write("dummy\n")  # no .sample_support.tsv companion

    sample_dir = tmpdir.join("output_dir", "f1").ensure(dir=True)
    class_file = sample_dir.join("s1_classification.txt")
    original = "isoform\tmin_sample_cov\niso1\t1\n"
    class_file.write(original)

    design_df = pd.DataFrame({
        "sampleID": ["s1"], "file_acc": ["f1"], "coverage": [str(coverage)],
    })
    annotate_with_sample_support(mock_args, design_df)

    assert class_file.read() == original  # unchanged, no min_sample_pct added


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
@patch('qc_pipeline.annotate_with_sample_support')
@patch('qc_pipeline.calculate_metrics_per_cell')
@patch('qc_pipeline.generate_report')
def test_pipeline_main_smoke(
        mock_generate_report,
        mock_calc,
        mock_sample_support,
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
@patch('qc_reports.reportAssetsPath', 'utilities')
def test_generate_report(mock_isfile, mock_run, mock_args, capsys):
    """Test report generation command."""
    mock_args.mode = "reads"
    df = pd.DataFrame({"sampleID": ["sample1"], "file_acc": ["file1"]})

    prefix = f"{mock_args.out_dir}/file1/sample1"
    class_file = f"{prefix}_classification.txt"
    junc_file = f"{prefix}_junctions.txt"
    cell_summary = f"{prefix}_SQANTI_cell_summary.txt.gz"
    clustering_file = os.path.join(os.path.dirname(prefix), "clustering", "umap_results.csv")

    # class, junc, and cell_summary exist; clustering does NOT
    def _isfile(path):
        return path != clustering_file
    mock_isfile.side_effect = _isfile

    generate_report(mock_args, df)

    captured = capsys.readouterr()
    assert "SQANTI3 report generated for file1" in captured.out

    expected_cmd = (
        f"Rscript utilities/SQANTI-sc_report.R "
        f"\"{class_file}\" \"{junc_file}\" {mock_args.report} \"{prefix}\" "
        f"{mock_args.mode} --cell_summary {cell_summary} "
        f"--refGTF \"{mock_args.refGTF}\""
    ).strip()

    actual_cmd = " ".join(mock_run.call_args[0][0].split())
    assert actual_cmd == expected_cmd


@patch('qc_reports.subprocess.run')
@patch('qc_reports.os.path.isfile')
@patch('qc_reports.reportAssetsPath', 'utilities')
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
        f"Rscript utilities/SQANTI-sc_report.R "
        f"\"{class_file}\" \"{junc_file}\" {mock_args.report} \"{prefix}\" "
        f"{mock_args.mode} --refGTF \"{mock_args.refGTF}\""
    ).strip()
    assert actual_cmd == expected_cmd


@patch('qc_reports.subprocess.run')
@patch('qc_reports.os.path.isfile')
@patch('qc_reports.reportAssetsPath', 'utilities')
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


# ==============================================================================
# Clustering Tests
# ==============================================================================

def test_prepare_anndata_reads_mode(mock_args, tmpdir):
    """Test prepare_anndata creates an AnnData with correct shape in reads mode."""
    mock_args.mode = 'reads'
    out_dir = str(tmpdir.join('output_dir'))
    os.makedirs(out_dir, exist_ok=True)
    mock_args.out_dir = out_dir
    file_acc = 'f1'
    sampleID = 's1'
    sample_dir = os.path.join(out_dir, file_acc)
    os.makedirs(sample_dir, exist_ok=True)
    prefix = os.path.join(sample_dir, sampleID)

    cls = pd.DataFrame({
        'isoform': ['i1', 'i2', 'i3'],
        'CB': ['CELL1', 'CELL1', 'CELL2'],
        'associated_gene': ['geneA', 'geneB', 'geneA'],
    })
    cls.to_csv(f"{prefix}_classification.txt", sep='\t', index=False)

    row = pd.Series({'sampleID': sampleID, 'file_acc': file_acc})
    adata = prepare_anndata(mock_args, row)

    assert adata is not None
    assert adata.shape[0] == 2   # 2 cells: CELL1, CELL2
    assert adata.shape[1] == 2   # 2 genes: geneA, geneB


def test_prepare_anndata_isoforms_mode(mock_args, tmpdir):
    """Test prepare_anndata correctly expands comma-separated CB and FL in isoforms mode."""
    mock_args.mode = 'isoforms'
    out_dir = str(tmpdir.join('output_dir'))
    os.makedirs(out_dir, exist_ok=True)
    mock_args.out_dir = out_dir
    file_acc = 'f1'
    sampleID = 's1'
    sample_dir = os.path.join(out_dir, file_acc)
    os.makedirs(sample_dir, exist_ok=True)
    prefix = os.path.join(sample_dir, sampleID)

    # Isoform i1 seen in two cells with 3 and 5 counts respectively
    cls = pd.DataFrame({
        'isoform': ['i1', 'i2'],
        'CB': ['CELL1,CELL2', 'CELL2'],
        'FL': ['3,5', '2'],
        'associated_gene': ['geneA', 'geneA'],
    })
    cls.to_csv(f"{prefix}_classification.txt", sep='\t', index=False)

    row = pd.Series({'sampleID': sampleID, 'file_acc': file_acc})
    adata = prepare_anndata(mock_args, row)

    assert adata is not None
    assert adata.shape[0] == 2   # 2 cells
    # CELL1: geneA=3, CELL2: geneA=5+2=7
    cell2_geneA = adata[adata.obs_names == 'CELL2', 'geneA'].X
    assert cell2_geneA[0, 0] == 7


def test_prepare_anndata_missing_file(mock_args, tmpdir, capsys):
    """Test prepare_anndata returns None when classification file is missing."""
    mock_args.mode = 'reads'
    mock_args.out_dir = str(tmpdir.join('output_dir'))
    os.makedirs(mock_args.out_dir, exist_ok=True)

    row = pd.Series({'sampleID': 'nonexist', 'file_acc': 'nonexist'})
    adata = prepare_anndata(mock_args, row)

    assert adata is None
    captured = capsys.readouterr()
    assert "Classification file not found" in captured.err


def test_prepare_anndata_no_valid_barcodes(mock_args, tmpdir, capsys):
    """Test prepare_anndata returns None when all CBs are empty/NA."""
    mock_args.mode = 'reads'
    out_dir = str(tmpdir.join('output_dir'))
    os.makedirs(out_dir, exist_ok=True)
    mock_args.out_dir = out_dir
    file_acc = 'f1'
    sampleID = 's1'
    sample_dir = os.path.join(out_dir, file_acc)
    os.makedirs(sample_dir, exist_ok=True)
    prefix = os.path.join(sample_dir, sampleID)

    cls = pd.DataFrame({
        'isoform': ['i1', 'i2'],
        'CB': ['', 'NA'],
        'associated_gene': ['geneA', 'geneB'],
    })
    cls.to_csv(f"{prefix}_classification.txt", sep='\t', index=False)

    row = pd.Series({'sampleID': sampleID, 'file_acc': file_acc})
    adata = prepare_anndata(mock_args, row)

    assert adata is None
    captured = capsys.readouterr()
    assert "No valid cell barcodes" in captured.out


@patch('sc_clustering.sc')
def test_run_clustering_analysis_smoke(mock_scanpy, mock_args, tmpdir):
    """Smoke test: run_clustering_analysis calls scanpy pipeline and writes CSV."""
    mock_args.mode = 'reads'
    out_dir = str(tmpdir.join('output_dir'))
    os.makedirs(out_dir, exist_ok=True)
    mock_args.out_dir = out_dir
    file_acc = 'f1'
    sampleID = 's1'
    sample_dir = os.path.join(out_dir, file_acc)
    os.makedirs(sample_dir, exist_ok=True)
    prefix = os.path.join(sample_dir, sampleID)

    # Create classification file with enough data
    cls = pd.DataFrame({
        'isoform': [f'i{i}' for i in range(20)],
        'CB': [f'CELL{i % 5}' for i in range(20)],
        'associated_gene': [f'gene{i % 3}' for i in range(20)],
    })
    cls.to_csv(f"{prefix}_classification.txt", sep='\t', index=False)

    # Mock scanpy to avoid requiring the full pipeline
    import numpy as np
    mock_adata = MagicMock()
    mock_adata.__getitem__.return_value = mock_adata
    mock_adata.n_vars = 3
    mock_adata.obs_names = [f'CELL{i}' for i in range(5)]
    mock_adata.obsm = {'X_umap': np.random.rand(5, 2), 'X_pca': np.random.rand(5, 3)}
    mock_adata.obs = pd.DataFrame({'leiden': ['0', '1', '0', '1', '0']}, index=mock_adata.obs_names)

    with patch('sc_clustering.prepare_anndata', return_value=mock_adata):
        row = pd.Series({'sampleID': sampleID, 'file_acc': file_acc})
        run_clustering_analysis(mock_args, row)

    # Check that the output CSV was written
    out_csv = os.path.join(sample_dir, 'clustering', 'umap_results.csv')
    assert os.path.isfile(out_csv)
    result_df = pd.read_csv(out_csv)
    assert 'Barcode' in result_df.columns
    assert 'Cluster' in result_df.columns
    assert 'UMAP_1' in result_df.columns
    assert 'UMAP_2' in result_df.columns


@patch('qc_reports.subprocess.run')
@patch('qc_reports.os.path.isfile')
@patch('qc_reports.reportAssetsPath', 'utilities')
def test_generate_report_with_clustering_flag(mock_isfile, mock_run, mock_args):
    """When clustering results exist, the report command should include --clustering."""
    mock_args.mode = "isoforms"
    df = pd.DataFrame({"sampleID": ["sample1"], "file_acc": ["file1"]})

    prefix = f"{mock_args.out_dir}/file1/sample1"
    clustering_file = os.path.join(os.path.dirname(prefix), "clustering", "umap_results.csv")

    def _isfile_side_effect(path):
        return True  # All files exist including clustering

    mock_isfile.side_effect = _isfile_side_effect

    generate_report(mock_args, df)

    actual_cmd = " ".join(mock_run.call_args[0][0].split())
    assert "--clustering" in actual_cmd
    assert "umap_results.csv" in actual_cmd


@patch('qc_reports.subprocess.run')
@patch('qc_reports.os.path.isfile')
@patch('qc_reports.reportAssetsPath', 'utilities')
def test_generate_report_with_optional_flags(mock_isfile, mock_run, mock_args):
    """Flags like --include_ORF, --CAGE_peak, --polyA_motif_list should appear in command."""
    mock_isfile.return_value = True
    mock_args.mode = "reads"
    mock_args.include_ORF = True
    mock_args.CAGE_peak = True
    mock_args.polyA_motif_list = True
    df = pd.DataFrame({"sampleID": ["sample1"], "file_acc": ["file1"]})

    generate_report(mock_args, df)

    actual_cmd = " ".join(mock_run.call_args[0][0].split())
    assert "--include_ORF" in actual_cmd
    assert "--CAGE_peak" in actual_cmd
    assert "--polyA_motif_list" in actual_cmd


@patch('qc_pipeline.build_parser')
@patch('qc_pipeline.run_sqanti3_qc')
@patch('qc_pipeline.annotate_with_ujc_hash')
@patch('qc_pipeline.annotate_with_cell_metadata')
@patch('qc_pipeline.annotate_with_sample_support')
@patch('qc_pipeline.calculate_metrics_per_cell')
@patch('qc_pipeline.generate_report')
@patch('qc_pipeline.run_clustering_analysis')
def test_pipeline_main_with_clustering(
        mock_clustering,
        mock_generate_report,
        mock_calc,
        mock_sample_support,
        mock_add,
        mock_hash,
        mock_get,
        mock_build_parser,
        tmpdir):
    """Pipeline main calls run_clustering_analysis when --run_clustering is set."""
    class DummyArgs:
        def __init__(self):
            self.inDESIGN = str(tmpdir.join('design.csv'))
            self.input_dir = str(tmpdir)
            self.out_dir = str(tmpdir)
            self.refGTF = __file__
            self.refFasta = __file__
            self.mode = 'isoforms'
            self.report = 'html'
            self.SKIPHASH = True
            self.write_per_cell_outputs = False
            self.log_level = 'INFO'
            self.run_clustering = True  # <-- enabled

    tmpdir.join('design.csv').write(
        'sampleID,file_acc,cell_association\ns1,f1,assoc.txt\n'
    )

    class DummyParser:
        def parse_args(self_inner):
            return DummyArgs()

    mock_build_parser.return_value = DummyParser()
    pipeline_main()

    assert mock_clustering.called
    assert mock_generate_report.called

# ==============================================================================
# FL Weighting Tests (isoforms mode)
# ==============================================================================

class TestFLWeightingIsoformsMode:
    """
    Tests verifying that in isoforms mode, transcript counts and junction counts
    are correctly weighted by the FL (full-length read count) values per cell.

    Core principle
    --------------
    Each row in the classification file is a unique isoform model.  The FL
    column contains comma-separated counts representing how many actual
    transcripts of that isoform were observed in each cell barcode listed in
    the CB column.  All per-cell metrics (structural-category percentages,
    junction counts, etc.) must be computed by *summing these FL counts*, NOT
    by simply counting isoform model rows.

    For junctions: each junction of a given isoform inherits the same FL counts
    as the isoform itself (one junction row per intron per isoform model, so
    3 exons → 2 junction rows, each with FL counts equal to the isoform's FL).
    """

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _cls_row(self, isoform, cb, fl, structural_category,
                 associated_gene="geneA", exons=2):
        """Return one classification row dict (isoforms mode)."""
        return {
            "isoform": isoform,
            "CB": cb,
            "FL": fl,
            "structural_category": structural_category,
            "associated_gene": associated_gene,
            "associated_transcript": "txA",
            "exons": exons,
            "length": 500,
            "ref_length": 600,
            "chrom": "chr1",
            "subcategory": "reference_match",
            "all_canonical": "True",
            "RTS_stage": "False",
            "predicted_NMD": "False",
            "within_CAGE_peak": "False",
            "polyA_motif_found": "False",
            "perc_A_downstream_TTS": "0",
            "diff_to_gene_TSS": "0",
            "coding": "coding",
            "min_cov": "0",
            "ratio_TSS": "0",
        }

    def _junc_row(self, isoform, junction_category="known", canonical="canonical",
                  rts=False):
        """Return one junction row dict (isoform key only; CB/FL come from cls join)."""
        return {
            "isoform": isoform,
            "junction_category": junction_category,
            "canonical": canonical,
            "RTS_junction": str(rts),
            "junction_number": "1",
            "chrom": "chr1",
            "strand": "+",
            "genomic_start_coord": "1000",
            "genomic_end_coord": "2000",
        }

    def _run(self, mock_args, tmpdir, cls_rows, junc_rows=None):
        """Write files, run calculate_metrics_per_cell, return cell summary DataFrame."""
        out_dir = str(tmpdir.join("output"))
        file_acc = "f1"
        sampleID = "s1"
        sample_dir = os.path.join(out_dir, file_acc)
        os.makedirs(sample_dir, exist_ok=True)
        prefix = os.path.join(sample_dir, sampleID)

        mock_args.mode = "isoforms"
        mock_args.out_dir = out_dir

        pd.DataFrame(cls_rows).to_csv(
            f"{prefix}_classification.txt", sep="\t", index=False
        )

        if junc_rows:
            pd.DataFrame(junc_rows).to_csv(
                f"{prefix}_junctions.txt", sep="\t", index=False
            )
        else:
            # Empty but valid junctions file
            pd.DataFrame(columns=["isoform"]).to_csv(
                f"{prefix}_junctions.txt", sep="\t", index=False
            )

        design_df = pd.DataFrame({"sampleID": [sampleID], "file_acc": [file_acc]})
        calculate_metrics_per_cell(mock_args, design_df)

        summary_path = f"{prefix}_SQANTI_cell_summary.txt.gz"
        assert os.path.isfile(summary_path), "Cell summary was not created"
        return pd.read_csv(summary_path, sep="\t", compression="gzip")

    # ------------------------------------------------------------------
    # Test 1 — total_reads reflects sum(FL), not number of isoform rows
    # ------------------------------------------------------------------

    def test_total_reads_is_fl_sum_not_row_count(self, mock_args, tmpdir):
        """
        A cell with 3 FSM isoforms each having FL=5 should have total_reads=15,
        not total_reads=3 (the number of isoform model rows).
        """
        cls_rows = [
            self._cls_row("iso1", "CB1", "5", "full-splice_match"),
            self._cls_row("iso2", "CB1", "5", "full-splice_match"),
            self._cls_row("iso3", "CB1", "5", "full-splice_match"),
        ]
        summary = self._run(mock_args, tmpdir, cls_rows)
        cb1 = summary[summary["CB"] == "CB1"].iloc[0]
        assert cb1["Transcripts_in_cell"] == 15, (
            f"Expected Transcripts_in_cell=15 (sum of FL counts), got {cb1['Transcripts_in_cell']}"
        )

    # ------------------------------------------------------------------
    # Test 2 — structural category % is FL-weighted
    # ------------------------------------------------------------------

    def test_structural_category_percentage_fl_weighted(self, mock_args, tmpdir):
        """
        CB1 has:
          iso1 (FSM, FL=1)  → 1 FSM transcript
          iso2 (NIC, FL=9)  → 9 NIC transcripts

        FL-weighted result  : FSM_prop = 10 %, NIC_prop = 90 %
        Naive row-count     : FSM_prop = 50 %, NIC_prop = 50 %  (wrong)
        """
        cls_rows = [
            self._cls_row("iso1", "CB1", "1", "full-splice_match"),
            self._cls_row("iso2", "CB1", "9", "novel_in_catalog"),
        ]
        summary = self._run(mock_args, tmpdir, cls_rows)
        cb1 = summary[summary["CB"] == "CB1"].iloc[0]

        assert abs(cb1["FSM_prop"] - 10.0) < 0.01, (
            f"Expected FSM_prop≈10 % (FL-weighted), got {cb1['FSM_prop']}"
        )
        assert abs(cb1["NIC_prop"] - 90.0) < 0.01, (
            f"Expected NIC_prop≈90 % (FL-weighted), got {cb1['NIC_prop']}"
        )

    # ------------------------------------------------------------------
    # Test 3 — junction counts are FL-weighted via isoform-key join
    # ------------------------------------------------------------------

    def test_junction_counts_fl_weighted_via_isoform_join(self, mock_args, tmpdir):
        """
        iso1 (FL=1 in CB1): 1 known_canonical junction
        iso2 (FL=9 in CB1): 1 novel_not_in_catalog junction

        Each junction inherits the isoform's FL count for that cell, so:
          Known_canonical_junctions      = 1   (1 × FL=1)
          Novel_non_canonical_junctions  = 9   (1 × FL=9)   ← novel_not_in_catalog + canonical=False → novel_non_canonical
          Known_canonical_junctions_prop = 1/(1+9) × 100 = 10 %

        In the junctions file there is no CB column — the CB/FL association is
        derived from the classification file via the isoform key, exactly as
        SQANTI-sc does it after the fix we applied.
        """
        cls_rows = [
            self._cls_row("iso1", "CB1", "1", "full-splice_match", exons=2),
            self._cls_row("iso2", "CB1", "9", "novel_not_in_catalog", exons=2),
        ]
        junc_rows = [
            self._junc_row("iso1", junction_category="known",   canonical="canonical"),
            self._junc_row("iso2", junction_category="novel",   canonical="non_canonical"),
        ]
        summary = self._run(mock_args, tmpdir, cls_rows, junc_rows)
        cb1 = summary[summary["CB"] == "CB1"].iloc[0]

        assert cb1["Known_canonical_junctions"] == 1, (
            f"Expected 1 known_canonical junction (FL=1), got {cb1['Known_canonical_junctions']}"
        )
        assert cb1["Novel_non_canonical_junctions"] == 9, (
            f"Expected 9 novel_non_canonical junctions (FL=9), got {cb1['Novel_non_canonical_junctions']}"
        )
        assert abs(cb1["Known_canonical_junctions_prop"] - 10.0) < 0.01, (
            f"Expected Known_canonical_junctions_prop≈10 %, got {cb1['Known_canonical_junctions_prop']}"
        )

    # ------------------------------------------------------------------
    # Test 4 — multiple junctions per isoform each inherit FL
    # ------------------------------------------------------------------

    def test_multi_junction_isoform_each_junction_gets_fl_count(self, mock_args, tmpdir):
        """
        iso1 has 3 exons → 2 junctions, FL=5 in CB1.
        Both junctions are known_canonical, so:
          Known_canonical_junctions = 2 × 5 = 10
        """
        cls_rows = [
            self._cls_row("iso1", "CB1", "5", "full-splice_match", exons=3),
        ]
        junc_rows = [
            self._junc_row("iso1", junction_category="known", canonical="canonical"),
            self._junc_row("iso1", junction_category="known", canonical="canonical"),
        ]
        summary = self._run(mock_args, tmpdir, cls_rows, junc_rows)
        cb1 = summary[summary["CB"] == "CB1"].iloc[0]

        assert cb1["Known_canonical_junctions"] == 10, (
            f"Expected Known_canonical_junctions=10 (2 junctions × FL=5), "
            f"got {cb1['Known_canonical_junctions']}"
        )

    # ------------------------------------------------------------------
    # Test 5 — same isoform appearing in multiple cells is weighted
    #           independently per cell
    # ------------------------------------------------------------------

    def test_fl_weighting_is_per_cell_independent(self, mock_args, tmpdir):
        """
        iso1 appears in CB1 (FL=3) and CB2 (FL=7) via comma-separated CB/FL.
        iso2 appears only in CB2 (FL=2).

        All isoforms are FSM; expected results:
          CB1: total_reads=3,  FSM_prop=100 %
          CB2: total_reads=9,  FSM_prop=100 %   (iso1 FL=7 + iso2 FL=2)
        """
        cls_rows = [
            self._cls_row("iso1", "CB1,CB2", "3,7", "full-splice_match"),
            self._cls_row("iso2", "CB2",     "2",   "full-splice_match"),
        ]
        summary = self._run(mock_args, tmpdir, cls_rows)

        cb1 = summary[summary["CB"] == "CB1"].iloc[0]
        cb2 = summary[summary["CB"] == "CB2"].iloc[0]

        assert cb1["Transcripts_in_cell"] == 3, (
            f"CB1 Transcripts_in_cell should be 3 (FL for CB1), got {cb1['Transcripts_in_cell']}"
        )
        assert cb2["Transcripts_in_cell"] == 9, (
            f"CB2 Transcripts_in_cell should be 9 (FL=7 + FL=2), got {cb2['Transcripts_in_cell']}"
        )
        assert abs(cb1["FSM_prop"] - 100.0) < 0.01
        assert abs(cb2["FSM_prop"] - 100.0) < 0.01

    # ------------------------------------------------------------------
    # Test 5b — annotated/novel READ counts are FL-weighted, split by gene
    #           annotation (associated_gene 'novel' prefix, not category)
    # ------------------------------------------------------------------

    def test_annotated_novel_reads_are_fl_weighted(self, mock_args, tmpdir):
        """
        Annotated_genes_transcripts / Novel_genes_transcripts are per-cell FL sums split by
        whether associated_gene starts with 'novel'.
          CB1: geneA FL=3 (annotated) + novelGene FL=1 (novel)
          CB2: geneA FL=10 + geneA FL=4 (annotated) + novelGene FL=6 (novel)
        """
        cls_rows = [
            self._cls_row("iso1", "CB1,CB2", "3,10", "full-splice_match", associated_gene="geneA"),
            self._cls_row("iso2", "CB2",     "4",    "full-splice_match", associated_gene="geneA"),
            self._cls_row("iso3", "CB1,CB2", "1,6",  "novel_in_catalog",  associated_gene="novelGene_1"),
        ]
        summary = self._run(mock_args, tmpdir, cls_rows)
        cb1 = summary[summary["CB"] == "CB1"].iloc[0]
        cb2 = summary[summary["CB"] == "CB2"].iloc[0]

        assert cb1["Annotated_genes_transcripts"] == 3
        assert cb1["Novel_genes_transcripts"] == 1
        assert cb2["Annotated_genes_transcripts"] == 14   # 10 + 4
        assert cb2["Novel_genes_transcripts"] == 6

    # ------------------------------------------------------------------
    # Test 5c — gene abundance bins store gene COUNTS in the 11 gapless
    #           half-open bins (0,1]..(9,10],(10,inf), FL-summed per gene
    # ------------------------------------------------------------------

    # All 11 abundance-bin column names, in order, for isoforms mode. In isoforms mode
    # the *_gene_reads_bin_* columns are renamed *_gene_transcripts_bin_* on output.
    ANNO_TX_BINS = [f"anno_gene_transcripts_bin_{i}" for i in range(1, 11)] + [
        "anno_gene_transcripts_bin_10plus"
    ]
    ANNO_ISO_BINS = [f"anno_gene_isoforms_bin_{i}" for i in range(1, 11)] + [
        "anno_gene_isoforms_bin_10plus"
    ]

    def test_gene_read_count_bins_are_fl_weighted_counts(self, mock_args, tmpdir):
        """
        Each bin column counts annotated genes whose per-cell FL sum falls in that
        bin. Bins are the 11 gapless half-open intervals (i-1, i] for i=1..10, then
        >10.
          CB1: geneA=3+2=5   -> bin_5
          CB2: geneA=10      -> bin_10 ;  geneB=4 -> bin_4
        Novel-gene bins are intentionally NOT computed (no plot consumes them).
        """
        cls_rows = [
            self._cls_row("iso1", "CB1,CB2", "3,10", "full-splice_match", associated_gene="geneA"),
            self._cls_row("iso2", "CB1",     "2",    "full-splice_match", associated_gene="geneA"),
            self._cls_row("iso3", "CB2",     "4",    "full-splice_match", associated_gene="geneB"),
            self._cls_row("iso4", "CB1,CB2", "1,6",  "novel_in_catalog",  associated_gene="novelGene_1"),
        ]
        summary = self._run(mock_args, tmpdir, cls_rows)
        cb1 = summary[summary["CB"] == "CB1"].iloc[0]
        cb2 = summary[summary["CB"] == "CB2"].iloc[0]

        # CB1 annotated: geneA = 3+2 = 5 -> lands in (4,5] and nowhere else
        assert cb1["anno_gene_transcripts_bin_5"] == 1
        assert sum(cb1[b] for b in self.ANNO_TX_BINS) == 1, "CB1 has exactly 1 annotated gene"

        # CB2 annotated: geneA = 10 -> (9,10] ; geneB = 4 -> (3,4]
        assert cb2["anno_gene_transcripts_bin_10"] == 1
        assert cb2["anno_gene_transcripts_bin_4"] == 1
        assert cb2["anno_gene_transcripts_bin_10plus"] == 0, "10 is <=10, must not spill into >10"
        assert sum(cb2[b] for b in self.ANNO_TX_BINS) == 2, "CB2 has exactly 2 annotated genes"

        # Novel-gene bin columns are not produced at all.
        assert not [c for c in summary.columns if c.startswith("novel_gene_transcripts_bin_")]

    # ------------------------------------------------------------------
    # Test 5d — bins are GAPLESS: float/EM abundances never fall through
    # ------------------------------------------------------------------

    def test_gene_abundance_bins_are_gapless_for_float_abundances(self, mock_args, tmpdir):
        """
        Regression guard for the old 4-bin scheme (1 / 2-5 / 6-9 / >=10), which had
        gaps at (1,2), (5,6) and (9,10) that silently DROPPED non-integer abundances.
        EM-based tools (Isosceles, bambu) and 1/n ambiguous splitting produce
        fractional FL, so every value must still land in exactly one bin.

        CB1: geneA=1.5 -> (1,2]=bin_2 ; geneB=5.5 -> (5,6]=bin_6
             geneC=9.5 -> (9,10]=bin_10 ; geneD=10.5 -> >10=bin_10plus
        """
        cls_rows = [
            self._cls_row("iso1", "CB1", "1.5",  "full-splice_match", associated_gene="geneA"),
            self._cls_row("iso2", "CB1", "5.5",  "full-splice_match", associated_gene="geneB"),
            self._cls_row("iso3", "CB1", "9.5",  "full-splice_match", associated_gene="geneC"),
            self._cls_row("iso4", "CB1", "10.5", "full-splice_match", associated_gene="geneD"),
        ]
        summary = self._run(mock_args, tmpdir, cls_rows)
        cb1 = summary[summary["CB"] == "CB1"].iloc[0]

        # Each fractional abundance lands in its half-open bin -- these are exactly
        # the values the old binning dropped.
        assert cb1["anno_gene_transcripts_bin_2"] == 1, "1.5 must land in (1,2]"
        assert cb1["anno_gene_transcripts_bin_6"] == 1, "5.5 must land in (5,6]"
        assert cb1["anno_gene_transcripts_bin_10"] == 1, "9.5 must land in (9,10]"
        assert cb1["anno_gene_transcripts_bin_10plus"] == 1, "10.5 must land in >10"

        # Nothing dropped: all 4 genes accounted for across the 11 bins.
        assert sum(cb1[b] for b in self.ANNO_TX_BINS) == 4

    # ------------------------------------------------------------------
    # Test 5e — unique-isoform bins count DISTINCT MODELS, not abundance
    # ------------------------------------------------------------------

    def test_unique_isoform_bins_count_distinct_models_not_abundance(self, mock_args, tmpdir):
        """
        anno_gene_isoforms_bin_* must count distinct transcript MODELS per gene per
        cell, which is a different quantity from the FL-sum abundance bins. This is
        the distinction the old report collapsed by re-plotting abundance under the
        "Unique Isoform Count" title.

        Both genes have abundance 2 in CB1, but differ in model count:
          geneA: ONE model with FL=2      -> abundance 2, 1 distinct isoform
          geneB: TWO models with FL=1 each -> abundance 2, 2 distinct isoforms
        """
        cls_rows = [
            self._cls_row("iso1", "CB1", "2", "full-splice_match", associated_gene="geneA"),
            self._cls_row("iso2", "CB1", "1", "full-splice_match", associated_gene="geneB"),
            self._cls_row("iso3", "CB1", "1", "full-splice_match", associated_gene="geneB"),
        ]
        summary = self._run(mock_args, tmpdir, cls_rows)
        cb1 = summary[summary["CB"] == "CB1"].iloc[0]

        # Abundance: both genes sum to 2 -> both in bin_2.
        assert cb1["anno_gene_transcripts_bin_2"] == 2
        assert sum(cb1[b] for b in self.ANNO_TX_BINS) == 2

        # Distinct models: geneA has 1, geneB has 2 -> the two families DIVERGE.
        assert cb1["anno_gene_isoforms_bin_1"] == 1, "geneA: one model carrying FL=2"
        assert cb1["anno_gene_isoforms_bin_2"] == 1, "geneB: two models of FL=1"
        assert sum(cb1[b] for b in self.ANNO_ISO_BINS) == 2

    # ------------------------------------------------------------------
    # Test 6 — reads mode uses 1 count per row regardless of FL column
    # ------------------------------------------------------------------

    def test_reads_mode_uses_one_count_per_row(self, mock_args, tmpdir):
        """
        In reads mode each classification row is a single read (CB is a single
        barcode, FL column is absent or irrelevant).  total_reads must equal
        the number of rows, not any FL value.
        """
        out_dir = str(tmpdir.join("output_reads"))
        file_acc = "f1"
        sampleID = "s1"
        sample_dir = os.path.join(out_dir, file_acc)
        os.makedirs(sample_dir, exist_ok=True)
        prefix = os.path.join(sample_dir, sampleID)

        mock_args.mode = "reads"
        mock_args.out_dir = out_dir

        # 3 reads all from CB1 — no FL column
        cls = pd.DataFrame([
            {"isoform": "r1", "CB": "CB1", "structural_category": "full-splice_match",
             "associated_gene": "geneA", "associated_transcript": "txA",
             "exons": 2, "length": 500, "ref_length": 600, "chrom": "chr1"},
            {"isoform": "r2", "CB": "CB1", "structural_category": "full-splice_match",
             "associated_gene": "geneA", "associated_transcript": "txA",
             "exons": 2, "length": 500, "ref_length": 600, "chrom": "chr1"},
            {"isoform": "r3", "CB": "CB1", "structural_category": "novel_in_catalog",
             "associated_gene": "geneA", "associated_transcript": "txA",
             "exons": 2, "length": 500, "ref_length": 600, "chrom": "chr1"},
        ])
        cls.to_csv(f"{prefix}_classification.txt", sep="\t", index=False)
        pd.DataFrame(columns=["isoform"]).to_csv(
            f"{prefix}_junctions.txt", sep="\t", index=False
        )

        design_df = pd.DataFrame({"sampleID": [sampleID], "file_acc": [file_acc]})
        calculate_metrics_per_cell(mock_args, design_df)

        summary = pd.read_csv(
            f"{prefix}_SQANTI_cell_summary.txt.gz", sep="\t", compression="gzip"
        )
        cb1 = summary[summary["CB"] == "CB1"].iloc[0]

        # 3 rows → 3 reads (count=1 per row, no FL weighting)
        assert cb1["Reads_in_cell"] == 3, (
            f"Reads mode: expected Reads_in_cell=3 (one per row), got {cb1['Reads_in_cell']}"
        )
        assert abs(cb1["FSM_prop"] - (2 / 3 * 100)) < 0.1, (
            f"Reads mode: expected FSM_prop≈66.7 %, got {cb1['FSM_prop']}"
        )

        # Gene-annotation columns are gene-based, not category-based: all 3 reads
        # belong to annotated geneA (even the novel-category read), 3 reads -> (2,3].
        assert cb1["Annotated_genes_reads"] == 3
        assert cb1["Novel_genes_reads"] == 0
        assert cb1["anno_gene_reads_bin_1"] == 0
        assert cb1["anno_gene_reads_bin_3"] == 1
        anno_read_bins = [f"anno_gene_reads_bin_{i}" for i in range(1, 11)] + [
            "anno_gene_reads_bin_10plus"
        ]
        assert sum(cb1[b] for b in anno_read_bins) == 1, "CB1 has exactly 1 annotated gene"
        # Novel-gene bin columns are not produced at all.
        assert not [c for c in summary.columns if c.startswith("novel_gene_reads_bin_")]


# ==============================================================================
# Export Scanpy / Seurat Tests
# ==============================================================================

# Add scripts path for export module imports
scripts_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../scripts")
)
if scripts_path not in sys.path:
    sys.path.insert(0, scripts_path)

from export_scanpy_seurat import (
    _prepare_classification,
    _build_gene_matrix,
    _build_isoform_matrix,
    _build_anndata,
)


def test_export_gene_matrix_reads_mode(tmpdir):
    """Gene matrix from reads-mode classification has correct shape and counts."""
    cls = pd.DataFrame({
        'isoform': ['r1', 'r2', 'r3', 'r4'],
        'CB': ['CELL1', 'CELL1', 'CELL2', 'CELL2'],
        'associated_gene': ['geneA', 'geneB', 'geneA', 'geneA'],
    })
    class_file = str(tmpdir.join("cls.txt"))
    cls.to_csv(class_file, sep='\t', index=False)

    prepared = _prepare_classification(class_file, "reads")
    mat, barcodes, genes = _build_gene_matrix(prepared, "reads")
    adata = _build_anndata(mat, barcodes, genes)

    assert adata.shape == (2, 2)  # 2 cells, 2 genes
    # CELL2 has 2 reads of geneA
    cell2_idx = list(adata.obs_names).index("CELL2")
    geneA_idx = list(adata.var_names).index("geneA")
    assert adata.X[cell2_idx, geneA_idx] == 2
    assert "isoform_counts" not in adata.obsm


def test_export_isoform_matrix_isoforms_mode(tmpdir):
    """Isoform matrix from isoforms-mode correctly sums FL weights and acts as layer."""
    cls = pd.DataFrame({
        'isoform': ['iso1', 'iso2'],
        'CB': ['CELL1,CELL2', 'CELL2'],
        'FL': ['3,5', '2'],
        'associated_gene': ['geneA', 'geneA'],
    })
    class_file = str(tmpdir.join("cls.txt"))
    cls.to_csv(class_file, sep='\t', index=False)

    prepared = _prepare_classification(class_file, "isoforms")
    gene_mat, barcodes, genes = _build_gene_matrix(prepared, "isoforms")
    iso_mat, isoforms, iso_genes = _build_isoform_matrix(prepared, barcodes)

    adata = _build_anndata(
        gene_mat, barcodes, genes,
        iso_mat=iso_mat, isoforms=isoforms, iso_genes=iso_genes
    )

    assert adata.shape[0] == 2  # 2 cells
    assert "isoform_counts" in adata.obsm
    assert adata.obsm["isoform_counts"].shape == (2, 2)
    assert "isoform_features" in adata.uns


def test_export_h5ad_with_cell_summary(tmpdir):
    """.h5ad .obs contains cell summary metrics when provided."""
    cls = pd.DataFrame({
        'isoform': ['r1', 'r2'],
        'CB': ['CELL1', 'CELL2'],
        'associated_gene': ['geneA', 'geneA'],
    })
    class_file = str(tmpdir.join("cls.txt"))
    cls.to_csv(class_file, sep='\t', index=False)

    summary = pd.DataFrame({
        'CB': ['CELL1', 'CELL2'],
        'Reads_in_cell': [10, 20],
        'Genes_in_cell': [5, 8],
    })
    summary_file = str(tmpdir.join("summary.txt.gz"))
    summary.to_csv(summary_file, sep='\t', index=False, compression='gzip')

    prepared = _prepare_classification(class_file, "reads")
    mat, barcodes, genes = _build_gene_matrix(prepared, "reads")

    adata = _build_anndata(mat, barcodes, genes,
                           cell_summary_file=summary_file)

    assert "Reads_in_cell" in adata.obs.columns
    assert "Genes_in_cell" in adata.obs.columns
    assert adata.obs.loc["CELL1", "Reads_in_cell"] == 10


def test_export_h5ad_with_clustering(tmpdir):
    """.h5ad includes UMAP in obsm and cluster labels in obs."""
    cls = pd.DataFrame({
        'isoform': ['r1', 'r2'],
        'CB': ['CELL1', 'CELL2'],
        'associated_gene': ['geneA', 'geneA'],
    })
    class_file = str(tmpdir.join("cls.txt"))
    cls.to_csv(class_file, sep='\t', index=False)

    clust = pd.DataFrame({
        'Barcode': ['CELL1', 'CELL2'],
        'Cluster': [0, 1],
        'UMAP_1': [1.0, 2.0],
        'UMAP_2': [3.0, 4.0],
    })
    clust_file = str(tmpdir.join("umap.csv"))
    clust.to_csv(clust_file, index=False)

    prepared = _prepare_classification(class_file, "reads")
    mat, barcodes, genes = _build_gene_matrix(prepared, "reads")

    adata = _build_anndata(mat, barcodes, genes,
                           clustering_file=clust_file)

    assert "cluster" in adata.obs.columns
    assert "X_umap" in adata.obsm
    assert adata.obsm["X_umap"].shape == (2, 2)


def test_export_h5ad_minimal(tmpdir):
    """.h5ad works with classification file only (no summary/clustering)."""
    cls = pd.DataFrame({
        'isoform': ['r1', 'r2', 'r3'],
        'CB': ['CELL1', 'CELL1', 'CELL2'],
        'associated_gene': ['geneA', 'geneB', 'geneA'],
    })
    class_file = str(tmpdir.join("cls.txt"))
    cls.to_csv(class_file, sep='\t', index=False)

    prepared = _prepare_classification(class_file, "reads")
    mat, barcodes, genes = _build_gene_matrix(prepared, "reads")

    adata = _build_anndata(mat, barcodes, genes)

    assert adata.shape == (2, 2)
    assert "X_umap" not in adata.obsm


# ==============================================================================
# Integrated h5ad Export Tests (src/sc_export.py)
# ==============================================================================

from sc_export import (
    _prepare_classification as sc_prepare_classification,
    _build_gene_matrix as sc_build_gene_matrix,
    _build_isoform_matrix as sc_build_isoform_matrix,
    _build_anndata as sc_build_anndata,
    export_h5ad,
)


def test_sc_export_h5ad_reads_mode(tmpdir):
    """h5ad from reads mode has correct gene x cell shape."""
    cls = pd.DataFrame({
        'isoform': ['r1', 'r2', 'r3', 'r4'],
        'CB': ['CELL1', 'CELL1', 'CELL2', 'CELL2'],
        'associated_gene': ['geneA', 'geneB', 'geneA', 'geneA'],
    })
    class_file = str(tmpdir.join("cls.txt"))
    cls.to_csv(class_file, sep='\t', index=False)

    prepared = sc_prepare_classification(class_file, "reads")
    mat, barcodes, genes = sc_build_gene_matrix(prepared, "reads")

    adata = sc_build_anndata(mat, barcodes, genes)

    assert adata.shape == (2, 2)  # 2 cells, 2 genes
    cell2_idx = list(adata.obs_names).index("CELL2")
    geneA_idx = list(adata.var_names).index("geneA")
    assert adata.X[cell2_idx, geneA_idx] == 2
    assert "isoform_counts" not in adata.obsm


def test_sc_export_h5ad_isoforms_mode(tmpdir):
    """h5ad from isoforms mode has gene matrix in .X and isoform data in .obsm."""
    cls = pd.DataFrame({
        'isoform': ['iso1', 'iso2'],
        'CB': ['CELL1,CELL2', 'CELL2'],
        'FL': ['3,5', '2'],
        'associated_gene': ['geneA', 'geneA'],
    })
    class_file = str(tmpdir.join("cls.txt"))
    cls.to_csv(class_file, sep='\t', index=False)

    prepared = sc_prepare_classification(class_file, "isoforms")
    gene_mat, barcodes, genes = sc_build_gene_matrix(prepared, "isoforms")
    iso_mat, isoforms, iso_genes = sc_build_isoform_matrix(prepared, barcodes)

    adata = sc_build_anndata(
        gene_mat, barcodes, genes,
        iso_mat=iso_mat, isoforms=isoforms, iso_genes=iso_genes,
    )

    assert adata.shape[0] == 2  # 2 cells
    cell2_idx = list(adata.obs_names).index("CELL2")
    geneA_idx = list(adata.var_names).index("geneA")
    assert adata.X[cell2_idx, geneA_idx] == 7

    assert "isoform_counts" in adata.obsm
    assert adata.obsm["isoform_counts"].shape == (2, 2)
    assert "isoform_features" in adata.uns


def test_sc_export_h5ad_with_cell_summary(tmpdir):
    """.obs contains cell summary QC metrics when summary file exists."""
    cls = pd.DataFrame({
        'isoform': ['r1', 'r2'],
        'CB': ['CELL1', 'CELL2'],
        'associated_gene': ['geneA', 'geneA'],
    })
    class_file = str(tmpdir.join("cls.txt"))
    cls.to_csv(class_file, sep='\t', index=False)

    summary = pd.DataFrame({
        'CB': ['CELL1', 'CELL2'],
        'Reads_in_cell': [10, 20],
        'Genes_in_cell': [5, 8],
    })
    summary_file = str(tmpdir.join("summary.txt.gz"))
    summary.to_csv(summary_file, sep='\t', index=False, compression='gzip')

    prepared = sc_prepare_classification(class_file, "reads")
    mat, barcodes, genes = sc_build_gene_matrix(prepared, "reads")
    adata = sc_build_anndata(mat, barcodes, genes,
                             cell_summary_file=summary_file)

    assert "Reads_in_cell" in adata.obs.columns
    assert "Genes_in_cell" in adata.obs.columns
    assert adata.obs.loc["CELL1", "Reads_in_cell"] == 10


def test_sc_export_h5ad_with_clustering(tmpdir):
    """.obsm has UMAP and .obs has cluster labels when clustering exists."""
    cls = pd.DataFrame({
        'isoform': ['r1', 'r2'],
        'CB': ['CELL1', 'CELL2'],
        'associated_gene': ['geneA', 'geneA'],
    })
    class_file = str(tmpdir.join("cls.txt"))
    cls.to_csv(class_file, sep='\t', index=False)

    clust = pd.DataFrame({
        'Barcode': ['CELL1', 'CELL2'],
        'Cluster': [0, 1],
        'UMAP_1': [1.0, 2.0],
        'UMAP_2': [3.0, 4.0],
    })
    clust_file = str(tmpdir.join("umap.csv"))
    clust.to_csv(clust_file, index=False)

    prepared = sc_prepare_classification(class_file, "reads")
    mat, barcodes, genes = sc_build_gene_matrix(prepared, "reads")
    adata = sc_build_anndata(mat, barcodes, genes,
                             clustering_file=clust_file)

    assert "cluster" in adata.obs.columns
    assert "X_umap" in adata.obsm
    import numpy as np
    assert adata.obsm["X_umap"].shape == (2, 2)


@patch('qc_pipeline.build_parser')
@patch('qc_pipeline.run_sqanti3_qc')
@patch('qc_pipeline.annotate_with_ujc_hash')
@patch('qc_pipeline.annotate_with_cell_metadata')
@patch('qc_pipeline.calculate_metrics_per_cell')
@patch('qc_pipeline.generate_report')
@patch('qc_pipeline.export_h5ad')
def test_pipeline_main_with_export_h5ad(
        mock_export,
        mock_generate_report,
        mock_calc,
        mock_add,
        mock_hash,
        mock_get,
        mock_build_parser,
        tmpdir):
    """Pipeline main calls export_h5ad when --export_h5ad is set."""
    class DummyArgs:
        def __init__(self):
            self.inDESIGN = str(tmpdir.join('design.csv'))
            self.input_dir = str(tmpdir)
            self.out_dir = str(tmpdir)
            self.refGTF = __file__
            self.refFasta = __file__
            self.mode = 'reads'
            self.report = 'pdf'
            self.SKIPHASH = True
            self.write_per_cell_outputs = False
            self.log_level = 'INFO'
            self.run_clustering = False
            self.export_h5ad = True  # <-- enabled

    tmpdir.join('design.csv').write(
        'sampleID,file_acc,cell_association\ns1,f1,assoc.txt\n'
    )

    class DummyParser:
        def parse_args(self_inner):
            return DummyArgs()

    mock_build_parser.return_value = DummyParser()
    pipeline_main()
    assert mock_export.called
    assert mock_generate_report.called
