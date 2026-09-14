import os
import sys
from unittest.mock import MagicMock

# --- MOCK FOR WINDOWS TESTING ---
if sys.platform == 'win32':
    sys.platform = 'linux'
    import collections
    os.uname = lambda: collections.namedtuple(
        'Uname', ['sysname', 'nodename', 'release', 'version', 'machine']
    )('Linux', 'localhost', '5.10.0', '1', 'x86_64')
    sys.modules['pysam'] = MagicMock()
# --------------------------------

import numpy as np
import pandas as pd
import pytest

sqanti_sc_src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../src"))
if sqanti_sc_src_path not in sys.path:
    sys.path.insert(0, sqanti_sc_src_path)

import filter_io
import cell_filter
import filter_args
from cell_filter import decide_cells, RESULT_CELL, RESULT_ARTIFACT, PASS, FAIL, NOT_EVALUATED
from cell_metrics import calculate_metrics_per_cell


def _summary_frame(rows):
    """Cell summary shaped like cell_metrics.py's output: CB first, numerics after."""
    return pd.DataFrame(rows)


def _verdict(summary, rules, **kw):
    kw.setdefault('mode', 'isoforms')
    kw.setdefault('sampleID', 's1')
    kw.setdefault('log', lambda *a, **k: None)
    return decide_cells(summary, rules, **kw).set_index('CB')


class _Args:
    def __init__(self, mode, out_dir):
        self.mode = mode
        self.out_dir = out_dir
        self.min_cov = 1
        self.ratio_TSS_threshold = 2.0
        self.ref_cov_min_pct = 45.0
        self.include_ORF = False


def _write_tsv(path, df):
    df.to_csv(path, sep='\t', index=False)
    return str(path)


def _read_tsv(path):
    return pd.read_csv(path, sep='\t', dtype=str, na_filter=False)


@pytest.fixture
def reads_classification(tmp_path):
    df = pd.DataFrame({
        'isoform': ['r1', 'r2', 'r3', 'r4', 'r5'],
        'CB': ['bc1', 'bc2', 'bc3', 'NA', 'unassigned'],
        'structural_category': ['full-splice_match'] * 5,
        'RTS_stage': ['FALSE', 'TRUE', 'FALSE', 'FALSE', 'FALSE'],
        'min_cov': ['1', 'NA', '3', 'NA', '2'],
    })
    return _write_tsv(tmp_path / "reads_classification.txt", df)


@pytest.fixture
def isoforms_classification(tmp_path):
    df = pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1'],
        'CB': ['bc1,bc2,bc3', 'bc2', 'bc1,bc3'],
        'FL': ['5,3,2', '7', '4,6'],
        'structural_category': ['full-splice_match', 'novel_in_catalog', 'full-splice_match'],
        'min_cov': ['NA', '4', 'NA'],
    })
    return _write_tsv(tmp_path / "iso_classification.txt", df)


class TestRoundTripFidelity:
    """Nothing the filter passes through may be altered. Reading with pandas'
    default NA handling would turn the literal 'NA' that classification_enrichment
    writes into NaN and emit it as '', corrupting untouched rows."""

    def test_isoforms_keep_all_is_identity(self, tmp_path, isoforms_classification):
        dst = tmp_path / "out.txt"
        filter_io.subset_classification(
            isoforms_classification, str(dst), 'isoforms', {'bc1', 'bc2', 'bc3'})
        pd.testing.assert_frame_equal(
            _read_tsv(isoforms_classification), _read_tsv(dst))

    def test_reads_keep_all_preserves_literal_NA_in_other_columns(
            self, tmp_path, reads_classification):
        dst = tmp_path / "out.txt"
        filter_io.subset_classification(
            reads_classification, str(dst), 'reads', {'bc1', 'bc2', 'bc3'})
        out = _read_tsv(dst)
        assert out['min_cov'].tolist() == ['1', 'NA', '3']
        assert out['RTS_stage'].tolist() == ['FALSE', 'TRUE', 'FALSE']


class TestCellFilterFLInvariant:
    """The write-side twin of TestFLWeightingIsoformsMode in sqanti_sc_test.py.

    In isoforms mode the comma-separated CB/FL pair is a sparse (transcript x cell)
    matrix, so removing a cell edits those strings; it does not drop the row. A row
    is dropped only when no cell survives.
    """

    def test_removing_one_cell_edits_strings_and_keeps_the_row(
            self, tmp_path, isoforms_classification):
        dst = tmp_path / "out.txt"
        filter_io.subset_classification(
            isoforms_classification, str(dst), 'isoforms', {'bc1', 'bc3'})
        out = _read_tsv(dst).set_index('isoform')
        assert out.loc['PB.1.1', 'CB'] == 'bc1,bc3'
        assert out.loc['PB.1.1', 'FL'] == '5,2'

    def test_model_supported_only_by_a_filtered_cell_is_dropped(
            self, tmp_path, isoforms_classification):
        dst = tmp_path / "out.txt"
        _, rows_out, surviving, _ = filter_io.subset_classification(
            isoforms_classification, str(dst), 'isoforms', {'bc1', 'bc3'})
        assert 'PB.2.1' not in surviving
        assert rows_out == 2

    def test_cb_and_fl_stay_positionally_aligned(self, tmp_path):
        src = _write_tsv(tmp_path / "c.txt", pd.DataFrame({
            'isoform': ['PB.1.1'],
            'CB': ['bcA,bcB,bcC,bcD'],
            'FL': ['10,20,30,40'],
        }))
        dst = tmp_path / "out.txt"
        filter_io.subset_classification(src, str(dst), 'isoforms', {'bcB', 'bcD'})
        out = _read_tsv(dst)
        assert out.loc[0, 'CB'] == 'bcB,bcD'
        assert out.loc[0, 'FL'] == '20,40'

    def test_mismatched_cb_fl_lengths_raise_a_named_error(self, tmp_path):
        src = _write_tsv(tmp_path / "c.txt", pd.DataFrame({
            'isoform': ['PB.1.1'], 'CB': ['bc1,bc2,bc3'], 'FL': ['5,3'],
        }))
        with pytest.raises(ValueError, match="different element counts"):
            filter_io.subset_classification(
                src, str(tmp_path / "out.txt"), 'isoforms', {'bc1'})

    def test_classification_without_fl_still_filters_cb(self, tmp_path):
        src = _write_tsv(tmp_path / "c.txt", pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'], 'CB': ['bc1,bc2', 'bc2'],
        }))
        dst = tmp_path / "out.txt"
        filter_io.subset_classification(src, str(dst), 'isoforms', {'bc1'})
        out = _read_tsv(dst)
        assert out['CB'].tolist() == ['bc1']


class TestReadsModeSubset:
    def test_filtered_and_sentinel_barcodes_are_dropped(
            self, tmp_path, reads_classification):
        dst = tmp_path / "out.txt"
        rows_in, rows_out, surviving, no_barcode = filter_io.subset_classification(
            reads_classification, str(dst), 'reads', {'bc1', 'bc3'})
        assert rows_in == 5
        assert rows_out == 2
        assert no_barcode == 2
        assert surviving == {'r1', 'r3'}


class TestChunking:
    def test_chunked_output_matches_unchunked(self, tmp_path, isoforms_classification):
        one = tmp_path / "one.txt"
        many = tmp_path / "many.txt"
        filter_io.subset_classification(
            isoforms_classification, str(one), 'isoforms', {'bc1', 'bc3'})
        filter_io.subset_classification(
            isoforms_classification, str(many), 'isoforms', {'bc1', 'bc3'}, chunksize=1)
        pd.testing.assert_frame_equal(_read_tsv(one), _read_tsv(many))

    def test_chunked_reads_output_matches_unchunked(self, tmp_path, reads_classification):
        one = tmp_path / "one.txt"
        many = tmp_path / "many.txt"
        filter_io.subset_classification(
            reads_classification, str(one), 'reads', {'bc1', 'bc3'})
        filter_io.subset_classification(
            reads_classification, str(many), 'reads', {'bc1', 'bc3'}, chunksize=2)
        pd.testing.assert_frame_equal(_read_tsv(one), _read_tsv(many))


class TestJunctionsSubset:
    def test_isoforms_mode_drops_the_stale_cb_column(self, tmp_path):
        src = _write_tsv(tmp_path / "j.txt", pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'CB': ['bc1,bc2,bc3', 'bc2'],
            'junction_category': ['known', 'novel'],
        }))
        dst = tmp_path / "out.txt"
        _, rows_out, dropped_cb = filter_io.subset_junctions(
            src, str(dst), 'isoforms', {'PB.1.1'})
        out = _read_tsv(dst)
        assert dropped_cb is True
        assert 'CB' not in out.columns
        assert rows_out == 1

    def test_reads_mode_keeps_the_cb_column(self, tmp_path):
        src = _write_tsv(tmp_path / "j.txt", pd.DataFrame({
            'isoform': ['r1', 'r2'], 'CB': ['bc1', 'bc2'],
            'junction_category': ['known', 'novel'],
        }))
        dst = tmp_path / "out.txt"
        _, _, dropped_cb = filter_io.subset_junctions(src, str(dst), 'reads', {'r1'})
        out = _read_tsv(dst)
        assert dropped_cb is False
        assert out['CB'].tolist() == ['bc1']

    def test_missing_junctions_file_is_not_an_error(self, tmp_path):
        assert filter_io.subset_junctions(
            str(tmp_path / "absent.txt"), str(tmp_path / "o.txt"), 'reads', set()
        ) == (0, 0, False)


class TestBarcodeLists:
    def test_single_column_applies_to_every_sample(self, tmp_path):
        path = tmp_path / "bc.txt"
        path.write_text("# a comment\n\nbc1\nbc2\n")
        scoped = filter_io.read_barcode_list(str(path))
        assert filter_io.barcodes_for_sample(scoped, 'anySample') == {'bc1', 'bc2'}

    def test_two_columns_scope_to_one_sample(self, tmp_path):
        path = tmp_path / "bc.txt"
        path.write_text("sampleA\tbc1\nsampleB\tbc2\n")
        scoped = filter_io.read_barcode_list(str(path))
        assert filter_io.barcodes_for_sample(scoped, 'sampleA') == {'bc1'}
        assert filter_io.barcodes_for_sample(scoped, 'sampleB') == {'bc2'}

    def test_global_and_scoped_entries_combine(self, tmp_path):
        path = tmp_path / "bc.txt"
        path.write_text("bc_global\nsampleA\tbc_a\n")
        scoped = filter_io.read_barcode_list(str(path))
        assert filter_io.barcodes_for_sample(scoped, 'sampleA') == {'bc_global', 'bc_a'}
        assert filter_io.barcodes_for_sample(scoped, 'sampleB') == {'bc_global'}


class TestSampleTable:
    def _make_outputs(self, tmp_path, file_acc='rep1', sampleID='s1'):
        d = tmp_path / file_acc
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{sampleID}_classification.txt").write_text("isoform\tCB\n")
        pd.DataFrame({'CB': ['bc1']}).to_csv(
            d / f"{sampleID}_SQANTI_cell_summary.txt.gz", sep='\t',
            index=False, compression='gzip')

    def test_reads_a_valid_design_without_rewriting_it(self, tmp_path):
        self._make_outputs(tmp_path)
        design = tmp_path / "design.csv"
        design.write_text("sampleID,file_acc\ns1,rep1\n")
        before = design.read_text()
        df = filter_io.read_sample_table(str(design), str(tmp_path))
        assert len(df) == 1
        assert design.read_text() == before

    def test_missing_required_column_is_a_clear_error(self, tmp_path):
        design = tmp_path / "design.csv"
        design.write_text("sampleID\ns1\n")
        with pytest.raises(ValueError, match="missing required column"):
            filter_io.read_sample_table(str(design), str(tmp_path))

    def test_missing_qc_output_names_the_path(self, tmp_path):
        design = tmp_path / "design.csv"
        design.write_text("sampleID,file_acc\ns1,rep1\n")
        with pytest.raises(ValueError, match="_classification.txt"):
            filter_io.read_sample_table(str(design), str(tmp_path))


class TestFilteredDesign:
    def test_file_acc_nests_the_run_directory(self, tmp_path):
        df = pd.DataFrame({
            'sampleID': ['s1', 's2'], 'file_acc': ['rep1', 'rep2'],
            'input_file': ['/a.bam', '/b.bam'], 'color_group': ['x', 'y'],
        })
        out_path = tmp_path / "d.csv"
        filter_io.write_filtered_design(str(out_path), df, 'cell_filter', str(tmp_path))
        out = pd.read_csv(out_path)
        assert out['file_acc'].tolist() == [
            os.path.join('rep1', 'cell_filter'), os.path.join('rep2', 'cell_filter')]
        assert out['source_file_acc'].tolist() == ['rep1', 'rep2']
        assert 'color_group' in out.columns
        assert 'input_file' not in out.columns


class TestSubsetEqualsRecompute:
    """Subsetting the cell summary must equal recomputing it from the filtered
    classification. That equality is what licenses skipping a re-run of
    calculate_metrics_per_cell -- which would re-parse a junctions file reaching
    ~70GB -- and it fails the day a cross-cell column is added to cell_metrics.py.

    Note it holds only for cell REMOVAL. Step 3's requantification changes the
    content of retained cells, so it must recompute.
    """

    def _cls_row(self, isoform, cb, fl, category, exons=2):
        return {
            "isoform": isoform, "CB": cb, "FL": fl,
            "structural_category": category, "associated_gene": "geneA",
            "associated_transcript": "txA", "exons": exons, "length": 500,
            "ref_length": 600, "chrom": "chr1", "subcategory": "reference_match",
            "all_canonical": "True", "RTS_stage": "False", "predicted_NMD": "False",
            "within_CAGE_peak": "False", "polyA_motif_found": "False",
            "perc_A_downstream_TTS": "0", "diff_to_gene_TSS": "0",
            "coding": "coding", "min_cov": "0", "ratio_TSS": "0",
        }

    def _summary(self, tmp_path, name, cls_path):
        out_dir = tmp_path / name
        sample_dir = out_dir / "f1"
        sample_dir.mkdir(parents=True, exist_ok=True)
        prefix = sample_dir / "s1"
        pd.read_csv(cls_path, sep='\t', dtype=str, na_filter=False).to_csv(
            f"{prefix}_classification.txt", sep='\t', index=False)
        pd.DataFrame(columns=["isoform"]).to_csv(
            f"{prefix}_junctions.txt", sep='\t', index=False)
        design = pd.DataFrame({'sampleID': ['s1'], 'file_acc': ['f1']})
        calculate_metrics_per_cell(_Args('isoforms', str(out_dir)), design)
        return pd.read_csv(f"{prefix}_SQANTI_cell_summary.txt.gz", sep='\t')

    def test_isoforms_subset_matches_recompute(self, tmp_path):
        src = _write_tsv(tmp_path / "cls.txt", pd.DataFrame([
            self._cls_row("PB.1.1", "bc1,bc2,bc3", "5,3,2", "full-splice_match"),
            self._cls_row("PB.2.1", "bc2", "7", "novel_in_catalog"),
            self._cls_row("PB.3.1", "bc1,bc3", "4,6", "full-splice_match", exons=1),
        ]))
        keep = {'bc1', 'bc3'}

        full = self._summary(tmp_path, "full", src)
        subset = full[full['CB'].isin(keep)].reset_index(drop=True)

        filtered_cls = tmp_path / "cls_filtered.txt"
        filter_io.subset_classification(src, str(filtered_cls), 'isoforms', keep)
        recomputed = self._summary(tmp_path, "recomputed", str(filtered_cls))

        pd.testing.assert_frame_equal(subset, recomputed, check_dtype=False)


class TestRulesSemantics:
    """SQANTI3's rules vocabulary: a list of numbers is a [min, max] range, a bare
    number is a minimum, a string or list of strings is membership."""

    def test_bare_number_is_a_minimum(self):
        s = _summary_frame({'CB': ['a', 'b'], 'Transcripts_in_cell': [5, 50]})
        out = _verdict(s, {'Transcripts_in_cell': 10})
        assert out.loc['a', 'filter_result'] == RESULT_ARTIFACT
        assert out.loc['b', 'filter_result'] == RESULT_CELL

    def test_two_element_list_is_a_range(self):
        s = _summary_frame({'CB': ['a', 'b', 'c'], 'MT_perc': [0.0, 30.0, 80.0],
                            'Transcripts_in_cell': [500, 500, 500]})
        out = _verdict(s, {'MT_perc': [0, 50]})
        assert out['filter_result'].tolist() == [RESULT_CELL, RESULT_CELL, RESULT_ARTIFACT]

    def test_string_rule_is_case_insensitive_membership(self):
        s = _summary_frame({'CB': ['a', 'b'], 'flag': ['TRUE', 'FALSE']})
        out = _verdict(s, {'flag': 'true'})
        assert out['filter_result'].tolist() == [RESULT_CELL, RESULT_ARTIFACT]

    def test_string_list_rule_is_membership(self):
        s = _summary_frame({'CB': ['a', 'b', 'c'], 'flag': ['x', 'y', 'z']})
        out = _verdict(s, {'flag': ['x', 'y']})
        assert out['filter_result'].tolist() == [RESULT_CELL, RESULT_CELL, RESULT_ARTIFACT]

    def test_reason_names_the_rule_and_the_observed_value(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [3]})
        out = _verdict(s, {'Transcripts_in_cell': 10})
        assert 'Transcripts_in_cell >= 10' in out.loc['a', 'filter_reason']
        assert 'got 3' in out.loc['a', 'filter_reason']

    def test_a_cell_failing_two_rules_reports_both(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [3], 'Annotated_genes': [1]})
        out = _verdict(s, {'Transcripts_in_cell': 10, 'Annotated_genes': 10})
        assert out.loc['a', 'filter_reason'].count(';') == 1

    def test_unknown_rule_column_is_a_clear_error(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [5]})
        with pytest.raises(ValueError, match="not present in the cell summary"):
            _verdict(s, {'No_Such_Column': 1})

    def test_unsupported_rule_shape_is_a_clear_error(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [5]})
        with pytest.raises(ValueError, match="unsupported rule"):
            _verdict(s, {'Transcripts_in_cell': {'min': 1}})

    def test_every_cell_failing_leaves_no_passing_barcode(self):
        s = _summary_frame({'CB': ['a', 'b'], 'Transcripts_in_cell': [1, 2]})
        out = _verdict(s, {'Transcripts_in_cell': 100})
        assert (out['filter_result'] == RESULT_ARTIFACT).all()


class TestDenominatorGating:
    """safe_prop returns 0 when its denominator is 0, and 0 is the good end of every
    artifact scale -- so an ungated max rule reads a missing value as a clean one.
    That is a silent false PASS, which is why gating exists and why the status has
    three values instead of two."""

    def test_healthy_depth_but_no_multiexon_reads_is_not_evaluated(self):
        s = _summary_frame({
            'CB': ['deep_but_monoexon'],
            'Transcripts_in_cell': [5000],
            'total_transcripts_no_monoexon': [3],
            'Non_canonical_prop_in_cell': [0.0],
        })
        out = _verdict(s, {'Non_canonical_prop_in_cell': [0, 10]}, min_depth_for_props=100)
        assert out.loc['deep_but_monoexon', 'Non_canonical_prop_in_cell_status'] == NOT_EVALUATED

    def test_not_evaluated_does_not_discard_the_cell(self):
        s = _summary_frame({
            'CB': ['a'], 'Transcripts_in_cell': [5000],
            'total_transcripts_no_monoexon': [3], 'Non_canonical_prop_in_cell': [99.0],
        })
        out = _verdict(s, {'Non_canonical_prop_in_cell': [0, 10]}, min_depth_for_props=100)
        assert out.loc['a', 'filter_result'] == RESULT_CELL

    def test_junction_props_gate_on_total_junctions_not_depth(self):
        s = _summary_frame({
            'CB': ['a'], 'Transcripts_in_cell': [5000], 'total_junctions': [2],
            'Novel_non_canonical_junctions_prop': [50.0],
        })
        out = _verdict(s, {'Novel_non_canonical_junctions_prop': [0, 10]},
                       min_depth_for_props=100)
        assert out.loc['a', 'Novel_non_canonical_junctions_prop_status'] == NOT_EVALUATED

    def test_shallow_cell_is_not_evaluated_for_depth_denominated_props(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [20],
                            'RTS_prop_in_cell': [5.0]})
        out = _verdict(s, {'RTS_prop_in_cell': [0, 1]}, min_depth_for_props=100)
        assert out.loc['a', 'RTS_prop_in_cell_status'] == NOT_EVALUATED

    def test_a_deep_enough_cell_is_evaluated_normally(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [5000],
                            'RTS_prop_in_cell': [40.0]})
        out = _verdict(s, {'RTS_prop_in_cell': [0, 5]}, min_depth_for_props=100)
        assert out.loc['a', 'RTS_prop_in_cell_status'] == FAIL
        assert out.loc['a', 'filter_result'] == RESULT_ARTIFACT

    def test_raw_counts_are_never_gated(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [3]})
        out = _verdict(s, {'Transcripts_in_cell': 10}, min_depth_for_props=100)
        assert out.loc['a', 'Transcripts_in_cell_status'] == FAIL

    def test_the_depth_column_follows_the_mode(self):
        assert cell_filter.depth_column('reads') == 'Reads_in_cell'
        assert cell_filter.depth_column('isoforms') == 'Transcripts_in_cell'
        assert cell_filter.denominator_column(
            'Non_canonical_prop_in_cell', 'reads') == 'total_reads_no_monoexon'


class TestManualLists:
    def _summary(self):
        return _summary_frame({'CB': ['a', 'b', 'c'], 'Transcripts_in_cell': [5000, 5000, 5000]})

    def test_drop_list_discards_a_passing_cell(self):
        out = _verdict(self._summary(), {'Transcripts_in_cell': 10}, drop={'b'})
        assert out.loc['b', 'filter_result'] == RESULT_ARTIFACT
        assert out.loc['b', 'filter_source'] == 'manual_drop'

    def test_keep_list_bypasses_a_failing_rule(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [1]})
        out = _verdict(s, {'Transcripts_in_cell': 10}, keep={'a'})
        assert out.loc['a', 'filter_result'] == RESULT_CELL
        assert out.loc['a', 'filter_source'] == 'manual_keep'

    def test_a_barcode_in_both_lists_is_a_hard_error_naming_it(self):
        with pytest.raises(ValueError, match="both --keep_barcodes and --drop_barcodes"):
            _verdict(self._summary(), {'Transcripts_in_cell': 10},
                     keep={'b'}, drop={'b'})

    def test_universe_excludes_unlisted_barcodes(self):
        out = _verdict(self._summary(), {'Transcripts_in_cell': 10}, universe={'a', 'b'})
        assert out.loc['c', 'filter_source'] == 'not_in_universe'
        assert out.loc['a', 'filter_result'] == RESULT_CELL

    def test_universe_members_still_face_the_rules(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [1]})
        out = _verdict(s, {'Transcripts_in_cell': 10}, universe={'a'})
        assert out.loc['a', 'filter_source'] == 'rules'

    def test_drop_beats_keep_only_via_the_error_not_silently(self):
        out = _verdict(self._summary(), {'Transcripts_in_cell': 10},
                       keep={'a'}, drop={'b'})
        assert out.loc['a', 'filter_source'] == 'manual_keep'
        assert out.loc['b', 'filter_source'] == 'manual_drop'


class TestSentinelBarcodes:
    """cell_metrics.py filters only notna() and != '', so an 'unassigned' row -- which
    SQANTI-sc_report.R guards against in ten places -- can sit in the cell summary
    aggregating every unbarcoded read. It would be rank 1 on a barcode-rank curve."""

    def test_sentinel_rows_are_marked_artifact_and_reported(self):
        messages = []
        s = _summary_frame({'CB': ['unassigned', 'a'], 'Transcripts_in_cell': [999999, 5000]})
        out = decide_cells(s, {'Transcripts_in_cell': 10}, mode='isoforms', sampleID='s1',
                           log=messages.append).set_index('CB')
        assert out.loc['unassigned', 'filter_result'] == RESULT_ARTIFACT
        assert out.loc['unassigned', 'filter_source'] == 'sentinel'
        assert any('sentinel barcode' in m for m in messages)

    def test_every_barcode_survives_into_the_verdict_table(self):
        s = _summary_frame({'CB': ['unassigned', 'NA', 'a'],
                            'Transcripts_in_cell': [10, 10, 5000]})
        out = _verdict(s, {'Transcripts_in_cell': 10})
        assert set(out.index) == {'unassigned', 'NA', 'a'}


class TestConditionalColumnWarning:
    def test_a_constant_conditional_column_warns(self):
        messages = []
        s = _summary_frame({'CB': ['a', 'b'], 'CAGE_peak_support_prop': [0.0, 0.0],
                            'Transcripts_in_cell': [5000, 5000]})
        decide_cells(s, {'CAGE_peak_support_prop': 1}, mode='isoforms', sampleID='s1',
                     log=messages.append)
        assert any('does not vary' in m for m in messages)

    def test_a_varying_conditional_column_does_not_warn(self):
        messages = []
        s = _summary_frame({'CB': ['a', 'b'], 'CAGE_peak_support_prop': [10.0, 80.0],
                            'Transcripts_in_cell': [5000, 5000]})
        decide_cells(s, {'CAGE_peak_support_prop': 1}, mode='isoforms', sampleID='s1',
                     log=messages.append)
        assert not any('does not vary' in m for m in messages)


class TestRobustSpread:
    """Ported from qc_cell_spread() in SQANTI-sc_multisample_report.R. R's mad()
    applies the 1.4826 constant by default and scipy's does not, which is the
    likeliest porting bug."""

    def test_mad_branch_matches_r_scaling(self):
        from cell_filter_methods import robust_spread
        assert robust_spread([1, 2, 3, 4, 5]) == pytest.approx(1.4826)

    def test_zero_inflation_falls_back_to_the_quantile_range(self):
        from cell_filter_methods import robust_spread
        values = [0.0] * 80 + [10.0] * 20
        assert robust_spread(values) == pytest.approx((10.0 - 0.0) / 2.563)

    def test_an_all_zero_feature_returns_nan_not_zero(self):
        from cell_filter_methods import robust_spread
        assert np.isnan(robust_spread([0.0] * 100))


class _FilterArgs:
    def __init__(self, out_dir, mode='isoforms', **kw):
        self.out_dir = str(out_dir)
        self.mode = mode
        self.run_name = 'cell_filter'
        self.rules = filter_args.DEFAULT_CELL_RULES
        self.auto_method = 'none'
        self.min_depth_for_props = 100
        self.keep_barcodes = None
        self.drop_barcodes = None
        self.barcode_universe = None
        self.apply = False
        for k, v in kw.items():
            setattr(self, k, v)


@pytest.fixture
def qc_run(tmp_path):
    """A minimal finished QC run: one isoforms-mode sample on disk."""
    sample_dir = tmp_path / "rep1"
    sample_dir.mkdir()
    prefix = sample_dir / "s1"
    pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1'],
        'CB': ['bc_good,bc_shallow', 'bc_shallow'],
        'FL': ['500,2', '3'],
        'min_cov': ['NA', '4'],
    }).to_csv(f"{prefix}_classification.txt", sep='\t', index=False)
    pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1'], 'CB': ['bc_good,bc_shallow', 'bc_shallow'],
        'junction_category': ['known', 'novel'],
    }).to_csv(f"{prefix}_junctions.txt", sep='\t', index=False)
    pd.DataFrame({
        'CB': ['bc_good', 'bc_shallow'],
        'Transcripts_in_cell': [500, 5],
        'Annotated_genes': [120, 2],
        'MT_perc': [3.0, 0.0],
        'RTS_prop_in_cell': [0.5, 0.0],
        'Intrapriming_prop_in_cell': [1.0, 0.0],
        'Non_canonical_prop_in_cell': [2.0, 0.0],
        'total_transcripts_no_monoexon': [400, 4],
    }).to_csv(f"{prefix}_SQANTI_cell_summary.txt.gz", sep='\t', index=False,
              compression='gzip')
    design = tmp_path / "design.csv"
    design.write_text("sampleID,file_acc\ns1,rep1\n")
    return tmp_path, str(design)


class TestRunCellFilter:
    def test_verdict_artifacts_are_written_without_apply(self, qc_run):
        out_dir, design = qc_run
        df = filter_io.read_sample_table(design, str(out_dir))
        cell_filter.run_cell_filter(_FilterArgs(out_dir), df, log=lambda *a: None)
        run_dir = out_dir / "rep1" / "cell_filter"
        for name in ('s1_CellFilter_cell_summary.txt.gz', 's1_pass_cells.txt',
                     's1_cell_filtering_reasons.txt', 's1_cell_filter_params.txt'):
            assert (run_dir / name).exists()

    def test_defaults_discard_the_shallow_cell(self, qc_run):
        out_dir, design = qc_run
        df = filter_io.read_sample_table(design, str(out_dir))
        cell_filter.run_cell_filter(_FilterArgs(out_dir), df, log=lambda *a: None)
        passing = (out_dir / "rep1" / "cell_filter" / "s1_pass_cells.txt").read_text().split()
        assert passing == ['bc_good']

    def test_verdict_table_keeps_every_barcode(self, qc_run):
        out_dir, design = qc_run
        df = filter_io.read_sample_table(design, str(out_dir))
        cell_filter.run_cell_filter(_FilterArgs(out_dir), df, log=lambda *a: None)
        verdict = pd.read_csv(
            out_dir / "rep1" / "cell_filter" / "s1_CellFilter_cell_summary.txt.gz", sep='\t')
        assert sorted(verdict['CB']) == ['bc_good', 'bc_shallow']
        assert set(verdict['filter_result']) == {'Cell', 'Artifact'}

    def test_no_filtered_data_is_written_without_apply(self, qc_run):
        out_dir, design = qc_run
        df = filter_io.read_sample_table(design, str(out_dir))
        cell_filter.run_cell_filter(_FilterArgs(out_dir), df, log=lambda *a: None)
        assert not (out_dir / "rep1" / "cell_filter" / "s1_classification.txt").exists()

    def test_apply_materialises_the_filtered_dataset(self, qc_run):
        out_dir, design = qc_run
        df = filter_io.read_sample_table(design, str(out_dir))
        filtered_design = cell_filter.run_cell_filter(
            _FilterArgs(out_dir, apply=True), df, log=lambda *a: None)
        run_dir = out_dir / "rep1" / "cell_filter"
        cls = _read_tsv(run_dir / "s1_classification.txt")
        assert cls['CB'].tolist() == ['bc_good']
        assert cls['FL'].tolist() == ['500']
        assert filtered_design is not None
        assert pd.read_csv(filtered_design)['file_acc'].tolist() == [
            os.path.join('rep1', 'cell_filter')]

    def test_params_file_records_the_lost_transcript_models(self, qc_run):
        out_dir, design = qc_run
        df = filter_io.read_sample_table(design, str(out_dir))
        cell_filter.run_cell_filter(_FilterArgs(out_dir, apply=True), df, log=lambda *a: None)
        params = dict(
            line.split('\t', 1) for line in
            (out_dir / "rep1" / "cell_filter" / "s1_cell_filter_params.txt")
            .read_text().strip().split('\n'))
        assert params['TranscriptModelsLostAllSupport'] == '1'
        assert params['BarcodesPassing'] == '1'

    def test_the_run_directory_makes_trials_coexist(self, qc_run):
        out_dir, design = qc_run
        df = filter_io.read_sample_table(design, str(out_dir))
        cell_filter.run_cell_filter(
            _FilterArgs(out_dir, run_name='trial_a'), df, log=lambda *a: None)
        cell_filter.run_cell_filter(
            _FilterArgs(out_dir, run_name='trial_b'), df, log=lambda *a: None)
        assert (out_dir / "rep1" / "trial_a" / "s1_pass_cells.txt").exists()
        assert (out_dir / "rep1" / "trial_b" / "s1_pass_cells.txt").exists()


class TestNonDestructiveness:
    """Every existing output, and the user's design CSV, must be byte-identical after
    a full --apply run. qc_io.fill_design_table rewrites its input in place, which is
    exactly the behaviour the filter must not inherit."""

    def test_inputs_are_untouched_by_apply(self, qc_run):
        import hashlib
        out_dir, design = qc_run
        targets = [
            out_dir / "rep1" / "s1_classification.txt",
            out_dir / "rep1" / "s1_junctions.txt",
            out_dir / "rep1" / "s1_SQANTI_cell_summary.txt.gz",
            out_dir / "design.csv",
        ]
        before = {p: hashlib.md5(p.read_bytes()).hexdigest() for p in targets}
        df = filter_io.read_sample_table(design, str(out_dir))
        cell_filter.run_cell_filter(_FilterArgs(out_dir, apply=True), df, log=lambda *a: None)
        after = {p: hashlib.md5(p.read_bytes()).hexdigest() for p in targets}
        assert before == after


class TestImportHygiene:
    """SQANTI3's src/commands.py calls _get_gtftogenepred_binary() at module scope and
    raises FileNotFoundError when the platform binary is missing. The standalone filter
    must never pull it in, so thresholds can be retuned anywhere."""

    def test_filter_pipeline_does_not_import_sqanti3_commands(self):
        import subprocess
        code = (
            "import os, sys; "
            f"sys.path.insert(0, {repr(sqanti_sc_src_path)}); "
            "import filter_pipeline; "
            "assert 'src.commands' not in sys.modules, 'SQANTI3 src.commands was imported'; "
            "print('ok')"
        )
        result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert 'ok' in result.stdout


class TestRulesFile:
    def test_the_shipped_default_parses_in_both_modes(self):
        reads = cell_filter.load_rules(filter_args.DEFAULT_CELL_RULES, 'reads')
        iso = cell_filter.load_rules(filter_args.DEFAULT_CELL_RULES, 'isoforms')
        assert 'Reads_in_cell' in reads and 'Transcripts_in_cell' not in reads
        assert 'Transcripts_in_cell' in iso and 'Reads_in_cell' not in iso
        assert cell_filter.DEPTH_TOKEN not in reads

    def test_the_default_names_no_conditional_column(self):
        rules = cell_filter.load_rules(filter_args.DEFAULT_CELL_RULES, 'isoforms')
        assert not set(rules) & set(cell_filter.CONDITIONAL_COLUMNS)

    def test_a_rules_file_without_the_all_key_is_rejected(self, tmp_path):
        path = tmp_path / "r.json"
        path.write_text('{"Transcripts_in_cell": 10}')
        with pytest.raises(ValueError, match='top-level "all" key'):
            cell_filter.load_rules(str(path), 'isoforms')


class TestCellSummarySubset:
    def test_subset_preserves_values_verbatim(self, tmp_path):
        src = tmp_path / "s.txt.gz"
        pd.DataFrame({
            'CB': ['bc1', 'bc2', 'bc3'],
            'Reads_in_cell': ['100', '2', '50'],
            'MT_perc': ['1.5', '0.0', '12.25'],
        }).to_csv(src, sep='\t', index=False, compression='gzip')
        dst = tmp_path / "out.txt.gz"
        rows_in, rows_out = filter_io.subset_cell_summary(str(src), str(dst), {'bc1', 'bc3'})
        out = pd.read_csv(dst, sep='\t', dtype=str, na_filter=False)
        assert (rows_in, rows_out) == (3, 2)
        assert out['CB'].tolist() == ['bc1', 'bc3']
        assert out['MT_perc'].tolist() == ['1.5', '12.25']
