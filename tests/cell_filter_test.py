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
    """The per-row CB column is a copy of the parent transcript's barcode list, so it
    goes stale as soon as cells are removed. The per-sample report groups the junctions
    by CB, so it has to be rewritten -- dropping it made that report fail outright."""

    def test_isoforms_mode_rewrites_the_stale_cb_list(self, tmp_path):
        src = _write_tsv(tmp_path / "j.txt", pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'CB': ['bc1,bc2,bc3', 'bc2'],
            'junction_category': ['known', 'novel'],
        }))
        dst = tmp_path / "out.txt"
        _, rows_out, rewrote_cb = filter_io.subset_junctions(
            src, str(dst), 'isoforms', {'PB.1.1'}, {'bc1', 'bc3'})
        out = _read_tsv(dst)
        assert rewrote_cb is True
        assert out['CB'].tolist() == ['bc1,bc3']
        assert rows_out == 1

    def test_isoforms_mode_keeps_the_cb_column_present(self, tmp_path):
        """The report reads it, so its absence is a failure, not an optimisation."""
        src = _write_tsv(tmp_path / "j.txt", pd.DataFrame({
            'isoform': ['PB.1.1'], 'CB': ['bc1,bc2'], 'junction_category': ['known'],
        }))
        dst = tmp_path / "out.txt"
        filter_io.subset_junctions(src, str(dst), 'isoforms', {'PB.1.1'}, {'bc1'})
        assert 'CB' in _read_tsv(dst).columns

    def test_reads_mode_keeps_the_cb_column_verbatim(self, tmp_path):
        src = _write_tsv(tmp_path / "j.txt", pd.DataFrame({
            'isoform': ['r1', 'r2'], 'CB': ['bc1', 'bc2'],
            'junction_category': ['known', 'novel'],
        }))
        dst = tmp_path / "out.txt"
        _, _, rewrote_cb = filter_io.subset_junctions(
            src, str(dst), 'reads', {'r1'}, {'bc1'})
        out = _read_tsv(dst)
        assert rewrote_cb is False
        assert out['CB'].tolist() == ['bc1']

    def test_missing_junctions_file_is_not_an_error(self, tmp_path):
        assert filter_io.subset_junctions(
            str(tmp_path / "absent.txt"), str(tmp_path / "o.txt"), 'reads', set(), set()
        ) == (0, 0, False)


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
    """SQANTI3's numeric rules vocabulary: a list of numbers is a [min, max] range and
    a bare number is a minimum. SQANTI3's string forms are rejected rather than
    supported -- cell_metrics.py coerces every column but CB to numeric, so a string
    rule matches nothing and would silently discard every cell."""

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

    def test_a_string_rule_is_rejected_rather_than_matching_nothing(self):
        s = _summary_frame({'CB': ['a', 'b'], 'Transcripts_in_cell': [5, 50]})
        with pytest.raises(ValueError, match="must be a number"):
            _verdict(s, {'Transcripts_in_cell': 'canonical'})

    def test_a_list_of_strings_is_rejected(self):
        s = _summary_frame({'CB': ['a', 'b'], 'Transcripts_in_cell': [5, 50]})
        with pytest.raises(ValueError, match="must be a number"):
            _verdict(s, {'Transcripts_in_cell': ['TRUE', 'FALSE']})

    def test_a_boolean_rule_is_rejected(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [5]})
        with pytest.raises(ValueError, match="unsupported rule"):
            _verdict(s, {'Transcripts_in_cell': True})

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


class TestUnmeasuredValues:
    """Post-#53/#54 an undefined proportion or a never-measured attribute is NA, not 0.
    A comparison against NA is False, so judging a cell on one discards it with a
    confidently wrong reason -- the cell is skipped for that criterion instead."""

    def test_an_NA_value_is_not_evaluated_rather_than_failed(self):
        s = _summary_frame({
            'CB': ['no_FSM_reads'],
            'Transcripts_in_cell': [5000],
            'FSM_RTS_prop': [np.nan],
        })
        out = _verdict(s, {'FSM_RTS_prop': [0, 5]})
        assert out.loc['no_FSM_reads', 'FSM_RTS_prop_status'] == NOT_EVALUATED

    def test_an_NA_value_does_not_discard_the_cell(self):
        s = _summary_frame({
            'CB': ['a'], 'Transcripts_in_cell': [5000], 'FSM_RTS_prop': [np.nan],
        })
        out = _verdict(s, {'FSM_RTS_prop': [0, 5]})
        assert out.loc['a', 'filter_result'] == RESULT_CELL
        assert out.loc['a', 'filter_reason'] == ''

    def test_a_cell_skipped_on_one_rule_is_still_judged_on_the_others(self):
        s = _summary_frame({
            'CB': ['a'], 'Transcripts_in_cell': [3], 'FSM_RTS_prop': [np.nan],
        })
        out = _verdict(s, {'FSM_RTS_prop': [0, 5], 'Transcripts_in_cell': 500})
        assert out.loc['a', 'FSM_RTS_prop_status'] == NOT_EVALUATED
        assert out.loc['a', 'Transcripts_in_cell_status'] == FAIL
        assert out.loc['a', 'filter_result'] == RESULT_ARTIFACT

    def test_a_measured_zero_is_still_judged(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [5000],
                            'RTS_prop_in_cell': [0.0]})
        out = _verdict(s, {'RTS_prop_in_cell': [0, 5]})
        assert out.loc['a', 'RTS_prop_in_cell_status'] == PASS

    def test_a_real_value_outside_the_range_still_fails(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [5000],
                            'RTS_prop_in_cell': [40.0]})
        out = _verdict(s, {'RTS_prop_in_cell': [0, 5]})
        assert out.loc['a', 'RTS_prop_in_cell_status'] == FAIL
        assert out.loc['a', 'filter_result'] == RESULT_ARTIFACT

    def test_the_depth_column_follows_the_mode(self):
        assert cell_filter.depth_column('reads') == 'Reads_in_cell'
        assert cell_filter.depth_column('isoforms') == 'Transcripts_in_cell'


class TestModeDetection:
    """The two depth columns are mutually exclusive and mode-specific, so the cell
    summary states its own mode and --mode does not exist."""

    def test_isoforms_summary_detects_isoforms(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [5]})
        assert cell_filter.detect_mode(s, 's1') == 'isoforms'

    def test_reads_summary_detects_reads(self):
        s = _summary_frame({'CB': ['a'], 'Reads_in_cell': [5]})
        assert cell_filter.detect_mode(s, 's1') == 'reads'

    def test_neither_column_is_a_clear_error_naming_both(self):
        s = _summary_frame({'CB': ['a'], 'MT_perc': [1.0]})
        with pytest.raises(ValueError, match="Reads_in_cell.*Transcripts_in_cell"):
            cell_filter.detect_mode(s, 's1')

    def test_both_columns_is_a_clear_error(self):
        s = _summary_frame({'CB': ['a'], 'Reads_in_cell': [5], 'Transcripts_in_cell': [5]})
        with pytest.raises(ValueError, match="both"):
            cell_filter.detect_mode(s, 's1')


class TestMeasuredEvidence:
    """The report gates three sections on CLI booleans naming SQANTI3 input files. The
    filter never re-runs SQANTI3, so it reads the same answer off the QC output, where
    cell_metrics.py wrote the whole family as NA when the evidence was absent."""

    def _summary(self, **overrides):
        base = {
            'CB': ['a', 'b'],
            'Transcripts_in_cell': [500, 500],
            'CAGE_peak_support_prop': [np.nan, np.nan],
            'PolyA_motif_support_prop': [np.nan, np.nan],
            'FSM_coding_prop': [np.nan, np.nan],
            'FSM_non_coding_prop': [np.nan, np.nan],
        }
        base.update(overrides)
        return _summary_frame(base)

    def test_all_NA_evidence_reads_as_not_measured(self):
        assert cell_filter.detect_measured_evidence(self._summary()) == {
            'CAGE_peak': False, 'polyA_motif_list': False, 'include_ORF': False}

    def test_a_measured_family_reads_as_measured(self):
        s = self._summary(CAGE_peak_support_prop=[10.0, 80.0])
        assert cell_filter.detect_measured_evidence(s)['CAGE_peak'] is True

    def test_a_measured_all_zero_family_still_reads_as_measured(self):
        s = self._summary(PolyA_motif_support_prop=[0.0, 0.0])
        assert cell_filter.detect_measured_evidence(s)['polyA_motif_list'] is True

    def test_ORF_is_detected_from_the_coding_columns(self):
        s = self._summary(FSM_coding_prop=[60.0, 40.0])
        assert cell_filter.detect_measured_evidence(s)['include_ORF'] is True

    def test_an_absent_column_reads_as_not_measured(self):
        s = _summary_frame({'CB': ['a'], 'Transcripts_in_cell': [500]})
        assert cell_filter.detect_measured_evidence(s) == {
            'CAGE_peak': False, 'polyA_motif_list': False, 'include_ORF': False}


class TestBarcodelessRows:
    """classification_enrichment.py fills a missing barcode with the literal 'NA', so
    cell_metrics.py aggregates every unbarcoded read into one row that looks like a huge
    cell. Nothing about it can be judged, so it is dropped before any statistic."""

    def test_a_barcodeless_row_is_dropped_and_reported(self):
        messages = []
        s = _summary_frame({'CB': ['unassigned', 'a'], 'Transcripts_in_cell': [999999, 5000]})
        out = decide_cells(s, {'Transcripts_in_cell': 10}, mode='isoforms', sampleID='s1',
                           log=messages.append).set_index('CB')
        assert 'unassigned' not in out.index
        assert any('no cell barcode' in m for m in messages)

    def test_every_placeholder_form_is_dropped(self):
        s = _summary_frame({'CB': ['unassigned', 'NA', '-', '*', 'a'],
                            'Transcripts_in_cell': [10, 10, 10, 10, 5000]})
        out = _verdict(s, {'Transcripts_in_cell': 10})
        assert set(out.index) == {'a'}

    def test_a_barcodeless_row_does_not_reach_the_statistics(self):
        """It would be rank 1 on any depth curve, so it must not be judged or counted."""
        s = _summary_frame({'CB': ['NA', 'a', 'b'],
                            'Transcripts_in_cell': [999999, 5000, 5000]})
        out = _verdict(s, {'Transcripts_in_cell': 10})
        assert len(out) == 2


class TestNeverMeasuredColumnWarning:
    """A rule on an attribute the QC run never measured judges nothing, because every
    cell is NA and therefore skipped. #54 made that state exactly detectable, replacing
    the hardcoded column list and the does-not-vary heuristic this used to need."""

    def test_an_all_NA_rule_column_warns(self):
        messages = []
        s = _summary_frame({'CB': ['a', 'b'],
                            'CAGE_peak_support_prop': [np.nan, np.nan],
                            'Transcripts_in_cell': [5000, 5000]})
        decide_cells(s, {'CAGE_peak_support_prop': 1}, mode='isoforms', sampleID='s1',
                     log=messages.append)
        assert any('NA for every cell' in m for m in messages)

    def test_a_measured_column_does_not_warn(self):
        messages = []
        s = _summary_frame({'CB': ['a', 'b'],
                            'CAGE_peak_support_prop': [10.0, 80.0],
                            'Transcripts_in_cell': [5000, 5000]})
        decide_cells(s, {'CAGE_peak_support_prop': 1}, mode='isoforms', sampleID='s1',
                     log=messages.append)
        assert not any('NA for every cell' in m for m in messages)

    def test_a_measured_all_zero_column_does_not_warn(self):
        messages = []
        s = _summary_frame({'CB': ['a', 'b'],
                            'CAGE_peak_support_prop': [0.0, 0.0],
                            'Transcripts_in_cell': [5000, 5000]})
        decide_cells(s, {'CAGE_peak_support_prop': 1}, mode='isoforms', sampleID='s1',
                     log=messages.append)
        assert not any('NA for every cell' in m for m in messages)


class _FilterArgs:
    def __init__(self, qc_dir, out_dir, **kw):
        self.qc_dir = str(qc_dir)
        self.out_dir = str(out_dir)
        self.rules = filter_args.DEFAULT_CELL_RULES
        for k, v in kw.items():
            setattr(self, k, v)


@pytest.fixture
def qc_run(tmp_path):
    """A minimal finished QC run: one isoforms-mode sample on disk, in its own qc/
    directory so the filter's writes are visibly separate from what it reads."""
    qc_dir = tmp_path / "qc"
    sample_dir = qc_dir / "rep1"
    sample_dir.mkdir(parents=True)
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
        'Transcripts_in_cell': [5000, 5],
        'Annotated_genes': [1200, 2],
        'MT_perc': [3.0, 0.0],
        'RTS_prop_in_cell': [0.5, 0.0],
        'Intrapriming_prop_in_cell': [1.0, 0.0],
        'Non_canonical_prop_in_cell': [2.0, 0.0],
        'total_transcripts_no_monoexon': [4000, 4],
    }).to_csv(f"{prefix}_SQANTI_cell_summary.txt.gz", sep='\t', index=False,
              compression='gzip')
    with open(f"{prefix}_corrected.gtf", 'w') as fh:
        for iso in ('PB.1.1', 'PB.2.1'):
            for feat in ('transcript', 'exon'):
                fh.write('1\tPacBio\t%s\t1\t100\t.\t+\t.\ttranscript_id "%s"; gene_id "g";\n'
                         % (feat, iso))
    with open(f"{prefix}_corrected.fasta", 'w') as fh:
        for iso, seq in (('PB.1.1', 'ACGT'), ('PB.2.1', 'TTTT')):
            fh.write('>%s\n%s\n%s\n' % (iso, seq, seq))
    design = tmp_path / "design.csv"
    design.write_text("sampleID,file_acc\ns1,rep1\n")
    return tmp_path, qc_dir, str(design)


class TestRunCellFilter:
    def _run(self, tmp_path, qc_dir, design, out_name="filter", **kw):
        df = filter_io.read_sample_table(design, str(qc_dir))
        args = _FilterArgs(qc_dir, tmp_path / out_name, **kw)
        return cell_filter.run_cell_filter(args, df, log=lambda *a: None)

    def test_verdict_artifacts_are_written_by_default(self, qc_run):
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design)
        out = tmp_path / "filter" / "rep1"
        for name in ('s1_CellFilter_cell_summary.txt.gz', 's1_pass_cells.txt',
                     's1_cell_filtering_reasons.txt', 's1_cell_filter_params.txt'):
            assert (out / name).exists()

    def test_defaults_discard_the_shallow_cell(self, qc_run):
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design)
        passing = (tmp_path / "filter" / "rep1" / "s1_pass_cells.txt").read_text().split()
        assert passing == ['bc_good']

    def test_verdict_table_keeps_every_barcode(self, qc_run):
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design)
        verdict = pd.read_csv(
            tmp_path / "filter" / "rep1" / "s1_CellFilter_cell_summary.txt.gz", sep='\t')
        assert sorted(verdict['CB']) == ['bc_good', 'bc_shallow']
        assert set(verdict['filter_result']) == {'Cell', 'Artifact'}

    def test_the_filtered_dataset_is_always_materialised(self, qc_run):
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design)
        cls = _read_tsv(tmp_path / "filter" / "rep1" / "s1_classification.txt")
        assert cls['CB'].tolist() == ['bc_good']
        assert cls['FL'].tolist() == ['500']

    def test_the_data_files_keep_their_QC_names(self, qc_run):
        """What lets the existing report and the multisample script run on the filter's
        output with the original design and no path rewriting."""
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design)
        out_names = {p.name for p in (tmp_path / "filter" / "rep1").iterdir()}
        assert {'s1_classification.txt', 's1_junctions.txt'} <= out_names

    def test_the_cell_summary_is_labelled_not_subset(self, qc_run):
        """SQANTI3 labels the table it judges and never emits a subset of it. Writing
        both would be the same barcodes twice, once with the verdict and once without."""
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design)
        out = tmp_path / "filter" / "rep1"
        assert (out / "s1_CellFilter_cell_summary.txt.gz").exists()
        assert not (out / "s1_SQANTI_cell_summary.txt.gz").exists()
        labelled = pd.read_csv(out / "s1_CellFilter_cell_summary.txt.gz", sep='\t')
        assert sorted(labelled['CB']) == ['bc_good', 'bc_shallow']

    def test_the_corrected_gtf_is_subset_to_surviving_models(self, qc_run):
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design)
        gtf = (tmp_path / "filter" / "rep1" / "s1_corrected.gtf").read_text()
        assert 'PB.1.1' in gtf
        assert 'PB.2.1' not in gtf

    def test_the_corrected_fasta_is_subset_to_surviving_models(self, qc_run):
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design)
        fasta = (tmp_path / "filter" / "rep1" / "s1_corrected.fasta").read_text()
        assert '>PB.1.1' in fasta
        assert '>PB.2.1' not in fasta
        # the sequence lines of a kept record must survive with it
        assert fasta.count('ACGT') == 2

    def test_mode_is_detected_and_returned(self, qc_run):
        tmp_path, qc_dir, design = qc_run
        mode, _ = self._run(tmp_path, qc_dir, design)
        assert mode == 'isoforms'

    def test_params_file_records_the_lost_transcript_models(self, qc_run):
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design)
        params = dict(
            line.split('\t', 1) for line in
            (tmp_path / "filter" / "rep1" / "s1_cell_filter_params.txt")
            .read_text().strip().split('\n'))
        assert params['TranscriptModelsLostAllSupport'] == '1'
        assert params['BarcodesPassing'] == '1'
        assert params['Mode'] == 'isoforms'

    def test_a_different_out_dir_makes_trials_coexist(self, qc_run):
        tmp_path, qc_dir, design = qc_run
        self._run(tmp_path, qc_dir, design, out_name="trial_a")
        self._run(tmp_path, qc_dir, design, out_name="trial_b")
        assert (tmp_path / "trial_a" / "rep1" / "s1_pass_cells.txt").exists()
        assert (tmp_path / "trial_b" / "rep1" / "s1_pass_cells.txt").exists()


class TestNonDestructiveness:
    """Every existing output, and the user's design CSV, must be byte-identical after
    a full --filter_outputs run. qc_io.fill_design_table rewrites its input in place, which is
    exactly the behaviour the filter must not inherit."""

    def test_inputs_are_untouched(self, qc_run):
        import hashlib
        tmp_path, qc_dir, design = qc_run
        targets = [
            qc_dir / "rep1" / "s1_classification.txt",
            qc_dir / "rep1" / "s1_junctions.txt",
            qc_dir / "rep1" / "s1_SQANTI_cell_summary.txt.gz",
            tmp_path / "design.csv",
        ]
        before = {p: hashlib.md5(p.read_bytes()).hexdigest() for p in targets}
        df = filter_io.read_sample_table(design, str(qc_dir))
        cell_filter.run_cell_filter(
            _FilterArgs(qc_dir, tmp_path / "filter"), df, log=lambda *a: None)
        after = {p: hashlib.md5(p.read_bytes()).hexdigest() for p in targets}
        assert before == after


class TestDownstreamOrdering:
    """Clustering must run BEFORE the report. SQANTI-sc_report.R finds umap_results.csv
    by looking next to its own output rather than using the --clustering path it is
    given, so the file has to exist by the time the report runs. Reversing these two
    produces a filtered report with no UMAP section and no error -- which is exactly
    what the first reads-mode end-to-end run turned up.
    """

    def _args(self, **kw):
        a = _FilterArgs('qc', 'out')
        a.report = 'skip'
        a.multisample_report = False
        a.run_clustering = False
        for k, v in kw.items():
            setattr(a, k, v)
        return a

    def test_clustering_runs_before_the_report(self, monkeypatch):
        import filter_pipeline
        calls = []
        monkeypatch.setitem(
            sys.modules, 'sc_clustering',
            type(sys)('sc_clustering'))
        sys.modules['sc_clustering'].run_clustering_analysis = \
            lambda a, r: calls.append('clustering')
        monkeypatch.setitem(sys.modules, 'qc_reports', type(sys)('qc_reports'))
        sys.modules['qc_reports'].generate_report = lambda a, d: calls.append('report')
        sys.modules['qc_reports'].generate_multisample_report = \
            lambda a, d: calls.append('multisample')

        df = pd.DataFrame({'sampleID': ['s1'], 'file_acc': ['rep1']})
        filter_pipeline._run_downstream(
            self._args(run_clustering=True, report='html', multisample_report=True), df)
        assert calls == ['clustering', 'report', 'multisample']

    def test_clustering_is_skipped_unless_asked(self, monkeypatch):
        import filter_pipeline
        calls = []
        monkeypatch.setitem(sys.modules, 'sc_clustering', type(sys)('sc_clustering'))
        sys.modules['sc_clustering'].run_clustering_analysis = \
            lambda a, r: calls.append('clustering')
        monkeypatch.setitem(sys.modules, 'qc_reports', type(sys)('qc_reports'))
        sys.modules['qc_reports'].generate_report = lambda a, d: calls.append('report')

        df = pd.DataFrame({'sampleID': ['s1'], 'file_acc': ['rep1']})
        filter_pipeline._run_downstream(self._args(report='html'), df)
        assert calls == ['report']


class TestClusteringArgs:
    def test_the_filter_takes_the_full_clustering_set(self):
        """Not just the toggle: a different cell set has to be re-clustered with the
        same settings the QC run used, or the two are not comparable."""
        parser = filter_args.build_filter_parser()
        ns = parser.parse_args(['cells', '-de', 'd.csv', '-q', 'qc', '--run_clustering',
                                '--resolution', '0.9', '--n_neighbors', '30'])
        assert ns.run_clustering is True
        assert ns.resolution == 0.9
        assert ns.n_neighbors == 30
        for flag in ('normalization', 'n_pc', 'n_top_genes', 'clustering_method', 'n_clusters'):
            assert hasattr(ns, flag), flag

    def test_both_parsers_emit_identical_clustering_flags(self):
        """One factory, so the QC and filter parsers cannot drift apart."""
        from qc_args import build_parser
        def clustering_flags(p):
            return sorted(a.option_strings[0] for a in p._actions
                          if a.option_strings and a.option_strings[0] in (
                              '--run_clustering', '--normalization', '--n_neighbors',
                              '--n_pc', '--resolution', '--n_top_genes',
                              '--clustering_method', '--n_clusters'))
        qc = clustering_flags(build_parser())
        filt = clustering_flags(
            filter_args.build_filter_parser()._subparsers._group_actions[0].choices['cells'])
        assert qc == filt
        assert len(qc) == 8


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

    def test_every_default_rule_is_numeric(self):
        rules = cell_filter.load_rules(filter_args.DEFAULT_CELL_RULES, 'isoforms')
        for column, rule in rules.items():
            assert cell_filter.describe_rule(column, rule)

    def test_a_rules_file_without_the_all_key_is_rejected(self, tmp_path):
        path = tmp_path / "r.json"
        path.write_text('{"Transcripts_in_cell": 10}')
        with pytest.raises(ValueError, match='top-level "all" key'):
            cell_filter.load_rules(str(path), 'isoforms')
