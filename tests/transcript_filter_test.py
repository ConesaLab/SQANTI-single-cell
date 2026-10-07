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

import argparse
import hashlib
import subprocess

import pandas as pd
import pytest

sqanti_sc_src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../src"))
if sqanti_sc_src_path not in sys.path:
    sys.path.insert(0, sqanti_sc_src_path)

import filter_args
import filter_io
import transcript_filter


def _cls_row(isoform, cb, fl, category="full-splice_match", perc_a="0", min_cov="NA"):
    return {
        "isoform": isoform, "CB": cb, "FL": fl,
        "structural_category": category, "associated_gene": "geneA",
        "associated_transcript": "txA", "exons": "2", "length": "500",
        "ref_length": "600", "chrom": "chr1", "subcategory": "reference_match",
        "all_canonical": "canonical", "RTS_stage": "FALSE", "predicted_NMD": "NA",
        "within_CAGE_peak": "NA", "polyA_motif_found": "NA",
        "perc_A_downstream_TTS": perc_a, "diff_to_gene_TSS": "0",
        "coding": "non_coding", "min_cov": min_cov, "ratio_TSS": "NA",
    }


MODELS = ('PB.1.1', 'PB.2.1', 'PB.3.1')
# What the fake SQANTI3 keeps: PB.2.1 and PB.3.1 are intra-primed. bc3 is observed
# only in PB.3.1, so it keeps no model at all.
PASSING = ('PB.1.1',)


def _write_sample(sample_dir, labelled=False):
    sample_dir.mkdir(parents=True)
    prefix = sample_dir / "s1"
    pd.DataFrame([
        _cls_row('PB.1.1', 'bc1,bc2', '5,3'),
        _cls_row('PB.2.1', 'bc2', '7', category='novel_in_catalog', perc_a='80'),
        _cls_row('PB.3.1', 'bc3', '2', perc_a='90'),
    ]).to_csv(f"{prefix}_classification.txt", sep='\t', index=False)
    pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1'],
        'junction_category': ['known', 'novel', 'known'],
        'canonical': ['canonical'] * 3,
    }).to_csv(f"{prefix}_junctions.txt", sep='\t', index=False)
    summary = pd.DataFrame({
        'CB': ['bc1', 'bc2', 'bc3'],
        'Transcripts_in_cell': [5, 10, 2],
    })
    if labelled:
        summary['filter_result'] = ['Cell', 'Cell', 'Cell']
        summary.loc[len(summary)] = ['bc_dropped', 1, 'Artifact']
        summary.to_csv(f"{prefix}_CellFilter_cell_summary.txt.gz", sep='\t',
                       index=False, compression='gzip')
    else:
        summary.to_csv(f"{prefix}_SQANTI_cell_summary.txt.gz", sep='\t',
                       index=False, compression='gzip')
    with open(f"{prefix}_corrected.gtf", 'w') as fh:
        for iso in MODELS:
            for feat in ('transcript', 'exon'):
                fh.write('1\tPacBio\t%s\t1\t100\t.\t+\t.\ttranscript_id "%s"; gene_id "g";\n'
                         % (feat, iso))
    with open(f"{prefix}_corrected.fasta", 'w') as fh:
        for iso in MODELS:
            fh.write('>%s\nACGT\n' % iso)


@pytest.fixture
def qc_run(tmp_path):
    """A minimal finished isoforms-mode QC run on disk, in its own qc/ directory."""
    qc_dir = tmp_path / "qc"
    _write_sample(qc_dir / "rep1")
    design = tmp_path / "design.csv"
    design.write_text("sampleID,file_acc\ns1,rep1\n")
    return tmp_path, qc_dir, str(design)


@pytest.fixture
def fake_sqanti3(monkeypatch):
    """Stands in for sqanti3_filter.py: records the command, writes the inclusion list
    the wrapper reads, and writes its labelled classification the way SQANTI3 does --
    through pandas, so NA, TRUE and integers come back rewritten."""
    calls = []

    def run(cmd, check=False, **kw):
        calls.append(cmd)
        out_dir, prefix = cmd[cmd.index('-d') + 1], cmd[cmd.index('-o') + 1]
        with open(os.path.join(out_dir, f"{prefix}_pass_isoforms.txt"), 'w') as fh:
            fh.write(''.join(f"{iso}\n" for iso in PASSING))
        cls = pd.read_csv(cmd[cmd.index('--sqanti_class') + 1], sep='\t')
        cls['filter_result'] = ['Isoform' if i in PASSING else 'Artifact'
                                for i in cls['isoform']]
        cls.to_csv(os.path.join(out_dir, f"{prefix}_RulesFilter_classification.txt"),
                   sep='\t', index=False)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(transcript_filter.subprocess, 'run', run)
    return calls


def _args(qc_dir, out_dir, *extra):
    argv = ['transcripts', '-de', 'design.csv', '-q', str(qc_dir), '-d', str(out_dir), *extra]
    return filter_args.build_filter_parser(
        method=filter_args.transcript_method(argv)).parse_args(argv)


def _run(tmp_path, qc_dir, design, *extra, log=None):
    args = _args(qc_dir, tmp_path / "filter", *extra)
    df = filter_io.read_sample_table(design, str(qc_dir), labelled_summary_ok=True)
    return transcript_filter.run_transcript_filter(
        args, df, log=log if log is not None else (lambda *a: None))


def _read_tsv(path):
    return pd.read_csv(path, sep='\t', dtype=str, na_filter=False)


class TestSqanti3Command:
    def _cmd(self, *extra):
        args = _args('qc', 'out', *extra)
        return transcript_filter.build_sqanti3_filter_command(
            args, 'qc/rep1/s1_classification.txt', 'out/rep1', 's1')

    def test_runs_the_rules_filter_on_the_sample_classification(self):
        cmd = self._cmd()
        assert cmd[1].endswith('sqanti3_filter.py')
        assert cmd[2] == 'rules'
        assert cmd[cmd.index('--sqanti_class') + 1] == os.path.abspath(
            'qc/rep1/s1_classification.txt')

    def test_writes_into_the_sample_directory_with_the_sample_prefix(self):
        cmd = self._cmd()
        assert cmd[cmd.index('-d') + 1] == os.path.abspath('out/rep1')
        assert cmd[cmd.index('-o') + 1] == 's1'

    def test_sqanti3_default_rules_apply_unless_a_file_is_given(self):
        assert '-j' not in self._cmd()
        cmd = self._cmd('--rules', 'my_rules.json')
        assert cmd[cmd.index('-j') + 1] == os.path.abspath('my_rules.json')

    def test_mono_exonic_models_are_never_discarded_wholesale(self):
        assert '-e' not in self._cmd()
        with pytest.raises(SystemExit):
            self._cmd('--filter_mono_exonic')

    def test_sqanti3_report_is_always_skipped(self):
        assert '--skip_report' in self._cmd()
        assert '--skip_report' in self._cmd('--report', 'html')


class TestRunTranscriptFilter:
    def test_classification_is_labelled_not_subset(self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        _run(tmp_path, qc_dir, design)
        out = tmp_path / "filter" / "rep1"
        cls = _read_tsv(out / "s1_RulesFilter_classification.txt")
        assert dict(zip(cls['isoform'], cls['filter_result'])) == {
            'PB.1.1': 'Isoform', 'PB.2.1': 'Artifact', 'PB.3.1': 'Artifact'}
        assert not (out / "s1_classification.txt").exists()

    def test_labelled_rows_are_the_original_rows_verbatim(self, qc_run, fake_sqanti3):
        """SQANTI3 writes its copy back through pandas: NA comes back empty, TRUE as
        True, integers as floats. cell_metrics.py compares those flags as exact strings,
        so the wrapper replaces that copy with the original rows plus the verdict."""
        tmp_path, qc_dir, design = qc_run
        _run(tmp_path, qc_dir, design)
        original = _read_tsv(qc_dir / "rep1" / "s1_classification.txt")
        labelled = _read_tsv(tmp_path / "filter" / "rep1" / "s1_RulesFilter_classification.txt")
        assert labelled.columns[-1] == 'filter_result'
        pd.testing.assert_frame_equal(labelled.drop(columns='filter_result'), original)
        assert labelled.loc[0, 'RTS_stage'] == 'FALSE'
        assert labelled.loc[0, 'min_cov'] == 'NA'

    def test_an_already_labelled_classification_is_refused(self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        path = qc_dir / "rep1" / "s1_classification.txt"
        cls = _read_tsv(path)
        cls['filter_result'] = 'Isoform'
        cls.to_csv(path, sep='\t', index=False)
        with pytest.raises(ValueError, match='already been through the transcript filter'):
            _run(tmp_path, qc_dir, design)
        assert not fake_sqanti3

    def test_junctions_gtf_and_fasta_follow_the_kept_models(self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        _run(tmp_path, qc_dir, design)
        out = tmp_path / "filter" / "rep1"
        assert _read_tsv(out / "s1_junctions.txt")['isoform'].tolist() == ['PB.1.1']
        gtf = (out / "s1_corrected.gtf").read_text()
        fasta = (out / "s1_corrected.fasta").read_text()
        for dropped in ('PB.2.1', 'PB.3.1'):
            assert dropped not in gtf
            assert dropped not in fasta
        assert gtf.count('PB.1.1') == 2
        assert '>PB.1.1\nACGT\n' in fasta

    def test_cell_summary_is_recomputed_from_the_kept_models(self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        _run(tmp_path, qc_dir, design)
        summary = pd.read_csv(
            tmp_path / "filter" / "rep1" / "s1_SQANTI_cell_summary.txt.gz", sep='\t')
        depth = dict(zip(summary['CB'], summary['Transcripts_in_cell']))
        # bc2 carried 3 + 7 before; only PB.1.1's 3 survives.
        assert depth == {'bc1': 5, 'bc2': 3}

    def test_a_cell_left_with_no_model_is_reported(self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        messages = []
        _run(tmp_path, qc_dir, design, log=messages.append)
        assert any('1 cell(s) kept no transcript model' in m for m in messages)

    def test_mode_and_evidence_are_returned(self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        mode, evidence = _run(tmp_path, qc_dir, design)
        assert mode == 'isoforms'
        assert evidence == {'CAGE_peak': False, 'polyA_motif_list': False,
                            'include_ORF': False}

    def test_one_sqanti3_run_per_sample(self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        _write_sample(qc_dir / "rep2")
        (tmp_path / "design.csv").write_text("sampleID,file_acc\ns1,rep1\ns1,rep2\n")
        _run(tmp_path, qc_dir, design)
        dirs = [cmd[cmd.index('-d') + 1] for cmd in fake_sqanti3]
        assert dirs == [os.path.abspath(tmp_path / "filter" / d) for d in ('rep1', 'rep2')]

    def test_a_sqanti3_failure_names_the_sample(self, qc_run, monkeypatch):
        tmp_path, qc_dir, design = qc_run

        def fail(cmd, check=False, **kw):
            raise subprocess.CalledProcessError(1, cmd)

        monkeypatch.setattr(transcript_filter.subprocess, 'run', fail)
        with pytest.raises(ValueError, match='sample s1'):
            _run(tmp_path, qc_dir, design)

    def test_inputs_are_untouched(self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        targets = sorted((qc_dir / "rep1").iterdir()) + [tmp_path / "design.csv"]
        before = {p: hashlib.md5(p.read_bytes()).hexdigest() for p in targets}
        _run(tmp_path, qc_dir, design)
        after = {p: hashlib.md5(p.read_bytes()).hexdigest() for p in targets}
        assert before == after
        assert sorted((qc_dir / "rep1").iterdir()) == targets[:-1]


class TestAfterTheCellFilter:
    """The intended order is cells first, then transcripts, with the cell filter's
    output directory as the input. That directory has the labelled summary and no
    plain one."""

    @pytest.fixture
    def cells_run(self, tmp_path):
        cells_dir = tmp_path / "filter_cells"
        _write_sample(cells_dir / "rep1", labelled=True)
        design = tmp_path / "design.csv"
        design.write_text("sampleID,file_acc\ns1,rep1\n")
        return tmp_path, cells_dir, str(design)

    def test_a_cell_filter_output_is_accepted_as_input(self, cells_run, fake_sqanti3):
        tmp_path, cells_dir, design = cells_run
        mode, _ = _run(tmp_path, cells_dir, design)
        assert mode == 'isoforms'
        assert (tmp_path / "filter" / "rep1" / "s1_SQANTI_cell_summary.txt.gz").exists()

    def test_cells_the_cell_filter_discarded_are_not_counted_as_lost(
            self, cells_run, fake_sqanti3):
        tmp_path, cells_dir, design = cells_run
        messages = []
        _run(tmp_path, cells_dir, design, log=messages.append)
        assert any('1 cell(s) kept no transcript model' in m for m in messages)

    def test_the_cell_filter_still_requires_the_plain_summary(self, cells_run):
        tmp_path, cells_dir, design = cells_run
        with pytest.raises(ValueError, match='_SQANTI_cell_summary'):
            filter_io.read_sample_table(design, str(cells_dir))


class TestLabelledClassificationReaders:
    """No kept-only copy exists, so every reader must skip the Artifact rows itself."""

    def _labelled(self, tmp_path):
        sample_dir = tmp_path / "rep1"
        sample_dir.mkdir()
        prefix = sample_dir / "s1"
        pd.DataFrame([
            dict(_cls_row('PB.1.1', 'bc1', '5'), filter_result='Isoform'),
            dict(_cls_row('PB.2.1', 'bc2', '7'), filter_result='Artifact'),
        ]).to_csv(f"{prefix}_RulesFilter_classification.txt", sep='\t', index=False)
        return prefix

    def test_the_labelled_file_is_found_before_the_plain_name(self, tmp_path):
        prefix = self._labelled(tmp_path)
        assert filter_io.classification_path(str(prefix)).endswith(
            's1_RulesFilter_classification.txt')
        assert filter_io.classification_path(str(tmp_path / "other")).endswith(
            'other_classification.txt')

    def test_drop_artifacts_removes_the_rows_and_the_column(self):
        df = pd.DataFrame({'isoform': ['a', 'b'], 'filter_result': ['Isoform', 'Artifact']})
        out = filter_io.drop_artifacts(df)
        assert out['isoform'].tolist() == ['a']
        assert 'filter_result' not in out.columns
        plain = pd.DataFrame({'isoform': ['a']})
        assert filter_io.drop_artifacts(plain) is plain

    def test_clustering_skips_artifacts(self, tmp_path):
        from sc_clustering import prepare_anndata
        self._labelled(tmp_path)
        args = argparse.Namespace(out_dir=str(tmp_path), mode='isoforms')
        adata = prepare_anndata(args, {'file_acc': 'rep1', 'sampleID': 's1'})
        assert list(adata.obs_names) == ['bc1']

    def test_the_report_is_handed_the_labelled_classification(self, tmp_path, monkeypatch):
        import qc_reports
        prefix = self._labelled(tmp_path)
        calls = []
        monkeypatch.setattr(qc_reports.subprocess, 'run', lambda cmd, **kw: calls.append(cmd))
        args = argparse.Namespace(out_dir=str(tmp_path), mode='isoforms', report='html')
        qc_reports.generate_report(args, pd.DataFrame({'sampleID': ['s1'], 'file_acc': ['rep1']}))
        assert f'"{prefix}_RulesFilter_classification.txt"' in calls[0]

    def _report_cmd(self, tmp_path, monkeypatch, **extra):
        import qc_reports
        calls = []
        monkeypatch.setattr(qc_reports.subprocess, 'run', lambda cmd, **kw: calls.append(cmd))
        args = argparse.Namespace(out_dir=str(tmp_path), mode='isoforms', report='html', **extra)
        qc_reports.generate_report(args, pd.DataFrame({'sampleID': ['s1'], 'file_acc': ['rep1']}))
        return calls[0]

    def test_the_transcript_filter_hands_the_report_its_input_summary_and_reasons(
            self, tmp_path, monkeypatch):
        prefix = self._labelled(tmp_path)
        open(f"{prefix}_filtering_reasons.txt", 'w').close()
        qc_dir = tmp_path / "qc"
        cmd = self._report_cmd(tmp_path, monkeypatch, subcommand='transcripts',
                               method='rules', qc_dir=str(qc_dir))
        assert f'--input_cell_summary "{qc_dir}/rep1/s1_SQANTI_cell_summary.txt.gz"' in cmd
        assert f'--transcript_filter_reasons "{prefix}_filtering_reasons.txt"' in cmd

    def test_a_cell_filter_input_is_passed_as_its_labelled_summary(self, tmp_path, monkeypatch):
        self._labelled(tmp_path)
        qc_sample = tmp_path / "cells" / "rep1"
        qc_sample.mkdir(parents=True)
        (qc_sample / "s1_CellFilter_cell_summary.txt.gz").touch()
        cmd = self._report_cmd(tmp_path, monkeypatch, subcommand='transcripts',
                               method='rules', qc_dir=str(tmp_path / "cells"))
        assert f'--input_cell_summary "{qc_sample}/s1_CellFilter_cell_summary.txt.gz"' in cmd

    @pytest.mark.parametrize('extra', [{}, {'subcommand': 'cells', 'qc_dir': 'qc'}])
    def test_other_callers_get_no_transcript_filter_inputs(self, tmp_path, monkeypatch, extra):
        prefix = self._labelled(tmp_path)
        open(f"{prefix}_filtering_reasons.txt", 'w').close()
        cmd = self._report_cmd(tmp_path, monkeypatch, **extra)
        assert '--input_cell_summary' not in cmd
        assert '--transcript_filter_reasons' not in cmd

    def test_the_ml_filter_hands_the_report_its_directory_not_reasons(
            self, tmp_path, monkeypatch):
        prefix = self._labelled(tmp_path)
        open(f"{prefix}_filtering_reasons.txt", 'w').close()
        cmd = self._report_cmd(tmp_path, monkeypatch, subcommand='transcripts',
                               method='ml', qc_dir=str(tmp_path / "qc"))
        assert f'--ml_dir "{tmp_path}/rep1"' in cmd
        assert '--transcript_filter_reasons' not in cmd

    def test_h5ad_export_skips_artifacts(self, tmp_path):
        from sc_export import _prepare_classification
        prefix = self._labelled(tmp_path)
        cls = _prepare_classification(f"{prefix}_RulesFilter_classification.txt", 'isoforms')
        assert set(cls['CB']) == {'bc1'}


class TestCellFilterAfterTranscripts:
    """Transcripts first, then cells: the cell filter reads the labelled classification
    and must keep its verdict column, so the readers downstream still skip artifacts."""

    def test_the_verdict_survives_the_cell_filter(self, qc_run, fake_sqanti3):
        import cell_filter
        tmp_path, qc_dir, design = qc_run
        _run(tmp_path, qc_dir, design)
        rules = tmp_path / "rules.json"
        rules.write_text('{"all": [{"depth": 1}]}')
        args = filter_args.build_filter_parser().parse_args(
            ['cells', '-de', design, '-q', str(tmp_path / "filter"),
             '-d', str(tmp_path / "cells"), '-j', str(rules)])
        df = filter_io.read_sample_table(design, args.qc_dir)
        cell_filter.run_cell_filter(args, df, log=lambda *a: None)
        cls = _read_tsv(tmp_path / "cells" / "rep1" / "s1_RulesFilter_classification.txt")
        # bc3 was lost to the transcript filter, so PB.3.1 has no cell left; PB.2.1 is
        # an artifact of a kept cell and keeps its label.
        assert dict(zip(cls['isoform'], cls['filter_result'])) == {
            'PB.1.1': 'Isoform', 'PB.2.1': 'Artifact'}
        assert not (tmp_path / "cells" / "rep1" / "s1_classification.txt").exists()


def _write_optional_files(prefix):
    """The four per-model files SQANTI3 QC writes only when asked, in their real formats."""
    with open(f"{prefix}_corrected.faa", 'w') as fh:
        for iso in MODELS:
            fh.write(f">{iso}\t{iso}.p1|169_aa|+|1|507\nMKV\n")
    with open(f"{prefix}_corrected.cds.gff3", 'w') as fh:
        for iso in MODELS:
            fh.write(f'1\tPacBio\tCDS\t1\t90\t.\t+\t.\ttranscript_id "{iso}"; gene_id "g";\n')
    with open(f"{prefix}_corrected.sam", 'w') as fh:
        fh.write("@HD\tVN:1.6\n@SQ\tSN:1\tLN:1000\n")
        for iso in MODELS:
            fh.write(f"{iso}\t0\t1\t1\t60\t90M\t*\t0\t0\tACGT\t*\n")
    with open(f"{prefix}.gff3", 'w') as fh:
        for iso in MODELS:
            fh.write(f"{iso}\ttappAS\ttranscript\t1\t90\t.\t+\t.\tID={iso}\n")


def _models_in(path):
    text = open(path).read()
    return {iso for iso in MODELS if iso in text}


class TestOptionalModelFiles:
    SUFFIXES = ('_corrected.faa', '_corrected.cds.gff3', '_corrected.sam', '.gff3')

    def test_first_field_subset_keeps_headers(self, tmp_path):
        src, dst = tmp_path / "in.sam", tmp_path / "out.sam"
        src.write_text("@HD\tVN:1.6\nr1\t0\t1\nr2\t0\t1\n")
        filter_io.subset_by_first_field(str(src), str(dst), {'r2'})
        assert dst.read_text() == "@HD\tVN:1.6\nr2\t0\t1\n"

    def test_only_files_present_are_written(self, tmp_path):
        (tmp_path / "in_corrected.sam").write_text("@HD\nr1\t0\n")
        written = filter_io.subset_optional_model_files(
            str(tmp_path / "in"), str(tmp_path / "out"), {'r1'})
        assert written == ['_corrected.sam']
        assert sorted(p.name for p in tmp_path.iterdir()) == [
            'in_corrected.sam', 'out_corrected.sam']

    def test_the_transcript_filter_cuts_them_to_the_passing_models(self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        _write_optional_files(qc_dir / "rep1" / "s1")
        _run(tmp_path, qc_dir, design)
        for suffix in self.SUFFIXES:
            assert _models_in(tmp_path / "filter" / "rep1" / f"s1{suffix}") == set(PASSING), suffix
        assert (tmp_path / "filter" / "rep1" / "s1_corrected.sam").read_text().startswith("@HD")

    def test_the_cell_filter_cuts_them_to_the_surviving_models(self, qc_run):
        import cell_filter
        tmp_path, qc_dir, design = qc_run
        _write_optional_files(qc_dir / "rep1" / "s1")
        rules = tmp_path / "rules.json"
        rules.write_text('{"all": [{"depth": 3}]}')
        args = filter_args.build_filter_parser().parse_args(
            ['cells', '-de', design, '-q', str(qc_dir), '-d', str(tmp_path / "cells"),
             '-j', str(rules)])
        cell_filter.run_cell_filter(args, filter_io.read_sample_table(design, str(qc_dir)),
                                    log=lambda *a: None)
        # bc3 (depth 2) is discarded, and PB.3.1 was seen only in bc3.
        for suffix in self.SUFFIXES:
            assert _models_in(tmp_path / "cells" / "rep1" / f"s1{suffix}") == {
                'PB.1.1', 'PB.2.1'}, suffix


class TestTranscriptsParser:
    def test_defaults(self):
        ns = _args('qc', 'out')
        assert ns.subcommand == 'transcripts'
        assert ns.rules is None
        assert ns.report == 'skip'

    def test_both_parsers_emit_identical_cell_metric_flags(self):
        """One factory, so the recomputed summary can use the QC run's settings."""
        from qc_args import build_parser

        def flags(parser):
            return sorted((a.option_strings[0], a.dest, a.default) for a in parser._actions
                          if a.option_strings and a.option_strings[0] in (
                              '--min_cov', '--ratio_TSS', '--ref_cov_min_pct'))

        transcripts = filter_args.build_filter_parser()._subparsers._group_actions[0] \
            .choices['transcripts']
        assert flags(build_parser()) == flags(transcripts)
        assert len(flags(transcripts)) == 3

    def test_takes_the_clustering_and_report_options(self):
        ns = _args('qc', 'out', '--run_clustering', '--resolution', '0.9',
                   '--report', 'html', '--multisample_report')
        assert ns.run_clustering is True
        assert ns.resolution == 0.9
        assert ns.multisample_report is True


# What the fake ML filter decides, as SQANTI3's R script would: PB.2.1 is a classifier
# negative and intra-primed, PB.3.1 only intra-primed.
ML_VERDICTS = {
    'PB.1.1': ('0.91', '0.09', 'Positive', 'FALSE'),
    'PB.2.1': ('0.12', '0.88', 'Negative', 'TRUE'),
    'PB.3.1': ('0.80', '0.20', 'Positive', 'TRUE'),
}


@pytest.fixture
def fake_sqanti3_ml(monkeypatch):
    """Stands in for sqanti3_filter.py ml: records the command and the input it was
    handed, and writes its classification as the R script does -- the input it read plus
    the four ML columns and the verdict."""
    calls = []

    def run(cmd, check=False, **kw):
        out_dir, prefix = cmd[cmd.index('-d') + 1], cmd[cmd.index('-o') + 1]
        ml_input = _read_tsv(cmd[cmd.index('--sqanti_class') + 1])
        calls.append((cmd, ml_input))
        with open(os.path.join(out_dir, f"{prefix}_pass_isoforms.txt"), 'w') as fh:
            fh.write(''.join(f"{iso}\n" for iso in PASSING))
        out = ml_input.copy()
        for i, column in enumerate(transcript_filter.ML_COLUMNS):
            out[column] = [ML_VERDICTS[iso][i] for iso in out['isoform']]
        out['filter_result'] = ['Isoform' if i in PASSING else 'Artifact'
                                for i in out['isoform']]
        out.to_csv(os.path.join(out_dir, f"{prefix}_ML_classification.txt"),
                   sep='\t', index=False)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(transcript_filter.subprocess, 'run', run)
    return calls


class TestMethodOption:
    def test_the_method_is_read_before_the_parser_is_built(self):
        assert filter_args.transcript_method(['transcripts', '--method', 'ml']) == 'ml'
        assert filter_args.transcript_method(['transcripts']) == 'rules'

    def test_an_unknown_method_is_rejected_by_the_parser(self):
        assert filter_args.transcript_method(['transcripts', '--method', 'svm']) == 'rules'
        with pytest.raises(SystemExit):
            _args('qc', 'out', '--method', 'svm')

    def test_j_is_the_rules_json_under_rules(self):
        ns = _args('qc', 'out', '-j', 'rules.json')
        assert ns.method == 'rules'
        assert ns.rules == 'rules.json'

    def test_j_is_the_probability_threshold_under_ml(self):
        ns = _args('qc', 'out', '--method', 'ml', '-j', '0.8')
        assert ns.threshold == 0.8
        assert not hasattr(ns, 'rules')
        with pytest.raises(SystemExit):
            _args('qc', 'out', '--method', 'ml', '-j', 'rules.json')

    @pytest.mark.parametrize('option', [['-t', '0.5'], ['-f'], ['--TP', 'tp.txt'],
                                        ['-i', '70']])
    def test_ml_options_are_refused_under_rules(self, option):
        with pytest.raises(SystemExit):
            _args('qc', 'out', *option)


class TestMLCommand:
    def _cmd(self, *extra):
        args = _args('qc', 'out', '--method', 'ml', *extra)
        return transcript_filter.build_sqanti3_filter_command(
            args, 'out/rep1/s1_ML_input.tmp', 'out/rep1', 's1')

    def test_runs_the_ml_filter_with_the_report_skipped(self):
        cmd = self._cmd()
        assert cmd[1].endswith('sqanti3_filter.py')
        assert cmd[2] == 'ml'
        assert cmd[cmd.index('-o') + 1] == 's1'
        assert '--skip_report' in cmd

    def test_sqanti3_defaults_apply_unless_an_option_is_given(self):
        cmd = self._cmd()
        for flag in ('-j', '-t', '-p', '-n', '-f', '-r', '-z', '-i', '-e',
                     '--intermediate_files'):
            assert flag not in cmd, flag

    def test_options_are_passed_under_sqanti3_flags(self):
        cmd = self._cmd('-j', '0.8', '-t', '0.7', '-z', '2000', '-i', '70', '-f',
                        '--intermediate_files', '-p', 'tp.txt', '-n', 'tn.txt',
                        '-r', 'drop.txt')
        values = {flag: cmd[cmd.index(flag) + 1] for flag in
                  ('-j', '-t', '-z', '-i', '-p', '-n', '-r')}
        assert values == {'-j': '0.8', '-t': '0.7', '-z': '2000', '-i': '70.0',
                          '-p': os.path.abspath('tp.txt'), '-n': os.path.abspath('tn.txt'),
                          '-r': os.path.abspath('drop.txt')}
        assert '-f' in cmd and '--intermediate_files' in cmd


class TestMLInput:
    def test_fl_becomes_the_models_total_over_its_cells(self, tmp_path):
        src, dst = tmp_path / "in.txt", tmp_path / "out.txt"
        pd.DataFrame([_cls_row('PB.1.1', 'bc1,bc2', '5,3'), _cls_row('PB.2.1', 'bc2', '7'),
                      _cls_row('PB.3.1', 'bc1,bc3', '0.5,1.25'),
                      _cls_row('PB.4.1', 'bc1', 'NA')]).to_csv(src, sep='\t', index=False)
        transcript_filter.write_ml_input(str(src), str(dst), 'isoforms')
        assert _read_tsv(dst)['FL'].tolist() == ['8', '7', '1.75', 'NA']

    def test_cell_and_junction_chain_columns_are_left_out(self, tmp_path):
        src, dst = tmp_path / "in.txt", tmp_path / "out.txt"
        rows = pd.DataFrame([dict(_cls_row('r1', 'bc1', 'NA'), UMI='ACGT',
                                  jxn_string='chr1_+_10_20', jxnHash='abc')])
        rows.to_csv(src, sep='\t', index=False)
        transcript_filter.write_ml_input(str(src), str(dst), 'reads')
        out = _read_tsv(dst)
        expected = rows.drop(columns=list(transcript_filter.NOT_ML_FEATURES))
        pd.testing.assert_frame_equal(out, expected)


class TestRunMLFilter:
    def _run_ml(self, tmp_path, qc_dir, design, *extra):
        return _run(tmp_path, qc_dir, design, '--method', 'ml', *extra)

    def test_sqanti3_reads_the_reduced_input_which_is_then_removed(
            self, qc_run, fake_sqanti3_ml):
        tmp_path, qc_dir, design = qc_run
        self._run_ml(tmp_path, qc_dir, design)
        (cmd, ml_input), = fake_sqanti3_ml
        assert 'CB' not in ml_input.columns
        assert ml_input['FL'].tolist() == ['8', '7', '2']
        assert not os.path.exists(cmd[cmd.index('--sqanti_class') + 1])

    def test_classification_is_the_input_plus_sqanti3_ml_columns_and_verdict(
            self, qc_run, fake_sqanti3_ml):
        tmp_path, qc_dir, design = qc_run
        self._run_ml(tmp_path, qc_dir, design)
        out = tmp_path / "filter" / "rep1"
        original = _read_tsv(qc_dir / "rep1" / "s1_classification.txt")
        labelled = _read_tsv(out / "s1_ML_classification.txt")
        added = [*transcript_filter.ML_COLUMNS, 'filter_result']
        assert list(labelled.columns) == list(original.columns) + added
        pd.testing.assert_frame_equal(labelled.drop(columns=added), original)
        assert labelled.loc[2, list(transcript_filter.ML_COLUMNS)].tolist() == list(
            ML_VERDICTS['PB.3.1'])
        assert dict(zip(labelled['isoform'], labelled['filter_result'])) == {
            'PB.1.1': 'Isoform', 'PB.2.1': 'Artifact', 'PB.3.1': 'Artifact'}
        assert sorted(p.name for p in out.glob('*classification*')) == [
            's1_ML_classification.txt']

    def test_readers_find_the_ml_classification(self, qc_run, fake_sqanti3_ml):
        tmp_path, qc_dir, design = qc_run
        self._run_ml(tmp_path, qc_dir, design)
        prefix = str(tmp_path / "filter" / "rep1" / "s1")
        assert filter_io.classification_path(prefix) == prefix + "_ML_classification.txt"
        summary = pd.read_csv(f"{prefix}_SQANTI_cell_summary.txt.gz", sep='\t')
        assert dict(zip(summary['CB'], summary['Transcripts_in_cell'])) == {'bc1': 5, 'bc2': 3}

    def test_sqanti3_rows_out_of_order_are_refused(self, qc_run, monkeypatch):
        tmp_path, qc_dir, design = qc_run

        def run(cmd, check=False, **kw):
            out_dir, prefix = cmd[cmd.index('-d') + 1], cmd[cmd.index('-o') + 1]
            out = _read_tsv(cmd[cmd.index('--sqanti_class') + 1]).iloc[::-1]
            for i, column in enumerate(transcript_filter.ML_COLUMNS):
                out[column] = [ML_VERDICTS[iso][i] for iso in out['isoform']]
            out.to_csv(os.path.join(out_dir, f"{prefix}_ML_classification.txt"),
                       sep='\t', index=False)
            open(os.path.join(out_dir, f"{prefix}_pass_isoforms.txt"), 'w').close()
            return subprocess.CompletedProcess(cmd, 0)

        monkeypatch.setattr(transcript_filter.subprocess, 'run', run)
        with pytest.raises(ValueError, match='same order'):
            self._run_ml(tmp_path, qc_dir, design)

    def test_the_input_copy_is_removed_when_sqanti3_fails(self, qc_run, monkeypatch):
        tmp_path, qc_dir, design = qc_run

        def fail(cmd, check=False, **kw):
            raise subprocess.CalledProcessError(1, cmd)

        monkeypatch.setattr(transcript_filter.subprocess, 'run', fail)
        with pytest.raises(ValueError, match='sample s1'):
            self._run_ml(tmp_path, qc_dir, design)
        assert list((tmp_path / "filter" / "rep1").glob('*.tmp')) == []

    def test_an_ml_run_refuses_a_directory_holding_a_rules_run(
            self, qc_run, fake_sqanti3):
        tmp_path, qc_dir, design = qc_run
        _run(tmp_path, qc_dir, design)
        calls = len(fake_sqanti3)
        with pytest.raises(ValueError, match='written by --method rules'):
            self._run_ml(tmp_path, qc_dir, design)
        assert len(fake_sqanti3) == calls

    def test_a_rules_run_refuses_a_directory_holding_an_ml_run(
            self, qc_run, fake_sqanti3_ml):
        tmp_path, qc_dir, design = qc_run
        self._run_ml(tmp_path, qc_dir, design)
        with pytest.raises(ValueError, match='written by --method ml'):
            _run(tmp_path, qc_dir, design)

    def test_training_lists_are_refused_for_several_samples(self, qc_run, fake_sqanti3_ml):
        tmp_path, qc_dir, design = qc_run
        _write_sample(qc_dir / "rep2")
        (tmp_path / "design.csv").write_text("sampleID,file_acc\ns1,rep1\ns1,rep2\n")
        tp = tmp_path / "tp.txt"
        tp.write_text("PB.1.1\n")
        with pytest.raises(ValueError, match='one sample'):
            self._run_ml(tmp_path, qc_dir, design, '--TP', str(tp), '--TN', str(tp))
        assert not fake_sqanti3_ml
        self._run_ml(tmp_path, qc_dir, design)
        assert len(fake_sqanti3_ml) == 2

    def test_training_lists_reach_sqanti3_for_one_sample(self, qc_run, fake_sqanti3_ml):
        tmp_path, qc_dir, design = qc_run
        tp = tmp_path / "tp.txt"
        tp.write_text("PB.1.1\n")
        self._run_ml(tmp_path, qc_dir, design, '--TP', str(tp), '--TN', str(tp))
        (cmd, _), = fake_sqanti3_ml
        assert cmd[cmd.index('-p') + 1] == str(tp)

    def test_the_log_names_the_ml_filter(self, qc_run, fake_sqanti3_ml):
        tmp_path, qc_dir, design = qc_run
        messages = []
        _run(tmp_path, qc_dir, design, '--method', 'ml', log=messages.append)
        assert any('1/3 transcript models passed the SQANTI3 ML filter' in m
                   for m in messages)


class TestEntryPoint:
    def test_main_dispatches_to_the_transcript_filter(self, qc_run, fake_sqanti3,
                                                      monkeypatch):
        import filter_pipeline
        tmp_path, qc_dir, design = qc_run
        monkeypatch.setattr(sys, 'argv', [
            'sqanti_sc_filter.py', 'transcripts', '-de', design, '-q', str(qc_dir),
            '-d', str(tmp_path / "filter")])
        filter_pipeline.main()
        assert len(fake_sqanti3) == 1
        assert (tmp_path / "filter" / "rep1" / "s1_RulesFilter_classification.txt").exists()

    def test_main_gives_j_its_ml_meaning(self, qc_run, fake_sqanti3_ml, monkeypatch):
        import filter_pipeline
        tmp_path, qc_dir, design = qc_run
        monkeypatch.setattr(sys, 'argv', [
            'sqanti_sc_filter.py', 'transcripts', '-de', design, '-q', str(qc_dir),
            '-d', str(tmp_path / "filter"), '--method', 'ml', '-j', '0.9'])
        filter_pipeline.main()
        (cmd, _), = fake_sqanti3_ml
        assert cmd[2] == 'ml'
        assert cmd[cmd.index('-j') + 1] == '0.9'
        assert (tmp_path / "filter" / "rep1" / "s1_ML_classification.txt").exists()
