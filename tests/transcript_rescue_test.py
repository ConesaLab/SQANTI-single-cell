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

import subprocess

import pandas as pd
import pytest

sqanti_sc_src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../src"))
if sqanti_sc_src_path not in sys.path:
    sys.path.insert(0, sqanti_sc_src_path)

import filter_args
import filter_io
import qc_args
import transcript_rescue


def _row(isoform, cb, fl, verdict, category='novel_in_catalog', gene='G1', ref='novel'):
    return {
        'isoform': isoform, 'CB': cb, 'FL': fl, 'structural_category': category,
        'associated_gene': gene, 'associated_transcript': ref, 'exons': '2',
        'length': '500', 'ref_length': '600', 'chrom': 'chr1', 'subcategory': 'x',
        'all_canonical': 'canonical', 'RTS_stage': 'FALSE', 'predicted_NMD': 'NA',
        'within_CAGE_peak': 'NA', 'polyA_motif_found': 'NA', 'perc_A_downstream_TTS': '0',
        'diff_to_gene_TSS': '0', 'coding': 'non_coding', 'min_cov': 'NA', 'ratio_TSS': 'NA',
        'filter_result': verdict,
    }


# Y1, Y2 and ENST9 are kept. F1 and I1 lose reference ENST1 every FSM it had. X1 maps to
# Y1 alone, X2 ties Y1 (4 reads in the sample) with Y2 (12), X3 ties two references with
# no reads, X5 maps to reference ENST9, which a kept model already carries under that
# name. X4 maps to nothing.
CLASSIFICATION = [
    _row('Y1', 'c1,c2', '3,1', 'Isoform'),
    _row('Y2', 'c3', '12', 'Isoform'),
    _row('ENST9', 'c3', '1', 'Isoform', category='full-splice_match', gene='G9', ref='ENST9'),
    _row('F1', 'c1,c2', '2,1', 'Artifact', category='full-splice_match', gene='G2', ref='ENST1'),
    _row('I1', 'c2,c4', '1,1', 'Artifact', category='incomplete-splice_match', gene='G2',
         ref='ENST1'),
    _row('X1', 'c1,c4', '1,2', 'Artifact', category='novel_not_in_catalog'),
    _row('X2', 'c2', '4', 'Artifact'),
    _row('X3', 'c5', '1', 'Artifact', gene='G3'),
    _row('X4', 'c1', '5', 'Artifact'),
    _row('X5', 'c1', '2', 'Artifact', gene='G9'),
]
RESCUE_TABLE = [
    ('F1', 'ENST1', 'automatic', 'reference', 'yes'),
    ('I1', 'ENST1', 'automatic', 'reference', 'no'),
    ('X1', 'Y1', 'rules_mapping', 'lr_defined', 'no'),
    ('X2', 'Y1', 'rules_mapping', 'lr_defined', 'no'),
    ('X2', 'Y2', 'rules_mapping', 'lr_defined', 'no'),
    ('X3', 'ENST7', 'rules_mapping', 'reference', 'yes'),
    ('X3', 'ENST8', 'rules_mapping', 'reference', 'yes'),
    ('X5', 'ENST9', 'rules_mapping', 'reference', 'yes'),
]
INCLUSION = ('ENST1', 'ENST7', 'ENST8', 'ENST9')
REFERENCES = ('ENST1', 'ENST7', 'ENST8', 'ENST9', 'ENST_unused')


def _write_filter_output(sample_dir, rows=CLASSIFICATION):
    sample_dir.mkdir(parents=True)
    prefix = str(sample_dir / 's1')
    pd.DataFrame(rows).to_csv(prefix + '_RulesFilter_classification.txt', sep='\t',
                              index=False)
    kept = [r for r in rows if r['filter_result'] == 'Isoform']
    pd.DataFrame({'isoform': [r['isoform'] for r in kept],
                  'canonical': ['canonical'] * len(kept),
                  'CB': [r['CB'] for r in kept]}).to_csv(prefix + '_junctions.txt', sep='\t',
                                                         index=False)
    pd.DataFrame({'CB': ['c1', 'c2', 'c3', 'c4'], 'Transcripts_in_cell': [4, 1, 13, 0]}) \
        .to_csv(prefix + '_SQANTI_cell_summary.txt.gz', sep='\t', index=False,
                compression='gzip')
    with open(prefix + '_corrected.gtf', 'w') as fh:
        for r in kept:
            fh.write('1\tPacBio\ttranscript\t1\t9\t.\t+\t.\ttranscript_id "%s"; gene_id "g";\n'
                     % r['isoform'])
    with open(prefix + '_corrected.fasta', 'w') as fh:
        fh.writelines('>%s\nACGT\n' % r['isoform'] for r in kept)
    return prefix


def _write_reference_qc(out_dir, prefix='reference'):
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, prefix)
    rows = []
    for iso in REFERENCES:
        r = _row(iso, '', 'NA', 'x', category='full-splice_match', gene='R', ref=iso)
        del r['CB'], r['filter_result']
        r['psauron_score'] = '0.9'
        rows.append(r)
    pd.DataFrame(rows).to_csv(base + '_classification.txt', sep='\t', index=False)
    pd.DataFrame({'isoform': ['ENST1', 'ENST1', 'ENST7'], 'canonical': ['canonical'] * 3,
                  'extra': ['e'] * 3}).to_csv(base + '_junctions.txt', sep='\t', index=False)
    with open(base + '_corrected.gtf', 'w') as fh:
        fh.write('#!genome-build GRCh38\n')
        for iso in REFERENCES:
            fh.write('1\tEnsembl\ttranscript\t1\t9\t.\t+\t.\ttranscript_id "%s"; gene_id "R";\n'
                     % iso)
    with open(base + '_corrected.fasta', 'w') as fh:
        fh.writelines('>%s\nTTTT\n' % iso for iso in REFERENCES)
    return base


@pytest.fixture
def runs(tmp_path):
    """A transcript filter's output (tf/), the QC run it read (pre/), and a design."""
    tf = tmp_path / 'tf'
    _write_filter_output(tf / 'rep1')
    pre = tmp_path / 'pre' / 'rep1'
    pre.mkdir(parents=True)
    (pre / 's1_corrected.fasta').write_text(''.join(
        '>%s\nACGT\n' % r['isoform'] for r in CLASSIFICATION))
    design = tmp_path / 'design.csv'
    design.write_text('sampleID,file_acc\ns1,rep1\n')
    return tmp_path, tf, pre.parent, str(design)


@pytest.fixture
def fake_sqanti3(monkeypatch):
    """Stands in for sqanti3_qc.py on the annotation and for sqanti3_rescue.py."""
    calls = []

    def run(cmd, check=False, cwd=None, **kw):
        out_dir, prefix = cmd[cmd.index('-d') + 1], cmd[cmd.index('-o') + 1]
        script = os.path.basename(cmd[1])
        entry = {'cmd': cmd, 'cwd': cwd, 'script': script}
        if script == 'sqanti3_qc.py':
            _write_reference_qc(out_dir, prefix)
        else:
            entry['input_columns'] = list(pd.read_csv(
                cmd[cmd.index('--filter_class') + 1], sep='\t', nrows=0).columns)
            with open(cmd[cmd.index('--corrected_isoforms_fasta') + 1]) as fh:
                entry['fasta_ids'] = [l[1:].split()[0] for l in fh if l.startswith('>')]
            os.makedirs(os.path.join(out_dir, 'logs'), exist_ok=True)
            pd.DataFrame(RESCUE_TABLE, columns=transcript_rescue.TABLE_COLUMNS).to_csv(
                os.path.join(out_dir, prefix + '_rescue_table.tsv'), sep='\t', index=False)
            pd.DataFrame({'isoform': INCLUSION}).to_csv(
                os.path.join(out_dir, prefix + '_rescue_inclusion_list.tsv'), sep='\t',
                index=False)
        calls.append(entry)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(transcript_rescue.subprocess, 'run', run)
    return calls


def _args(tf, pre, out_dir, *extra):
    argv = ['rescue', '-de', 'design.csv', '-q', str(tf), '-d', str(out_dir),
            '--pre_filter_dir', str(pre), '--refGTF', 'ann.gtf', '--refFasta', 'genome.fa',
            *extra]
    return filter_args.build_filter_parser().parse_args(argv)


def _rescue(tmp_path, tf, pre, design, *extra, log=None):
    args = _args(tf, pre, tmp_path / 'rescue', *extra)
    df = filter_io.read_sample_table(design, str(tf), labelled_summary_ok=True)
    messages = []
    transcript_rescue.run_transcript_rescue(args, df, log=messages.append if log is None else log)
    return str(tmp_path / 'rescue' / 'rep1' / 's1'), messages


def _read(path):
    return pd.read_csv(path, sep='\t', dtype=str, na_filter=False).set_index('isoform')


def _cells(row):
    return dict(zip(row['CB'].split(','), map(float, row['FL'].split(','))))


class TestReadsMoveCellByCell:
    def test_each_target_gains_reads_in_the_removed_models_cells(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, _ = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        cls = _read(out + '_rescued_classification.txt')
        # Automatic: ENST1 holds F1's and I1's reads, cell by cell.
        assert _cells(cls.loc['ENST1']) == {'c1': 2, 'c2': 2, 'c4': 1}
        # X1 alone, plus a quarter of X2 (Y1 has 4 of the two targets' 16 reads).
        assert _cells(cls.loc['Y1']) == {'c1': 4, 'c2': 2, 'c4': 2}
        assert _cells(cls.loc['Y2']) == {'c3': 12, 'c2': 3}
        # Tied references hold no reads, so X3 is split evenly and kept fractional.
        assert _cells(cls.loc['ENST7']) == {'c5': 0.5}
        assert _cells(cls.loc['ENST8']) == {'c5': 0.5}

    def test_untouched_counts_keep_their_text(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, _ = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        assert _read(out + '_rescued_classification.txt').loc['Y2', 'FL'] == '12,3'

    def test_removed_rows_and_the_verdict_are_gone(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, _ = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        cls = _read(out + '_rescued_classification.txt')
        assert set(cls.index) == {'Y1', 'Y2', 'ENST9', 'ENST1', 'ENST7', 'ENST8'}
        assert 'filter_result' not in cls.columns
        assert 'psauron_score' not in cls.columns

    def test_reads_are_conserved_except_the_untargeted(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, messages = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        cls = _read(out + '_rescued_classification.txt')
        total = sum(sum(_cells(r).values()) for _, r in cls.iterrows())
        assert total == pytest.approx(17 + 15)  # kept 17, removed 20 of which X4's 5 dropped
        assert any('1 had no target' in m for m in messages)
        assert any('reads moved 15 of 20 removed (5 dropped)' in m for m in messages)

    def test_a_reference_already_kept_under_its_own_name_gains_the_reads_in_place(
            self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, _ = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        cls = pd.read_csv(out + '_rescued_classification.txt', sep='\t', dtype=str)
        assert (cls['isoform'] == 'ENST9').sum() == 1
        assert _cells(cls.set_index('isoform').loc['ENST9']) == {'c3': 1, 'c1': 2}

    def test_cells_detected_follows_the_moved_reads(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, _ = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        cls = _read(out + '_rescued_classification.txt')
        assert cls.loc['Y1', 'cells_detected'] == '3'
        assert cls.loc['ENST1', 'cells_detected'] == '3'
        assert cls.loc['ENST1', 'max_cluster_cells'] == 'NA'


class TestOtherOutputs:
    def test_junction_barcodes_follow_their_model(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, _ = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        junc = pd.read_csv(out + '_junctions.txt', sep='\t', dtype=str, na_filter=False)
        cls = _read(out + '_rescued_classification.txt')
        assert list(junc.columns) == ['isoform', 'canonical', 'CB']
        assert junc.loc[junc['isoform'] == 'Y1', 'CB'].item() == cls.loc['Y1', 'CB']
        enst1 = junc[junc['isoform'] == 'ENST1']
        assert len(enst1) == 2 and (enst1['CB'] == cls.loc['ENST1', 'CB']).all()

    def test_gtf_and_fasta_gain_the_reintroduced_references(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, _ = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        gtf = open(out + '_corrected.gtf').read()
        for iso in ('Y1', 'ENST1', 'ENST7', 'ENST8'):
            assert f'transcript_id "{iso}"' in gtf
        assert gtf.count('transcript_id "ENST9"') == 1
        assert 'ENST_unused' not in gtf and '#!genome-build' not in gtf
        fasta = open(out + '_corrected.fasta').read()
        assert '>ENST1\nTTTT' in fasta and '>Y1\nACGT' in fasta

    def test_sqanti3_records_are_kept_and_its_work_dir_removed(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, _ = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        assert len(pd.read_csv(out + '_rescue_table.tsv', sep='\t')) == len(RESCUE_TABLE)
        assert os.path.isfile(out + '_rescue_inclusion_list.tsv')
        sample_dir = os.path.dirname(out)
        assert not os.path.exists(os.path.join(sample_dir, 'sqanti3_rescue'))
        assert os.path.isdir(os.path.join(sample_dir, 'logs', 'rescue'))
        assert os.path.isfile(out + '_SQANTI_cell_summary.txt.gz')

    def test_readers_find_the_rescued_classification(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        out, _ = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        assert filter_io.classification_path(out) == out + '_rescued_classification.txt'


class TestSqanti3Calls:
    def test_rescue_command_and_reduced_input(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        _rescue(tmp_path, tf, pre, design, '--mode', 'full', '-j', 'my_rules.json')
        call = [c for c in fake_sqanti3 if c['script'] == 'sqanti3_rescue.py'][0]
        cmd = call['cmd']
        assert cmd[cmd.index('-s') + 1] == 'rules'
        assert cmd[cmd.index('--mode') + 1] == 'full'
        assert cmd[cmd.index('-j') + 1] == os.path.abspath('my_rules.json')
        assert '-q' not in cmd and '--requant' not in cmd and '-e' not in cmd
        assert cmd[cmd.index('--corrected_isoforms_fasta') + 1] == \
            os.path.abspath(str(pre / 'rep1' / 's1_corrected.fasta'))
        assert os.path.isabs(cmd[cmd.index('-k') + 1])
        assert call['cwd'] == cmd[cmd.index('-d') + 1]
        assert call['input_columns'] == list(transcript_rescue.RESCUE_INPUT_COLUMNS)

    def test_one_reference_qc_per_short_read_set(self, tmp_path, fake_sqanti3):
        args = _args('tf', 'pre', tmp_path / 'out', '--CAGE_peak', 'cage.bed', '--include_ORF')
        df = pd.DataFrame({'sampleID': ['a', 'b', 'c'], 'file_acc': ['a', 'b', 'c'],
                           'coverage': ['sj1.tab', 'sj1.tab', 'sj2.tab']})
        prefixes = transcript_rescue.reference_classifications(args, df, lambda *a: None)
        assert len(fake_sqanti3) == 2
        assert prefixes['a'] == prefixes['b'] != prefixes['c']
        cmd = fake_sqanti3[0]['cmd']
        assert cmd[cmd.index('--isoforms') + 1] == cmd[cmd.index('--refGTF') + 1]
        assert cmd[cmd.index('--coverage') + 1] == 'sj1.tab'
        assert cmd[cmd.index('--CAGE_peak') + 1] == os.path.abspath('cage.bed')
        assert '--include_ORF' in cmd

    def test_coverage_file_of_file_names_is_expanded(self, tmp_path, fake_sqanti3):
        fofn = tmp_path / 'cov.fofn'
        fofn.write_text('/x/a.SJ.out.tab\n/x/b.SJ.out.tab\n')
        args = _args('tf', 'pre', tmp_path / 'out')
        df = pd.DataFrame({'sampleID': ['a'], 'file_acc': ['a'], 'coverage': [str(fofn)]})
        transcript_rescue.reference_classifications(args, df, lambda *a: None)
        cmd = fake_sqanti3[0]['cmd']
        assert cmd[cmd.index('--coverage') + 1] == '/x/a.SJ.out.tab,/x/b.SJ.out.tab'

    @pytest.mark.parametrize('mode, strategy', [('full', 'rules'), ('automatic', 'ml')])
    def test_ml_output_full_rescue_checks_references_by_rules(self, runs, fake_sqanti3,
                                                              mode, strategy):
        tmp_path, tf, pre, design = runs
        prefix = str(tf / 'rep1' / 's1')
        cls = pd.read_csv(prefix + '_RulesFilter_classification.txt', sep='\t', dtype=str)
        cls['POS_MLprob'] = '0.9'
        cls.to_csv(prefix + '_ML_classification.txt', sep='\t', index=False)
        os.remove(prefix + '_RulesFilter_classification.txt')
        _, messages = _rescue(tmp_path, tf, pre, design, '--mode', mode)
        call = [c for c in fake_sqanti3 if c['script'] == 'sqanti3_rescue.py'][0]
        assert call['cmd'][call['cmd'].index('-s') + 1] == strategy
        assert ('POS_MLprob' in call['input_columns']) == (strategy == 'ml')
        assert any('checks reference targets with the rules' in m
                   for m in messages) == (mode == 'full')

    def test_a_model_written_twice_is_left_out_of_the_alignment(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        fasta = pre / 'rep1' / 's1_corrected.fasta'
        fasta.write_text(fasta.read_text() + '>X1\nAC\n')
        _, messages = _rescue(tmp_path, tf, pre, design, '--mode', 'full')
        ids = [c for c in fake_sqanti3 if c['script'] == 'sqanti3_rescue.py'][0]['fasta_ids']
        assert 'X1' not in ids and 'X2' in ids and len(ids) == len(set(ids))
        assert any('1 model id(s) appear more than once' in m and '(e.g. X1)' in m
                   for m in messages)

    def test_automatic_rescue_passes_the_fasta_untouched(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        fasta = pre / 'rep1' / 's1_corrected.fasta'
        fasta.write_text(fasta.read_text() + '>X1\nAC\n')
        _rescue(tmp_path, tf, pre, design)
        cmd = [c for c in fake_sqanti3 if c['script'] == 'sqanti3_rescue.py'][0]['cmd']
        assert cmd[cmd.index('--corrected_isoforms_fasta') + 1] == os.path.abspath(str(fasta))

    def test_k_replaces_the_reference_qc(self, runs, fake_sqanti3):
        tmp_path, tf, pre, design = runs
        base = _write_reference_qc(str(tmp_path / 'myref'), 'gencode')
        _rescue(tmp_path, tf, pre, design, '-k', base + '_classification.txt')
        assert all(c['script'] != 'sqanti3_qc.py' for c in fake_sqanti3)

    def test_nothing_lost_skips_automatic_rescue(self, tmp_path, fake_sqanti3):
        rows = [r for r in CLASSIFICATION if r['isoform'] not in ('F1', 'I1')]
        tf = tmp_path / 'tf'
        _write_filter_output(tf / 'rep1', rows)
        pre = tmp_path / 'pre' / 'rep1'
        pre.mkdir(parents=True)
        (pre / 's1_corrected.fasta').write_text('>X1\nACGT\n')
        design = tmp_path / 'design.csv'
        design.write_text('sampleID,file_acc\ns1,rep1\n')
        out, messages = _rescue(tmp_path, tf, pre.parent, str(design))
        assert all(c['script'] != 'sqanti3_rescue.py' for c in fake_sqanti3)
        assert any('nothing to rescue' in m for m in messages)
        assert set(_read(out + '_rescued_classification.txt').index) == {'Y1', 'Y2', 'ENST9'}
        assert len(pd.read_csv(out + '_rescue_table.tsv', sep='\t')) == 0


class TestRefusals:
    def _check(self, runs, *extra, mode='isoforms'):
        tmp_path, tf, pre, design = runs
        args = _args(tf, pre, tmp_path / 'rescue', *extra)
        df = filter_io.read_sample_table(design, str(tf), labelled_summary_ok=True)
        transcript_rescue.check_inputs(args, df, mode)

    def test_reads_mode(self, runs):
        with pytest.raises(ValueError, match='isoforms-mode'):
            self._check(runs, mode='reads')

    def test_not_a_transcript_filter_output(self, runs):
        tmp_path, tf, pre, design = runs
        prefix = str(tf / 'rep1' / 's1')
        os.replace(prefix + '_RulesFilter_classification.txt', prefix + '_classification.txt')
        with pytest.raises(ValueError, match="output directory of 'transcripts'"):
            self._check(runs)

    def test_missing_pre_filter_fasta(self, runs):
        os.remove(str(runs[2] / 'rep1' / 's1_corrected.fasta'))
        with pytest.raises(ValueError, match='s1_corrected.fasta'):
            self._check(runs)

    def test_out_dir_not_empty(self, runs):
        tmp_path = runs[0]
        (tmp_path / 'rescue' / 'rep1').mkdir(parents=True)
        (tmp_path / 'rescue' / 'rep1' / 'old_mapping.sam').write_text('')
        with pytest.raises(ValueError, match='not empty'):
            self._check(runs)

    def test_k_without_its_siblings(self, runs, tmp_path):
        lone = tmp_path / 'lone_classification.txt'
        lone.write_text('isoform\n')
        with pytest.raises(ValueError, match='_junctions.txt'):
            self._check(runs, '-k', str(lone))

    def test_ref_gtf_required(self, runs):
        tmp_path, tf, pre, design = runs
        argv = ['rescue', '-de', design, '-q', str(tf), '--pre_filter_dir', str(pre),
                '--refFasta', 'g.fa']
        args = filter_args.build_filter_parser().parse_args(argv)
        df = filter_io.read_sample_table(design, str(tf), labelled_summary_ok=True)
        with pytest.raises(ValueError, match='--refGTF'):
            transcript_rescue.check_inputs(args, df, 'isoforms')


class TestParser:
    def test_reference_qc_defaults_are_the_qc_runs(self):
        rescue = _args('tf', 'pre', 'out')
        qc = qc_args.build_parser().parse_args(
            ['--refFasta', 'g', '--refGTF', 'a', '-de', 'd', '-m', 'isoforms'])
        for option in ('sites', 'window', 'ratio_TSS_metric', 'min_ref_len', 'chunks',
                       'include_ORF', 'CAGE_peak', 'polyA_motif_list'):
            assert getattr(rescue, option) == getattr(qc, option), option

    def test_mode_defaults_to_automatic_and_is_not_the_run_mode(self):
        args = _args('tf', 'pre', 'out')
        assert args.rescue_mode == 'automatic'
        assert not hasattr(args, 'mode')
