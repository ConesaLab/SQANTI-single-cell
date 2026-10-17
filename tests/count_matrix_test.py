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

import pandas as pd
import pytest
import scipy.io

sqanti_sc_src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../src"))
if sqanti_sc_src_path not in sys.path:
    sys.path.insert(0, sqanti_sc_src_path)

import filter_args
import qc_args
import rescue_args
from count_matrix import write_count_matrices
from filter_io import LABELLED_CLASSIFICATIONS, RESCUED_CLASSIFICATION


def _write(tmp_path, rows, suffix='_classification.txt'):
    sample_dir = tmp_path / "f1"
    sample_dir.mkdir()
    pd.DataFrame(rows).to_csv(sample_dir / f"s1{suffix}", sep='\t', index=False)
    return str(tmp_path), str(sample_dir / "s1_counts")


def _run(out_dir, mode):
    write_count_matrices(MagicMock(out_dir=out_dir, mode=mode),
                         pd.DataFrame({'file_acc': ['f1'], 'sampleID': ['s1']}))


def _read(mex):
    with open(os.path.join(mex, 'matrix.mtx')) as fh:
        header = fh.readline()
    matrix = scipy.io.mmread(os.path.join(mex, 'matrix.mtx')).tocsr()
    features = pd.read_csv(os.path.join(mex, 'features.tsv'), sep='\t', header=None)
    barcodes = pd.read_csv(os.path.join(mex, 'barcodes.tsv'), header=None)[0].tolist()
    return header, matrix, features, barcodes


def _iso(isoform, cb, fl, gene='G1', **extra):
    return {'isoform': isoform, 'associated_gene': gene, 'CB': cb, 'FL': fl, **extra}


class TestIsoformsMode:
    def test_reproduces_the_counts_in_the_input_layout(self, tmp_path):
        out_dir, mex = _write(tmp_path, [_iso('iso1', 'B,A', '5,2'), _iso('iso2', 'A', '3', gene='G2')])
        _run(out_dir, 'isoforms')
        header, m, features, barcodes = _read(mex)
        assert header.strip() == '%%MatrixMarket matrix coordinate integer general'
        assert barcodes == ['A', 'B']
        assert features.values.tolist() == [['iso1', 'G1', 'Gene Expression'],
                                            ['iso2', 'G2', 'Gene Expression']]
        assert m.toarray().tolist() == [[2, 5], [3, 0]]

    def test_fractional_counts_make_a_real_matrix(self, tmp_path):
        out_dir, mex = _write(tmp_path, [_iso('iso1', 'A', '0.5'), _iso('iso2', 'A', '1.25')])
        _run(out_dir, 'isoforms')
        header, m, _, _ = _read(mex)
        assert header.strip() == '%%MatrixMarket matrix coordinate real general'
        assert m.toarray().ravel().tolist() == [0.5, 1.25]

    def test_models_without_cells_keep_their_row(self, tmp_path):
        out_dir, mex = _write(tmp_path, [_iso('iso1', 'A', '4'), _iso('iso2', 'NA', 'NA')])
        _run(out_dir, 'isoforms')
        _, m, features, barcodes = _read(mex)
        assert features[0].tolist() == ['iso1', 'iso2'] and barcodes == ['A']
        assert m.toarray().tolist() == [[4], [0]]

    @pytest.mark.parametrize('suffix', list(LABELLED_CLASSIFICATIONS.values()))
    def test_isoform_filter_artifacts_are_left_out(self, tmp_path, suffix):
        rows = [_iso('iso1', 'A', '4', filter_result='Isoform'),
                _iso('iso2', 'A', '9', filter_result='Artifact')]
        out_dir, mex = _write(tmp_path, rows, suffix)
        _run(out_dir, 'isoforms')
        _, m, features, _ = _read(mex)
        assert features[0].tolist() == ['iso1'] and m.toarray().tolist() == [[4]]

    def test_reads_the_rescued_classification(self, tmp_path):
        out_dir, mex = _write(tmp_path, [_iso('iso1', 'A', '1.5'), _iso('ENST1', 'A', '2.5')],
                              RESCUED_CLASSIFICATION)
        _run(out_dir, 'isoforms')
        _, m, features, _ = _read(mex)
        assert features[0].tolist() == ['iso1', 'ENST1'] and m.toarray().ravel().tolist() == [1.5, 2.5]


class TestReadsMode:
    def test_counts_reads_per_gene_and_cell(self, tmp_path):
        rows = [{'isoform': f'r{i}', 'associated_gene': g, 'CB': cb}
                for i, (g, cb) in enumerate([('G2', 'A'), ('G1', 'A'), ('G1', 'A'), ('G1', 'B'),
                                             ('G1', 'NA')])]
        out_dir, mex = _write(tmp_path, rows)
        _run(out_dir, 'reads')
        header, m, features, barcodes = _read(mex)
        assert 'integer' in header
        assert features[0].tolist() == ['G1', 'G2'] and barcodes == ['A', 'B']
        assert m.toarray().tolist() == [[2, 1], [1, 0]]


@pytest.mark.parametrize('parser', [
    lambda: qc_args.build_parser(),
    lambda: filter_args.build_filter_parser()._subparsers._group_actions[0].choices['cells'],
    lambda: filter_args.build_filter_parser()._subparsers._group_actions[0].choices['isoforms'],
    lambda: rescue_args.build_rescue_parser(),
])
def test_every_entry_point_takes_the_output_options(parser):
    dests = {a.dest for a in parser()._actions}
    assert {'count_matrix', 'export_h5ad', 'sqanti_reads_tables'} <= dests
    assert 'write_per_cell_outputs' not in dests
