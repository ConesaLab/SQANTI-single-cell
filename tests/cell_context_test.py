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

sqanti_sc_src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../src"))
if sqanti_sc_src_path not in sys.path:
    sys.path.insert(0, sqanti_sc_src_path)

import cell_context
from cell_context import CELLS_DETECTED, MAX_CLUSTER_CELLS, MAX_CLUSTER_PCT


def _read(path):
    return pd.read_csv(path, sep='\t', dtype=str, na_filter=False)


def _clusters(mapping):
    return pd.Series(mapping, dtype=str)


class TestIsoformsMode:
    def test_counts_listed_cells_with_reads(self):
        df = pd.DataFrame({'isoform': ['A', 'B'], 'CB': ['c1,c2,c3', 'NA'],
                           'FL': ['3,0,1', 'NA']})
        out = cell_context.annotate_frame(df, 'isoforms')
        assert out[CELLS_DETECTED].tolist() == ['2', '0']

    def test_cluster_columns_are_not_measured_without_clusters(self):
        df = pd.DataFrame({'isoform': ['A'], 'CB': ['c1'], 'FL': ['1']})
        out = cell_context.annotate_frame(df, 'isoforms')
        assert out[MAX_CLUSTER_CELLS].tolist() == ['NA']
        assert out[MAX_CLUSTER_PCT].tolist() == ['NA']

    def test_count_and_share_are_maximised_separately(self):
        """A rare-cell cluster is not held to a big cluster's absolute count: the most
        cells come from the big cluster, the highest share from the small one."""
        clusters = _clusters({**{f'big{i}': '0' for i in range(10)}, 's1': '1', 's2': '1'})
        df = pd.DataFrame({'isoform': ['A'], 'CB': ['big0,big1,big2,s1'],
                           'FL': ['1,1,1,1']})
        out = cell_context.annotate_frame(df, 'isoforms', clusters)
        assert out.loc[0, MAX_CLUSTER_CELLS] == '3'
        assert out.loc[0, MAX_CLUSTER_PCT] == '50.0000'

    def test_cells_outside_the_umap_count_only_as_detected(self):
        df = pd.DataFrame({'isoform': ['A'], 'CB': ['novel_only'], 'FL': ['4']})
        out = cell_context.annotate_frame(df, 'isoforms', _clusters({'c1': '0'}))
        assert out.loc[0, CELLS_DETECTED] == '1'
        assert out.loc[0, MAX_CLUSTER_CELLS] == '0'
        assert out.loc[0, MAX_CLUSTER_PCT] == '0.0000'

    def test_an_association_without_counts_counts_every_listed_cell(self):
        df = pd.DataFrame({'isoform': ['A'], 'CB': ['c1,c2'], 'FL': ['NA']})
        out = cell_context.annotate_frame(df, 'isoforms')
        assert out.loc[0, CELLS_DETECTED] == '2'


class TestReadsMode:
    def _reads(self):
        return pd.DataFrame({
            'isoform': ['r1', 'r2', 'r3', 'r4'],
            'jxn_string': ['u1', 'u1', 'u1', 'u2'],
            'CB': ['c1', 'c1', 'c2', 'NA'],
        })

    def test_reads_of_one_junction_chain_share_its_cell_count(self):
        out = cell_context.annotate_frame(self._reads(), 'reads')
        assert out[CELLS_DETECTED].tolist() == ['2', '2', '2', '0']

    def test_reads_cluster_columns_use_the_chain_cells(self):
        out = cell_context.annotate_frame(
            self._reads(), 'reads', _clusters({'c1': '0', 'c2': '0', 'c3': '1'}))
        assert out[MAX_CLUSTER_CELLS].tolist() == ['2', '2', '2', '0']
        assert out[MAX_CLUSTER_PCT].tolist()[0] == '100.0000'

    def test_novel_mono_exons_pool_into_one_key_per_strand(self):
        """classification_enrichment keys a mono-exon as chr_strand_monoexon_<transcript>,
        so every novel mono-exon on a strand shares one key. Documented, not fixed."""
        df = pd.DataFrame({'isoform': ['r1', 'r2'],
                           'jxn_string': ['chr1_+_monoexon_novel'] * 2,
                           'CB': ['c1', 'c2']})
        out = cell_context.annotate_frame(df, 'reads')
        assert out[CELLS_DETECTED].tolist() == ['2', '2']

    def test_without_junction_chains_nothing_is_measured(self):
        df = pd.DataFrame({'isoform': ['r1'], 'CB': ['c1']})
        out = cell_context.annotate_frame(df, 'reads')
        assert out.loc[0, list(cell_context.COLUMNS)].tolist() == ['NA', 'NA', 'NA']


class TestRefresh:
    def test_rewrites_in_place_and_keeps_other_fields_verbatim(self, tmp_path):
        path = tmp_path / "s_classification.txt"
        pd.DataFrame({'isoform': ['A'], 'CB': ['c1'], 'FL': ['2'],
                      'min_cov': ['NA'], 'RTS_stage': ['FALSE']}).to_csv(
            path, sep='\t', index=False)
        cell_context.refresh(str(path), 'isoforms', _clusters({'c1': '0'}))
        out = _read(path)
        assert out.loc[0, 'min_cov'] == 'NA'
        assert out.loc[0, 'RTS_stage'] == 'FALSE'
        assert out.loc[0, MAX_CLUSTER_CELLS] == '1'
        assert not os.path.exists(f"{path}.tmp")

    def test_a_refresh_without_clusters_clears_the_cluster_columns(self, tmp_path):
        path = tmp_path / "s_classification.txt"
        pd.DataFrame({'isoform': ['A'], 'CB': ['c1'], 'FL': ['2'],
                      CELLS_DETECTED: ['5'], MAX_CLUSTER_CELLS: ['4'],
                      MAX_CLUSTER_PCT: ['10.0000']}).to_csv(path, sep='\t', index=False)
        cell_context.refresh(str(path), 'isoforms')
        out = _read(path)
        assert out.loc[0, [CELLS_DETECTED, MAX_CLUSTER_CELLS, MAX_CLUSTER_PCT]].tolist() == \
            ['1', 'NA', 'NA']

    def test_artifact_rows_keep_the_values_they_were_judged_on(self, tmp_path):
        path = tmp_path / "s_RulesFilter_classification.txt"
        pd.DataFrame({
            'isoform': ['r1', 'r2'], 'jxn_string': ['u1', 'u1'], 'CB': ['c1', 'c2'],
            CELLS_DETECTED: ['2', '2'], MAX_CLUSTER_CELLS: ['2', '2'],
            MAX_CLUSTER_PCT: ['50.0000', '50.0000'],
            'filter_result': ['Isoform', 'Artifact'],
        }).to_csv(path, sep='\t', index=False)
        cell_context.refresh(str(path), 'reads', _clusters({'c1': '0', 'c2': '0'}))
        out = _read(path).set_index('isoform')
        # The kept read's chain is now seen in one cell: the artifact read's cell is gone.
        assert out.loc['r1', CELLS_DETECTED] == '1'
        assert out.loc['r2', [CELLS_DETECTED, MAX_CLUSTER_PCT]].tolist() == ['2', '50.0000']

    def test_the_verdict_stays_the_last_column(self, tmp_path):
        path = tmp_path / "s_RulesFilter_classification.txt"
        pd.DataFrame({'isoform': ['A'], 'CB': ['c1'], 'FL': ['1'],
                      'filter_result': ['Isoform']}).to_csv(path, sep='\t', index=False)
        cell_context.refresh(str(path), 'isoforms')
        assert _read(path).columns[-1] == 'filter_result'

    @pytest.mark.parametrize('chunksize', [1, 2, 100])
    def test_reads_counts_do_not_depend_on_chunking(self, tmp_path, chunksize):
        path = tmp_path / "s_classification.txt"
        pd.DataFrame({'isoform': ['r1', 'r2', 'r3'], 'jxn_string': ['u1', 'u2', 'u1'],
                      'CB': ['c1', 'c2', 'c3']}).to_csv(path, sep='\t', index=False)
        cell_context.refresh(str(path), 'reads', chunksize=chunksize)
        assert _read(path)[CELLS_DETECTED].tolist() == ['2', '1', '2']

    def test_an_empty_classification_gets_the_columns(self, tmp_path):
        path = tmp_path / "s_classification.txt"
        path.write_text("isoform\tCB\tFL\n")
        cell_context.refresh(str(path), 'isoforms')
        assert list(_read(path).columns) == ['isoform', 'CB', 'FL', *cell_context.COLUMNS]


class TestClusteringHook:
    def test_clustering_writes_the_cluster_columns(self, tmp_path, monkeypatch):
        """Whenever clustering runs, the classification next to it gets the columns
        computed from those clusters."""
        import anndata
        import numpy as np
        import sc_clustering

        prefix = tmp_path / "rep1" / "s1"
        prefix.parent.mkdir()
        pd.DataFrame({'isoform': ['A', 'B'], 'associated_gene': ['g1', 'g2'],
                      'CB': ['c1,c2', 'c3'], 'FL': ['1,1', '2']}).to_csv(
            f"{prefix}_classification.txt", sep='\t', index=False)

        def fake_anndata(args, row):
            adata = anndata.AnnData(np.ones((3, 2)))
            adata.obs_names = ['c1', 'c2', 'c3']
            return adata

        def fake_umap(adata):
            adata.obsm['X_umap'] = np.zeros((adata.n_obs, 2))

        def fake_leiden(adata, **kw):
            adata.obs['leiden'] = pd.Categorical(['0', '0', '1'])

        monkeypatch.setattr(sc_clustering, 'prepare_anndata', fake_anndata)
        for name in ('normalize_total', 'log1p', 'scale', 'neighbors'):
            monkeypatch.setattr(sc_clustering.sc.pp, name, lambda *a, **k: None)
        monkeypatch.setattr(sc_clustering.sc.pp, 'highly_variable_genes',
                            lambda adata, **k: adata.var.__setitem__('highly_variable', True))
        monkeypatch.setattr(sc_clustering.sc.tl, 'pca', lambda *a, **k: None)
        monkeypatch.setattr(sc_clustering.sc.tl, 'umap', fake_umap)
        monkeypatch.setattr(sc_clustering.sc.tl, 'leiden', fake_leiden)

        args = MagicMock(out_dir=str(tmp_path), mode='isoforms', normalization='log1p',
                         n_top_genes=2, n_pc=2, n_neighbors=2, resolution=1.0,
                         clustering_method='leiden')
        sc_clustering.run_clustering_analysis(args, {'file_acc': 'rep1', 'sampleID': 's1'})

        out = _read(f"{prefix}_classification.txt").set_index('isoform')
        assert out.loc['A', MAX_CLUSTER_CELLS] == '2'
        assert out.loc['A', MAX_CLUSTER_PCT] == '100.0000'
        assert out.loc['B', MAX_CLUSTER_CELLS] == '1'
