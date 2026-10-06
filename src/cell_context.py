import os

import pandas as pd

import filter_io
from filter_io import (CHUNKSIZE, READ_KW, RESULT_ARTIFACT, RESULT_COLUMN,
                       SENTINEL_BARCODES)

CELLS_DETECTED = 'cells_detected'
MAX_CLUSTER_CELLS = 'max_cluster_cells'
MAX_CLUSTER_PCT = 'max_cluster_pct'
COLUMNS = (CELLS_DETECTED, MAX_CLUSTER_CELLS, MAX_CLUSTER_PCT)
NOT_MEASURED = 'NA'


def read_clusters(umap_csv):
    """Barcode -> cluster label, from sc_clustering's umap_results.csv."""
    clusters = pd.read_csv(umap_csv, dtype=str, usecols=['Barcode', 'Cluster'])
    return clusters.set_index('Barcode')['Cluster']


def _isoform_pairs(df):
    """(row, cell) for every cell a model is detected in: listed with FL > 0."""
    cb_lists = df['CB'].astype(str).str.split(',')
    fl_lists = (df['FL'] if 'FL' in df.columns else pd.Series('NA', index=df.index))
    fl_lists = fl_lists.astype(str).str.split(',')
    # An association file without counts leaves FL as a single NA: every listed cell counts.
    unpaired = cb_lists.str.len() != fl_lists.str.len()
    fl_lists[unpaired] = cb_lists[unpaired].map(lambda cells: ['NA'] * len(cells))
    flat = pd.DataFrame({'CB': cb_lists, 'FL': fl_lists}).explode(['CB', 'FL'])
    # A non-numeric FL counts as one read, as cell_metrics.py counts it.
    flat = flat[pd.to_numeric(flat['FL'], errors='coerce').fillna(1) > 0]
    flat = flat[~flat['CB'].isin(SENTINEL_BARCODES)]
    return pd.DataFrame({'key': flat.index, 'CB': flat['CB'].to_numpy()})


def _read_pairs(df):
    """In reads mode a row is one read, so the model is its junction chain."""
    pairs = df.loc[~df['CB'].isin(SENTINEL_BARCODES), ['jxn_string', 'CB']]
    return pairs.rename(columns={'jxn_string': 'key'})


def _stats(pairs, clusters):
    pairs = pairs.drop_duplicates()
    stats = pairs.groupby('key', sort=False).size().rename(CELLS_DETECTED).to_frame()
    if clusters is None:
        return stats
    sizes = clusters.value_counts()
    clustered = pairs.assign(cluster=pairs['CB'].map(clusters)).dropna(subset=['cluster'])
    per_cluster = clustered.groupby(['key', 'cluster'], sort=False).size().rename('n')
    per_cluster = per_cluster.reset_index()
    per_cluster['pct'] = 100 * per_cluster['n'] / per_cluster['cluster'].map(sizes)
    best = per_cluster.groupby('key', sort=False).agg(
        **{MAX_CLUSTER_CELLS: ('n', 'max'), MAX_CLUSTER_PCT: ('pct', 'max')})
    return stats.join(best)


def _format(stats, keys, clustered):
    values = stats.reindex(keys)
    out = pd.DataFrame(index=keys.index)
    out[CELLS_DETECTED] = values[CELLS_DETECTED].fillna(0).astype(int).astype(str).to_numpy()
    if clustered:
        out[MAX_CLUSTER_CELLS] = (values[MAX_CLUSTER_CELLS].fillna(0)
                                  .astype(int).astype(str).to_numpy())
        out[MAX_CLUSTER_PCT] = values[MAX_CLUSTER_PCT].fillna(0).map('{:.4f}'.format).to_numpy()
    else:
        out[MAX_CLUSTER_CELLS] = NOT_MEASURED
        out[MAX_CLUSTER_PCT] = NOT_MEASURED
    return out


def _values(df, mode, clusters, read_stats=None):
    """The three columns for the rows of df, as strings, aligned to its index."""
    if 'CB' not in df.columns or (mode == 'reads' and 'jxn_string' not in df.columns):
        return pd.DataFrame(NOT_MEASURED, index=df.index, columns=list(COLUMNS))
    if mode == 'reads':
        stats = read_stats if read_stats is not None else _stats(_read_pairs(df), clusters)
        return _format(stats, df['jxn_string'], clusters is not None)
    keys = pd.Series(df.index, index=df.index)
    return _format(_stats(_isoform_pairs(df), clusters), keys, clusters is not None)


def _place(df):
    """New columns go before the verdict, which stays last as SQANTI3 writes it."""
    if RESULT_COLUMN in df.columns and df.columns[-1] != RESULT_COLUMN:
        df = df[[c for c in df.columns if c != RESULT_COLUMN] + [RESULT_COLUMN]]
    return df


def annotate_frame(df, mode, clusters=None):
    """Add or replace the three columns on an in-memory classification."""
    values = _values(df, mode, clusters)
    for column in COLUMNS:
        df[column] = values[column]
    return _place(df)


def _kept(df):
    if RESULT_COLUMN not in df.columns:
        return pd.Series(True, index=df.index)
    return df[RESULT_COLUMN] != RESULT_ARTIFACT


def _reads_stats(class_file, header, clusters):
    usecols = [c for c in ('jxn_string', 'CB', RESULT_COLUMN) if c in header]
    narrow = pd.read_csv(class_file, sep='\t', usecols=usecols, dtype='category',
                         na_filter=False)
    pairs = _read_pairs(narrow[_kept(narrow)]).drop_duplicates()
    return _stats(pairs.astype(str), clusters)


def refresh(class_file, mode, clusters=None, chunksize=CHUNKSIZE):
    """Rewrite the classification with the three columns recomputed for the cells it holds
    now. Artifact rows keep the values they were judged on; every reader skips them."""
    header = filter_io.classification_header(class_file)
    read_stats = None
    if mode == 'reads' and {'jxn_string', 'CB'} <= set(header):
        read_stats = _reads_stats(class_file, header, clusters)

    tmp = f"{class_file}.tmp"
    wrote_header = False
    with open(tmp, 'w') as out_fh:
        for chunk in pd.read_csv(class_file, chunksize=chunksize, **READ_KW):
            kept = _kept(chunk)
            values = _values(chunk[kept], mode, clusters, read_stats)
            for column in COLUMNS:
                if column not in chunk.columns:
                    chunk[column] = NOT_MEASURED
                chunk.loc[kept, column] = values[column]
            _place(chunk).to_csv(out_fh, sep='\t', index=False, header=not wrote_header)
            wrote_header = True
        if not wrote_header:
            columns = [c for c in header if c != RESULT_COLUMN]
            columns += [c for c in COLUMNS if c not in columns]
            columns += [RESULT_COLUMN] if RESULT_COLUMN in header else []
            out_fh.write('\t'.join(columns) + '\n')
    os.replace(tmp, class_file)
