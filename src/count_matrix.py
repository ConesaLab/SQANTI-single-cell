"""Per-sample count matrices in MEX format, the layout isoforms mode reads as `abundance`."""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

from cell_metrics import cb_fl_coo
from filter_io import (RESULT_COLUMN, SENTINEL_BARCODES, classification_path, drop_artifacts,
                       sample_prefix)

FEATURE_TYPE = 'Gene Expression'


def _read_classification(class_file, columns):
    cls = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False,
                      usecols=lambda c: c in (*columns, RESULT_COLUMN))
    for col in columns:
        if col not in cls.columns:
            cls[col] = np.nan
    return drop_artifacts(cls).reset_index(drop=True)


def isoform_matrix(cls):
    """Isoforms x cells: every kept model in classification order, with the barcodes that
    hold a count. Returns (matrix, features, barcodes)."""
    row, col, data, cells = cb_fl_coo(cls)
    real = ~np.isin(cells, SENTINEL_BARCODES)
    keep = real[col]
    col = (np.cumsum(real) - 1)[col[keep]]
    matrix = sp.coo_matrix((data[keep], (row[keep], col)), shape=(len(cls), int(real.sum())))
    features = pd.DataFrame({'id': cls['isoform'],
                             'name': cls['associated_gene'].fillna(cls['isoform']),
                             'type': FEATURE_TYPE})
    return matrix.tocsr(), features, cells[real]


def gene_matrix(cls):
    """Genes x cells from a reads-mode classification, one read per row. Returns
    (matrix, features, barcodes)."""
    cls = cls[cls['CB'].notna() & ~cls['CB'].isin(SENTINEL_BARCODES) & cls['associated_gene'].notna()]
    gene_codes, genes = pd.factorize(cls['associated_gene'], sort=True)
    cell_codes, cells = pd.factorize(cls['CB'], sort=True)
    matrix = sp.coo_matrix((np.ones(len(cls)), (gene_codes, cell_codes)),
                           shape=(len(genes), len(cells)))
    features = pd.DataFrame({'id': genes, 'name': genes, 'type': FEATURE_TYPE})
    return matrix.tocsr(), features, np.asarray(cells)


def write_mex(out_dir, matrix, features, barcodes):
    """matrix.mtx (features x barcodes, 1-based), features.tsv and barcodes.tsv, uncompressed.
    The field is integer unless a count is fractional."""
    os.makedirs(out_dir, exist_ok=True)
    m = matrix.tocoo()
    nonzero = m.data != 0
    rows, cols, values = m.row[nonzero] + 1, m.col[nonzero] + 1, m.data[nonzero]
    integral = bool(np.all(values == np.round(values)))
    entries = pd.DataFrame({'row': rows, 'col': cols,
                            'value': values.astype(np.int64) if integral else values})
    entries = entries.sort_values(['col', 'row'], kind='stable')
    with open(os.path.join(out_dir, 'matrix.mtx'), 'w', newline='\n') as fh:
        fh.write(f"%%MatrixMarket matrix coordinate {'integer' if integral else 'real'} general\n")
        fh.write(f"{matrix.shape[0]} {matrix.shape[1]} {len(entries)}\n")
        entries.to_csv(fh, sep=' ', header=False, index=False, lineterminator='\n')
    features.to_csv(os.path.join(out_dir, 'features.tsv'), sep='\t', header=False, index=False,
                    lineterminator='\n')
    with open(os.path.join(out_dir, 'barcodes.tsv'), 'w', newline='\n') as fh:
        fh.writelines(f"{b}\n" for b in barcodes)


def write_count_matrices(args, df):
    for _, row in df.iterrows():
        prefix = sample_prefix(args.out_dir, row['file_acc'], row['sampleID'])
        class_file = classification_path(prefix)
        if not os.path.isfile(class_file):
            print(f"[INFO] Classification file not found for {row['file_acc']}. "
                  f"Skipping the count matrix.", file=sys.stdout)
            continue
        if args.mode == 'isoforms':
            cls = _read_classification(class_file, ('isoform', 'associated_gene', 'CB', 'FL'))
            matrix, features, barcodes = isoform_matrix(cls)
        else:
            cls = _read_classification(class_file, ('associated_gene', 'CB'))
            matrix, features, barcodes = gene_matrix(cls)
        out_dir = f"{prefix}_counts"
        write_mex(out_dir, matrix, features, barcodes)
        print(f"**** Count matrix written: {out_dir} ({matrix.shape[0]} "
              f"{'isoforms' if args.mode == 'isoforms' else 'genes'} x {len(barcodes)} cells)",
              file=sys.stdout)
