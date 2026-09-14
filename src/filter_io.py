import os
import sys

import pandas as pd

CHUNKSIZE = 500000

# classification_enrichment.py fills every missing field with the literal string
# 'NA', so pandas' default NA handling would read those back as NaN and write
# them out as '' -- corrupting rows the filter never touched. na_filter=False
# keeps every field an exact string on the way through.
_READ_KW = dict(sep='\t', dtype=str, na_filter=False, low_memory=False)

SENTINEL_BARCODES = ('', 'NA', 'unassigned', '-', '*')


def sample_prefix(out_dir, file_acc, sampleID):
    return os.path.join(out_dir, str(file_acc), str(sampleID))


def read_sample_table(design_path, out_dir):
    """Design CSV -> validated frame. Unlike qc_io.fill_design_table this resolves
    nothing, writes nothing, and needs only the two columns that locate an output."""
    if not os.path.isfile(design_path):
        raise ValueError(f"ERROR: design file not found: {design_path}")
    df = pd.read_csv(design_path, sep=',')
    missing = [c for c in ('sampleID', 'file_acc') if c not in df.columns]
    if missing:
        raise ValueError(
            f"ERROR: design file {design_path} is missing required column(s): {', '.join(missing)}"
        )
    for _, row in df.iterrows():
        prefix = sample_prefix(out_dir, row['file_acc'], row['sampleID'])
        for suffix in ('_classification.txt', '_SQANTI_cell_summary.txt.gz'):
            if not os.path.isfile(prefix + suffix):
                raise ValueError(
                    f"ERROR: expected SQANTI-sc output not found for sample "
                    f"{row['sampleID']}: {prefix + suffix}"
                )
    return df


def read_barcode_list(path):
    """Barcode list -> {sampleID or None: set(barcodes)}. One column applies to every
    sample; two tab-separated columns scope each entry to one sample, which matters
    because the same 16bp barcode appears in every 10x sample of a design."""
    scoped = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n').rstrip('\r')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) == 1:
                key, barcode = None, parts[0].strip()
            else:
                key, barcode = parts[0].strip(), parts[1].strip()
            if barcode:
                scoped.setdefault(key, set()).add(barcode)
    return scoped


def barcodes_for_sample(scoped, sampleID):
    if not scoped:
        return set()
    return set(scoped.get(None, set())) | set(scoped.get(str(sampleID), set()))


def read_cell_summary(path):
    return pd.read_csv(path, **_READ_KW)


def subset_cell_summary(src, dst, keep_cells):
    """Row subset, values preserved verbatim. Exact because every column of the
    summary is a per-cell aggregation with a per-cell denominator -- there is no
    cross-cell normalisation anywhere in cell_metrics.py, so subsetting and
    recomputing give the same table. Stops being true once step 3 redistributes
    counts into retained cells."""
    df = pd.read_csv(src, **_READ_KW)
    cb_col = df.columns[0]
    out = df[df[cb_col].isin(keep_cells)]
    out.to_csv(dst, sep='\t', index=False, compression='gzip')
    return len(df), len(out)


def _split_cb_fl(chunk):
    cb_lists = chunk['CB'].astype(str).str.split(',')
    if 'FL' not in chunk.columns:
        return cb_lists, None
    fl_lists = chunk['FL'].astype(str).str.split(',')
    bad = cb_lists.str.len() != fl_lists.str.len()
    if bad.any():
        offending = chunk.loc[bad, 'isoform'].head(3).tolist() if 'isoform' in chunk.columns else []
        raise ValueError(
            "ERROR: CB and FL have different element counts on the same row "
            f"(first offenders: {offending}). CB[i] must pair with FL[i]."
        )
    return cb_lists, fl_lists


def _filter_isoform_row_cells(chunk, keep_cells):
    """Edit the comma-separated CB/FL strings in place of dropping rows.

    In isoforms mode one row is one transcript model spanning many cells, and the
    CB/FL pair IS the sparse (transcript x cell) matrix cell_metrics.py parses into
    COO arrays. Dropping a row would delete the model from every cell, so cells are
    removed from the lists and only a row left with no surviving cell is dropped.
    """
    cb_lists, fl_lists = _split_cb_fl(chunk)
    if fl_lists is None:
        flat = pd.DataFrame({'CB': cb_lists}).explode('CB')
    else:
        flat = pd.DataFrame({'CB': cb_lists, 'FL': fl_lists}).explode(['CB', 'FL'])
    flat = flat[flat['CB'].isin(keep_cells)]

    out = chunk.copy()
    grouped = flat.groupby(level=0, sort=False)
    out['CB'] = grouped['CB'].agg(','.join).reindex(chunk.index)
    if fl_lists is not None:
        out['FL'] = grouped['FL'].agg(','.join).reindex(chunk.index)
    return out[out['CB'].notna()]


def subset_classification(src, dst, mode, keep_cells, chunksize=CHUNKSIZE):
    """Write the classification restricted to keep_cells. Returns (rows_in, rows_out,
    surviving isoform ids, rows dropped for having no barcode at all)."""
    keep_cells = set(keep_cells)
    rows_in = rows_out = no_barcode = 0
    surviving = set()
    header = True
    with open(dst, 'w') as out_fh:
        for chunk in pd.read_csv(src, chunksize=chunksize, **_READ_KW):
            rows_in += len(chunk)
            if 'CB' not in chunk.columns:
                raise ValueError(f"ERROR: {src} has no CB column; run the QC pipeline first.")
            if mode == 'isoforms':
                kept = _filter_isoform_row_cells(chunk, keep_cells)
            else:
                blank = chunk['CB'].isin(SENTINEL_BARCODES)
                no_barcode += int(blank.sum())
                kept = chunk[~blank & chunk['CB'].isin(keep_cells)]
            rows_out += len(kept)
            if 'isoform' in kept.columns:
                surviving.update(kept['isoform'].tolist())
            kept.to_csv(out_fh, sep='\t', index=False, header=header)
            header = False
    return rows_in, rows_out, surviving, no_barcode


def subset_junctions(src, dst, mode, keep_isoforms, chunksize=CHUNKSIZE):
    """Subset by surviving isoform. In isoforms mode the per-row CB column added by
    annotate_with_cell_metadata is a stale copy of the transcript's full barcode list
    once cells are removed; it is dropped rather than rewritten because nothing reads
    it in that mode (cell_metrics excludes it from _junc_usecols, write_cv_by_cell
    sources CB from the classification) and rewriting costs a full pass over a file
    that reaches ~70GB."""
    if not os.path.isfile(src):
        return 0, 0, False
    keep_isoforms = set(keep_isoforms)
    rows_in = rows_out = 0
    dropped_cb = False
    header = True
    with open(dst, 'w') as out_fh:
        for chunk in pd.read_csv(src, chunksize=chunksize, **_READ_KW):
            rows_in += len(chunk)
            if mode == 'isoforms' and 'CB' in chunk.columns:
                chunk = chunk.drop(columns=['CB'])
                dropped_cb = True
            kept = chunk[chunk['isoform'].isin(keep_isoforms)] if 'isoform' in chunk.columns else chunk
            rows_out += len(kept)
            kept.to_csv(out_fh, sep='\t', index=False, header=header)
            header = False
    return rows_in, rows_out, dropped_cb


def write_barcode_list(path, barcodes):
    with open(path, 'w') as fh:
        for bc in barcodes:
            fh.write(f"{bc}\n")


def write_filtered_design(path, df, run_name, out_dir):
    """Emit a design CSV pointing at the filtered outputs. file_acc carries the run
    directory so os.path.join(out_dir, file_acc, sampleID) -- rebuilt verbatim in
    every downstream module -- resolves without any of them changing. The input-file
    columns are deliberately dropped so feeding this back to sqanti_sc.py fails loudly
    instead of silently re-running QC."""
    keep_cols = [c for c in ('sampleID', 'color_group', 'shape_group', 'shade_group')
                 if c in df.columns]
    out = df[keep_cols].copy()
    out['file_acc'] = [os.path.join(str(fa), run_name) for fa in df['file_acc']]
    out['source_file_acc'] = df['file_acc'].astype(str).values
    ordered = ['sampleID', 'file_acc', 'source_file_acc'] + [
        c for c in keep_cols if c != 'sampleID']
    out[ordered].to_csv(path, sep=',', index=False)
    return path
