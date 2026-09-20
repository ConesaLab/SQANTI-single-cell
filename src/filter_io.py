import os
import re

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


def read_sample_table(design_path, qc_dir):
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
        prefix = sample_prefix(qc_dir, row['file_acc'], row['sampleID'])
        for suffix in ('_classification.txt', '_SQANTI_cell_summary.txt.gz'):
            if not os.path.isfile(prefix + suffix):
                raise ValueError(
                    f"ERROR: expected SQANTI-sc output not found for sample "
                    f"{row['sampleID']}: {prefix + suffix}"
                )
    return df


def read_cell_summary(path):
    return pd.read_csv(path, **_READ_KW)


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


def subset_junctions(src, dst, mode, keep_isoforms, keep_cells, chunksize=CHUNKSIZE):
    """Subset by surviving isoform, and in isoforms mode rewrite the per-row CB list.

    That column is a copy of the parent transcript's barcode list, so it goes stale the
    moment cells are removed. It must be rewritten rather than dropped: the per-sample
    report groups the junctions by CB. Rewriting is free here because the subset already
    streams the whole file.
    """
    if not os.path.isfile(src):
        return 0, 0, False
    keep_isoforms = set(keep_isoforms)
    keep_cells = set(keep_cells)
    rows_in = rows_out = 0
    rewrote_cb = False
    header = True
    with open(dst, 'w') as out_fh:
        for chunk in pd.read_csv(src, chunksize=chunksize, **_READ_KW):
            rows_in += len(chunk)
            kept = chunk[chunk['isoform'].isin(keep_isoforms)] if 'isoform' in chunk.columns else chunk
            if mode == 'isoforms' and 'CB' in kept.columns and len(kept):
                cb_lists = kept['CB'].astype(str).str.split(',')
                flat = pd.DataFrame({'CB': cb_lists}).explode('CB')
                flat = flat[flat['CB'].isin(keep_cells)]
                kept = kept.copy()
                kept['CB'] = flat.groupby(level=0, sort=False)['CB'].agg(','.join).reindex(kept.index)
                rewrote_cb = True
            rows_out += len(kept)
            kept.to_csv(out_fh, sep='\t', index=False, header=header)
            header = False
    return rows_in, rows_out, rewrote_cb


_GTF_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')


def subset_gtf(src, dst, keep_isoforms):
    """Subset the corrected GTF to the surviving models, by transcript_id.

    Line-based rather than parsed: the corrected GTF reaches many GB, and SQANTI3's own
    GTF reader would pull `src.commands` into the filter, which raises at import time
    when a platform binary is missing. A feature line with no transcript_id cannot be
    attributed to a model, so it is dropped and counted.
    """
    if not os.path.isfile(src):
        return 0, 0, 0
    keep_isoforms = set(keep_isoforms)
    lines_in = lines_out = unattributed = 0
    with open(src) as in_fh, open(dst, 'w') as out_fh:
        for line in in_fh:
            if line.startswith('#'):
                out_fh.write(line)
                continue
            lines_in += 1
            match = _GTF_TRANSCRIPT_ID.search(line)
            if match is None:
                unattributed += 1
                continue
            if match.group(1) in keep_isoforms:
                out_fh.write(line)
                lines_out += 1
    return lines_in, lines_out, unattributed


def subset_fasta(src, dst, keep_isoforms):
    """Subset the corrected FASTA to the surviving models.

    Streamed record by record because the corrected FASTA reaches tens of GB. The id is
    the first whitespace-delimited token of the header, which is how SQANTI3 writes it
    and how the classification names the model.
    """
    if not os.path.isfile(src):
        return 0, 0
    keep_isoforms = set(keep_isoforms)
    records_in = records_out = 0
    keeping = False
    with open(src) as in_fh, open(dst, 'w') as out_fh:
        for line in in_fh:
            if line.startswith('>'):
                records_in += 1
                parts = line[1:].split()
                keeping = bool(parts) and parts[0] in keep_isoforms
                if keeping:
                    records_out += 1
            if keeping:
                out_fh.write(line)
    return records_in, records_out


def write_barcode_list(path, barcodes):
    with open(path, 'w') as fh:
        for bc in barcodes:
            fh.write(f"{bc}\n")
