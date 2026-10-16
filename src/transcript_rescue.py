import argparse
import os
import shutil
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

import cell_context
import filter_io
from cell_metrics import calculate_metrics_per_cell
from filter_io import (CHUNKSIZE, LABELLED_CLASSIFICATIONS, READ_KW, RESCUED_CLASSIFICATION,
                       RESULT_ARTIFACT, RESULT_COLUMN, SENTINEL_BARCODES)
from paths import sqantiqcPath
from qc_io import coverage_argument
from isoform_filter import _cells_in, detect_run_mode

FSM = 'full-splice_match'
# What sqanti3_rescue.py reads from the classification. Its copy leaves out the CB/FL
# lists, which SQANTI3 would otherwise load whole.
RESCUE_INPUT_COLUMNS = ('isoform', 'structural_category', 'associated_gene',
                        'associated_transcript', 'exons', RESULT_COLUMN)
RESCUE_TABLE = '_rescue_table.tsv'
INCLUSION_LIST = '_rescue_inclusion_list.tsv'
RESIDUAL = '_rescue_residual.tsv'
TABLE_COLUMNS = ('artifact', 'assigned_transcript', 'rescue_mode', 'origin', 'reintroduced')
REFERENCE_PREFIX = 'reference'
REFERENCE_SIBLINGS = ('_junctions.txt', '_corrected.gtf', '_corrected.fasta')
# Per-model files beside the classification, and how to cut each to a set of models.
MODEL_FILES = (
    ('_corrected.gtf', filter_io.subset_gtf),
    ('_corrected.fasta', filter_io.subset_fasta),
    ('_corrected.faa', filter_io.subset_fasta),
    ('_corrected.cds.gff3', filter_io.subset_gtf),
)


def detect_strategy(prefix):
    """The method of the isoform filter that removed the models."""
    for method, suffix in LABELLED_CLASSIFICATIONS.items():
        if os.path.isfile(prefix + suffix):
            return method
    raise ValueError(
        f"ERROR: {prefix} has no isoform filter output. Rescue reads the output "
        "directory of sqanti_sc_filter.py isoforms, which labels each model kept or "
        "removed."
    )


def rescue_strategy(method, rescue_mode):
    """Full rescue on ML output judges reference targets by the rules: SQANTI3 would apply the
    sample's random forest to them, which was trained on count and cell columns a reference
    transcript does not have. Targets are then chosen by alignment alone."""
    return 'rules' if rescue_mode == 'full' else method


def _value(row, column):
    value = row.get(column)
    if value is None or pd.isna(value) or not str(value).strip():
        return None
    return str(value)


def short_reads(row):
    return _value(row, 'coverage'), _value(row, 'SR_bam')


def check_inputs(args, df, mode):
    """Everything is checked before anything is written, so a refusal leaves no trace."""
    if mode != 'isoforms':
        raise ValueError(
            "ERROR: rescue works on isoforms-mode runs. In reads mode each row is one read, "
            "so a removed read has no count to move."
        )
    if args.refClassif:
        if not args.refClassif.endswith('_classification.txt'):
            raise ValueError("ERROR: -k must be a SQANTI3 QC *_classification.txt.")
        missing = [args.refClassif] if not os.path.isfile(args.refClassif) else []
        base = reference_base(args.refClassif)
        missing += [base + s for s in REFERENCE_SIBLINGS if not os.path.isfile(base + s)]
        if missing:
            raise ValueError(
                "ERROR: -k must be the classification of a SQANTI3 QC run on the reference "
                f"annotation, with its other outputs beside it. Missing: {', '.join(missing)}"
            )
    for _, row in df.iterrows():
        sampleID, file_acc = row['sampleID'], row['file_acc']
        prefix = filter_io.sample_prefix(args.qc_dir, file_acc, sampleID)
        detect_strategy(prefix)
        pre = filter_io.sample_prefix(args.pre_filter_dir, file_acc, sampleID)
        for path in (prefix + '_corrected.gtf', pre + '_corrected.fasta'):
            if not os.path.isfile(path):
                raise ValueError(f"ERROR: expected file not found for sample {sampleID}: "
                                 f"{path}")
        out_dir = os.path.join(args.out_dir, str(file_acc))
        if os.path.isdir(out_dir) and os.listdir(out_dir):
            raise ValueError(
                f"ERROR: {out_dir} is not empty. Write the rescue to a new --out_dir: "
                "SQANTI3 rescue silently reuses an earlier mapping it finds."
            )


def _fasta_id(header):
    parts = header[1:].split()
    return parts[0] if parts else ''


def repeated_fasta_ids(path):
    seen, repeated = set(), set()
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                iso = _fasta_id(line)
                (repeated if iso in seen else seen).add(iso)
    return sorted(repeated)


def write_without(src, dst, isoforms):
    isoforms = set(isoforms)
    keeping = False
    with open(src) as in_fh, open(dst, 'w') as out_fh:
        for line in in_fh:
            if line.startswith('>'):
                keeping = _fasta_id(line) not in isoforms
            if keeping:
                out_fh.write(line)


def alignment_fasta(pre_prefix, work_prefix, sampleID, log):
    """The models' sequences for full rescue to align. SQANTI3 QC can write a model twice
    when its input GTF does not keep a transcript's lines together, and SQANTI3 rescue cannot
    read an alignment whose targets repeat a name, so a repeated model is left out."""
    fasta = pre_prefix + '_corrected.fasta'
    repeated = repeated_fasta_ids(fasta)
    if not repeated:
        return fasta
    log(f"**** {sampleID}: {len(repeated)} model id(s) appear more than once in {fasta} "
        f"(e.g. {repeated[0]}); left out of the alignment")
    unique = work_prefix + '_corrected_unique.fasta'
    write_without(fasta, unique, repeated)
    return unique


def reference_base(class_file):
    suffix = '_classification.txt'
    return class_file[:-len(suffix)] if class_file.endswith(suffix) else class_file


def reference_qc_command(args, coverage, sr_bam, out_dir):
    """SQANTI3 QC on the annotation itself, with the options of the QC run: its rows join
    the sample's classification and are read by the same cell metrics."""
    cmd = [sys.executable, os.path.join(sqantiqcPath, 'sqanti3_qc.py'),
           '--isoforms', os.path.abspath(args.refGTF),
           '--refGTF', os.path.abspath(args.refGTF),
           '--refFasta', os.path.abspath(args.refFasta),
           '--min_ref_len', str(args.min_ref_len), '-w', str(args.window), '-s', args.sites,
           '--ratio_TSS_metric', args.ratio_TSS_metric, '-n', str(args.chunks),
           '-d', os.path.abspath(out_dir), '-o', REFERENCE_PREFIX, '--report', 'skip']
    if args.include_ORF:
        cmd.append('--include_ORF')
    for flag, value in (('--CAGE_peak', args.CAGE_peak),
                        ('--polyA_motif_list', args.polyA_motif_list),
                        ('--polyA_peak', args.polyA_peak), ('--phyloP_bed', args.phyloP_bed)):
        if value:
            cmd += [flag, os.path.abspath(value)]
    if coverage:
        cmd += ['--coverage', coverage_argument(coverage)]
    if sr_bam:
        cmd += ['--SR_bam', sr_bam]
    return cmd


def _run(cmd, what, log_dir, cwd=None):
    print(' '.join(cmd), file=sys.stdout)
    try:
        subprocess.run(cmd, check=True, cwd=cwd)
    except subprocess.CalledProcessError as exc:
        raise ValueError(f"ERROR: {what} failed (exit status {exc.returncode}). Its log is "
                         f"in {log_dir}.")


def reference_classifications(args, df, log):
    """Sample -> prefix of its reference QC run. One run per distinct short-read set, since
    the reference rows must carry the evidence the sample's own models were judged on."""
    if args.refClassif:
        base = reference_base(os.path.abspath(args.refClassif))
        return {row['sampleID']: base for _, row in df.iterrows()}
    groups = {}
    for _, row in df.iterrows():
        groups.setdefault(short_reads(row), []).append(row['sampleID'])
    prefixes = {}
    for n, (key, samples) in enumerate(groups.items(), start=1):
        out_dir = os.path.abspath(os.path.join(args.out_dir, 'reference_qc', str(n)))
        os.makedirs(out_dir, exist_ok=True)
        log(f"**** SQANTI3 QC on the reference annotation for sample(s) {', '.join(samples)}")
        _run(reference_qc_command(args, *key, out_dir),
             'SQANTI3 QC on the reference annotation', os.path.join(out_dir, 'logs'))
        base = os.path.join(out_dir, REFERENCE_PREFIX)
        if not os.path.isfile(base + '_classification.txt'):
            raise ValueError(f"ERROR: SQANTI3 QC wrote no {base}_classification.txt")
        for sampleID in samples:
            prefixes[sampleID] = base
    return prefixes


def lost_references(class_file, mono_exonic=None, chunksize=CHUNKSIZE):
    """References with an FSM and no kept FSM: what automatic rescue puts back. SQANTI3
    counts mono-exonic FSMs too, unless -e none."""
    seen, kept = set(), set()
    for chunk in pd.read_csv(class_file, chunksize=chunksize, **READ_KW,
                             usecols=['structural_category', 'associated_transcript', 'exons',
                                      RESULT_COLUMN]):
        fsm = chunk[chunk['structural_category'] == FSM]
        if mono_exonic == 'none':
            fsm = fsm[pd.to_numeric(fsm['exons'], errors='coerce') > 1]
        seen.update(fsm['associated_transcript'])
        kept.update(fsm.loc[fsm[RESULT_COLUMN] != RESULT_ARTIFACT, 'associated_transcript'])
    return seen - kept


def write_rescue_input(src, dst, strategy, chunksize=CHUNKSIZE):
    columns = list(RESCUE_INPUT_COLUMNS) + (['POS_MLprob'] if strategy == 'ml' else [])
    header = True
    with open(dst, 'w') as out_fh:
        for chunk in pd.read_csv(src, chunksize=chunksize, usecols=columns, **READ_KW):
            chunk[columns].to_csv(out_fh, sep='\t', index=False, header=header)
            header = False


def build_rescue_command(args, strategy, rescue_input, prefix, fasta, reference,
                         work_dir, sampleID):
    cmd = [sys.executable, os.path.join(sqantiqcPath, 'sqanti3_rescue.py'),
           '--filter_class', os.path.abspath(rescue_input),
           '-rg', os.path.abspath(args.refGTF), '-rf', os.path.abspath(args.refFasta),
           '--corrected_isoforms_fasta', os.path.abspath(fasta),
           '--filtered_isoforms_gtf', os.path.abspath(prefix + '_corrected.gtf'),
           '-k', os.path.abspath(reference + '_classification.txt'),
           '--mode', args.rescue_mode, '-s', strategy,
           '-o', str(sampleID), '-d', os.path.abspath(work_dir), '-l', args.log_level]
    if strategy == 'rules' and args.rules:
        cmd += ['-j', os.path.abspath(args.rules)]
    if args.rescue_mono_exonic:
        cmd += ['-e', args.rescue_mono_exonic]
    return cmd


def _cell_counts(cb, fl):
    """A row's cells and counts. A count that is not a number counts as one, as the cell
    metrics count it; an association file without counts lists the cells only."""
    cells = [] if cb in SENTINEL_BARCODES else cb.split(',')
    counts = fl.split(',') if fl else []
    if len(counts) != len(cells):
        counts = ['NA'] * len(cells)
    try:
        values = np.array(counts, dtype=float)
    except ValueError:
        values = pd.to_numeric(pd.Series(counts, dtype=object), errors='coerce').to_numpy(float)
    values = np.where(np.isnan(values), 1.0, values)
    keep = np.array([c not in SENTINEL_BARCODES for c in cells], dtype=bool)
    return [c for c, k in zip(cells, keep) if k], values[keep]


def _fmt(value):
    return f"{value:.15g}"


def _number(text):
    try:
        value = float(text)
    except ValueError:
        return 1.0
    return 1.0 if np.isnan(value) else value


def add_counts(cb, fl, added):
    """CB/FL of a row with `added` ({cell: reads}) summed in. Untouched cells keep their
    original text; new cells go at the end."""
    cells = [] if cb in SENTINEL_BARCODES else cb.split(',')
    counts = fl.split(',') if fl else []
    if len(counts) != len(cells):
        counts = ['NA'] * len(cells)
    position = {cell: i for i, cell in enumerate(cells)}
    for cell in sorted(added):
        if cell in position:
            i = position[cell]
            counts[i] = _fmt(_number(counts[i]) + added[cell])
        else:
            cells.append(cell)
            counts.append(_fmt(added[cell]))
    return ','.join(cells), ','.join(counts)


def residual_name(gene):
    """SQANTI3's requantification row for the counts of removed models without a target."""
    return 'general_residual' if str(gene).startswith('novel') else f'{gene}_residual'


def plan_moves(class_file, table, chunksize=CHUNKSIZE):
    """{target: {cell: counts}} from the rescue table, and {residual row: {cell: counts}}. A
    removed model's counts stay in their cells; tied targets share them in proportion to
    their counts across the sample, evenly when none has any (a reference transcript has
    none), as SQANTI3's requantification splits a count -- but without rounding to whole
    counts. A removed model without a target gives its counts to its gene's residual row."""
    targets_of = table.groupby('artifact', sort=False)['assigned_transcript'].apply(list).to_dict()
    target_ids = set(table['assigned_transcript'])
    mode_of = table.groupby('artifact', sort=False)['rescue_mode'].first().to_dict()
    artifact_cells, sample_total = {}, defaultdict(float)
    residual = defaultdict(lambda: defaultdict(float))
    stats = defaultdict(float)
    for chunk in pd.read_csv(class_file, chunksize=chunksize, **READ_KW,
                             usecols=['isoform', 'associated_gene', 'CB', 'FL', RESULT_COLUMN]):
        for iso, gene, cb, fl, verdict in zip(chunk['isoform'], chunk['associated_gene'],
                                              chunk['CB'], chunk['FL'], chunk[RESULT_COLUMN]):
            removed = verdict == RESULT_ARTIFACT
            if not (removed or iso in target_ids):
                continue
            cells, values = _cell_counts(cb, fl)
            if removed:
                stats['removed_models'] += 1
                stats['removed_counts'] += values.sum()
                if iso in targets_of:
                    artifact_cells[iso] = (cells, values)
                else:
                    row = residual[residual_name(gene)]
                    for cell, value in zip(cells, values):
                        row[cell] += value
            elif iso in target_ids:
                sample_total[iso] = values.sum()
    added = defaultdict(lambda: defaultdict(float))
    for artifact, targets in targets_of.items():
        cells, values = artifact_cells.get(artifact, ([], np.array([])))
        weights = np.array([sample_total.get(t, 0.0) for t in targets])
        shares = (weights / weights.sum() if weights.sum() > 0
                  else np.full(len(targets), 1.0 / len(targets)))
        for target, share in zip(targets, shares):
            for cell, value in zip(cells, values):
                added[target][cell] += value * share
        kind = 'automatic' if mode_of[artifact] == 'automatic' else 'mapping'
        stats[f'{kind}_models'] += 1
        stats[f'{kind}_counts'] += values.sum()
        stats['tied_models'] += len(targets) > 1
    stats['residual_models'] = (stats['removed_models'] - stats['automatic_models']
                                - stats['mapping_models'])
    stats['residual_counts'] = (stats['removed_counts'] - stats['automatic_counts']
                               - stats['mapping_counts'])
    return ({t: dict(cells) for t, cells in added.items()},
            {r: dict(cells) for r, cells in residual.items()}, dict(stats))


def write_residual(path, residual):
    with open(path, 'w') as fh:
        fh.write('isoform\tCB\tFL\n')
        for name in sorted(residual):
            cb, fl = add_counts('', '', residual[name])
            fh.write(f'{name}\t{cb}\t{fl}\n')


def _columns(path):
    """Header as pandas reads it: SQANTI3 repeats a coverage column name when two coverage
    files share one, and pandas suffixes the repeat."""
    return list(pd.read_csv(path, nrows=0, **READ_KW).columns)


def _rows_of(path, isoforms, chunksize=CHUNKSIZE):
    isoforms = set(isoforms)
    parts = [chunk[chunk['isoform'].isin(isoforms)]
             for chunk in pd.read_csv(path, chunksize=chunksize, **READ_KW)]
    return pd.concat(parts) if parts else pd.DataFrame(columns=['isoform'])


def write_rescued_classification(src, dst, added, reference_class, references,
                                 chunksize=CHUNKSIZE):
    """Kept models with their moved counts, then the reintroduced references from the
    reference QC. Returns (references appended, {isoform: new CB} for every changed row).
    A reference already present as a kept model (Bambu names its models after the
    reference) gains the counts in place."""
    columns = [c for c in _columns(src) if c != RESULT_COLUMN]
    kept_ids, new_cb = set(), {}
    hit_ids = set(added)
    with open(dst, 'w') as out_fh:
        out_fh.write('\t'.join(columns) + '\n')
        for chunk in pd.read_csv(src, chunksize=chunksize, **READ_KW):
            kept = chunk[chunk[RESULT_COLUMN] != RESULT_ARTIFACT][columns].copy()
            kept_ids.update(kept['isoform'])
            for i in kept.index[kept['isoform'].isin(hit_ids)]:
                iso = kept.at[i, 'isoform']
                kept.at[i, 'CB'], kept.at[i, 'FL'] = add_counts(kept.at[i, 'CB'],
                                                               kept.at[i, 'FL'], added[iso])
                new_cb[iso] = kept.at[i, 'CB']
            kept.to_csv(out_fh, sep='\t', index=False, header=False)
        appended = sorted(set(references) - kept_ids)
        if appended:
            rows = _rows_of(reference_class, appended)
            missing = set(appended) - set(rows['isoform'])
            if missing:
                raise ValueError(
                    f"ERROR: {len(missing)} reintroduced reference transcript(s) are not in "
                    f"{reference_class} (e.g. {sorted(missing)[0]}): was the reference QC "
                    "run on --refGTF?")
            rows = rows.reindex(columns=columns, fill_value='NA').reset_index(drop=True)
            for i, iso in enumerate(rows['isoform']):
                rows.at[i, 'CB'], rows.at[i, 'FL'] = add_counts('', '', added.get(iso, {}))
                new_cb[iso] = rows.at[i, 'CB']
            rows.to_csv(out_fh, sep='\t', index=False, header=False)
    return appended, new_cb


def write_rescued_junctions(src, reference_junctions, dst, new_cb, appended,
                            chunksize=CHUNKSIZE):
    """In isoforms mode a junction row's CB copies its model's, so it follows the moved
    counts; the reintroduced references bring their junctions from the reference QC."""
    if not os.path.isfile(src):
        return
    columns = _columns(src)
    with open(dst, 'w') as out_fh:
        out_fh.write('\t'.join(columns) + '\n')
        for chunk in pd.read_csv(src, chunksize=chunksize, **READ_KW):
            if 'CB' in chunk.columns:
                hit = chunk['isoform'].isin(set(new_cb))
                chunk.loc[hit, 'CB'] = chunk.loc[hit, 'isoform'].map(new_cb)
            chunk.to_csv(out_fh, sep='\t', index=False, header=False)
        if appended:
            rows = _rows_of(reference_junctions, appended)
            rows = rows.reindex(columns=columns, fill_value='NA')
            if 'CB' in rows.columns:
                rows['CB'] = rows['isoform'].map(new_cb).fillna('')
            rows.to_csv(out_fh, sep='\t', index=False, header=False)


def write_model_files(prefix, reference, out_prefix, appended):
    """The kept models' GTF, FASTA and ORF files, with the reintroduced references' entries
    from the reference QC appended."""
    written = []
    for suffix, subset in MODEL_FILES:
        if not os.path.isfile(prefix + suffix):
            continue
        shutil.copyfile(prefix + suffix, out_prefix + suffix)
        if appended and os.path.isfile(reference + suffix):
            tmp = out_prefix + suffix + '.tmp'
            subset(reference + suffix, tmp, appended)
            with open(tmp) as add_fh, open(out_prefix + suffix, 'a') as out_fh:
                for line in add_fh:
                    if not line.startswith('#'):
                        out_fh.write(line)
            os.remove(tmp)
        written.append(suffix)
    return written


def _empty_records(work_prefix):
    pd.DataFrame(columns=list(TABLE_COLUMNS)).to_csv(work_prefix + RESCUE_TABLE, sep='\t',
                                                     index=False)
    pd.DataFrame(columns=['isoform']).to_csv(work_prefix + INCLUSION_LIST, sep='\t',
                                             index=False)


def rescue_sample(args, row, reference, log):
    sampleID, file_acc = row['sampleID'], row['file_acc']
    prefix = filter_io.sample_prefix(args.qc_dir, file_acc, sampleID)
    pre_prefix = filter_io.sample_prefix(args.pre_filter_dir, file_acc, sampleID)
    in_class = filter_io.classification_path(prefix)
    method = detect_strategy(prefix)
    strategy = rescue_strategy(method, args.rescue_mode)
    out_dir = os.path.join(args.out_dir, str(file_acc))
    out_prefix = os.path.join(out_dir, str(sampleID))
    work_dir = os.path.join(out_dir, 'sqanti3_rescue')
    work_prefix = os.path.join(work_dir, str(sampleID))
    os.makedirs(work_dir, exist_ok=True)

    if args.rescue_mode == 'automatic' and not lost_references(in_class,
                                                               args.rescue_mono_exonic):
        # SQANTI3's automatic rescue raises when it finds nothing to rescue.
        log(f"**** {sampleID}: no reference transcript lost all its FSMs; nothing to rescue")
        _empty_records(work_prefix)
    else:
        if strategy != method:
            log(f"**** {sampleID}: ML filter output: full rescue checks reference targets with "
                "the rules and chooses targets by alignment alone")
        rescue_input = work_prefix + '_rescue_input.tmp'
        write_rescue_input(in_class, rescue_input, strategy)
        fasta = (alignment_fasta(pre_prefix, work_prefix, sampleID, log)
                 if args.rescue_mode == 'full' else pre_prefix + '_corrected.fasta')
        _run(build_rescue_command(args, strategy, rescue_input, prefix, fasta, reference,
                                  work_dir, sampleID),
             f"SQANTI3 rescue for sample {sampleID}", os.path.join(work_dir, 'logs'),
             cwd=work_dir)
        os.remove(rescue_input)

    table = pd.read_csv(work_prefix + RESCUE_TABLE, sep='\t', dtype=str, na_filter=False)
    inclusion = pd.read_csv(work_prefix + INCLUSION_LIST, sep='\t', dtype=str,
                            na_filter=False)['isoform']
    references = set(inclusion) | set(table.loc[table['origin'] == 'reference',
                                                'assigned_transcript'])
    added, residual, stats = plan_moves(in_class, table)
    write_residual(out_prefix + RESIDUAL, residual)
    out_class = out_prefix + RESCUED_CLASSIFICATION
    appended, new_cb = write_rescued_classification(
        in_class, out_class, added, reference + '_classification.txt', references)
    write_rescued_junctions(prefix + '_junctions.txt', reference + '_junctions.txt',
                            out_prefix + '_junctions.txt', new_cb, appended)
    write_model_files(prefix, reference, out_prefix, appended)
    cell_context.refresh(out_class)
    for suffix in (RESCUE_TABLE, INCLUSION_LIST):
        shutil.move(work_prefix + suffix, out_prefix + suffix)
    if os.path.isdir(os.path.join(work_dir, 'logs')):
        os.makedirs(os.path.join(out_dir, 'logs'), exist_ok=True)
        shutil.move(os.path.join(work_dir, 'logs'), os.path.join(out_dir, 'logs', 'rescue'))
    shutil.rmtree(work_dir)

    log(f"**** {sampleID}: {len(appended)} reference transcript(s) reintroduced; of "
        f"{int(stats.get('removed_models', 0))} removed models, "
        f"{int(stats.get('automatic_models', 0))} gave their counts to a reference by "
        f"automatic rescue, {int(stats.get('mapping_models', 0))} to a target by mapping "
        f"({int(stats.get('tied_models', 0))} split between tied targets), and "
        f"{int(stats.get('residual_models', 0))} had no target")
    removed_counts = stats.get('removed_counts', 0.0)
    if removed_counts:
        log(f"**** {sampleID}: counts moved {_fmt(removed_counts - stats['residual_counts'])} "
            f"of {_fmt(removed_counts)} removed; {_fmt(stats['residual_counts'])} without a "
            f"target are in {os.path.basename(out_prefix + RESIDUAL)}")
    return stats


def run_transcript_rescue(args, df, log=print):
    mode, evidence = detect_run_mode(args, df)
    check_inputs(args, df, mode)
    references = reference_classifications(args, df, log)
    for _, row in df.iterrows():
        rescue_sample(args, row, references[row['sampleID']], log)

    # Rescue changes what the kept cells hold, so the summary is recomputed rather than
    # carried over; a cell the filter emptied can come back.
    metrics_args = argparse.Namespace(**{**vars(args), 'mode': mode, **evidence})
    calculate_metrics_per_cell(metrics_args, df)
    for _, row in df.iterrows():
        out_summary = filter_io.sample_prefix(
            args.out_dir, row['file_acc'], row['sampleID']) + '_SQANTI_cell_summary.txt.gz'
        if not os.path.isfile(out_summary):
            raise ValueError(f"ERROR: the cell summary was not recomputed for sample "
                             f"{row['sampleID']}: {out_summary}")
        before = _cells_in(filter_io.cell_summary_path(
            filter_io.sample_prefix(args.qc_dir, row['file_acc'], row['sampleID'])))
        after = _cells_in(out_summary)
        if after - before:
            log(f"**** {row['sampleID']}: {len(after - before)} cell(s) the isoform filter "
                "had emptied hold rescued counts again")
        if before - after:
            log(f"**** {row['sampleID']}: {len(before - after)} cell(s) absent after rescue")
    return mode, evidence
