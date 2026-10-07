import argparse
import os
import subprocess
import sys

import pandas as pd

import filter_io
from cell_filter import detect_measured_evidence, detect_mode
from filter_io import (CHUNKSIZE, LABELLED_CLASSIFICATIONS, READ_KW, RESULT_ARTIFACT,
                       RESULT_COLUMN)
from cell_metrics import calculate_metrics_per_cell
from paths import sqantiqcPath

# What SQANTI3's ML filter adds to its classification, in its order, before the verdict.
ML_COLUMNS = ('POS_MLprob', 'NEG_MLprob', 'ML_classifier', 'intra_priming')
# Columns SQANTI-sc adds that describe cells or key a junction chain, not a model's
# quality. A text column that reaches the random forest is used as its rank.
NOT_ML_FEATURES = ('CB', 'UMI', 'jxn_string', 'jxnHash')


def detect_run_mode(args, df):
    """Mode and measured evidence of the run being filtered, read off its cell summaries
    as the cell filter reads them. A cell filter's labelled summary serves equally."""
    modes = set()
    evidence = {}
    for _, row in df.iterrows():
        prefix = filter_io.sample_prefix(args.qc_dir, row['file_acc'], row['sampleID'])
        summary = filter_io.read_cell_summary(filter_io.cell_summary_path(prefix))
        mode = detect_mode(summary, row['sampleID'])
        if modes and mode not in modes:
            raise ValueError(
                f"ERROR: sample {row['sampleID']} is {mode} mode but an earlier sample in "
                f"the design is {modes.pop()} mode; a QC run is one mode throughout."
            )
        modes.add(mode)
        for flag, was_measured in detect_measured_evidence(summary).items():
            evidence[flag] = evidence.get(flag, False) or was_measured
    return modes.pop(), evidence


def build_sqanti3_filter_command(args, class_file, out_dir, sampleID):
    # SQANTI3's own report is skipped, as the QC run skips SQANTI3 QC's: the reports are
    # SQANTI-sc's. -e is never passed: mono-exonic models are judged by the rules like any
    # other, and a rule on the exons column expresses the same thing when wanted.
    cmd = [sys.executable, os.path.join(sqantiqcPath, 'sqanti3_filter.py'), args.method,
           '--sqanti_class', os.path.abspath(class_file),
           '-d', os.path.abspath(out_dir), '-o', str(sampleID),
           '-l', args.log_level, '--skip_report']
    if args.method == 'rules':
        if args.rules:
            cmd += ['-j', os.path.abspath(args.rules)]
        return cmd

    for flag, value in (('-j', args.threshold), ('-t', args.percent_training),
                        ('-z', args.max_class_size), ('-i', args.intrapriming)):
        if value is not None:
            cmd += [flag, str(value)]
    for flag, path in (('-p', args.TP), ('-n', args.TN), ('-r', args.remove_columns)):
        if path:
            cmd += [flag, os.path.abspath(path)]
    if args.force_fsm_in:
        cmd.append('-f')
    if args.intermediate_files:
        cmd.append('--intermediate_files')
    return cmd


def _fl_totals(fl):
    counts = pd.to_numeric(fl.astype(str).str.split(',').explode(), errors='coerce')
    totals = counts.groupby(level=0).sum(min_count=1)
    return totals.map(lambda v: 'NA' if pd.isna(v) else f"{v:.15g}")


def write_ml_input(src, dst, mode, chunksize=CHUNKSIZE):
    """The classification as SQANTI3's ML filter reads it. In isoforms mode FL becomes the
    model's total over its cells, as SQANTI3 adds up per-sample FL columns."""
    header = True
    with open(dst, 'w') as out_fh:
        for chunk in pd.read_csv(src, chunksize=chunksize, **READ_KW):
            chunk = chunk.drop(columns=[c for c in NOT_ML_FEATURES if c in chunk.columns])
            if mode == 'isoforms' and 'FL' in chunk.columns:
                chunk['FL'] = _fl_totals(chunk['FL'])
            chunk.to_csv(out_fh, sep='\t', index=False, header=header)
            header = False


def run_sqanti3_filter(args, class_file, out_dir, sampleID):
    cmd = build_sqanti3_filter_command(args, class_file, out_dir, sampleID)
    print(' '.join(cmd), file=sys.stdout)
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise ValueError(
            f"ERROR: the SQANTI3 filter failed for sample {sampleID} (exit status "
            f"{exc.returncode}). Its log is in {os.path.join(out_dir, 'logs')}."
        )


def _cells_in(summary_path):
    summary = pd.read_csv(summary_path, sep='\t', dtype=str, na_filter=False,
                          usecols=lambda c: c in ('CB', RESULT_COLUMN))
    keep = ~summary['CB'].isin(filter_io.SENTINEL_BARCODES)
    if RESULT_COLUMN in summary.columns:
        keep &= summary[RESULT_COLUMN] != RESULT_ARTIFACT
    return set(summary.loc[keep, 'CB'])


def check_not_already_filtered(args, df):
    for _, row in df.iterrows():
        in_class = filter_io.classification_path(
            filter_io.sample_prefix(args.qc_dir, row['file_acc'], row['sampleID']))
        if RESULT_COLUMN in filter_io.classification_header(in_class):
            raise ValueError(
                f"ERROR: {in_class} has already been through the transcript filter. "
                "Filtering it again could keep models whose junctions, GTF and FASTA "
                "entries were already removed; run on the QC run or a cell filter's "
                "output instead."
            )


def check_out_dir_method(args, df):
    """Readers pick up whichever labelled classification they find, and the report the
    rules filter's reasons file, so one directory holds one method's output."""
    for _, row in df.iterrows():
        prefix = filter_io.sample_prefix(args.out_dir, row['file_acc'], row['sampleID'])
        for method, suffix in LABELLED_CLASSIFICATIONS.items():
            if method != args.method and os.path.isfile(prefix + suffix):
                raise ValueError(
                    f"ERROR: {prefix + suffix} was written by --method {method}. Write "
                    f"--method {args.method} to a different --out_dir."
                )


def check_training_lists(args, df):
    """Model IDs are numbered per sample, so one TP/TN list cannot serve several samples."""
    if args.method == 'ml' and (args.TP or args.TN) and len(df) > 1:
        raise ValueError(
            "ERROR: --TP and --TN list the models of one sample, and model IDs are not "
            "shared between samples. Run each sample with a one-row design, or leave them "
            "out so SQANTI3 builds the lists from each sample."
        )


def run_transcript_filter(args, df, log=print):
    check_not_already_filtered(args, df)
    check_out_dir_method(args, df)
    check_training_lists(args, df)
    mode, evidence = detect_run_mode(args, df)

    for _, row in df.iterrows():
        sampleID, file_acc = row['sampleID'], row['file_acc']
        in_prefix = filter_io.sample_prefix(args.qc_dir, file_acc, sampleID)
        in_class = filter_io.classification_path(in_prefix)
        out_dir = os.path.join(args.out_dir, str(file_acc))
        os.makedirs(out_dir, exist_ok=True)
        out_prefix = os.path.join(out_dir, str(sampleID))

        labelled = out_prefix + LABELLED_CLASSIFICATIONS[args.method]
        if args.method == 'ml':
            ml_input = f"{out_prefix}_ML_input.tmp"
            write_ml_input(in_class, ml_input, mode)
            try:
                run_sqanti3_filter(args, ml_input, out_dir, sampleID)
            finally:
                os.remove(ml_input)
        else:
            run_sqanti3_filter(args, in_class, out_dir, sampleID)
        passing = filter_io.read_id_list(f"{out_prefix}_pass_isoforms.txt")

        # Replaces SQANTI3's own copy: the rules filter writes it back through pandas ('NA'
        # becomes empty, TRUE True, -1 -1.0, and cell_metrics.py compares those flags as
        # exact strings), the ML filter from the reduced input above. Same rows and
        # verdicts, original values.
        if args.method == 'ml':
            sqanti3_copy = f"{labelled}.sqanti3"
            os.replace(labelled, sqanti3_copy)
            rows_in, rows_out = filter_io.write_labelled_classification(
                in_class, labelled, passing, columns_from=sqanti3_copy, columns=ML_COLUMNS)
            os.remove(sqanti3_copy)
        else:
            rows_in, rows_out = filter_io.write_labelled_classification(
                in_class, labelled, passing)
        if os.path.isfile(f"{in_prefix}_junctions.txt"):
            filter_io.subset_by_isoform(
                f"{in_prefix}_junctions.txt", f"{out_prefix}_junctions.txt", passing)
        filter_io.subset_gtf(
            f"{in_prefix}_corrected.gtf", f"{out_prefix}_corrected.gtf", passing)
        filter_io.subset_fasta(
            f"{in_prefix}_corrected.fasta", f"{out_prefix}_corrected.fasta", passing)
        extra = filter_io.subset_optional_model_files(in_prefix, out_prefix, passing)
        log(f"**** {sampleID}: {rows_out}/{rows_in} transcript models passed the "
            f"SQANTI3 {'ML' if args.method == 'ml' else 'rules'} filter")
        if extra:
            log(f"**** {sampleID}: also cut to the passing models: {', '.join(extra)}")

    # Removing models changes every cell's metrics, so the summary is recomputed from
    # the kept models rather than carried over.
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
        lost = len(before - _cells_in(out_summary))
        if lost:
            log(f"**** {row['sampleID']}: {lost} cell(s) kept no transcript model and are "
                f"absent from the filtered cell summary")

    return mode, evidence
