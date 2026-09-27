import argparse
import os
import subprocess
import sys

import pandas as pd

import filter_io
from cell_filter import detect_measured_evidence, detect_mode
from filter_io import RESULT_ARTIFACT, RESULT_COLUMN
from cell_metrics import calculate_metrics_per_cell
from paths import sqantiqcPath


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
    cmd = [sys.executable, os.path.join(sqantiqcPath, 'sqanti3_filter.py'), 'rules',
           '--sqanti_class', os.path.abspath(class_file),
           '-d', os.path.abspath(out_dir), '-o', str(sampleID),
           '-l', args.log_level, '--skip_report']
    if args.rules:
        cmd += ['-j', os.path.abspath(args.rules)]
    return cmd


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


def run_transcript_filter(args, df, log=print):
    check_not_already_filtered(args, df)
    mode, evidence = detect_run_mode(args, df)

    for _, row in df.iterrows():
        sampleID, file_acc = row['sampleID'], row['file_acc']
        in_prefix = filter_io.sample_prefix(args.qc_dir, file_acc, sampleID)
        in_class = filter_io.classification_path(in_prefix)
        out_dir = os.path.join(args.out_dir, str(file_acc))
        os.makedirs(out_dir, exist_ok=True)
        out_prefix = os.path.join(out_dir, str(sampleID))

        run_sqanti3_filter(args, in_class, out_dir, sampleID)
        passing = filter_io.read_id_list(f"{out_prefix}_pass_isoforms.txt")

        # Replaces SQANTI3's own copy, which it writes back through pandas: 'NA' becomes
        # empty, TRUE becomes True and -1 becomes -1.0, and cell_metrics.py compares
        # those flags as exact strings. Same rows, columns and verdicts, original values.
        rows_in, rows_out = filter_io.write_labelled_classification(
            in_class, f"{out_prefix}_RulesFilter_classification.txt", passing)
        if os.path.isfile(f"{in_prefix}_junctions.txt"):
            filter_io.subset_by_isoform(
                f"{in_prefix}_junctions.txt", f"{out_prefix}_junctions.txt", passing)
        filter_io.subset_gtf(
            f"{in_prefix}_corrected.gtf", f"{out_prefix}_corrected.gtf", passing)
        filter_io.subset_fasta(
            f"{in_prefix}_corrected.fasta", f"{out_prefix}_corrected.fasta", passing)
        extra = filter_io.subset_optional_model_files(in_prefix, out_prefix, passing)
        log(f"**** {sampleID}: {rows_out}/{rows_in} transcript models passed the "
            f"SQANTI3 rules filter")
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
