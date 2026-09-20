import json
import os

import numpy as np
import pandas as pd

import filter_io
from filter_io import SENTINEL_BARCODES

PASS = 'pass'
FAIL = 'fail'
NOT_EVALUATED = 'not_evaluated'

RESULT_CELL = 'Cell'
RESULT_ARTIFACT = 'Artifact'

DEPTH_TOKEN = 'depth'

_MODE_DEPTH_COLUMNS = {'isoforms': 'Transcripts_in_cell', 'reads': 'Reads_in_cell'}


def depth_column(mode):
    return _MODE_DEPTH_COLUMNS[mode]


def detect_mode(summary, sampleID):
    """The two depth columns are mutually exclusive and mode-specific, so the cell
    summary states its own mode and the user does not have to."""
    found = [mode for mode, column in _MODE_DEPTH_COLUMNS.items()
             if column in summary.columns]
    if len(found) == 1:
        return found[0]
    if not found:
        raise ValueError(
            f"ERROR: cannot tell which mode sample {sampleID} was run in: its cell "
            f"summary has neither '{_MODE_DEPTH_COLUMNS['reads']}' (reads mode) nor "
            f"'{_MODE_DEPTH_COLUMNS['isoforms']}' (isoforms mode)."
        )
    raise ValueError(
        f"ERROR: cell summary for sample {sampleID} has both "
        f"'{_MODE_DEPTH_COLUMNS['reads']}' and '{_MODE_DEPTH_COLUMNS['isoforms']}'; "
        "cannot tell which mode it was run in."
    )


# The report gates these sections on a CLI boolean, which the filter cannot supply:
# those flags name files for SQANTI3, and the filter never re-runs it. cell_metrics.py
# writes the whole family as NA when the evidence was absent, so the QC output answers
# the same question directly.
_EVIDENCE_COLUMNS = {
    'CAGE_peak': ('CAGE_peak_support_prop',),
    'polyA_motif_list': ('PolyA_motif_support_prop',),
    'include_ORF': ('FSM_coding_prop', 'FSM_non_coding_prop'),
}


def detect_measured_evidence(summary):
    """Which optional SQANTI3 evidence the QC run behind this summary actually measured."""
    measured = {}
    for flag, columns in _EVIDENCE_COLUMNS.items():
        present = [c for c in columns if c in summary.columns]
        measured[flag] = bool(present) and any(
            pd.to_numeric(summary[c], errors='coerce').notna().any() for c in present)
    return measured


def load_rules(path, mode):
    with open(path) as fh:
        raw = json.load(fh)
    if not isinstance(raw, dict):
        raise ValueError(f"ERROR: {path} must contain a JSON object.")
    # Keyed rather than flat so step 2's per-cluster thresholds slot in without a
    # format change.
    if 'all' not in raw:
        raise ValueError(
            f"ERROR: {path} must have a top-level \"all\" key holding the rules.")
    unknown = [k for k in raw if k != 'all']
    if unknown:
        raise ValueError(
            f"ERROR: {path} has unsupported top-level key(s): {', '.join(unknown)}. "
            "Only \"all\" is supported for cell rules."
        )
    # The depth column is named per mode. A token keeps one shipped default working
    # in both without silently ignoring a rule that names the other mode's column.
    # Rebuilt in place rather than popped so the JSON's order survives into the
    # reason strings.
    return {(depth_column(mode) if k == DEPTH_TOKEN else k): v
            for k, v in raw['all'].items()}


def _reject_rule(column, rule):
    # cell_metrics.py coerces every column but CB to numeric, so SQANTI3's string
    # rule forms cannot match anything here -- and would silently discard every cell.
    if isinstance(rule, (str, bool)) or (
            isinstance(rule, list) and any(isinstance(x, (str, bool)) for x in rule)):
        raise ValueError(
            f"ERROR: unsupported rule for column '{column}': {rule!r}. Cell metrics are "
            "numeric, so a rule must be a number (minimum) or a [min, max] pair."
        )
    raise ValueError(f"ERROR: unsupported rule for column '{column}': {rule!r}")


def describe_rule(column, rule):
    if isinstance(rule, bool):
        _reject_rule(column, rule)
    if isinstance(rule, list) and rule and all(
            isinstance(x, (int, float)) and not isinstance(x, bool) for x in rule):
        return f"{column} in [{min(rule)}, {max(rule)}]"
    if isinstance(rule, (int, float)):
        return f"{column} >= {rule}"
    _reject_rule(column, rule)


def _evaluate_rule(values, column, rule):
    """Vectorised. Returns a boolean Series (True = passes). SQANTI3's numeric rules
    semantics: a list of numbers is a [min, max] range, a bare number is a minimum."""
    if isinstance(rule, bool):
        _reject_rule(column, rule)
    if isinstance(rule, list) and rule and all(
            isinstance(x, (int, float)) and not isinstance(x, bool) for x in rule):
        return (values >= min(rule)) & (values <= max(rule))
    if isinstance(rule, (int, float)):
        return values >= rule
    _reject_rule(column, rule)


def validate_rules(rules, summary, mode, sampleID, log=print):
    unknown = [c for c in rules if c not in summary.columns]
    if unknown:
        raise ValueError(
            f"ERROR: rule column(s) not present in the cell summary: "
            f"{', '.join(sorted(unknown))}. Available depth column for mode "
            f"'{mode}' is '{depth_column(mode)}'."
        )
    for column in rules:
        if pd.to_numeric(summary[column], errors='coerce').isna().all():
            log(f"[WARNING] {sampleID}: rule column '{column}' is NA for every cell, so "
                f"this rule judges nothing. The attribute was never measured in the QC "
                f"run that produced this summary.")


def apply_rules(summary, rules):
    """Returns (status_frame, reasons Series). Status is per criterion and has three
    values because a cell whose value is NA is one we have no evidence about: a
    comparison against NA is False, so judging it anyway discards a healthy cell with
    a confidently wrong reason."""
    status = pd.DataFrame(index=summary.index)
    reasons = pd.Series([[] for _ in range(len(summary))], index=summary.index)
    for column, rule in rules.items():
        values = pd.to_numeric(summary[column], errors='coerce')
        evaluable = values.notna()
        passed = _evaluate_rule(values, column, rule)

        col_status = np.where(~evaluable, NOT_EVALUATED, np.where(passed, PASS, FAIL))
        status[f"{column}_status"] = col_status

        failing = col_status == FAIL
        if failing.any():
            text = describe_rule(column, rule)
            shown = summary[column].astype(str)
            for idx in summary.index[failing]:
                reasons.at[idx] = reasons.at[idx] + [f"{text} (got {shown.at[idx]})"]
    return status, reasons


def decide_cells(summary, rules, mode, sampleID, log=print):
    """Per-barcode verdict from the rules alone."""
    cb_col = summary.columns[0]
    cb = summary[cb_col].astype(str)

    # classification_enrichment.py fills a missing barcode with the literal 'NA', so
    # cell_metrics.py aggregates every unbarcoded read into one row that looks like a
    # huge cell. Nothing about it can be judged, so it is dropped rather than reported.
    sentinel = cb.isin(SENTINEL_BARCODES)
    n_sentinel = int(sentinel.sum())
    if n_sentinel:
        log(f"[WARNING] {sampleID}: dropped {n_sentinel} row(s) with no cell barcode "
            f"({', '.join(sorted(set(cb[sentinel])))}).")

    real = summary[~sentinel.values]
    validate_rules(rules, real, mode, sampleID, log=log)
    status, reasons = apply_rules(real, rules)

    result = pd.Series(RESULT_CELL, index=real.index)
    source = pd.Series('pass', index=real.index)
    reason_text = pd.Series('', index=real.index)

    failed_rules = reasons.apply(len) > 0
    result[failed_rules] = RESULT_ARTIFACT
    source[failed_rules] = 'rules'
    reason_text[failed_rules] = reasons[failed_rules].apply('; '.join)

    verdict = pd.concat([real.reset_index(drop=True), status.reset_index(drop=True)], axis=1)
    verdict['filter_result'] = result.values
    verdict['filter_source'] = source.values
    verdict['filter_reason'] = reason_text.values

    return verdict


def write_verdict_artifacts(verdict, run_prefix, params):
    """The three files SQANTI3's rules filter emits, in cell vocabulary, plus the
    parameters. Like SQANTI3's classification the cell summary is labelled in place and
    never subset, so this table carries every judged barcode."""
    cb_col = verdict.columns[0]
    verdict.to_csv(f"{run_prefix}_CellFilter_cell_summary.txt.gz",
                   sep='\t', index=False, compression='gzip')

    passing = verdict.loc[verdict['filter_result'] == RESULT_CELL, cb_col]
    filter_io.write_barcode_list(f"{run_prefix}_pass_cells.txt", passing)

    artifacts = verdict[verdict['filter_result'] == RESULT_ARTIFACT]
    artifacts[[cb_col, 'filter_source', 'filter_reason']].to_csv(
        f"{run_prefix}_cell_filtering_reasons.txt", sep='\t', index=False)

    with open(f"{run_prefix}_cell_filter_params.txt", 'w') as fh:
        for key, value in params.items():
            fh.write(f"{key}\t{value}\n")
    return set(passing)


def run_cell_filter(args, df, log=print):
    modes = set()
    evidence = {flag: False for flag in _EVIDENCE_COLUMNS}
    for _, row in df.iterrows():
        sampleID, file_acc = row['sampleID'], row['file_acc']
        qc_prefix = filter_io.sample_prefix(args.qc_dir, file_acc, sampleID)
        out_dir = os.path.join(args.out_dir, str(file_acc))
        os.makedirs(out_dir, exist_ok=True)
        out_prefix = os.path.join(out_dir, str(sampleID))

        summary = filter_io.read_cell_summary(f"{qc_prefix}_SQANTI_cell_summary.txt.gz")
        mode = detect_mode(summary, sampleID)
        if modes and mode not in modes:
            raise ValueError(
                f"ERROR: sample {sampleID} is {mode} mode but an earlier sample in the "
                f"design is {modes.pop()} mode; a QC run is one mode throughout."
            )
        modes.add(mode)
        for flag, was_measured in detect_measured_evidence(summary).items():
            evidence[flag] = evidence[flag] or was_measured
        rules = load_rules(args.rules, mode)

        verdict = decide_cells(summary, rules, mode, sampleID, log=log)

        params = {
            'SampleID': sampleID,
            'Mode': mode,
            'QCDir': os.path.abspath(args.qc_dir),
            'RulesFile': os.path.abspath(args.rules),
            'BarcodesIn': len(verdict),
            'BarcodesPassing': int((verdict['filter_result'] == RESULT_CELL).sum()),
        }
        for source, n in verdict.loc[
                verdict['filter_result'] == RESULT_ARTIFACT, 'filter_source'
        ].value_counts().items():
            params[f"ArtifactsBy_{source}"] = int(n)

        keep_cells = write_verdict_artifacts(verdict, out_prefix, params)
        log(f"**** {sampleID}: {params['BarcodesPassing']}/{params['BarcodesIn']} "
            f"barcodes passed the cell filter")

        rows_in, rows_out, surviving, no_barcode = filter_io.subset_classification(
            f"{qc_prefix}_classification.txt", f"{out_prefix}_classification.txt",
            mode, keep_cells)
        j_in, j_out, rewrote_cb = filter_io.subset_junctions(
            f"{qc_prefix}_junctions.txt", f"{out_prefix}_junctions.txt",
            mode, surviving, keep_cells)
        g_in, g_out, g_unattributed = filter_io.subset_gtf(
            f"{qc_prefix}_corrected.gtf", f"{out_prefix}_corrected.gtf", surviving)
        f_in, f_out = filter_io.subset_fasta(
            f"{qc_prefix}_corrected.fasta", f"{out_prefix}_corrected.fasta", surviving)

        params.update({
            'ClassificationRowsIn': rows_in,
            'ClassificationRowsOut': rows_out,
            # A cell-level filter removes transcript models: a model observed in no
            # retained cell is not observed. Counted because it surprises people.
            'TranscriptModelsLostAllSupport': rows_in - rows_out,
            'ClassificationRowsWithoutBarcode': no_barcode,
            'JunctionRowsIn': j_in,
            'JunctionRowsOut': j_out,
            'JunctionsCBColumnRewritten': rewrote_cb,
            'GTFLinesIn': g_in,
            'GTFLinesOut': g_out,
            'GTFLinesWithoutTranscriptId': g_unattributed,
            'FastaRecordsIn': f_in,
            'FastaRecordsOut': f_out,
        })
        write_verdict_artifacts(verdict, out_prefix, params)

    return (modes.pop() if modes else None), evidence
