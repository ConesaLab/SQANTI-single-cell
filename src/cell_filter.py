import json
import os

import numpy as np
import pandas as pd

import filter_io
from filter_io import SENTINEL_BARCODES


RESULT_CELL = 'Cell'
RESULT_ARTIFACT = 'Artifact'
# The one column added to the cell summary, named as SQANTI3 names its own.
RESULT_COLUMN = 'filter_result'

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
    """Read the rules file into a list of rule-sets, SQANTI3's shape.

    Rules within a set are ANDed; the sets are ORed, so a set is an alternative way
    for a cell to be acceptable. SQANTI3 uses that to let one kind of evidence stand
    in for another -- it waives the canonical-junction requirement when short reads
    back the junction up -- and the same applies to cells: "deep enough", or
    "shallower but still detecting plenty of genes".
    """
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
    rulesets = raw['all']
    if not isinstance(rulesets, list) or not rulesets:
        raise ValueError(
            f"ERROR: \"all\" in {path} must be a non-empty list of rule-sets, as in "
            "SQANTI3's filter rules: [{\"depth\": 500}]. A single rule-set is a "
            "one-element list."
        )
    for rs in rulesets:
        if not isinstance(rs, dict) or not rs:
            raise ValueError(
                f"ERROR: every entry of \"all\" in {path} must be a non-empty object "
                f"of column -> threshold; got {rs!r}."
            )
    # The depth column is named per mode. A token keeps one shipped default working
    # in both without silently ignoring a rule that names the other mode's column.
    # Rebuilt in place rather than popped so the JSON's order survives into the
    # reason strings.
    return [{(depth_column(mode) if k == DEPTH_TOKEN else k): v for k, v in rs.items()}
            for rs in rulesets]


def rule_columns(rulesets):
    """Every column named by any rule-set, in first-seen order."""
    seen = []
    for rs in rulesets:
        for column in rs:
            if column not in seen:
                seen.append(column)
    return seen


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
    """Human-readable form of the rule itself, for logs and the params record."""
    if isinstance(rule, bool):
        _reject_rule(column, rule)
    if isinstance(rule, list) and rule and all(
            isinstance(x, (int, float)) and not isinstance(x, bool) for x in rule):
        return f"{column} in [{min(rule)}, {max(rule)}]"
    if isinstance(rule, (int, float)):
        return f"{column} >= {rule}"
    _reject_rule(column, rule)


def _failure_reasons(column, rule, values, shown):
    """SQANTI3's wording for a failed comparison: '{column}: {value} < {threshold}'.

    SQANTI3 stores a range as separate Min_Threshold and Max_Threshold rules, so a
    range violation reports only the bound that was crossed. Matching that keeps one
    vocabulary across both filters' reasons files.
    """
    if isinstance(rule, list):
        low, high = min(rule), max(rule)
        return np.where(values < low,
                        column + ": " + shown + " < " + str(low),
                        column + ": " + shown + " > " + str(high))
    return column + ": " + shown + " < " + str(rule)


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


def validate_rules(rulesets, summary, mode, sampleID, log=print):
    columns = rule_columns(rulesets)
    unknown = [c for c in columns if c not in summary.columns]
    if unknown:
        raise ValueError(
            f"ERROR: rule column(s) not present in the cell summary: "
            f"{', '.join(sorted(unknown))}. Available depth column for mode "
            f"'{mode}' is '{depth_column(mode)}'."
        )
    for column in columns:
        if pd.to_numeric(summary[column], errors='coerce').isna().all():
            log(f"[WARNING] {sampleID}: rule column '{column}' is NA for every cell, so "
                f"this rule judges nothing. The attribute was never measured in the QC "
                f"run that produced this summary.")


def apply_rules(summary, rulesets):
    """Evaluate the rule-sets. Returns (status_frame, reasons Series, passed Series).

    Rule-sets are ORed and the rules inside one are ANDed, as in SQANTI3. Two
    consequences worth stating because they shape what the reasons mean:

    A cell fails only when EVERY rule-set fails, so no single rule "killed" it.
    Reasons are therefore collected from every rule-set, as SQANTI3's get_reasons
    does, and a cell commonly lists several.

    NA means the rule is SKIPPED for that cell, where SQANTI3 fails it. The divergence
    is deliberate: SQANTI3's NA is usually a missing input, uniform across all rows,
    and its alternatives route around it; ours is a per-cell fact -- a cell with no
    multi-exonic reads has no denominator for Non_canonical_prop_in_cell -- so failing
    it would judge the cell on its own composition rather than its quality.
    """
    index = summary.index
    passed_any = pd.Series(False, index=index)
    reason_masks = []

    for rules in rulesets:
        this_set = pd.Series(True, index=index)
        for column, rule in rules.items():
            values = pd.to_numeric(summary[column], errors='coerce')
            evaluable = values.notna()
            ok = _evaluate_rule(values, column, rule)
            this_set &= ok | ~evaluable
            failing = evaluable & ~ok
            if failing.any():
                reason_masks.append((column, rule, failing))
        passed_any |= this_set

    # Reasons are recorded for DISCARDED cells only, as in SQANTI3, which builds them
    # from the artifact rows alone. A kept cell may well have failed rules in an
    # alternative it did not need, and saying so would read as a contradiction next to
    # its own Cell verdict.
    discarded = ~passed_any
    reasons = pd.Series([[] for _ in range(len(summary))], index=index)
    seen = set()
    for column, rule, failing in reason_masks:
        key = (column, repr(rule))
        if key in seen:
            continue
        seen.add(key)
        relevant = failing & discarded
        if not relevant.any():
            continue
        shown = summary[column].astype(str)
        values = pd.to_numeric(summary[column], errors='coerce')
        texts = _failure_reasons(column, rule, values, shown)
        for pos in np.flatnonzero(relevant.to_numpy()):
            reasons.at[index[pos]] = reasons.at[index[pos]] + [texts[pos]]
    return reasons, passed_any


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
    reasons, passed_any = apply_rules(real, rules)

    result = pd.Series(RESULT_CELL, index=real.index)
    reason_text = pd.Series('', index=real.index)

    # The verdict comes from the OR over rule-sets, NOT from whether any reason was
    # collected: with alternatives a cell can fail rules in one set and still be kept
    # because another set accepted it.
    discarded = ~passed_any
    result[discarded] = RESULT_ARTIFACT
    reason_text[discarded] = reasons[discarded].apply('; '.join)

    # filter_result is the ONLY column added to the summary, as SQANTI3 adds only
    # filter_result to its classification. The reasons travel in their own file.
    labelled = real.reset_index(drop=True).copy()
    labelled[RESULT_COLUMN] = result.values
    return labelled, reason_text.reset_index(drop=True)


def write_filter_artifacts(labelled, reasons, run_prefix, params):
    """The three files SQANTI3's rules filter emits, in cell vocabulary, plus the
    parameters. Like SQANTI3's classification the cell summary is labelled in place and
    never subset, so it carries every judged barcode and exactly one added column."""
    cb_col = labelled.columns[0]
    labelled.to_csv(f"{run_prefix}_CellFilter_cell_summary.txt.gz",
                    sep='\t', index=False, compression='gzip')

    passing = labelled.loc[labelled[RESULT_COLUMN] == RESULT_CELL, cb_col]
    filter_io.write_barcode_list(f"{run_prefix}_pass_cells.txt", passing)

    # Artifacts only, as SQANTI3 builds its reasons from the artifact rows alone.
    # SQANTI3 carries structural_category here because its rules are keyed by it; our
    # rules have one key, and no cell-summary column is categorical, so there is
    # nothing equivalent to report.
    discarded = labelled[RESULT_COLUMN] == RESULT_ARTIFACT
    pd.DataFrame({cb_col: labelled.loc[discarded, cb_col],
                  'filter_reason': reasons[discarded.values]}).to_csv(
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

        labelled, reasons = decide_cells(summary, rules, mode, sampleID, log=log)

        params = {
            'SampleID': sampleID,
            'Mode': mode,
            'QCDir': os.path.abspath(args.qc_dir),
            'RulesFile': os.path.abspath(args.rules),
            'BarcodesIn': len(labelled),
            'BarcodesPassing': int((labelled[RESULT_COLUMN] == RESULT_CELL).sum()),
        }

        keep_cells = write_filter_artifacts(labelled, reasons, out_prefix, params)
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
        write_filter_artifacts(labelled, reasons, out_prefix, params)

    return (modes.pop() if modes else None), evidence
