import json
import os
import sys

import numpy as np
import pandas as pd

import filter_io
from cell_filter_methods import AUTO_METHODS
from filter_io import SENTINEL_BARCODES

PASS = 'pass'
FAIL = 'fail'
NOT_EVALUATED = 'not_evaluated'

RESULT_CELL = 'Cell'
RESULT_ARTIFACT = 'Artifact'

# Columns cell_metrics.py writes as a constant 0 when the matching SQANTI3 run flag
# was absent, which is indistinguishable from genuine zero support. A rule on one of
# these fails every cell in the sample because of a command-line flag rather than the
# data, so they are kept out of the defaults and warned about if a user adds them.
CONDITIONAL_COLUMNS = (
    'PolyA_motif_support_prop', 'CAGE_peak_support_prop',
    'TSS_ratio_validated_prop', 'srjunctions_support_prop', 'NMD_prop_in_cell',
)

_JUNCTION_PROPS = (
    'Known_canonical_junctions_prop', 'Known_non_canonical_junctions_prop',
    'Novel_canonical_junctions_prop', 'Novel_non_canonical_junctions_prop',
)


def depth_column(mode):
    return 'Transcripts_in_cell' if mode == 'isoforms' else 'Reads_in_cell'


def denominator_column(column, mode):
    """The column that must be non-trivial for `column` to carry information.

    safe_prop() returns 0 when its denominator is 0, and 0 is the good end of every
    artifact scale -- so an ungated `max` rule reads a missing value as a clean one.
    The denominator is per proportion, not per cell: a cell with 5,000 reads of which
    3 are multi-exonic has a healthy depth and a meaningless non-canonical share.
    """
    if column in _JUNCTION_PROPS:
        return 'total_junctions'
    if column in ('Non_canonical_prop_in_cell', 'Canonical_prop_in_cell'):
        return ('total_transcripts_no_monoexon' if mode == 'isoforms'
                else 'total_reads_no_monoexon')
    # '_prop' anywhere, not just as a suffix: the per-cell artifact loads are named
    # RTS_prop_in_cell / Intrapriming_prop_in_cell, so a suffix test misses them.
    if '_prop' in column or column.endswith('_perc'):
        # Per-category shares divide by their own category count, which is narrower
        # still; depth is the conservative stand-in rather than an exact gate.
        return depth_column(mode)
    return None


DEPTH_TOKEN = '@depth'


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


def describe_rule(column, rule):
    if isinstance(rule, list) and all(isinstance(x, (int, float)) for x in rule) and rule:
        return f"{column} in [{min(rule)}, {max(rule)}]"
    if isinstance(rule, (int, float)):
        return f"{column} >= {rule}"
    if isinstance(rule, str):
        return f"{column} == {rule}"
    if isinstance(rule, list):
        return f"{column} in {{{', '.join(map(str, rule))}}}"
    raise ValueError(f"ERROR: unsupported rule for column '{column}': {rule!r}")


def _evaluate_rule(values, column, rule):
    """Vectorised. Returns a boolean Series (True = passes) and the reason text used
    when it does not. SQANTI3's rules semantics: a list of numbers is a [min, max]
    range, a bare number is a minimum, a string or list of strings is membership."""
    if isinstance(rule, bool):
        raise ValueError(f"ERROR: unsupported rule for column '{column}': {rule!r}")
    if isinstance(rule, list) and rule and all(isinstance(x, (int, float)) for x in rule):
        numeric = pd.to_numeric(values, errors='coerce')
        return (numeric >= min(rule)) & (numeric <= max(rule))
    if isinstance(rule, (int, float)):
        return pd.to_numeric(values, errors='coerce') >= rule
    if isinstance(rule, str):
        return values.astype(str).str.lower() == rule.lower()
    if isinstance(rule, list):
        allowed = {str(x).lower() for x in rule}
        return values.astype(str).str.lower().isin(allowed)
    raise ValueError(f"ERROR: unsupported rule for column '{column}': {rule!r}")


def validate_rules(rules, summary, mode, sampleID, log=print):
    unknown = [c for c in rules if c not in summary.columns]
    if unknown:
        raise ValueError(
            f"ERROR: rule column(s) not present in the cell summary: "
            f"{', '.join(sorted(unknown))}. Available depth column for mode "
            f"'{mode}' is '{depth_column(mode)}'."
        )
    for column in rules:
        if column in CONDITIONAL_COLUMNS:
            values = pd.to_numeric(summary[column], errors='coerce')
            if values.nunique(dropna=True) <= 1:
                log(f"[WARNING] {sampleID}: rule column '{column}' does not vary in this "
                    f"sample. cell_metrics.py writes it as a constant when the matching "
                    f"SQANTI3 flag was not used, so this rule reflects the command line "
                    f"rather than the data.")


def apply_rules(summary, rules, mode, min_depth_for_props):
    """Returns (status_frame, reasons Series). Status is per criterion and has three
    values, because a cell too shallow for a proportion to mean anything is a cell we
    have no evidence about -- which must not read as evidence of cleanliness."""
    status = pd.DataFrame(index=summary.index)
    reasons = pd.Series([[] for _ in range(len(summary))], index=summary.index)
    for column, rule in rules.items():
        values = summary[column]
        passed = _evaluate_rule(values, column, rule)

        denom = denominator_column(column, mode)
        if denom is not None and denom in summary.columns:
            evaluable = pd.to_numeric(summary[denom], errors='coerce').fillna(0) >= min_depth_for_props
        else:
            evaluable = pd.Series(True, index=summary.index)

        col_status = np.where(~evaluable, NOT_EVALUATED, np.where(passed, PASS, FAIL))
        status[f"{column}_status"] = col_status

        failing = col_status == FAIL
        if failing.any():
            text = describe_rule(column, rule)
            shown = values.astype(str)
            for idx in summary.index[failing]:
                reasons.at[idx] = reasons.at[idx] + [f"{text} (got {shown.at[idx]})"]
    return status, reasons


def decide_cells(summary, rules, mode, sampleID, keep=None, drop=None, universe=None,
                 auto_method='none', auto_params=None, min_depth_for_props=100,
                 log=print):
    """Per-barcode verdict. Composition order is fixed and documented: drop, then
    keep (which bypasses everything after it), then the universe restriction, then
    the rules, then the automatic method.
    """
    keep = set(keep or ())
    drop = set(drop or ())
    both = keep & drop
    if both:
        listed = ', '.join(sorted(both)[:5])
        raise ValueError(
            f"ERROR: {len(both)} barcode(s) appear in both --keep_barcodes and "
            f"--drop_barcodes for sample {sampleID}: {listed}"
            f"{' ...' if len(both) > 5 else ''}. Remove them from one of the lists."
        )
    if auto_method not in AUTO_METHODS:
        raise ValueError(
            f"ERROR: unknown automatic method '{auto_method}'. "
            f"Available: {', '.join(sorted(AUTO_METHODS))}."
        )

    cb_col = summary.columns[0]
    cb = summary[cb_col].astype(str)

    sentinel = cb.isin(SENTINEL_BARCODES)
    n_sentinel = int(sentinel.sum())
    if n_sentinel:
        log(f"[WARNING] {sampleID}: {n_sentinel} sentinel barcode(s) "
            f"({', '.join(sorted(set(cb[sentinel])))}) excluded before any statistic. "
            f"These aggregate reads that were never assigned to a cell.")

    real = summary[~sentinel.values]
    validate_rules(rules, real, mode, sampleID, log=log)
    status, reasons = apply_rules(real, rules, mode, min_depth_for_props)

    auto_reasons = AUTO_METHODS[auto_method](real, auto_params or {})

    real_cb = real[cb_col].astype(str)
    source = pd.Series('pass', index=real.index)
    result = pd.Series(RESULT_CELL, index=real.index)
    reason_text = pd.Series('', index=real.index)

    in_drop = real_cb.isin(drop)
    in_keep = real_cb.isin(keep) & ~in_drop
    out_of_universe = (~real_cb.isin(universe)) if universe is not None else pd.Series(
        False, index=real.index)

    failed_rules = reasons.apply(len) > 0
    auto_failed = real_cb.isin(auto_reasons)

    decided = in_drop | in_keep
    result[in_drop] = RESULT_ARTIFACT
    source[in_drop] = 'manual_drop'
    reason_text[in_drop] = 'listed in --drop_barcodes'
    source[in_keep] = 'manual_keep'

    sel = out_of_universe & ~decided
    result[sel] = RESULT_ARTIFACT
    source[sel] = 'not_in_universe'
    reason_text[sel] = 'not listed in --barcode_universe'
    decided = decided | sel

    sel = failed_rules & ~decided
    result[sel] = RESULT_ARTIFACT
    source[sel] = 'rules'
    reason_text[sel] = reasons[sel].apply('; '.join)
    decided = decided | sel

    sel = auto_failed & ~decided
    result[sel] = RESULT_ARTIFACT
    source[sel] = f"auto:{auto_method}"
    reason_text[sel] = real_cb[sel].map(auto_reasons)

    verdict = pd.concat([real.reset_index(drop=True), status.reset_index(drop=True)], axis=1)
    verdict['filter_result'] = result.values
    verdict['filter_source'] = source.values
    verdict['filter_reason'] = reason_text.values

    if n_sentinel:
        sentinel_rows = summary[sentinel.values].copy()
        for col in status.columns:
            sentinel_rows[col] = NOT_EVALUATED
        sentinel_rows['filter_result'] = RESULT_ARTIFACT
        sentinel_rows['filter_source'] = 'sentinel'
        sentinel_rows['filter_reason'] = 'not a cell barcode'
        verdict = pd.concat([verdict, sentinel_rows], ignore_index=True)

    return verdict


def write_verdict_artifacts(verdict, run_prefix, params):
    """The four files SQANTI3's rules filter emits, in cell vocabulary. The verdict
    table keeps every barcode and adds a column; it is not a subset."""
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
    run_name = getattr(args, 'run_name', 'cell_filter')
    rules = load_rules(args.rules, args.mode)
    keep_scoped = filter_io.read_barcode_list(args.keep_barcodes) if args.keep_barcodes else {}
    drop_scoped = filter_io.read_barcode_list(args.drop_barcodes) if args.drop_barcodes else {}
    universe_scoped = (filter_io.read_barcode_list(args.barcode_universe)
                       if args.barcode_universe else None)

    for _, row in df.iterrows():
        sampleID, file_acc = row['sampleID'], row['file_acc']
        prefix = filter_io.sample_prefix(args.out_dir, file_acc, sampleID)
        run_dir = os.path.join(args.out_dir, str(file_acc), run_name)
        os.makedirs(run_dir, exist_ok=True)
        run_prefix = os.path.join(run_dir, str(sampleID))

        summary = filter_io.read_cell_summary(f"{prefix}_SQANTI_cell_summary.txt.gz")
        universe = (filter_io.barcodes_for_sample(universe_scoped, sampleID)
                    if universe_scoped is not None else None)
        verdict = decide_cells(
            summary, rules, args.mode, sampleID,
            keep=filter_io.barcodes_for_sample(keep_scoped, sampleID),
            drop=filter_io.barcodes_for_sample(drop_scoped, sampleID),
            universe=universe,
            auto_method=args.auto_method,
            min_depth_for_props=args.min_depth_for_props,
            log=log,
        )

        params = {
            'SampleID': sampleID,
            'Mode': args.mode,
            'RunName': run_name,
            'RulesFile': os.path.abspath(args.rules),
            'AutoMethod': args.auto_method,
            'MinDepthForProps': args.min_depth_for_props,
            'KeepBarcodesFile': args.keep_barcodes or 'NA',
            'DropBarcodesFile': args.drop_barcodes or 'NA',
            'BarcodeUniverseFile': args.barcode_universe or 'NA',
            'Apply': bool(args.apply),
            'BarcodesIn': len(verdict),
            'BarcodesPassing': int((verdict['filter_result'] == RESULT_CELL).sum()),
        }
        for source, n in verdict.loc[
                verdict['filter_result'] == RESULT_ARTIFACT, 'filter_source'
        ].value_counts().items():
            params[f"ArtifactsBy_{source}"] = int(n)

        keep_cells = write_verdict_artifacts(verdict, run_prefix, params)
        log(f"**** {sampleID}: {params['BarcodesPassing']}/{params['BarcodesIn']} "
            f"barcodes passed the cell filter")

        if not args.apply:
            continue

        rows_in, rows_out, surviving, no_barcode = filter_io.subset_classification(
            f"{prefix}_classification.txt", f"{run_prefix}_classification.txt",
            args.mode, keep_cells)
        j_in, j_out, dropped_cb = filter_io.subset_junctions(
            f"{prefix}_junctions.txt", f"{run_prefix}_junctions.txt",
            args.mode, surviving)
        cells_in, cells_out = filter_io.subset_cell_summary(
            f"{prefix}_SQANTI_cell_summary.txt.gz",
            f"{run_prefix}_SQANTI_cell_summary.txt.gz", keep_cells)

        params.update({
            'ClassificationRowsIn': rows_in,
            'ClassificationRowsOut': rows_out,
            # A cell-level filter removes transcript models: a model observed in no
            # retained cell is not observed. Counted because it surprises people.
            'TranscriptModelsLostAllSupport': rows_in - rows_out,
            'ClassificationRowsWithoutBarcode': no_barcode,
            'JunctionRowsIn': j_in,
            'JunctionRowsOut': j_out,
            'JunctionsCBColumnDropped': dropped_cb,
            'CellSummaryRowsIn': cells_in,
            'CellSummaryRowsOut': cells_out,
        })
        write_verdict_artifacts(verdict, run_prefix, params)

    if not args.apply:
        return None

    filtered_design = os.path.join(args.out_dir, f"{run_name}_design.csv")
    filter_io.write_filtered_design(filtered_design, df, run_name, args.out_dir)
    return filtered_design
