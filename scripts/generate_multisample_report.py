#!/usr/bin/env python3
"""
generate_multisample_report.py

This script generates SQANTI-sc multisample reports from existing per-sample outputs,
without rerunning the full SQANTI-sc pipeline.

It reads a tab-separated input file (FOFN) listing one sample per line. Columns:
  Col 1 (required): path to *_SQANTI_cell_summary.txt.gz
  Col 2 (optional): path to *_classification.txt (enables length distribution plots)
  Col 3 (optional): color_group — condition used to color PCA dots by hue (e.g. technology)
  Col 4 (optional): shape_group — condition used to encode PCA dot shape (e.g. tissue)
  Col 5 (optional): shade_group — condition used to lighten/darken the color hue (e.g. cell count)

Lines starting with '#' are ignored. If col 2 is absent but col 3+ are present, leave col 2
empty (i.e. use a tab with nothing between it and the next tab).

Usage:
    python generate_multisample_report.py --mode reads --fofn samples.txt [-o ./results]

    python generate_multisample_report.py --mode isoforms --fofn samples.txt \
        --report both -o ./multisample -p my_comparison

FOFN Format (tab-separated text file):
    # col1: cell_summary                              col2: classification (opt)               col3: color  col4: shape  col5: shade
    /results/SRX123/pb_brain_5k_SQANTI_cell_summary.txt.gz  /results/SRX123/pb_brain_5k_classification.txt  PacBio  brain  5k
    /results/SRX456/ont_lung_500_SQANTI_cell_summary.txt.gz /results/SRX456/ont_lung_500_classification.txt ONT     lung   500
    /results/SRX789/pb_heart_SQANTI_cell_summary.txt.gz
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate SQANTI-sc multisample report from existing per-sample outputs "
            "without rerunning sqanti_sc.py."
        )
    )

    parser.add_argument(
        "--fofn", required=True,
        help=(
            "Tab-separated file listing samples. Column 1: path to "
            "*_SQANTI_cell_summary.txt.gz (required). Column 2: path to "
            "*_classification.txt (optional). Lines starting with '#' are ignored."
        ),
    )
    parser.add_argument(
        "--mode", required=True, choices=["reads", "isoforms"],
        help="Mode used to generate the original outputs.",
    )
    parser.add_argument(
        "-o", "--output_dir", default=".",
        help="Directory to save outputs (default: current directory).",
    )
    parser.add_argument(
        "-p", "--prefix", default="SQANTI_sc_multisample_report",
        help="Output prefix for the report files (default: SQANTI_sc_multisample_report).",
    )
    parser.add_argument(
        "--report", default="both", choices=["pdf", "html", "both"],
        help="Report format to generate (default: both).",
    )
    parser.add_argument(
        "--rscript", default="Rscript",
        help="Path to the Rscript executable (default: Rscript).",
    )

    return parser.parse_args()


def parse_fofn(fofn_path):
    """Read the FOFN and return validated cell summary, classification, and group lists.

    Returns a tuple of five parallel lists:
        (cell_summaries, class_files, color_groups, shape_groups, shade_groups)
    Group lists are empty strings for rows where the corresponding column is absent.
    """
    cell_summaries = []
    class_files = []
    color_groups = []
    shape_groups = []
    shade_groups = []

    try:
        with open(fofn_path, "r") as f:
            for line_num, line in enumerate(f, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                cell_summary = fields[0].strip()

                if not cell_summary:
                    print(f"[WARNING] Empty cell summary path at line {line_num}; skipping.", file=sys.stderr)
                    continue

                if not os.path.isfile(cell_summary):
                    print(f"[WARNING] Cell summary not found at line {line_num}: {cell_summary}; skipping.", file=sys.stderr)
                    continue

                cell_summaries.append(os.path.abspath(cell_summary))

                # Col 2: classification file (optional)
                if len(fields) > 1:
                    class_file = fields[1].strip()
                    if class_file and os.path.isfile(class_file):
                        class_files.append(os.path.abspath(class_file))
                    elif class_file:
                        print(
                            f"[WARNING] Classification file not found at line {line_num}: {class_file}; skipping.",
                            file=sys.stderr,
                        )

                # Cols 3–5: group annotation (optional)
                color_groups.append(fields[2].strip() if len(fields) > 2 else "")
                shape_groups.append(fields[3].strip() if len(fields) > 3 else "")
                shade_groups.append(fields[4].strip() if len(fields) > 4 else "")

    except Exception as e:
        print(f"[ERROR] Failed to read FOFN: {e}", file=sys.stderr)
        sys.exit(1)

    return cell_summaries, class_files, color_groups, shape_groups, shade_groups


def main():
    args = parse_args()

    cell_summaries, class_files, color_groups, shape_groups, shade_groups = parse_fofn(args.fofn)

    if len(cell_summaries) < 2:
        print(
            "[ERROR] Need at least 2 valid cell summary files to generate a multisample report.",
            file=sys.stderr,
        )
        sys.exit(1)

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    report_script = (
        Path(__file__).resolve().parents[1] / "src" / "report_assets" / "SQANTI-sc_multisample_report.R"
    )
    if not report_script.is_file():
        print(f"[ERROR] Report script not found: {report_script}", file=sys.stderr)
        sys.exit(1)

    cmd = [
        args.rscript,
        str(report_script),
        "--files", ",".join(cell_summaries),
        "--out_dir", output_dir,
        "--mode", args.mode,
        "--report", args.report,
        "--prefix", args.prefix,
    ]
    if len(class_files) >= 2:
        cmd.extend(["--class_files", ",".join(class_files)])

    # Pass group annotation vectors if any non-empty values were provided
    if any(v for v in color_groups):
        cmd.extend(["--color_group", ",".join(color_groups)])
    if any(v for v in shape_groups):
        cmd.extend(["--shape_group", ",".join(shape_groups)])
    if any(v for v in shade_groups):
        cmd.extend(["--shade_group", ",".join(shade_groups)])

    print("**** Generating multisample SQANTI-sc report...", file=sys.stdout)
    print(f"[INFO] Samples included: {len(cell_summaries)}", file=sys.stdout)
    if len(class_files) >= 2:
        print(f"[INFO] Classification files included for length plots: {len(class_files)}", file=sys.stdout)
    if any(v for v in color_groups):
        print(f"[INFO] Group annotation: color_group provided.", file=sys.stdout)
    if any(v for v in shape_groups):
        print(f"[INFO] Group annotation: shape_group provided.", file=sys.stdout)
    if any(v for v in shade_groups):
        print(f"[INFO] Group annotation: shade_group provided.", file=sys.stdout)

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"[ERROR] Multisample report generation failed: {exc}", file=sys.stderr)
        sys.exit(exc.returncode if exc.returncode is not None else 1)

    print("**** Multisample SQANTI-sc report generated.", file=sys.stdout)
    if args.report in ("pdf", "both"):
        print(f"[INFO] PDF: {os.path.join(output_dir, args.prefix + '.pdf')}", file=sys.stdout)
    if args.report in ("html", "both"):
        print(f"[INFO] HTML: {os.path.join(output_dir, args.prefix + '.html')}", file=sys.stdout)


if __name__ == "__main__":
    main()
