import os
import sys

import pandas as pd

import filter_io
from cell_filter import run_cell_filter
from filter_args import build_filter_parser


def _run_downstream(args, filtered_design):
    """Re-run the optional stages on the filtered data by calling the existing QC
    functions with the filtered design. Imported lazily so the filter stays usable
    without scanpy, anndata or R installed."""
    df = pd.read_csv(filtered_design, sep=',')

    if args.run_clustering:
        from sc_clustering import run_clustering_analysis
        for _, row in df.iterrows():
            run_clustering_analysis(args, row)

    if args.export_h5ad:
        from sc_export import export_h5ad
        export_h5ad(args, df)

    if args.report != 'skip':
        from qc_reports import generate_report
        generate_report(args, df)

    if args.multisample_report:
        from qc_reports import generate_multisample_report
        generate_multisample_report(args, df)


def main():
    args = build_filter_parser().parse_args()
    if args.subcommand is None:
        build_filter_parser().print_help()
        sys.exit(1)

    try:
        df = filter_io.read_sample_table(args.inDESIGN, args.out_dir)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)

    try:
        filtered_design = run_cell_filter(args, df)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)

    if filtered_design is not None:
        print(f"**** Filtered design written: {filtered_design}")
        _run_downstream(args, filtered_design)


if __name__ == "__main__":
    main()
