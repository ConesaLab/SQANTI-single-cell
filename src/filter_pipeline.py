import sys

import filter_io
from cell_filter import run_cell_filter
from filter_args import build_filter_parser


def _run_downstream(args, df):
    """Re-run the reports on the filtered data. The filter's out_dir has the same
    <dir>/<file_acc>/<sampleID>_* layout as a QC run, so the existing functions need
    the original design and no path rewriting. Imported lazily so the filter stays
    usable without R installed."""
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
        df = filter_io.read_sample_table(args.inDESIGN, args.qc_dir)
        mode, evidence = run_cell_filter(args, df)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)

    if mode is not None and (args.report != 'skip' or args.multisample_report):
        args.mode = mode
        for flag, was_measured in evidence.items():
            setattr(args, flag, was_measured)
        _run_downstream(args, df)


if __name__ == "__main__":
    main()
