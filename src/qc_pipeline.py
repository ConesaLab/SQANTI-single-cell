import sys
import paths

from qc_args import build_parser, warn_ignored_options
from qc_io import fill_design_table
from sqanti3_qc_runner import run_sqanti3_qc
from classification_enrichment import annotate_with_ujc_hash, annotate_with_cell_metadata, annotate_with_sample_support
from cell_metrics import calculate_metrics_per_cell
from qc_outputs import write_outputs
from qc_reports import generate_report, generate_multisample_report
from sc_clustering import run_clustering_analysis

def main():
    ap = build_parser()
    args = ap.parse_args()

    warn_ignored_options(args)

    try:
        df = fill_design_table(args)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    run_sqanti3_qc(args, df)

    if not args.SKIPHASH and args.mode != 'isoforms':
        annotate_with_ujc_hash(args, df)

    annotate_with_cell_metadata(args, df)

    annotate_with_sample_support(args, df)

    calculate_metrics_per_cell(args, df)

    if args.run_clustering:
        for _, row in df.iterrows():
            run_clustering_analysis(args, row)

    write_outputs(args, df)

    generate_report(args, df)

    if getattr(args, 'multisample_report', False):
        generate_multisample_report(args, df)

if __name__ == "__main__":
    main()


