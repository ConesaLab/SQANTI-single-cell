import sys
import paths

from qc_args import build_parser
from qc_io import fill_design_table
from sqanti3_qc_runner import run_sqanti3_qc
from classification_enrichment import annotate_with_ujc_hash, annotate_with_cell_metadata
from cell_metrics import calculate_metrics_per_cell
from qc_outputs import write_gene_counts_by_cell, write_ujc_counts_by_cell, write_cv_by_cell
from qc_reports import generate_report, generate_multisample_report
from sc_clustering import run_clustering_analysis

def main():
    ap = build_parser(version_str='1.0.0')
    args = ap.parse_args()

    try:
        df = fill_design_table(args)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    run_sqanti3_qc(args, df)

    if not args.SKIPHASH and args.mode != 'isoforms':
        annotate_with_ujc_hash(args, df)

    annotate_with_cell_metadata(args, df)

    calculate_metrics_per_cell(args, df)

    if getattr(args, 'write_per_cell_outputs', False):
        write_gene_counts_by_cell(args, df)
        write_ujc_counts_by_cell(args, df)
        write_cv_by_cell(args, df)
    else:
        print("[INFO] Skipping per-cell outputs (enable with --write_per_cell_outputs).", file=sys.stdout)

    if args.run_clustering:
        for _, row in df.iterrows():
            run_clustering_analysis(args, row)

    generate_report(args, df)

    if getattr(args, 'multisample_report', False):
        generate_multisample_report(args, df)

if __name__ == "__main__":
    main()


