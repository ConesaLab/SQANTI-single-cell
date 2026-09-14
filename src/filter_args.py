import argparse
import os

from qc_args import add_clustering_args

DEFAULT_CELL_RULES = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'filter_assets', 'cell_filter_default.json')


def add_cell_filter_args(group):
    """Also the hook for a future --run_cell_filter stage in qc_pipeline, which is
    why the arguments live in a factory rather than inline in the parser."""
    group.add_argument('-j', '--rules', default=DEFAULT_CELL_RULES,
                       help='JSON of per-metric cell rules. Default: the bundled '
                            'cell_filter_default.json.')
    group.add_argument('--auto_method', default='none',
                       help='Automatic data-driven method. Default: none.')
    group.add_argument('--min_depth_for_props', type=int, default=100,
                       help='A proportion is only evaluated when its own denominator '
                            'reaches this many reads/transcripts; below it the criterion '
                            'is recorded as not_evaluated rather than passed. Default: 100.')
    group.add_argument('--keep_barcodes',
                       help='Barcodes to retain unconditionally, bypassing every rule.')
    group.add_argument('--drop_barcodes',
                       help='Barcodes to discard unconditionally.')
    group.add_argument('--barcode_universe',
                       help='Restrict to these barcodes; listed barcodes still face the rules.')
    group.add_argument('--apply', action='store_true', default=False,
                       help='Write the filtered classification, junctions and cell summary '
                            'as well as the verdict. Default: False (verdict only).')
    return group


def build_filter_parser(version_str: str = '1.2.0'):
    ap = argparse.ArgumentParser(
        description="SQANTI-sc filter: quality filtering of cell barcodes and transcript "
                    "models from an existing SQANTI-sc QC run"
    )
    ap.add_argument('-v', '--version', action='version', version='sqanti-sc-filter ' + version_str)

    common = argparse.ArgumentParser(add_help=False)
    apr = common.add_argument_group("Required arguments")
    apr.add_argument('-de', '--design', dest="inDESIGN", required=True,
                     help='Design file with sampleID and file_acc, as used for the QC run.')
    apr.add_argument('-m', '--mode', choices=["isoforms", "reads"], required=True,
                     help='Input data type.')
    apc = common.add_argument_group("Common options")
    apc.add_argument('-d', '--out_dir', default=".",
                     help='Output directory of the QC run. Default: current dir.')
    apc.add_argument('--run_name', default='cell_filter',
                     help='Subdirectory the filter writes into, so threshold trials can '
                          'coexist. Default: cell_filter.')
    apc.add_argument('-l', '--log_level', default='INFO',
                     choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'],
                     help='Set logging level. Default: INFO.')

    sub = ap.add_subparsers(dest='subcommand')
    cells = sub.add_parser('cells', parents=[common],
                           help='Filter cell barcodes on per-cell QC metrics.')
    add_cell_filter_args(cells.add_argument_group("Cell filter options"))

    apd = cells.add_argument_group("Downstream re-run options (require --apply)")
    apd.add_argument('--report', choices=["pdf", "html", "both", "skip"], default="skip",
                     help="Re-render the per-sample report on filtered data. Default: skip.")
    apd.add_argument('--refGTF', help='Reference annotation (GTF), for the report.')
    apd.add_argument('--CAGE_peak', help="FANTOM5 CAGE Peak (BED), for the report.")
    apd.add_argument('--polyA_motif_list', help="Ranked list of polyA motifs, for the report.")
    apd.add_argument('--include_ORF', action="store_true", default=False,
                     help="Report ORF sections. Default: False.")
    apd.add_argument('--ignore_cell_summary', action="store_true", default=False,
                     help="Don't save cell summary table in report. Default: False.")
    apd.add_argument('--export_h5ad', action='store_true', default=False,
                     help='Export an AnnData .h5ad per sample from the filtered data.')
    apd.add_argument('--multisample_report', action='store_true', default=False,
                     help='Generate a multisample report from the filtered cell summaries.')
    apd.add_argument('--multisample_report_prefix', default='SQANTI_sc_multisample_report',
                     help='Output prefix for the multisample report. '
                          'Default: SQANTI_sc_multisample_report.')
    apd.add_argument('--pca_features', default=None,
                     help='Optional file with one cell-summary column name per line, '
                          'replacing the curated feature set in the multisample report.')
    add_clustering_args(
        cells.add_argument_group("Clustering and UMAP options (require --apply)"))

    return ap
