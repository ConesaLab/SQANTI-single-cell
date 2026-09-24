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
    apr.add_argument('-q', '--qc_dir', required=True,
                     help='Output directory of the QC run, read but never written.')
    apc = common.add_argument_group("Common options")
    apc.add_argument('-d', '--out_dir', default="filter",
                     help='Directory the filter writes into, laid out exactly like a QC '
                          'run so every downstream tool works on it unchanged. Threshold '
                          'trials are just different directories. Default: filter.')
    apc.add_argument('-l', '--log_level', default='INFO',
                     choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'],
                     help='Set logging level. Default: INFO.')

    sub = ap.add_subparsers(dest='subcommand')
    cells = sub.add_parser('cells', parents=[common],
                           help='Filter cell barcodes on per-cell QC metrics.')
    add_cell_filter_args(cells.add_argument_group("Cell filter options"))

    apd = cells.add_argument_group("Report options")
    apd.add_argument('--report', choices=["pdf", "html", "both", "skip"], default="skip",
                     help="Re-render the per-sample report on filtered data. Default: skip.")
    # The only report input that is not already in the QC outputs: the report reads the
    # annotation itself to compare reference and sample transcript lengths. The optional
    # evidence flags have no counterpart here because SQANTI3 is not re-run.
    apd.add_argument('--refGTF', help='Reference annotation (GTF), for the report.')
    apd.add_argument('--multisample_report', action='store_true', default=False,
                     help='Generate a multisample report from the filtered cell summaries.')
    apd.add_argument('--multisample_report_prefix', default='SQANTI_sc_multisample_report',
                     help='Output prefix for the multisample report. '
                          'Default: SQANTI_sc_multisample_report.')

    # Clustering is a pipeline stage the report consumes, not a report setting: with a
    # different set of cells the embedding has to be refitted, and reusing the QC run's
    # coordinates would place the kept cells in a space the discarded ones shaped.
    add_clustering_args(
        cells.add_argument_group("Clustering and UMAP options (re-run on the kept cells)"))

    return ap
