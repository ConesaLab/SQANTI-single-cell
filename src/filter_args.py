import argparse
import os

from qc_args import add_cell_metrics_args, add_clustering_args

DEFAULT_CELL_RULES = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'filter_assets', 'cell_filter_default.json')


def add_cell_filter_args(group):
    """Also the hook for a future --run_cell_filter stage in qc_pipeline, which is
    why the arguments live in a factory rather than inline in the parser."""
    group.add_argument('-j', '--rules', default=DEFAULT_CELL_RULES,
                       help='JSON of per-metric cell rules. Default: the bundled '
                            'cell_filter_default.json.')
    return group


TRANSCRIPT_METHODS = ('rules', 'ml')


def transcript_method(argv):
    """-j and the method options depend on --method, as they depend on the subcommand in
    SQANTI3, so the method is read before the parser is built. An unknown value is left
    for the real parser to reject."""
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument('--method')
    method = pre.parse_known_args(argv)[0].method
    return method if method in TRANSCRIPT_METHODS else 'rules'


def add_transcript_filter_args(group, method='rules'):
    group.add_argument('--method', choices=TRANSCRIPT_METHODS, default='rules',
                       help="SQANTI3 filter to run: 'rules' or 'ml' (machine learning). "
                            "The options below, -j included, are those of the chosen "
                            "method; 'transcripts --method ml -h' lists the ML ones. "
                            "Default: rules.")
    if method == 'rules':
        group.add_argument('-j', '--rules', default=None,
                           help='SQANTI3 rules filter JSON, keyed by structural category. '
                                'Default: SQANTI3\'s own filter_default.json.')
        return group

    # Unset values are not passed, so SQANTI3's own defaults apply.
    group.add_argument('-j', '--threshold', type=float, default=None,
                       help='Probability threshold to classify a model as an isoform. '
                            'Default: SQANTI3\'s.')
    group.add_argument('-t', '--percent_training', type=float, default=None,
                       help='Proportion of the training set used to train; the rest tests '
                            'the classifier. Default: SQANTI3\'s.')
    group.add_argument('-p', '--TP',
                       help='File of true-positive model IDs, one per line, no header. '
                            'Without --TP and --TN SQANTI3 builds both lists from the data.')
    group.add_argument('-n', '--TN',
                       help='File of true-negative model IDs, one per line, no header.')
    group.add_argument('-f', '--force_fsm_in', action='store_true', default=False,
                       help='Keep every FSM model whatever the classifier says.')
    group.add_argument('--intermediate_files', action='store_true', default=False,
                       help='Also write SQANTI3\'s ML input table.')
    group.add_argument('-r', '--remove_columns',
                       help='File of classification columns, one per line, no header, to '
                            'leave out of training. Cell barcodes, UMIs and junction '
                            'chains are always left out.')
    group.add_argument('-z', '--max_class_size', type=int, default=None,
                       help='Largest number of models in each of the TP and TN lists. '
                            'Default: SQANTI3\'s.')
    group.add_argument('-i', '--intrapriming', type=float, default=None,
                       help='Adenine percentage at the genomic 3\' end that flags a model '
                            'as intra-priming. Default: SQANTI3\'s.')
    return group


RESCUE_MODES = ('automatic', 'full')


def add_rescue_args(group):
    group.add_argument('--pre_filter_dir', required=True,
                       help="The directory the transcript filter read (its -q). Full rescue "
                            "aligns the removed models, whose sequences only that directory "
                            "still has.")
    group.add_argument('--refFasta', required=True,
                       help='Reference genome (FASTA) the QC run used. --refGTF is required '
                            'too.')
    group.add_argument('--mode', dest='rescue_mode', choices=RESCUE_MODES, default='automatic',
                       help="SQANTI3 rescue mode. 'automatic' puts back the reference "
                            "transcripts whose every FSM was removed; 'full' also aligns the "
                            "removed ISM, NIC and NNC models to the kept models and reference "
                            "transcripts of their genes. Default: automatic.")
    group.add_argument('-k', '--refClassif',
                       help='Classification of a SQANTI3 QC run on --refGTF with the QC run\'s '
                            'orthogonal data, used instead of running one. Its '
                            '_junctions.txt, _corrected.gtf and _corrected.fasta must sit '
                            'beside it.')
    group.add_argument('-j', '--rules', default=None,
                       help='Rules JSON full rescue checks reference targets with: the one the '
                            'transcript filter used, or any for ML filter output, whose '
                            'reference targets are checked by rules too. Default: SQANTI3\'s '
                            'filter_default.json.')
    return group


def add_reference_qc_args(group):
    """The QC options that shape a classification, with the QC parser's defaults: the
    reference QC's rows join the sample's classification."""
    group.add_argument('--CAGE_peak', help="FANTOM5 CAGE Peak (BED).")
    group.add_argument('--polyA_motif_list', help="Ranked list of polyA motifs (text).")
    group.add_argument('--polyA_peak', help="PolyA Peak (BED).")
    group.add_argument('--phyloP_bed', help="PhyloP conservation scores (BED).")
    group.add_argument('--include_ORF', action="store_true", default=False,
                       help="Include ORF prediction.")
    group.add_argument('-s', '--sites', default="ATAC,GCAG,GTAG",
                       help='Canonical splice sites. Default: ATAC,GCAG,GTAG.')
    group.add_argument('-w', '--window', default=20, type=int,
                       help='Window size for Adenine content screen. Default: 20.')
    group.add_argument('--ratio_TSS_metric', default='max',
                       choices=['max', 'mean', 'median', '3quartile'],
                       help='Metric for ratio_TSS column. Default: max.')
    group.add_argument('--min_ref_len', type=int, default=0,
                       help="Min reference transcript length. Default: 0.")
    group.add_argument('-n', '--chunks', default=10, type=int,
                       help='Chunks to split the reference QC into. Default: 10.')
    return group


def add_downstream_args(parser):
    apd = parser.add_argument_group("Report options")
    apd.add_argument('--report', choices=["pdf", "html", "both", "skip"], default="skip",
                     help="Re-render the per-sample report on filtered data. Default: skip.")
    # The only report input that is not already in the QC outputs: the report reads the
    # annotation itself to compare reference and sample transcript lengths. The optional
    # evidence flags have no counterpart here because SQANTI3 QC is not re-run.
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
        parser.add_argument_group("Clustering and UMAP options (re-run on the filtered data)"))


def build_filter_parser(version_str: str = '1.2.0', method: str = 'rules'):
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
                     help='Output directory of the QC run, or of an earlier filter step. '
                          'Read but never written.')
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
    add_downstream_args(cells)

    transcripts = sub.add_parser(
        'transcripts', parents=[common],
        help='Filter transcript models with the SQANTI3 rules or ML filter.')
    add_transcript_filter_args(transcripts.add_argument_group("Transcript filter options"),
                               method)
    add_cell_metrics_args(transcripts.add_argument_group(
        "Cell summary options (pass the values the QC run used)"))
    add_downstream_args(transcripts)

    rescue = sub.add_parser(
        'rescue', parents=[common],
        help="Put back what the transcript filter removed, as SQANTI3 rescue does, moving "
             "each removed model's reads to its target cell by cell.")
    add_rescue_args(rescue.add_argument_group("Rescue options"))
    add_reference_qc_args(rescue.add_argument_group(
        "Reference QC options, unless -k is given (pass the values the QC run used)"))
    add_cell_metrics_args(rescue.add_argument_group(
        "Cell summary options (pass the values the QC run used)"))
    add_downstream_args(rescue)

    return ap
