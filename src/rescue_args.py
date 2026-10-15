import argparse

from filter_args import add_downstream_args
from qc_args import add_cell_metrics_args

RESCUE_MODES = ('automatic', 'full')


def add_rescue_args(group):
    group.add_argument('--mode', dest='rescue_mode', choices=RESCUE_MODES, default='automatic',
                       help="SQANTI3 rescue mode. 'automatic' puts back the reference "
                            "transcripts whose every FSM was removed; 'full' also aligns the "
                            "removed ISM, NIC and NNC models to the kept models and reference "
                            "transcripts of their genes. Default: automatic.")
    group.add_argument('-e', '--rescue_mono_exonic', choices=('all', 'fsm', 'none'),
                       default=None,
                       help="Which removed mono-exonic models can be rescued, as SQANTI3's -e: "
                            "'all', 'fsm' (mono-exonic FSMs only, by automatic rescue) or "
                            "'none'. Default: SQANTI3's (all).")
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


def build_rescue_parser(version_str: str = '1.2.0'):
    ap = argparse.ArgumentParser(
        description="SQANTI-sc rescue: put back what the transcript filter removed, as "
                    "SQANTI3 rescue does, moving each removed model's reads to its target "
                    "cell by cell"
    )
    ap.add_argument('-v', '--version', action='version', version='sqanti-sc-rescue ' + version_str)

    apr = ap.add_argument_group("Required arguments")
    apr.add_argument('-de', '--design', dest="inDESIGN", required=True,
                     help='Design file with sampleID and file_acc, as used for the QC run.')
    apr.add_argument('-q', '--qc_dir', required=True,
                     help="Output directory of the transcript filter (sqanti_sc_filter.py "
                          "transcripts). Read but never written.")
    apr.add_argument('--pre_filter_dir', required=True,
                     help="The directory the transcript filter read (its -q). Full rescue "
                          "aligns the removed models, whose sequences only that directory "
                          "still has.")
    apr.add_argument('--refGTF', required=True,
                     help='Reference annotation (GTF) the QC run used.')
    apr.add_argument('--refFasta', required=True,
                     help='Reference genome (FASTA) the QC run used.')
    apc = ap.add_argument_group("Common options")
    apc.add_argument('-d', '--out_dir', default="rescue",
                     help='Directory the rescue writes into, laid out exactly like a QC run '
                          'so every downstream tool works on it unchanged. Must be new or '
                          'empty. Default: rescue.')
    apc.add_argument('-l', '--log_level', default='INFO',
                     choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'],
                     help='Set logging level. Default: INFO.')

    add_rescue_args(ap.add_argument_group("Rescue options"))
    add_reference_qc_args(ap.add_argument_group(
        "Reference QC options, unless -k is given (pass the values the QC run used)"))
    add_cell_metrics_args(ap.add_argument_group(
        "Cell summary options (pass the values the QC run used)"))
    add_downstream_args(ap, ref_gtf=False)
    return ap
