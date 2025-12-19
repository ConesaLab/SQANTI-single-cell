import argparse

def build_parser(version_str: str = '0.2.1'):
    ap = argparse.ArgumentParser(
        description="SQANTI-sc: Structural and Quality Annotation of Novel Transcript Isoforms at the Single-Cell level"
    )
    # Required arguments
    apr = ap.add_argument_group("Required arguments")
    apr.add_argument('--refFasta', required=True,
                     help='Reference genome (Fasta).')
    apr.add_argument('--refGTF', required=True,
                     help='Reference annotation (GTF).')
    apr.add_argument('-de', '--design', dest="inDESIGN", required=True,
                     help='Design file with sampleID and file_acc.')
    apr.add_argument('-m', '--mode', choices=["isoforms", "reads"],
                     required=True, help='Input data type.')

    # SQANTI-sc specific options
    apsc = ap.add_argument_group("SQANTI-sc specific options")
    apsc.add_argument('-d', '--out_dir', default=".",
                      help='Output directory. Default: current dir.')
    apsc.add_argument('-i', '--input_dir', default='.',
                      help='Input directory for fastq/gff/bam files.')
    apsc.add_argument('--report', choices=["pdf", "html", "both", "skip"],
                      default="pdf", help="Report format. Default: pdf.")
    apsc.add_argument('-@', '--samtools_cpus', default=10, type=int,
                      help='Threads for samtools. Default: 10.')
    apsc.add_argument('--verbose', action="store_true", default=False,
                      help='Print all steps.')

    apsc.add_argument('--skip_hash', dest="SKIPHASH", action='store_true',
                      help='Skip UJC hashing step.')
    apsc.add_argument('--ignore_cell_summary', action="store_true",
                      default=False,
                      help="Don't save cell summary table in report.")
    apsc.add_argument('--multisample_report', action='store_true', default=False,
                     help='Generate a multisample cohort report from per-sample cell summaries.')
    apsc.add_argument('--multisample_report_prefix', default='SQANTI_sc_multisample_report',
                     help='Output prefix for the multisample report (default: SQANTI_sc_multisample_report).')
    apsc.add_argument('--write_per_cell_outputs', action='store_true', default=False,
                       help='When set, writes per-cell gene_counts.csv, ujc_counts.csv, and cv.csv for each sample.')

    # SQANTI3 Customization and filtering
    apc = ap.add_argument_group("SQANTI3 customization and filtering")
    apc.add_argument('--min_ref_len', type=int, default=0,
                     help="Min reference transcript length. Default: 0.")
    apc.add_argument('--genename', action='store_true',
                     help='Use gene_name tag instead of gene_id.')
    apc.add_argument('--novel_gene_prefix',
                     help='Prefix for novel gene IDs.')
    apc.add_argument('--ref_cov_min_pct', type=float, default=45.0,
                     help='Minimum %% of reference transcript length a read must cover to count towards coverage plots (default: 45.0).')
    apc.add_argument('-s', '--sites', default="ATAC,GCAG,GTAG",
                     help='Canonical splice sites. Default: ATAC,GCAG,GTAG.')
    apc.add_argument('-w', '--window', default=20, type=int,
                     help='Window size for Adenine content screen. Default: 20.')

    # SQANTI3 Aligner and mapping options
    apa = ap.add_argument_group("SQANTI3 aligner and mapping options")
    apa.add_argument('--aligner_choice', default='minimap2',
                     choices=['minimap2', "uLTRA", "gmap", "deSALT"],
                     help="Aligner to use. Default: minimap2.")
    apa.add_argument('-x', '--gmap_index',
                     help='Path to gmap_build index prefix.')

    # SQANTI3 ORF prediction
    apo = ap.add_argument_group("SQANTI3 ORF prediction")
    apo.add_argument('--skipORF', action="store_true", default=False,
                     help="Skip ORF prediction.")
    apo.add_argument('--orf_input',
                     help="Input fasta for ORF prediction.")

    # SQANTI3 Orthogonal data inputs
    apod = ap.add_argument_group("SQANTI3 orthogonal data inputs")
    apod.add_argument('--short_reads',
                     help='FOFN of short-read RNA-Seq files (FASTA/FASTQ).')
    apod.add_argument('--SR_bam',
                     help='Directory or FOFN of short-read BAM files.')
    apod.add_argument('--CAGE_peak', help="FANTOM5 CAGE Peak (BED).")
    apod.add_argument('--polyA_motif_list',
                     help="Ranked list of polyA motifs (text).")
    apod.add_argument('--polyA_peak', help="PolyA Peak (BED).")
    apod.add_argument('--phyloP_bed', help="PhyloP conservation scores (BED).")
    apod.add_argument('-e', '--expression',
                     help='Expression matrix (e.g., Kallisto tsv).')
    apod.add_argument('-c', '--coverage',
                     help='Junction coverage files (single, comma-separated, or pattern).')

    # SQANTI3 Functional annotation
    apf = ap.add_argument_group("SQANTI3 functional annotation")
    apf.add_argument('--isoAnnotLite', action='store_true',
                     help='Run isoAnnot Lite for tappAS compatibility.')
    apf.add_argument('--gff3',
                     help='Precomputed tappAS GFF3 for functional transfer.')

    # SQANTI3 Output options
    apout = ap.add_argument_group("SQANTI3 output options")
    apout.add_argument('--isoform_hits', action='store_true',
                       help='Report all FSM/ISM hits in a separate file.')
    apout.add_argument('--ratio_TSS_metric', default='max',
                       choices=['max', 'mean', 'median', '3quartile'],
                       help='Metric for ratio_TSS column. Default: max.')

    # SQANTI3 Performance options
    app = ap.add_argument_group("SQANTI3 performance options")
    app.add_argument('-t', '--mapping_cpus', default=10, type=int,
                     help='Threads for alignment. Default: 10.')
    app.add_argument('-n', '--chunks', default=10, type=int,
                     help='Number of chunks for analysis. Default: 10.')
    app.add_argument('-l', '--log_level', default='INFO',
                     choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'],
                     help='Set logging level. Default: INFO.')

    # SQANTI3 Optional arguments
    apm = ap.add_argument_group("SQANTI3 optional arguments")
    apm.add_argument('--is_fusion', action="store_true",
                     help="Input are fusion isoforms (GTF required).")
    apm.add_argument('-v', '--version', action='version',
                     version='sqanti-sc ' + version_str)

    # SQANTI-sc Clustering and UMAP options
    apcl = ap.add_argument_group("SQANTI-sc clustering and UMAP options")
    apcl.add_argument('--run_clustering', action='store_true', default=False,
                      help='Run cell clustering and UMAP analysis.')
    apcl.add_argument('--normalization', choices=['log1p', 'sqrt', 'pearson'], default='log1p',
                      help='Normalization method. Default: log1p.')
    apcl.add_argument('--n_neighbors', type=int, default=15,
                      help='Number of neighbors for UMAP. Default: 15.')
    apcl.add_argument('--n_pc', type=int, default=30,
                      help='Number of principal components. Default: 30.')
    apcl.add_argument('--resolution', type=float, default=0.5,
                      help='Resolution for Leiden clustering. Default: 0.5.')
    apcl.add_argument('--n_top_genes', type=int, default=2000,
                      help='Number of highly variable genes. Default: 2000.')
    apcl.add_argument('--clustering_method', choices=['leiden', 'louvain', 'kmeans'], default='leiden',
                      help='Clustering method. Default: leiden.')
    apcl.add_argument('--n_clusters', type=int, default=10,
                      help='Number of clusters for K-means. Default: 10.')

    return ap


