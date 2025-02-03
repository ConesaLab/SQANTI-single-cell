#!/usr/bin/env python3
import subprocess, os, re, sys, glob
import argparse
import pandas as pd
import numpy as np
import shutil
import hashlib
import pysam

#!/usr/bin/env python3
# SQANTI_Single_Cell: Structural and Quality Annotation of Novel Transcripts in reads at the single cell level
# Author: Carlos Blanco

__author__  = "carlos.blanco@csic.es"
__version__ = '1.0'  # Python 3.11

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
sys.path.insert(0, utilitiesPath)

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
sqantiqcPath = os.path.join(os.path.dirname(os.path.realpath(__file__)))


FIELDS_JUNC = ['isoform', 'chrom', 'strand', 'junction_number', 'genomic_start_coord',
                   'genomic_end_coord', 'junction_category',
                   'diff_to_Ref_start_site', 'diff_to_Ref_end_site', 'canonical']

FIELDS_CLASS = ['isoform', 'chrom', 'strand', 'length',  'exons',  'structural_category',
                'associated_gene', 'associated_transcript',  'ref_length', 'ref_exons',
                'subcategory', 'RTS_stage', 'all_canonical',
                'predicted_NMD', 'perc_A_downstream_TTS', "jxn_string"]

RSCRIPTPATH = shutil.which('Rscript')


def fill_design_table(args):
    df = pd.read_csv(args.inDESIGN, sep = ",")
    # If number of columns is less than 2, probably wrongly formatted
    if df.shape[1] < 2:
        print("ERROR: is incorrectly formatted, is it not separated by commas?".format(args.inDESIGN), file=sys.stderr)
        sys.exit(-1)
    
    # Create the new columns
    df['classification_file'] = args.input_dir + '/' + df['file_acc'] + '/' + df['sampleID'] + '_classification.txt'
    df['junction_file'] = args.input_dir + '/' + df['file_acc'] + '/' + df['sampleID'] + '_junctions.txt'
    df.to_csv(args.inDESIGN, sep = ',', index = False)
    return(df)


def convert_and_runSQANTI3(args, df):
    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']

        bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
        bam_files = glob.glob(bam_pattern)

        if not bam_files:
            print(
                f"ERROR: The file_acc {file_acc} you included in your design file does not correspond to any .bam file in the {args.input_dir} directory",
                file=sys.stderr,
            )
            continue

        bam_file = bam_files[0]

        if not os.path.isfile(bam_file):
            print(f"[ERROR] BAM file {bam_file} is not accessible. Skipping...", file=sys.stderr)
            continue

        # Check required files
        if not os.path.isfile(args.genome):
            print("[ERROR] You inputted bam files to map but no reference genome FASTA", file=sys.stderr)
            sys.exit(-1)

        if not os.path.isfile(args.annotation):
            print("[ERROR] You inputted bam files to map but no reference annotation GTF", file=sys.stderr)
            sys.exit(-1)

        # Convert unmapped BAM to FASTQ
        cmd_tofastq = f"samtools fastq -@ {args.samtools_cpus} {bam_file} -n > {args.input_dir}/{file_acc}_tmp.fastq"
        try:
            subprocess.run(cmd_tofastq, shell=True, check=True)
        except subprocess.CalledProcessError:
            print(f"ERROR running command: {0}\n Missing SAMTOOLS".format(cmd_tofastq), file=sys.stderr)
            sys.exit(-1)

        # Locate generated FASTQ
        fastq_pattern = os.path.join(args.input_dir, f"{file_acc}*.fastq")
        fastq_files = glob.glob(fastq_pattern)
        if not fastq_files:
            print(f"[ERROR] FASTQ file not generated for {file_acc}", file=sys.stderr)
            continue

        fastq_file = fastq_files[0]

        # Run SQANTI3
        cmd_sqanti = (
            f"python {sqantiqcPath}/sqanti3_qc.py {fastq_file} {args.annotation} {args.genome} "
            f"--skipORF --min_ref_len {args.min_ref_len} --aligner_choice {args.aligner_choice} "
            f"-t {args.mapping_cpus} -n {args.chunks} -d {args.input_dir}/{file_acc} "
            f"-o {sampleID} -s {args.sites} --fasta --report {args.report}"
        )
        if args.force_id_ignore:
            cmd_sqanti += " --force_id_ignore"

        try:
            subprocess.run(cmd_sqanti, shell=True, check=True)
        except subprocess.CalledProcessError:
            print(f"[ERROR] Failed to run SQANTI3: {cmd_sqanti}", file=sys.stderr)
            sys.exit(-1)

        # Cleanup
        try:
            os.remove(f"{fastq_file}")
            os.remove(f"{file_acc}_tmp.renamed.fasta")
        except FileNotFoundError:
            print(f"[WARNING] Failed to remove temporary file: {fastq_file}", file=sys.stderr)


def make_UJC_hash(args, df):

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        # Input dir, sqanti3 dir, samplename
        outputPathPrefix = os.path.join(args.input_dir, file_acc, sampleID)

        print(f"**** Calculating UJCs for {file_acc}...", file=sys.stdout)

        ## Take the corrected GTF
        introns_cmd = f"""gtftools -i {outputPathPrefix}tmp_introns.bed -c "$(cut -f 1 {outputPathPrefix}_corrected.gtf.cds.gff | sort | uniq | paste -sd ',' - | sed 's/chr//g')" {outputPathPrefix}_corrected.gtf.cds.gff"""
        ujc_cmd = f"""awk -F'\t' -v OFS="\t" '{{print $5,"chr"$1,$4,$2+1"_"$3}}' {outputPathPrefix}tmp_introns.bed | bedtools groupby -g 1 -c 2,3,4 -o distinct,distinct,collapse | sed 's/,/_/g' | awk -F'\t' -v OFS="\t" '{{print $1,$2"_"$3"_"$4}}' > {outputPathPrefix}tmp_UJC.txt"""

        if subprocess.check_call(introns_cmd, shell=True) != 0:
            print("ERROR running command: {0}\n Missing GTFTOOLS".format(introns_cmd), file=sys.stderr)
            sys.exit(-1)

        if os.path.exists(f"{outputPathPrefix}_corrected.gtf.ensembl"):
            os.remove(f"{outputPathPrefix}_corrected.gtf.ensembl")

        if subprocess.check_call(ujc_cmd, shell=True) != 0:
            print("ERROR running command: {0}\n Missing BEDTOOLS".format(introns_cmd), file=sys.stderr)
            sys.exit(-1)
        os.remove(f"{outputPathPrefix}tmp_introns.bed")

        ## Pandas merge to the left
        classfile = f"{outputPathPrefix}_classification.txt"
        clas_df = pd.read_csv(classfile, sep="\t", usecols=[0, 1, 2, 7], dtype="str")
        clas_df.columns = ["isoform", "chr", "strand", "associated_transcript"]
        ujc_df = pd.read_csv(f"{outputPathPrefix}tmp_UJC.txt", sep="\t", names=["isoform", "jxn_string"], dtype="str")

        merged_df = pd.merge(clas_df, ujc_df, on="isoform", how="left")
        # Fill missing values in UJC column using the transcript ID
        merged_df["jxn_string"] = merged_df.apply(
            lambda row: row["chr"] + "_" + row["strand"] + "_" + "monoexon" + "_" + row["associated_transcript"]
            if pd.isna(row["jxn_string"])
            else row["jxn_string"],
            axis=1,
        )

        merged_df['jxnHash'] = merged_df['jxn_string'].apply(
            lambda x: hashlib.sha256(x.encode('utf-8')).hexdigest()
        )

        merged_df.to_csv(f"{outputPathPrefix}_temp.txt", index=False, sep="\t")

        cmd_paste = f"""bash -c 'paste <(cat {classfile} | tr -d '\r') <(cut -f 5,6 {outputPathPrefix}_temp.txt | tr -d '\r') > {outputPathPrefix}_classification_tmp.txt'"""
        subprocess.call(cmd_paste, shell=True)

        print(f"**** UJCs successfully calculated and added to {file_acc} classification", file=sys.stdout)

        os.remove(f"{outputPathPrefix}tmp_UJC.txt")
        os.remove(f"{outputPathPrefix}_temp.txt")


def make_metadata(args, df):
    """
    Parse the BAM file to extract isoform, UMI, and cell barcode information 
    and generate a df to add the UMI and cell barcode information to the classification.
    """

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.input_dir, file_acc, sampleID)

        bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
        bam_files = glob.glob(bam_pattern)
        if not bam_files:
            print(f"[ERROR] BAM file for {file_acc} not found in {args.input_dir}. Skipping...", file=sys.stderr)
            continue

        bam_file = bam_files[0]

        if not os.path.isfile(bam_file):
            print(f"[ERROR] BAM file {bam_file} is not accessible. Skipping...", file=sys.stderr)
            continue

        metadata_list = []

        try:
            with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
                for read in bam:
                    if read.has_tag("XM") and read.has_tag("CB"):
                        umi = read.get_tag("XM")
                        cell_barcode = read.get_tag("CB")
                        isoform = read.query_name
                        metadata_list.append({
                            "isoform": isoform,
                            "XM": umi,
                            "CB": cell_barcode
                        })
        except Exception as e:
            print(f"[ERROR] Failed to parse BAM file {bam_file}: {str(e)}", file=sys.stderr)
            continue

        # Create a DataFrame from the metadata list
        metadata_df = pd.DataFrame(metadata_list)

        if metadata_df.empty:
            print(f"[INFO] No valid reads found in the BAM file for {file_acc} for metadata creation.", file=sys.stdout)
            return None

        # Merge the metadata_df with the classification file
        classification_path = f"{outputPathPrefix}_classification_tmp.txt"
        if os.path.isfile(classification_path):
            try:
                classification_df = pd.read_csv(classification_path, sep='\t')
                merged_df = pd.merge(classification_df, metadata_df, on="isoform", how="left")
                merged_df.to_csv(f"{outputPathPrefix}_classification.txt", index = False, sep = "\t")

            except Exception as e:
                print(f"[ERROR] Could not merge metadata with classification table: {str(e)}", file=sys.stderr)
        else:
            print(f"[INFO] Classification file for {file_acc} not found.", file=sys.stdout)

        print(f"**** UMI and cell barcode information successfully added to {file_acc} classification", file=sys.stdout)

        os.remove(f"{outputPathPrefix}_classification_tmp.txt")


def make_cell_summary(args, df):
    """
    Create a classification-like file with the information at the cell level from the SQANTI3 classification file.

    Args:
        args: A namespace object containing input_dir (directory with the classification files).
        df: A pandas DataFrame with metadata containing 'file_acc' and 'sampleID' columns.

    Output:
        A new cell-level summary table saved as a .txt file.
    """
    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.input_dir, file_acc, sampleID)

        classification_path = f"{outputPathPrefix}_classification.txt"

        if os.path.isfile(classification_path):
            try:
                # Read the classification file
                classification_df = pd.read_csv(classification_path, sep="\t")

                # Group by cell barcode
                grouped = classification_df.groupby('CB')

                # Aggregating metrics for each cell
                cell_summary = grouped.agg(
                    total_reads=('CB', 'size'),
                    total_umi=('XM', lambda x: x.nunique()),
                    total_genes=('associated_gene', lambda x: x[classification_df.loc[x.index, 'structural_category'].isin(['genic_intron', 'intergenic']) == False].nunique()),
                    unique_junction_chains=('jxn_string', lambda x: x[classification_df.loc[x.index, 'exons'] > 1].nunique()),
                    median_read_length=('length', 'median'),
                    median_exons=('exons', 'median'),
                    mitochondrial_proportion=('chrom', lambda x: (x == 'MT').mean()),
                    fsm_proportion=('structural_category', lambda x: (x == 'full-splice_match').mean()),
                    ism_proportion=('structural_category', lambda x: (x == 'incomplete-splice_match').mean()),
                    nic_proportion=('structural_category', lambda x: (x == 'novel_in_catalog').mean()),
                    nnc_proportion=('structural_category', lambda x: (x == 'novel_not_in_catalog').mean()),
                    genic_proportion=('structural_category', lambda x: (x == 'genic').mean()),
                    genic_intron_proportion=('structural_category', lambda x: (x == 'genic_intron').mean()),
                    intergenic_proportion=('structural_category', lambda x: (x == 'intergenic').mean()),
                    antisense_proportion=('structural_category', lambda x: (x == 'antisense').mean()),
                    fusion_proportion=('structural_category', lambda x: (x == 'fusion').mean()),
                    rts_artifact_proportion = ('RTS_stage', lambda x: (x.astype(str) == 'True').mean()),
                    intrapriming_proportion=('perc_A_downstream_TTS', lambda x: (x >= 60).mean()),
                    all_canonical_proportion=(
                        'all_canonical', 
                        lambda x: (
                            x[
                                classification_df.loc[x.index, 'subcategory']
                                .isin(['mono-exon', 'mono-exon_by_intron_retention']) == False
                            ]
                            .eq('non_canonical')
                            .mean()
                        )
                    ),
                    monoexon_proportion=('exons', lambda x: (x == 1).mean())
                ).reset_index()

                # Save to output file
                output_file = f"{outputPathPrefix}_cell_summary.txt"
                cell_summary.to_csv(output_file, sep='\t', index=False)
                print(f"Cell-level summary table saved to {output_file}")

            except Exception as e:
                print(f"Error processing {classification_path}: {e}")


def main():
    global utilitiesPath
    global sqantiqcPath

    #arguments
    parser = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")
    parser.add_argument('--genome', type=str, help='\t\tReference genome (Fasta format).', default = False, required = False)
    parser.add_argument('--annotation', type=str, help='\t\tReference annotation file (GTF format).', default = False, required = True)
    parser.add_argument('-de', '--design', type=str, dest="inDESIGN" ,required=True, help='Path to design file, must have sampleID and file_acc column.')
    parser.add_argument('-i', '--input_dir', type=str, default = './', help = '\t\tPath to directory where fastq/GTF files are stored. Or path to parent directory with children directories of SQANTI3 runs. Default: Directory where the script was run.')
    parser.add_argument('-f', '--factor', type=str, dest="inFACTOR" ,required=False, help='This is the column name that plots are to be faceted by. Default: None')
    parser.add_argument('-p','--prefix', type=str, dest="PREFIX", required=False, help='SQANTI-sc output filename prefix. Default: sqantiSingleCell')
    parser.add_argument('-d','--dir', type=str, help='\t\tDirectory for output sqanti_sc files. Default: Directory where the script was run.', default = "./", required=False)
    parser.add_argument('--min_ref_len', type=int, default=0, help="\t\tMinimum reference transcript length. Default: 0 bp")
    parser.add_argument('--force_id_ignore', action="store_true", default=False, help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")
    parser.add_argument('--aligner_choice', type=str, choices=['minimap2', "uLTRA"], default='minimap2', help="\t\tDefault: minimap2")
    parser.add_argument('-@', '--samtools_cpus', default=10, type=int, help='\t\tNumber of threads used during conversion from bam to fastq. Default: 10')    
    parser.add_argument('-t', '--mapping_cpus', default=10, type=int, help='\t\tNumber of threads used during alignment by aligners. Default: 10')
    parser.add_argument('-n', '--chunks', default=10, type=int, help='\t\tNumber of chunks to split SQANTI3 analysis in for speed up. Default: 10')
    parser.add_argument('-s','--sites', type=str, default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)
    parser.add_argument('-ge','--gene_expression', type=int, dest="ANNOTEXP", required=False, help='Expression cut off level for determining underannotated genes. Default = 100', default = 100)
    parser.add_argument('-je','--jxn_expression', type=int, dest="JXNEXP", required=False, help='Coverage threshold for detected reference donors and acceptor. Default = 10', default = 10)
    parser.add_argument('-pc','--perc_coverage', type=int, dest="PERCCOV", required=False, help='Percent gene coverage of UJC for determining well-covered unannotated transcripts. Default = 20', default = 20)
    parser.add_argument('-pj','--perc_junctions', type=int, dest="PERCMAXJXN", required=False, help='Percent of the max junctions in gene for determining near full-length putative novel transcripts. Default = 80', default = 80)
    parser.add_argument('-fl','--factor_level', type=str, dest="FACTORLVL", required=False, help='Factor level to evaluate for underannotation', default = None)
    parser.add_argument('--all_tables', dest="ALLTABLES", action='store_true', help='Export all output tables. Default tables are gene counts, ujc counts, length_summary, cv and and underannotated gene tables')
    parser.add_argument('--pca_tables', dest="PCATABLES", action='store_true', help='Export table for making PCA plots')
    parser.add_argument('--report', type=str, choices = ["pdf", "html", "both", "skip"], default = "skip", help = "\t\tDefault: skip")
    parser.add_argument('--verbose', help = 'If verbose is run, it will print all steps, by default it is FALSE', action="store_true")
    parser.add_argument('-v', '--version', help="Display program version number.", action='version', version='sqanti-sc '+str(__version__))

    args = parser.parse_args()

    # Check and read design file
    df = fill_design_table(args)

    # Check method and run SQANTI3
    convert_and_runSQANTI3(args, df)

    # Make UJC and hash
    make_UJC_hash(args, df)

    # Make metadata table
    make_metadata(args, df)

    # Make cell summary
    make_cell_summary(args, df)

    # Run plotting script
    # plotting_script_path = os.path.join(os.path.dirname(__file__), 'utilities', 'sqanti_reads_tables_and_plots_02ndk.py')

    #cmd_plotting = f"python {plotting_script_path} --ref {args.annotation} --design {args.inDESIGN} -o {args.dir} --gene-expression {args.ANNOTEXP} --jxn-expression {args.JXNEXP} --perc-coverage {args.PERCCOV} --perc-junctions {args.PERCMAXJXN} --report {args.report}"
    #if args.inFACTOR:
    #    cmd_plotting = cmd_plotting + f" --factor {args.inFACTOR}"
    #if args.FACTORLVL != None:
    #    cmd_plotting = cmd_plotting + f" --factor-level {args.FACTORLVL}"
    #if args.PREFIX:
    #    cmd_plotting = cmd_plotting + f" --prefix {args.PREFIX}"
    #else:
    #    cmd_plotting = cmd_plotting + " --prefix sqantiSingleCell"
    #if args.ALLTABLES:
    #    cmd_plotting = cmd_plotting + " --all-tables"
    #if args.PCATABLES:
    #    cmd_plotting = cmd_plotting + "--pca-tables"
    
    #subprocess.call(cmd_plotting, shell = True)


if __name__ == "__main__":
    main()
