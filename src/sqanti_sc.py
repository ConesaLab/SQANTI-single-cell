#!/usr/bin/env python3
import subprocess, os, re, sys, glob
import argparse
import pandas as pd
import numpy as np
import shutil
import hashlib
import pysam

sqanti3_src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../SQANTI3/src"))
if sqanti3_src_path not in sys.path:
    sys.path.insert(0, sqanti3_src_path)

from commands import (
    RSCRIPTPATH, run_command
)

#!/usr/bin/env python3
# SQANTI_Single_Cell: Structural and Quality Annotation of Novel Transcripts in reads at the single cell level
# Author: Carlos Blanco

__author__  = "carlos.blanco@csic.es"
__version__ = '1.0'  # Python 3.11

utilitiesPath = os.path.abspath(os.path.join(os.path.dirname(__file__), "../SQANTI3/src/utilities"))
sys.path.insert(0, utilitiesPath)

sqantiqcPath = os.path.abspath(os.path.join(os.path.dirname(__file__), "../SQANTI3"))


FIELDS_JUNC = ['isoform', 'chrom', 'strand', 'junction_number', 'genomic_start_coord',
                   'genomic_end_coord', 'junction_category',
                   'diff_to_Ref_start_site', 'diff_to_Ref_end_site', 'canonical']

FIELDS_CLASS = ['isoform', 'chrom', 'strand', 'length',  'exons',  'structural_category',
                'associated_gene', 'associated_transcript',  'ref_length', 'ref_exons',
                'subcategory', 'RTS_stage', 'all_canonical',
                'predicted_NMD', 'perc_A_downstream_TTS', "jxn_string"]

RSCRIPTPATH = shutil.which('Rscript')


def fill_design_table(args):
    df = pd.read_csv(args.inDESIGN, sep=",")

    # If number of columns is less than 2, probably wrongly formatted
    if df.shape[1] < 2:
        print(f"ERROR: {args.inDESIGN} is incorrectly formatted, is it not separated by commas?", file=sys.stderr)
        sys.exit(1)

    # Check for required columns
    required_columns = {'sampleID', 'file_acc'}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        print(f"ERROR: Missing required columns: {', '.join(missing_columns)}", file=sys.stderr)
        sys.exit(1)

    # Create the new columns
    df['classification_file'] = args.input_dir + '/' + df['file_acc'] + '/' + df['sampleID'] + '_classification.txt'
    df['junction_file'] = args.input_dir + '/' + df['file_acc'] + '/' + df['sampleID'] + '_junctions.txt'
    df.to_csv(args.inDESIGN, sep=',', index=False)
    return df


def get_files_runSQANTI3(args, df):
    def check_files_exist(*files):
        for file in files:
            if not os.path.isfile(file):
                print(f'[ERROR] Missing file: {file}', file=sys.stdout)
                sys.exit(-1)

    def build_sqanti_command(input_file):
        cmd = (
            f"python {sqantiqcPath}/sqanti3_qc.py {input_file} {args.annotation} {args.genome} "
            f"--min_ref_len {args.min_ref_len} --aligner_choice {args.aligner_choice} "
            f"-t {args.cpus} -n {args.chunks} -d {args.out_dir}/{file_acc} -o {sampleID} -s {args.sites} --report skip"
        )
        if args.force_id_ignore:
            cmd += " --force_id_ignore"
        if args.skipORF:
            cmd += " --skipORF"
        if is_fastq:
            cmd += " --fasta"
        return cmd

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']

        # Check for .gtf or .gff file
        gtf_pattern = os.path.join(args.input_dir, f"{file_acc}*.g*f")
        gtf_files = glob.glob(gtf_pattern)
        if gtf_files:
            gtf_file = gtf_files[0]
            check_files_exist(args.genome, args.annotation)
            if args.verbose:
                print(f'[INFO] Running SQANTI-sc qc for sample {gtf_file}', file=sys.stdout)
            cmd_sqanti = build_sqanti_command(gtf_file)
            print(cmd_sqanti, file=sys.stdout)
            subprocess.run(cmd_sqanti, shell=True, check=True)
            continue

        # Check for .fastq or .fasta files
        fastq_pattern = os.path.join(args.input_dir, f"{file_acc}*.fast*")
        fastq_files = glob.glob(fastq_pattern)
        if fastq_files:
            fastq_file = fastq_files[0]
            check_files_exist(args.genome, args.annotation)
            if args.verbose:
                print(f'[INFO] Running SQANTI-sc qc for sample {fastq_file}', file=sys.stdout)
            cmd_sqanti = build_sqanti_command(fastq_file, is_fastq=True)
            print(cmd_sqanti, file=sys.stdout)
            subprocess.run(cmd_sqanti, shell=True, check=True)
            continue

        # Check for BAM files if mode is "reads"
        if args.mode == "reads":
            bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
            bam_files = glob.glob(bam_pattern)
            if bam_files:
                bam_file = bam_files[0]
                check_files_exist(args.genome, args.annotation)
                if args.verbose:
                    print(f'[INFO] Running SQANTI-sc qc for sample {bam_file}', file=sys.stdout)

                # Convert unmapped BAM to FASTQ
                cmd_tofastq = f"samtools fastq -@ {args.samtools_cpus} {bam_file} -n > {args.input_dir}/{file_acc}_tmp.fastq"
                try:
                    subprocess.run(cmd_tofastq, shell=True, check=True)
                except subprocess.CalledProcessError:
                    print(f"ERROR running command: {cmd_tofastq}\n Missing SAMTOOLS", file=sys.stderr)
                    sys.exit(-1)

                # Locate generated FASTQ
                fastq_pattern = os.path.join(args.input_dir, f"{file_acc}*.fastq")
                fastq_files = glob.glob(fastq_pattern)
                if not fastq_files:
                    print(f"[ERROR] FASTQ file not generated for {file_acc}", file=sys.stderr)
                    continue

                fastq_file = fastq_files[0]

                # Run SQANTI3
                cmd_sqanti = build_sqanti_command(fastq_file, is_fastq=True)
                print(cmd_sqanti, file=sys.stdout)
                subprocess.run(cmd_sqanti, shell=True, check=True)
                continue

        # If none of the conditions are met, raise an error
        print(f"ERROR: The file_acc you included in your design file does not correspond to .fastq, .gtf or directories with junctions and classification files in the {args.input_dir} directory", file=sys.stdout)

    # Cleanup
    try:
        os.remove(f"{args.input_dir}/{file_acc}_tmp.fastq")
        os.remove(f"{file_acc}_tmp.renamed.fasta")
    except FileNotFoundError:
        print(f"[WARNING] Failed to remove temporary file: {file_acc}_tmp.fastq", file=sys.stderr)


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


def add_cell_data(args, df):
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

        metadata_dict = {}

        try:
            with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
                for read in bam:
                    if read.has_tag("XM") and read.has_tag("CB"):
                        umi = read.get_tag("XM")
                        cell_barcode = read.get_tag("CB")
                        isoform = read.query_name
                        metadata_dict[isoform] = {"UMI": umi, "CB": cell_barcode}
        except Exception as e:
            print(f"[ERROR] Failed to parse BAM file {bam_file}: {str(e)}", file=sys.stderr)
            continue

        if not metadata_dict:
            print(f"[INFO] No valid reads found in the BAM file for {file_acc} for metadata creation.", file=sys.stdout)
            return None

        # Load the classification file
        classification_path = f"{outputPathPrefix}_classification_tmp.txt"
        if os.path.isfile(classification_path):
            try:
                classification_df = pd.read_csv(classification_path, sep='\t', low_memory=False)
                
                # Add UMI and CB information
                classification_df['UMI'] = classification_df['isoform'].map(lambda x: metadata_dict.get(x, {}).get("UMI", None))
                classification_df['CB'] = classification_df['isoform'].map(lambda x: metadata_dict.get(x, {}).get("CB", None))
                
                # Fill _dup2 supplementary alignments with their respective primary alignment's UMI and CB
                for index, row in classification_df.iterrows():
                    match = re.match(r"(.+)_dup\d+$", row['isoform'])  # Match _dup<number>
                    if match:
                        primary_isoform = match.group(1)  # Extract the original molecule ID
                        if primary_isoform in metadata_dict:
                            classification_df.at[index, 'UMI'] = metadata_dict[primary_isoform]['UMI']
                            classification_df.at[index, 'CB'] = metadata_dict[primary_isoform]['CB']
                            
                classification_df.to_csv(f"{outputPathPrefix}_classification.txt", index=False, sep="\t")

            except Exception as e:
                print(f"[ERROR] Could not merge metadata with classification table: {str(e)}", file=sys.stderr)
        else:
            print(f"[INFO] Classification file for {file_acc} not found.", file=sys.stdout)

        print(f"**** UMI and cell barcode information successfully added to {file_acc} reads classification", file=sys.stdout)

        os.remove(f"{outputPathPrefix}_classification_tmp.txt")


def generate_report(args, df):
    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.input_dir, file_acc, sampleID)

        classification_path = f"{outputPathPrefix}_classification.txt"
        if os.path.isfile(classification_path):
            try:
                print("**** Generating SQANTI3 report....", file=sys.stderr)
                cmd = f"{RSCRIPTPATH} {utilitiesPath}/SQANTI-sc_reads.R {classification_path} {utilitiesPath} {args.report}"
                run_command(cmd,"SQANTI3 report")

            except Exception as e:
                print(f"Error generating report for {classification_path}: {e}")


def main():
    global utilitiesPath
    global sqantiqcPath

    #arguments
    parser = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")
    # Positional arguments
    apr = ap.add_argument_group("Required arguments")
    parser.add_argument('--genome', type=str, required=True, help='\t\tReference genome (Fasta format).')
    parser.add_argument('--annotation', type=str, required=True, help='\t\tReference annotation file (GTF format).')
    parser.add_argument('-de', '--design', type=str, dest="inDESIGN", required=True, help='Path to design file, must have sampleID and file_acc column.')
    parser.add_argument('-m', '--mode', type=str, choices = ["isoforms", "reads"], required=True, help = '\t\tType of data to run SQANTI3 on (reads or isoforms)')

    # Optional arguments
    apr = ap.add_argument_group("Optional arguments")
    parser.add_argument('-i', '--input_dir', type=str, default = './', help = '\t\tPath to directory where fastq files are stored. Or path to parent directory with children directories of SQANTI3 runs. Default: Directory where the script was run.')
    parser.add_argument('-f', '--factor', type=str, dest="inFACTOR" ,required=False, help='This is the column name that plots are to be faceted by. Default: None')
    parser.add_argument('-p','--prefix', type=str, dest="PREFIX", required=False, help='SQANTI-sc output filename prefix. Default: sqantiSingleCell')
    parser.add_argument('-o','--out_dir', type=str, help='\t\tDirectory for output sqanti_sc files. Default: Directory where the script was run.', default = "./", required=False)
    parser.add_argument('--min_ref_len', type=int, default=0, help="\t\tMinimum reference transcript length. Default: 0 bp")
    parser.add_argument('--force_id_ignore', action="store_true", default=False, help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")
    parser.add_argument('--skipORF', action="store_true", default=False, help="\t\t Skip ORF prediction to save time.")
    parser.add_argument('--aligner_choice', type=str, choices=['minimap2', "uLTRA"], default='minimap2', help="\t\tDefault: minimap2")
    parser.add_argument('-@', '--samtools_cpus', default=10, type=int, help='\t\tNumber of threads used during conversion from bam to fastq. Default: 10')    
    parser.add_argument('-t', '--mapping_cpus', default=10, type=int, help='\t\tNumber of threads used during alignment by aligners. Default: 10')
    parser.add_argument('-n', '--chunks', default=10, type=int, help='\t\tNumber of chunks to split SQANTI3 analysis in for speed up. Default: 10')
    parser.add_argument('-s','--sites', type=str, default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)
    parser.add_argument('--skip_hash', dest="SKIPHASH", action='store_true', help='Skip the hashing step')
    parser.add_argument('--report', type=str, choices = ["pdf", "html", "both", "skip"], default = "pdf", help = "\t\tDefault: pdf")
    parser.add_argument('-c' '--ignore_cell_summary', action="store_true", default=False, help="\t\t Add this flag to not save the cell summary table generated during the report to save space. Do not add if running the cell filter module afterwards, as this table is used to inform the cell filtering.")
    parser.add_argument('--verbose', help = 'If verbose is run, it will print all steps, by default it is FALSE', action="store_true", default=False)
    parser.add_argument('-v', '--version', help="Display program version number.", action='version', version='sqanti-sc '+str(__version__))

    args = parser.parse_args()

    # Check and read design file
    df = fill_design_table(args)

    # Run SQANTI3 for reads
    if args.reads:
        get_files_runSQANTI3(args, df)

    if args.isoforms:
        get_files_runSQANTI3(args, df)

    # Make UJC and hash
    if not args.SKIPHASH:
        make_UJC_hash(args, df)

    # Make metadata table
    add_cell_data(args, df)

    # Run plotting script
    generate_report(args, df)


if __name__ == "__main__":
    main()
