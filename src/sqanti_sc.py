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

RSCRIPTPATH = shutil.which('Rscript')


def fill_design_table(args):
    # Read the design CSV file into a DataFrame
    df = pd.read_csv(args.inDESIGN, sep=",")

    # Check for required columns
    required_columns = {'sampleID', 'file_acc'}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(f"ERROR: Missing required columns in design table: {missing_columns}")

    # Create the new columns based on the existing ones
    df['classification_file'] = args.input_dir + '/' + df['file_acc'] + '/' + df['sampleID'] + '_classification.txt'
    df['junction_file'] = args.input_dir + '/' + df['file_acc'] + '/' + df['sampleID'] + '_junctions.txt'

    # Write the DataFrame back to the CSV file
    df.to_csv(args.inDESIGN, sep=',', index=False)

    return df


def get_files_runSQANTI3(args, df):
    def check_files_exist(*files):
        """Check if each file exists, return False if missing instead of exiting."""
        missing_files = [file for file in files if not os.path.isfile(file)]
        if missing_files:
            print(f'[ERROR] Missing file(s): {", ".join(missing_files)}', file=sys.stderr)
            return False
        return True

    def build_sqanti_command(input_file, is_fastq=False):
        """Build the SQANTI command to run."""
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

    def run_command(cmd, error_message):
        """Run a shell command and handle errors."""
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError:
            print(f"[ERROR] {error_message}: {cmd}", file=sys.stderr)
            sys.exit(-1)

    # Validate genome and annotation before processing files
    if not check_files_exist(args.genome, args.annotation):
        print("[ERROR] Genome or annotation file is missing.", file=sys.stderr)
        sys.exit(-1)

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']

        # Check for .gtf or .gff file
        gtf_pattern = os.path.join(args.input_dir, f"{file_acc}*.g*f")
        gtf_files = glob.glob(gtf_pattern)
        if gtf_files:
            gtf_file = gtf_files[0]
            if args.verbose:
                print(f'[INFO] Running SQANTI-sc qc for sample {gtf_file}', file=sys.stdout)
            cmd_sqanti = build_sqanti_command(gtf_file)
            print(cmd_sqanti, file=sys.stdout)
            run_command(cmd_sqanti, "SQANTI3 failed to execute")
            continue

        # Check for .fastq or .fasta files
        fastq_pattern = os.path.join(args.input_dir, f"{file_acc}*.fast*")
        fastq_files = glob.glob(fastq_pattern)
        if fastq_files:
            fastq_file = fastq_files[0]
            if args.verbose:
                print(f'[INFO] Running SQANTI-sc qc for sample {fastq_file}', file=sys.stdout)
            cmd_sqanti = build_sqanti_command(fastq_file, is_fastq=True)
            print(cmd_sqanti, file=sys.stdout)
            run_command(cmd_sqanti, "SQANTI3 failed to execute")
            continue

        # Check for BAM files if mode is "reads"
        if args.mode == "reads":
            bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
            bam_files = glob.glob(bam_pattern)
            if bam_files:
                bam_file = bam_files[0]
                if args.verbose:
                    print(f'[INFO] Running SQANTI-sc qc for sample {bam_file}', file=sys.stdout)

                # Convert unmapped BAM to FASTQ
                tmp_fastq = f"{args.input_dir}/{file_acc}_tmp.fastq"
                cmd_tofastq = f"samtools fastq -@ {args.samtools_cpus} {bam_file} -n > {tmp_fastq}"
                run_command(cmd_tofastq, "Samtools failed to execute")

                # Ensure FASTQ was generated
                if not os.path.isfile(tmp_fastq):
                    print(f"[ERROR] FASTQ file not generated for {file_acc}", file=sys.stderr)
                    continue

                # Run SQANTI3 in fasta mode
                cmd_sqanti = build_sqanti_command(tmp_fastq, is_fastq=True)
                print(cmd_sqanti, file=sys.stdout)
                run_command(cmd_sqanti, "SQANTI3 failed to execute")

                # Cleanup temporary files
                try:
                    os.remove(tmp_fastq)
                except FileNotFoundError:
                    print(f"[WARNING] Failed to remove temporary file: {tmp_fastq}", file=sys.stderr)

                continue

        # No valid files found
        print(f"[ERROR] No valid .fastq, .gtf, or .bam files found for {file_acc} in {args.input_dir}", file=sys.stderr)


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
    Extract isoform, UMI, and cell barcode information from BAM file or a cell association file (depending on user's input)
    and add it to classification.
    """

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.input_dir, file_acc, sampleID)

        bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
        bam_files = glob.glob(bam_pattern)
        cell_association_file = os.path.join(args.input_dir, f"{file_acc}_cell_association.txt")

        cell_association_dict = {}

        # 1. Try to read cell data from BAM file (if exists)
        if bam_files:
            bam_file = bam_files[0]
            if not os.path.isfile(bam_file):
                print(f"[ERROR] BAM file {bam_file} is not accessible. Skipping...", file=sys.stderr)
                continue

            try:
                with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
                    for read in bam:
                        if read.has_tag("XM") and read.has_tag("CB"):
                            umi = read.get_tag("XM")
                            cell_barcode = read.get_tag("CB")
                            isoform = read.query_name
                            cell_association_dict[isoform] = {"UMI": umi, "CB": cell_barcode}
            except Exception as e:
                print(f"[ERROR] Failed to parse BAM file {bam_file}: {str(e)}", file=sys.stderr)
                # Fallback to cell association file (if available)
        # 2. If not bam file, try to read from cell association file
        elif os.path.isfile(cell_association_file):
            try:
                cell_association_df = pd.read_csv(cell_association_file, sep='\t', dtype=str) # Added: load cell association
                # Checks if `id` and `cell_barcode` columns exist, umi is not mandatory
                if 'id' not in cell_association_df.columns or 'cell_barcode' not in cell_association_df.columns:
                     print(f"[ERROR] Cell association file {cell_association_file} does not contain the required columns ('id' and 'cell_barcode'). Skipping...", file=sys.stderr)
                     continue
                for _, row in cell_association_df.iterrows():
                    isoform = row['id']
                    cell_barcode = row['cell_barcode']
                    umi = row.get('umi')
                    cell_association_dict[isoform] = {"CB": cell_barcode}
                    if umi:
                        cell_association_dict[isoform]["UMI"] = umi
            except Exception as e:
                print(f"[ERROR] Failed to parse cell_association file {cell_association_file}: {str(e)}", file=sys.stderr)
                continue

        # 3. Add cell information to the classification file
        if not cell_association_dict:
            print(f"[INFO] No valid reads found in the BAM or cell association file for {file_acc} for cell association creation.", file=sys.stdout)
            return None

        # Load the classification file
        classification_path = f"{outputPathPrefix}_classification_tmp.txt"
        if os.path.isfile(classification_path):
            try:
                classification_df = pd.read_csv(classification_path, sep='\t', low_memory=False)
                
                # Add UMI and CB information
                classification_df['UMI'] = classification_df['isoform'].map(lambda x: cell_association_dict.get(x, {}).get("UMI", None))
                classification_df['CB'] = classification_df['isoform'].map(lambda x: cell_association_dict.get(x, {}).get("CB", None))
                
                # Fill _dup2 supplementary alignments with their respective primary alignment's UMI and CB
                for index, row in classification_df.iterrows():
                    match = re.match(r"(.+)_dup\d+$", row['isoform'])
                    if match:
                        primary_isoform = match.group(1)
                        if primary_isoform in cell_association_dict:
                            classification_df.at[index, 'UMI'] = cell_association_dict[primary_isoform].get('UMI', None)
                            classification_df.at[index, 'CB'] = cell_association_dict[primary_isoform]['CB']
                            
                classification_df.to_csv(f"{outputPathPrefix}_classification.txt", index=False, sep="\t")

            except Exception as e:
                print(f"[ERROR] Could not merge cell association with classification table: {str(e)}", file=sys.stderr)
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
    apr = parser.add_argument_group("Required arguments")
    apr.add_argument('--genome', type=str, required=True, help='\t\tReference genome (Fasta format).')
    apr.add_argument('--annotation', type=str, required=True, help='\t\tReference annotation file (GTF format).')
    apr.add_argument('-de', '--design', type=str, dest="inDESIGN", required=True, help='Path to design file, must have sampleID and file_acc column.')
    apr.add_argument('-m', '--mode', type=str, choices = ["isoforms", "reads"], required=True, help = '\t\tType of data to run SQANTI3 on (reads or isoforms)')

    # Optional arguments
    apo = parser.add_argument_group("Optional arguments")
    apo.add_argument('-i', '--input_dir', type=str, default = './', help = '\t\tPath to directory where fastq files are stored. Or path to parent directory with children directories of SQANTI3 runs. Default: Directory where the script was run.')
    apo.add_argument('-f', '--factor', type=str, dest="inFACTOR" ,required=False, help='This is the column name that plots are to be faceted by. Default: None')
    apo.add_argument('-p','--prefix', type=str, dest="PREFIX", required=False, help='SQANTI-sc output filename prefix. Default: sqantiSingleCell')
    apo.add_argument('-o','--out_dir', type=str, help='\t\tDirectory for output sqanti_sc files. Default: Directory where the script was run.', default = "./", required=False)
    apo.add_argument('--min_ref_len', type=int, default=0, help="\t\tMinimum reference transcript length. Default: 0 bp")
    apo.add_argument('--force_id_ignore', action="store_true", default=False, help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")
    apo.add_argument('--skipORF', action="store_true", default=False, help="\t\t Skip ORF prediction to save time.")
    apo.add_argument('--aligner_choice', type=str, choices=['minimap2', "uLTRA"], default='minimap2', help="\t\tDefault: minimap2")
    apo.add_argument('-@', '--samtools_cpus', default=10, type=int, help='\t\tNumber of threads used during conversion from bam to fastq. Default: 10')    
    apo.add_argument('-t', '--mapping_cpus', default=10, type=int, help='\t\tNumber of threads used during alignment by aligners. Default: 10')
    apo.add_argument('-n', '--chunks', default=10, type=int, help='\t\tNumber of chunks to split SQANTI3 analysis in for speed up. Default: 10')
    apo.add_argument('-s','--sites', type=str, default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)
    apo.add_argument('--skip_hash', dest="SKIPHASH", action='store_true', help='Skip the hashing step')
    apo.add_argument('--report', type=str, choices = ["pdf", "html", "both", "skip"], default = "pdf", help = "\t\tDefault: pdf")
    apo.add_argument('-c' '--ignore_cell_summary', action="store_true", default=False, help="\t\t Add this flag to not save the cell summary table generated during the report to save space. Do not add if running the cell filter module afterwards, as this table is used to inform the cell filtering.")
    apo.add_argument('--verbose', help = 'If verbose is run, it will print all steps, by default it is FALSE', action="store_true", default=False)
    apo.add_argument('-v', '--version', help="Display program version number.", action='version', version='sqanti-sc '+str(__version__))

    args = parser.parse_args()

    # Check and read design file
    try:
        df = fill_design_table(args)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Run SQANTI3 for reads
    if args.reads:
        get_files_runSQANTI3(args, df)

    if args.isoforms:
        get_files_runSQANTI3(args, df)

    # Make UJC and hash
    if not args.SKIPHASH:
        make_UJC_hash(args, df)

    # Add cell information to classification
    add_cell_data(args, df)

    # Run plotting script
    generate_report(args, df)


if __name__ == "__main__":
    main()
