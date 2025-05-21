#!/usr/bin/env python3
import subprocess, os, re, sys, glob
import argparse
import pandas as pd
import numpy as np
import hashlib
import pysam

sqanti3_src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../SQANTI3"))
if sqanti3_src_path not in sys.path:
    sys.path.insert(0, sqanti3_src_path)

from src.commands import run_command
from src.module_logging import qc_logger, update_logger

#!/usr/bin/env python3
# SQANTI_Single_Cell: Structural and Quality Annotation of transcripts and reads at the single cell level
# Author: Carlos Blanco

__author__  = "carlos.blanco@csic.es"
__version__ = '0.1.1'  # Python 3.11

utilitiesPath = os.path.abspath(os.path.join(os.path.dirname(__file__), "../utilities"))
sys.path.insert(0, utilitiesPath)

sqantiqcPath = os.path.abspath(os.path.join(os.path.dirname(__file__), "../SQANTI3"))


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
            f"python {sqantiqcPath}/sqanti3_qc.py --isoforms {input_file} --refGTF {args.refGTF} --refFasta {args.refFasta} "
            f"--min_ref_len {args.min_ref_len} --aligner_choice {args.aligner_choice} -w {args.window} "
            f"-t {args.mapping_cpus} -n {args.chunks} -d {args.out_dir}/{file_acc} -o {sampleID} "
            f"-s {args.sites} --ratio_TSS_metric {args.ratio_TSS_metric} --novel_gene_prefix {args.novel_gene_prefix} --report skip"
        )

        if is_fastq:
            cmd += " --fasta"
        if args.aligner_choice == "gmap":
            cmd += f" -x {args.gmap_index}"
        
        flag_args = [
            ('force_id_ignore', '--force_id_ignore'),
            ('skipORF', '--skipORF'),
            ('genename', '--genename'),
            ('saturation', '--saturation'),
            ('isoAnnotLite', '--isoAnnotLite'),
            ('isoform_hits', '--isoform_hits')
        ]

        for arg, flag in flag_args:
            if getattr(args, arg):
                cmd += f" {flag}"

        optional_files = [
            ("CAGE_peak", "--CAGE_peak"),
            ("polyA_motif_list", "--polyA_motif_list"),
            ("polyA_peak", "--polyA_peak"),
            ("phyloP_bed", "--phyloP_bed"),
            ("orf_input", "--orf_input"),
            ("expression", "--expression"),
            ("coverage", "--coverage"),
            ("gff3", "--gff3"),
            ("short_reads", "--short_reads"),
            ("SR_bam", "--SR_bam"),
            ]

        for arg_name, flag in optional_files:
            if getattr(args, arg_name):
                cmd += f" {flag} {getattr(args, arg_name)}"

        return cmd

    # Validate genome and annotation before processing files
    if not check_files_exist(args.refFasta, args.refGTF):
        print("[ERROR] Genome or annotation file is missing.", file=sys.stderr)
        sys.exit(-1)

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']

        # Create the logs directory for this file_acc
        file_acc_dir = os.path.join(args.out_dir, file_acc)
        logs_dir = os.path.join(file_acc_dir, "logs")
        os.makedirs(logs_dir, exist_ok=True)

        # Update the logger to write to the correct log file
        log_file = os.path.join(logs_dir, "sqanti3.log")
        update_logger(qc_logger, file_acc_dir, args.log_level)

        # Check for .gtf or .gff file
        gtf_pattern = os.path.join(args.input_dir, f"{file_acc}*.g*f")
        gtf_files = glob.glob(gtf_pattern)
        if gtf_files:
            gtf_file = gtf_files[0]
            if args.verbose:
                print(f'[INFO] Running SQANTI-sc qc for sample {gtf_file}', file=sys.stdout)
            
            cmd_sqanti = build_sqanti_command(gtf_file)

            if args.is_fusion:
                cmd_sqanti += " --is_fusion"
                if not args.skipORF:
                    if not args.orf_input:
                        raise ValueError("[ERROR] --orf_input must be provided when using --is_fusion with a GTF input file, unless --skipORF is specified.")
                    cmd_sqanti += f" --orf_input {args.orf_input}"

            print(cmd_sqanti, file=sys.stdout)
            run_command(cmd_sqanti, qc_logger, out_file=log_file, description="SQANTI3 failed to execute")
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
            run_command(cmd_sqanti, qc_logger, out_file=log_file, description="SQANTI3 failed to execute")
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
                run_command(cmd_tofastq, qc_logger, out_file=log_file, description="Samtools failed to execute")

                # Ensure FASTQ was generated
                if not os.path.isfile(tmp_fastq):
                    print(f"[ERROR] FASTQ file not generated for {file_acc}", file=sys.stderr)
                    continue

                # Run SQANTI3 in fasta mode
                cmd_sqanti = build_sqanti_command(tmp_fastq, is_fastq=True)
                print(cmd_sqanti, file=sys.stdout)
                run_command(cmd_sqanti, qc_logger, out_file=log_file, description="SQANTI3 failed to execute")

                # Cleanup temporary files
                try:
                    os.remove(tmp_fastq)
                except FileNotFoundError:
                    print(f"[WARNING] Failed to remove temporary file: {tmp_fastq}", file=sys.stderr)

                continue

        # No valid files found
        print(f"[ERROR] No valid .fastq, .gtf, or .bam files found for {file_acc} in {args.input_dir}", file=sys.stderr)


def make_UJC_hash(args, df):
    def format_chr(chr_value):
        return f"chr{chr_value}" if not str(chr_value).startswith("chr") else str(chr_value)

    def create_jxn_string(row):
        base_string = format_chr(row["chr"]) + "_" + row["strand"] + "_"
        
        if pd.isna(row["jxn_string"]):
            # Mono-exonic case
            return f"{base_string}monoexon_{row['associated_transcript']}"
        else:
            # Multi-exonic case
            junctions = row["jxn_string"].split("_")
            return f"{base_string}" + "_".join(junctions[2:-1])

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)

        print(f"**** Calculating UJCs for {file_acc}...", file=sys.stdout)

        introns_cmd = f"""gtftools -i {outputPathPrefix}tmp_introns.bed -c "$(cut -f 1 {outputPathPrefix}_corrected.gtf.cds.gff | sort | uniq | paste -sd ',' - | sed 's/chr//g')" {outputPathPrefix}_corrected.gtf.cds.gff"""
        ujc_cmd = f"""awk -F'\t' -v OFS="\t" '{{print $5,$1,$4,$2+1"_"$3}}' {outputPathPrefix}tmp_introns.bed | bedtools groupby -g 1 -c 2,3,4 -o distinct,distinct,collapse | sed 's/,/_/g' | awk -F'\t' -v OFS="\t" '{{print $1,$2"_"$3"_"$4}}' > {outputPathPrefix}tmp_UJC.txt"""

        if subprocess.check_call(introns_cmd, shell=True) != 0:
            print("ERROR running command: {0}\n Missing GTFTOOLS".format(introns_cmd), file=sys.stderr)
            sys.exit(-1)

        if os.path.exists(f"{outputPathPrefix}_corrected.gtf.ensembl"):
            os.remove(f"{outputPathPrefix}_corrected.gtf.ensembl")

        if subprocess.check_call(ujc_cmd, shell=True) != 0:
            print("ERROR running command: {0}\n Missing BEDTOOLS".format(introns_cmd), file=sys.stderr)
            sys.exit(-1)
        os.remove(f"{outputPathPrefix}tmp_introns.bed")

        # Add "chr" to the chromosome if missing
        with open(f"{outputPathPrefix}tmp_UJC.txt", 'r') as f:
            lines = f.readlines()

        with open(f"{outputPathPrefix}tmp_UJC.txt", 'w') as f:
            for line in lines:
                parts = line.strip().split('\t')
                chr_part = parts[1].split('_')[0]
                parts[1] = '_'.join([format_chr(chr_part)] + parts[1].split('_')[1:])
                f.write('\t'.join(parts) + '\n')

        classfile = f"{outputPathPrefix}_classification.txt"
        clas_df = pd.read_csv(classfile, sep="\t", usecols=[0, 1, 2, 7], dtype="str")
        clas_df.columns = ["isoform", "chr", "strand", "associated_transcript"]
        ujc_df = pd.read_csv(f"{outputPathPrefix}tmp_UJC.txt", sep="\t", names=["isoform", "jxn_string"], dtype="str")

        merged_df = pd.merge(clas_df, ujc_df, on="isoform", how="left")
        
        merged_df["jxn_string"] = merged_df.apply(create_jxn_string, axis=1)

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
    and adds them to classification.
    """

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)

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

        # 2. If not bam file, try to read from cell association file
        elif os.path.isfile(cell_association_file):
            try:
                cell_association_df = pd.read_csv(cell_association_file, sep='\t', dtype=str)
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

                # Fix boolean values
                classification_df['RTS_stage'] = classification_df['RTS_stage'].map({True: 'TRUE', False: 'FALSE'})
                classification_df['predicted_NMD'] = classification_df['predicted_NMD'].map({True: 'TRUE', False: 'FALSE'})

                # Replace empty strings with NaN, then NaN with 'NA'
                classification_df = classification_df.fillna('NA')

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
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        
        classification_file = f"{outputPathPrefix}_classification.txt"
        
        print(f"**** Generating SQANTI3 report for {file_acc}...", file=sys.stdout)

        if os.path.isfile(classification_file):
            try:
                ignore_flag  = " --ignore_cell_summary" if args.ignore_cell_summary else ""
                cmd = (
                    f"Rscript {utilitiesPath}/SQANTI-sc_reads.R "
                    f"{classification_file} "
                    f"{args.report} "
                    f"{outputPathPrefix}"
                    f"{ignore_flag}"
                )
                subprocess.run(cmd, shell=True, check=True)
                print(f"**** SQANTI3 report successfully generated for {file_acc}", file=sys.stdout)
            except subprocess.CalledProcessError as e:
                print(f"Error generating report for {classification_file}: {str(e)}", file=sys.stdout)
        else:
            print(f"Classification file for {file_acc} not found. Skipping report generation.", file=sys.stdout)


def main():
    global utilitiesPath
    global sqantiqcPath

    #arguments
    ap = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")
    # Required arguments
    apr = ap.add_argument_group("Required arguments")
    apr.add_argument('--refFasta', type=str, required=True, help='\t\tReference genome (Fasta format).')
    apr.add_argument('--refGTF', type=str, required=True, help='\t\tReference annotation file (GTF format).')
    apr.add_argument('-de', '--design', type=str, dest="inDESIGN", required=True, help='Path to design file, must have sampleID and file_acc column.')
    apr.add_argument('-m', '--mode', type=str, choices = ["isoforms", "reads"], required=True, help = '\t\tType of data to run SQANTI3 on (reads or isoforms)')

    # Optional arguments
    apsc = ap.add_argument_group("Optional arguments")
    apsc.add_argument('-d','--out_dir', type=str, help='\t\tDirectory for output sqanti_sc files. Default: Directory where the script was run.', default = "./", required=False)
    apsc.add_argument('-i', '--input_dir', type=str, default = './', help = '\t\tPath to directory where fastq files are stored. Or path to parent directory with children directories of SQANTI3 runs. Default: Directory where the script was run.')
    apsc.add_argument('--report', type=str, choices = ["pdf", "html", "both", "skip"], default = "pdf", help = "\t\tDefault: pdf")
    apsc.add_argument('-@', '--samtools_cpus', default=10, type=int, help='\t\tNumber of threads used during conversion from bam to fastq. Default: 10')    
    apsc.add_argument('--verbose', help = 'If verbose is run, it will print all steps, by default it is FALSE', action="store_true", default=False)
    apsc.add_argument('-f', '--factor', type=str, dest="inFACTOR" ,required=False, help='This is the column name that plots are to be faceted by. Default: None')
    apsc.add_argument('--skip_hash', dest="SKIPHASH", action='store_true', help='Skip the hashing step')
    apsc.add_argument('--ignore_cell_summary', action="store_true", default=False, help="\t\t Add this flag to not save the cell summary table generated during the report to save space. Do not add if running the cell filter module afterwards, as this table is used to inform the cell filtering.")
    apsc.add_argument('-cm', '--count_matrix', help='Cellxisoform count matrix')

    # SQANTI3 arguments
    # Customization and filtering args
    apc = ap.add_argument_group("SQANTI3 customization and filtering")
    apc.add_argument('--min_ref_len', type=int, default=0, help="\t\tMinimum reference transcript length. Default: 0 bp")
    apc.add_argument('--force_id_ignore', action="store_true", default=False, help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")
    apc.add_argument('--genename', action='store_true' ,help='Use gene_name tag from GTF to define genes. Default: gene_id used to define genes')
    apc.add_argument('--short_reads', type=str, help='File Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.')
    apc.add_argument('--SR_bam', type=str, help=' Directory or fofn file with the sorted bam files of Short Reads RNA-Seq mapped against the genome')
    apc.add_argument('--novel_gene_prefix', default=None, help='Prefix for novel isoforms (default: None)')

    # Aligner and mapping options
    apa = ap.add_argument_group("SQANTI3 aligner and mapping options")
    apa.add_argument('--aligner_choice', type=str, choices=['minimap2', "uLTRA", "gmap", "deSALT"], default='minimap2', help="\t\tDefault: minimap2")
    apa.add_argument('-x','--gmap_index', help='Path and prefix of the reference index created by gmap_build. Mandatory if using GMAP unless -g option is specified.')
    apa.add_argument('-s','--sites', type=str, default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)

    # ORF prediction
    apo = ap.add_argument_group("SQANTI3 ORF prediction")
    apo.add_argument('--skipORF', action="store_true", default=False, help="\t\t Skip ORF prediction to save time.")
    apo.add_argument('--orf_input', type=str, help="Input fasta to run ORF on. By default, ORF is run on genome-corrected fasta - this overrides it. If input is fusion (--is_fusion), this must be provided for ORF prediction.")

    # Functional annotation
    apf = ap.add_argument_group("SQANTI3 functional annotation")
    apf.add_argument('--CAGE_peak', type=str, help="FANTOM5 Cage Peak (BED format, optional)")
    apf.add_argument('--polyA_motif_list', type=str, help="Ranked list of polyA motifs (text, optional)")
    apf.add_argument('--polyA_peak', type=str, help="PolyA Peak (BED format, optional)")
    apf.add_argument('--phyloP_bed', type=str, help="PhyloP BED for conservation score (BED, optional)")

    # Output options
    apout = ap.add_argument_group("SQANTI3 output options")
    apout.add_argument('--saturation', action="store_true", default=False, help='Include saturation curves into report')
    apout.add_argument('--isoform_hits' , action='store_true', help=' Report all FSM/ISM isoform hits in a separate file')
    apout.add_argument('--ratio_TSS_metric' , choices=['max', 'mean', 'median', '3quartile'], default='max', help=' Define which statistic metric should be reported in the ratio_TSS column (default: %(default)s)')

    # Performance options
    app = ap.add_argument_group("SQANTI3 performance options")
    app.add_argument('-t', '--mapping_cpus', default=10, type=int, help='\t\tNumber of threads used during alignment by SQANTI3 aligners. Default: 10')
    app.add_argument('-n', '--chunks', default=10, type=int, help='\t\tNumber of chunks to split SQANTI3 analysis in for speed up. Default: 10')
    app.add_argument("-l","--log_level", default="INFO",choices=["ERROR","WARNING","INFO","DEBUG"],
                        help="Set the logging level %(default)s")

    # Optional arguments
    apm = ap.add_argument_group("SQANTI3 optional arguments")
    apm.add_argument('--is_fusion', action="store_true", help="Input are fusion isoforms, must supply GTF as input")
    apm.add_argument('-e','--expression', type=str, help='Expression matrix (supported: Kallisto tsv)')
    apm.add_argument('-c','--coverage', help='Junction coverage files (provide a single file, comma-delmited filenames, or a file pattern, ex: "mydir/*.junctions").')
    apm.add_argument('-w','--window', default=20, type=int, help='Size of the window in the genomic DNA screened for Adenine content downstream of TTS (default: %(default)s)')
    apm.add_argument('-v', '--version', help="Display program version number.", action='version', version='sqanti-sc '+str(__version__))
    apm.add_argument('--isoAnnotLite', action='store_true', help='Run isoAnnot Lite to output a tappAS-compatible gff3 file')
    apm.add_argument('--gff3' ,type=str, help='Precomputed tappAS species specific GFF3 file. It will serve as reference to transfer functional attributes')

    args = ap.parse_args()

    # Check and read design file
    try:
        df = fill_design_table(args)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Run SQANTI3 
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
