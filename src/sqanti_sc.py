#!/usr/bin/env python3
import subprocess
import os
import re
import sys
import glob
import argparse
import pandas as pd
import hashlib
import pysam

sqanti3_src_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../SQANTI3")
)
if sqanti3_src_path not in sys.path:
    sys.path.insert(0, sqanti3_src_path)

from src.commands import run_command
from src.module_logging import qc_logger, update_logger

#!/usr/bin/env python3
# SQANTI_Single_Cell: Structural and Quality Annotation of transcripts and
# reads at the single cell level
# Author: Carlos Blanco

__author__ = "carlos.blanco@csic.es"
__version__ = '0.2.0'  # Python 3.11

utilitiesPath = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../utilities")
)
sys.path.insert(0, utilitiesPath)

sqantiqcPath = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../SQANTI3")
)


def fill_design_table(args):
    # Read the design CSV file into a DataFrame
    df = pd.read_csv(args.inDESIGN, sep=",")

    # Check for required columns
    required_columns = {'sampleID', 'file_acc'}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(f"ERROR: Missing required columns: {missing_columns}")

    # Create the new columns based on the existing ones
    df['classification_file'] = df.apply(
        lambda row: os.path.join(
            args.out_dir, row['file_acc'], f"{row['sampleID']}_classification.txt"
        ),
        axis=1
    )
    df['junction_file'] = df.apply(
        lambda row: os.path.join(
            args.out_dir, row['file_acc'], f"{row['sampleID']}_junctions.txt"
        ),
        axis=1
    )

    if 'cell_association_file' not in df.columns:
        df['cell_association_file'] = ''

    for index, row in df.iterrows():
        if pd.isna(row['cell_association_file']) or not row['cell_association_file']:
            file_acc = row['file_acc']
            bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
            bam_files = glob.glob(bam_pattern)
            assoc_file_path = os.path.join(
                args.input_dir, f"{file_acc}_cell_association.txt"
            )

            if bam_files:
                df.at[index, 'cell_association_file'] = os.path.abspath(bam_files[0])
            elif os.path.isfile(assoc_file_path):
                df.at[index, 'cell_association_file'] = os.path.abspath(
                    assoc_file_path
                )

    # Write the DataFrame back to the CSV file
    df.to_csv(args.inDESIGN, sep=',', index=False)

    return df


def get_files_runSQANTI3(args, df):
    def check_files_exist(*files):
        """Check if each file exists, return False if missing instead of exiting."""
        missing_files = [f for f in files if not os.path.isfile(f)]
        if missing_files:
            print(f'[ERROR] Missing file(s): {", ".join(missing_files)}',
                  file=sys.stderr)
            return False
        return True

    def build_sqanti_command(input_file, file_acc, sampleID, is_fastq=False):
        """Build the SQANTI command to run."""
        cmd_parts = [
            "python", f"{sqantiqcPath}/sqanti3_qc.py",
            "--isoforms", input_file,
            "--refGTF", args.refGTF,
            "--refFasta", args.refFasta,
            "--min_ref_len", str(args.min_ref_len),
            "--aligner_choice", args.aligner_choice,
            "-w", str(args.window),
            "-t", str(args.mapping_cpus),
            "-n", str(args.chunks),
            "-d", os.path.join(args.out_dir, file_acc),
            "-o", sampleID,
            "-s", args.sites,
            "--ratio_TSS_metric", args.ratio_TSS_metric,
            "--report", "skip"
        ]

        if args.novel_gene_prefix:
            cmd_parts.extend(["--novel_gene_prefix", args.novel_gene_prefix])

        if is_fastq:
            cmd_parts.append("--fasta")

        if args.aligner_choice == "gmap" and args.gmap_index:
            cmd_parts.extend(["-x", args.gmap_index])

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
                cmd_parts.append(flag)

        optional_files = [
            ("CAGE_peak", "--CAGE_peak"),
            ("polyA_motif_list", "--polyA_motif_list"),
            ("polyA_peak", "--polyA_peak"),
            ("phyloP_bed", "--phyloP_bed"),
            ("orf_input", "--orf_input"),
            ("expression", "--expression"),
            ("coverage", "--coverage"),
            ("fl_count", "--fl_count"),
            ("gff3", "--gff3"),
            ("short_reads", "--short_reads"),
            ("SR_bam", "--SR_bam"),
        ]
        for arg_name, flag in optional_files:
            if getattr(args, arg_name):
                cmd_parts.extend([flag, getattr(args, arg_name)])

        return " ".join(cmd_parts)

    # Validate genome and annotation before processing files
    if not check_files_exist(args.refFasta, args.refGTF):
        print("[ERROR] Genome or annotation file is missing.", file=sys.stderr)
        sys.exit(-1)

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        file_acc_dir = os.path.join(args.out_dir, file_acc)
        logs_dir = os.path.join(file_acc_dir, "logs")
        os.makedirs(logs_dir, exist_ok=True)
        log_file = os.path.join(logs_dir, "sqanti3.log")
        update_logger(qc_logger, file_acc_dir, args.log_level)

        gtf_pattern = os.path.join(args.input_dir, f"{file_acc}*.g*f")
        gtf_files = glob.glob(gtf_pattern)
        if gtf_files:
            gtf_file = gtf_files[0]
            if args.verbose:
                print(f'[INFO] Running SQANTI-sc qc for {gtf_file}', file=sys.stdout)
            cmd = build_sqanti_command(gtf_file, file_acc, sampleID)
            if args.is_fusion:
                cmd += " --is_fusion"
                if not args.skipORF:
                    if not args.orf_input:
                        raise ValueError(
                            "[ERROR] --orf_input must be provided for fusion "
                            "unless --skipORF is specified."
                        )
                    cmd += f" --orf_input {args.orf_input}"
            print(cmd, file=sys.stdout)
            run_command(cmd, qc_logger, out_file=log_file,
                        description="SQANTI3 failed to execute")
            continue

        fastq_pattern = os.path.join(args.input_dir, f"{file_acc}*.fast*")
        fastq_files = glob.glob(fastq_pattern)
        if fastq_files:
            fastq_file = fastq_files[0]
            if args.verbose:
                print(f'[INFO] Running SQANTI-sc qc for {fastq_file}', file=sys.stdout)
            cmd = build_sqanti_command(fastq_file, file_acc, sampleID, is_fastq=True)
            print(cmd, file=sys.stdout)
            run_command(cmd, qc_logger, out_file=log_file,
                        description="SQANTI3 failed to execute")
            continue

        if args.mode == "reads":
            bam_pattern = os.path.join(args.input_dir, f"{file_acc}*.bam")
            bam_files = glob.glob(bam_pattern)
            if bam_files:
                bam_file = bam_files[0]
                if args.verbose:
                    print(f'[INFO] Running SQANTI-sc qc for {bam_file}', file=sys.stdout)

                tmp_fastq = f"{args.input_dir}/{file_acc}_tmp.fastq"
                cmd_tofastq = (
                    f"samtools fastq -@ {args.samtools_cpus} {bam_file} "
                    f"-n > {tmp_fastq}"
                )
                run_command(cmd_tofastq, qc_logger, out_file=log_file,
                            description="Samtools failed to execute")

                if not os.path.isfile(tmp_fastq):
                    print(f"[ERROR] FASTQ not generated for {file_acc}",
                          file=sys.stderr)
                    continue

                cmd = build_sqanti_command(
                    tmp_fastq, file_acc, sampleID, is_fastq=True
                )
                print(cmd, file=sys.stdout)
                run_command(cmd, qc_logger, out_file=log_file,
                            description="SQANTI3 failed to execute")

                try:
                    os.remove(tmp_fastq)
                except FileNotFoundError:
                    print(f"[WARNING] Failed to remove temp file: {tmp_fastq}",
                          file=sys.stderr)
                continue

        print(f"[ERROR] No valid input files found for {file_acc}",
              file=sys.stderr)


def make_UJC_hash(args, df):
    def format_chr(chr_val):
        return f"chr{chr_val}" if not str(chr_val).startswith("chr") else str(chr_val)

    def create_jxn_string(row):
        base = f"{format_chr(row['chr'])}_{row['strand']}_"
        if pd.isna(row["jxn_string"]):
            return f"{base}monoexon_{row['associated_transcript']}"
        junctions = row["jxn_string"].split("_")
        return f"{base}" + "_".join(junctions[2:-1])

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        print(f"**** Calculating UJCs for {file_acc}...", file=sys.stdout)

        # Ensure the corrected GTF contains gene_id attributes on every exon/CDS line so that
        # downstream `gtftools` does not fail with `IndexError: list index out of range`.
        # `gffread` will rewrite the file adding the missing attributes. We write to a
        # temporary file and then move it back so the final filename remains unchanged.
        gffread_tmp = f"{outputPathPrefix}_corrected.gtf.gffread_tmp"
        gffread_cmd = f"gffread {outputPathPrefix}_corrected.gtf -T -o {gffread_tmp} && mv {gffread_tmp} {outputPathPrefix}_corrected.gtf"

        try:
            subprocess.check_call(gffread_cmd, shell=True)
        except subprocess.CalledProcessError:
            print(f"ERROR running command: {gffread_cmd}\n Missing or failed gffread", file=sys.stderr)
            sys.exit(-1)

        introns_cmd = f"""gtftools -i {outputPathPrefix}tmp_introns.bed -c "$(cut -f 1 {outputPathPrefix}_corrected.gtf | sort | uniq | paste -sd ',' - | sed 's/chr//g')" {outputPathPrefix}_corrected.gtf"""
        ujc_cmd = f"""awk -F'\t' -v OFS="\t" '{{print $5,"chr"$1,$4,$2+1"_"$3}}' {outputPathPrefix}tmp_introns.bed | bedtools groupby -g 1 -c 2,3,4 -o distinct,distinct,collapse | sed 's/,/_/g' | awk -F'\t' -v OFS="\t" '{{print $1,$2"_"$3"_"$4}}' > {outputPathPrefix}tmp_UJC.txt"""

        if subprocess.check_call(introns_cmd, shell=True) != 0:
            print(f"ERROR running command: {introns_cmd}\nMissing GTFTOOLS",
                  file=sys.stderr)
            sys.exit(-1)

        if os.path.exists(f"{outputPathPrefix}_corrected.gtf.ensembl"):
            os.remove(f"{outputPathPrefix}_corrected.gtf.ensembl")

        if subprocess.check_call(ujc_cmd, shell=True) != 0:
            print(f"ERROR running command: {ujc_cmd}\nMissing BEDTOOLS",
                  file=sys.stderr)
            sys.exit(-1)
        os.remove(f"{outputPathPrefix}tmp_introns.bed")

        with open(f"{outputPathPrefix}tmp_UJC.txt", 'r') as f_in, \
             open(f"{outputPathPrefix}tmp_UJC_chr.txt", 'w') as f_out:
            for line in f_in:
                parts = line.strip().split('\t')
                chr_part = parts[1].split('_')[0]
                parts[1] = '_'.join([format_chr(chr_part)] + parts[1].split('_')[1:])
                f_out.write('\t'.join(parts) + '\n')
        os.rename(f"{outputPathPrefix}tmp_UJC_chr.txt", f"{outputPathPrefix}tmp_UJC.txt")

        classfile = f"{outputPathPrefix}_classification.txt"
        clas_df = pd.read_csv(
            classfile, sep="\t", usecols=[0, 1, 2, 7], dtype="str"
        )
        clas_df.columns = ["isoform", "chr", "strand", "associated_transcript"]
        ujc_df = pd.read_csv(
            f"{outputPathPrefix}tmp_UJC.txt", sep="\t",
            names=["isoform", "jxn_string"], dtype="str"
        )
        merged_df = pd.merge(clas_df, ujc_df, on="isoform", how="left")
        merged_df["jxn_string"] = merged_df.apply(create_jxn_string, axis=1)
        merged_df['jxnHash'] = merged_df['jxn_string'].apply(
            lambda x: hashlib.sha256(x.encode('utf-8')).hexdigest()
        )

        original_class_df = pd.read_csv(classfile, sep='\t', dtype=str, low_memory=False)
        new_cols_df = merged_df[['isoform', 'jxn_string', 'jxnHash']]

        # Merge new columns with the original classification file dataframe.
        final_df = pd.merge(original_class_df, new_cols_df, on='isoform', how='left')

        # Save the result to the temporary file that add_cell_data expects.
        final_df.to_csv(f"{outputPathPrefix}_classification_tmp.txt", index=False, sep='\t')

        print(f"**** UJCs successfully calculated and added to {file_acc} classification", file=sys.stdout)

        os.remove(f"{outputPathPrefix}tmp_UJC.txt")


def add_cell_data(args, df):
    """
    Extract isoform, UMI, and cell barcode information from BAM file or a
    cell association file and adds them to classification.
    """
    def _parse_isoform_tsv_for_cell_association(tsv_file):
        """Parse a TSV file for isoforms and return a dict with cell association data."""
        cell_dict = {}
        try:
            df_assoc = pd.read_csv(tsv_file, sep='\t', dtype=str, comment='#')
            required = ['pbid', 'cell_barcodes']
            if not all(col in df_assoc.columns for col in required):
                print(f"[ERROR] {tsv_file} missing required columns: 'pbid' and 'cell_barcodes'.",
                      file=sys.stderr)
                return cell_dict
            for _, row in df_assoc.iterrows():
                pbid, cbs = row['pbid'], row['cell_barcodes']
                cell_dict[pbid] = {"CB": cbs}
        except Exception as e:
            print(f"[ERROR] Failed to parse {tsv_file}: {e}", file=sys.stderr)
        return cell_dict

    def _parse_bam_for_cell_association(bam_file):
        """Parse a BAM file and return a dict with cell association data."""
        cell_dict = {}
        print(f"Extracting CBs and UMIs from {bam_file} ...", file=sys.stderr)
        if not os.path.isfile(bam_file):
            print(f"[ERROR] BAM file {bam_file} not accessible. Skipping...",
                  file=sys.stderr)
            return cell_dict
        try:
            with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
                for read in bam:
                    if read.has_tag("XM") and read.has_tag("CB"):
                        umi = read.get_tag("XM")
                        cb = read.get_tag("CB")
                        isoform = read.query_name
                        cell_dict[isoform] = {"UMI": umi, "CB": cb}
        except Exception as e:
            print(f"[ERROR] Failed to parse BAM file {bam_file}: {e}",
                  file=sys.stderr)
        return cell_dict

    def _parse_tsv_for_cell_association(tsv_file):
        """Parse a TSV file and return a dict with cell association data."""
        cell_dict = {}
        try:
            df_assoc = pd.read_csv(tsv_file, sep='\t', dtype=str, comment='#')
            required = ['id', 'cell_barcode']
            if not all(col in df_assoc.columns for col in required):
                print(f"[ERROR] {tsv_file} missing required columns.",
                      file=sys.stderr)
                return cell_dict
            for _, row in df_assoc.iterrows():
                iso, cb = row['id'], row['cell_barcode']
                umi = row.get('umi')
                cell_dict[iso] = {"CB": cb}
                if umi:
                    cell_dict[iso]["UMI"] = umi
        except Exception as e:
            print(f"[ERROR] Failed to parse {tsv_file}: {e}", file=sys.stderr)
        return cell_dict

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        assoc_path = row.get('cell_association_file')

        cell_data = {}
        if pd.notna(assoc_path) and assoc_path:
            print(f"INFO: Using cell association file: {assoc_path}",
                  file=sys.stdout)
            if args.mode == "isoforms":
                    cell_data = _parse_isoform_tsv_for_cell_association(assoc_path)
            elif args.mode == "reads":
                if assoc_path.endswith(".bam"):
                    cell_data = _parse_bam_for_cell_association(assoc_path)
                else:
                    cell_data = _parse_tsv_for_cell_association(assoc_path)

        tmp_class = f"{outputPathPrefix}_classification_tmp.txt"
        final_class = f"{outputPathPrefix}_classification.txt"

        if not cell_data:
            print(f"[INFO] No cell data for {file_acc}. Skipping.",
                  file=sys.stdout)
            if os.path.exists(tmp_class):
                os.rename(tmp_class, final_class)
            continue

        if os.path.isfile(tmp_class):
            try:
                class_df = pd.read_csv(tmp_class, sep='\t',
                                       low_memory=False, dtype=str)

                if args.mode == "reads":
                    class_df['UMI'] = class_df['isoform'].map(
                        lambda x: cell_data.get(x, {}).get("UMI")
                    )
                    class_df['CB'] = class_df['isoform'].map(
                        lambda x: cell_data.get(x, {}).get("CB")
                    )

                    for i, r in class_df.iterrows():
                        match = re.match(r"(.+)_dup\d+$", r['isoform'])
                        if match:
                            primary = match.group(1)
                            if primary in cell_data:
                                class_df.at[i, 'UMI'] = cell_data[primary].get('UMI')
                                class_df.at[i, 'CB'] = cell_data[primary].get('CB')
                elif args.mode == "isoforms":
                    class_df['CB'] = class_df['isoform'].map(
                        lambda x: cell_data.get(x, {}).get("CB")
                    )

                for col in class_df.select_dtypes(include=['bool']).columns:
                    class_df[col] = class_df[col].map({True: 'TRUE', False: 'FALSE'})

                class_df = class_df.fillna('NA')
                class_df.to_csv(final_class, index=False, sep="\t")
            except Exception as e:
                print(f"[ERROR] Could not merge cell data: {e}",
                      file=sys.stderr)
        else:
            print(f"[INFO] Classification file for {file_acc} not found.",
                  file=sys.stdout)

        junctions_path = f"{outputPathPrefix}_junctions.txt"
        if os.path.isfile(junctions_path):
            try:
                junc_df = pd.read_csv(junctions_path, sep='\t', dtype=str)
                if 'isoform' in junc_df.columns:
                    junc_df['CB'] = junc_df['isoform'].map(
                        lambda x: cell_data.get(x, {}).get("CB")
                    )
                    junc_df.to_csv(junctions_path, sep='\t', index=False)
                    print(f"**** CB column added to {file_acc} junctions.",
                          file=sys.stdout)
                else:
                    print(f"[ERROR] 'isoform' col not in {junctions_path}",
                          file=sys.stderr)
            except Exception as e:
                print(f"[ERROR] Could not add CB to junctions: {e}",
                      file=sys.stderr)
        else:
            print(f"[INFO] Junctions file for {file_acc} not found.",
                  file=sys.stdout)

        print(f"**** Cell data added to {file_acc}", file=sys.stdout)
        if os.path.exists(tmp_class):
            os.remove(tmp_class)


def generate_report(args, df):
    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{outputPathPrefix}_classification.txt"
        junc_file = f"{outputPathPrefix}_junctions.txt"
        print(f"**** Generating SQANTI3 report for {file_acc}...",
              file=sys.stdout)

        if os.path.isfile(class_file):
            try:
                flags = []
                if args.ignore_cell_summary:
                    flags.append("--ignore_cell_summary")
                if args.skipORF:
                    flags.append("--skipORF")
                if args.CAGE_peak:
                    flags.append("--CAGE_peak")
                if args.polyA_motif_list:
                    flags.append("--polyA_motif_list")

                cmd = (
                    f"Rscript {utilitiesPath}/SQANTI-sc_reads.R "
                    f"{class_file} {junc_file} {args.report} {outputPathPrefix} "
                    f"{args.mode} {' '.join(flags)}"
                )
                subprocess.run(cmd, shell=True, check=True)
                print(f"**** SQANTI3 report generated for {file_acc}",
                      file=sys.stdout)
            except subprocess.CalledProcessError as e:
                print(f"Error generating report for {class_file}: {e}",
                      file=sys.stdout)
        else:
            print(f"Classification file for {file_acc} not found.",
                  file=sys.stdout)


def main():
    global utilitiesPath, sqantiqcPath

    ap = argparse.ArgumentParser(
        description="SQANTI Single-Cell: Structural and Quality Annotation"
    )
    apr = ap.add_argument_group("Required arguments")
    apr.add_argument('--refFasta', required=True,
                     help='Reference genome (Fasta).')
    apr.add_argument('--refGTF', required=True,
                     help='Reference annotation (GTF).')
    apr.add_argument('-de', '--design', dest="inDESIGN", required=True,
                     help='Design file with sampleID and file_acc.')
    apr.add_argument('-m', '--mode', choices=["isoforms", "reads"],
                     required=True, help='Input data type.')

    apsc = ap.add_argument_group("Optional arguments")
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
    apsc.add_argument('-f', '--factor', dest="inFACTOR",
                      help='Column name for faceting plots.')
    apsc.add_argument('--skip_hash', dest="SKIPHASH", action='store_true',
                      help='Skip UJC hashing step.')
    apsc.add_argument('--ignore_cell_summary', action="store_true",
                      default=False,
                      help="Don't save cell summary table in report.")
    apsc.add_argument('-cm', '--count_matrix',
                      help='Cellxisoform count matrix.')

    apc = ap.add_argument_group("SQANTI3 customization and filtering")
    apc.add_argument('--min_ref_len', type=int, default=0,
                     help="Min reference transcript length. Default: 0.")
    apc.add_argument('--force_id_ignore', action="store_true", default=False,
                     help="Allow non-PacBio transcript IDs.")
    apc.add_argument('--genename', action='store_true',
                     help='Use gene_name tag instead of gene_id.')
    apc.add_argument('--short_reads',
                     help='FOFN of short-read RNA-Seq files (FASTA/FASTQ).')
    apc.add_argument('--SR_bam',
                     help='Directory or FOFN of short-read BAM files.')
    apc.add_argument('--novel_gene_prefix',
                     help='Prefix for novel gene IDs.')

    apa = ap.add_argument_group("SQANTI3 aligner and mapping options")
    apa.add_argument('--aligner_choice', default='minimap2',
                     choices=['minimap2', "uLTRA", "gmap", "deSALT"],
                     help="Aligner to use. Default: minimap2.")
    apa.add_argument('-x', '--gmap_index',
                     help='Path to gmap_build index prefix.')
    apa.add_argument('-s', '--sites', default="ATAC,GCAG,GTAG",
                     help='Canonical splice sites. Default: ATAC,GCAG,GTAG.')

    apo = ap.add_argument_group("SQANTI3 ORF prediction")
    apo.add_argument('--skipORF', action="store_true", default=False,
                     help="Skip ORF prediction.")
    apo.add_argument('--orf_input',
                     help="Input fasta for ORF prediction.")

    apf = ap.add_argument_group("SQANTI3 functional annotation")
    apf.add_argument('--CAGE_peak', help="FANTOM5 CAGE Peak (BED).")
    apf.add_argument('--polyA_motif_list',
                     help="Ranked list of polyA motifs (text).")
    apf.add_argument('--polyA_peak', help="PolyA Peak (BED).")
    apf.add_argument('--phyloP_bed', help="PhyloP conservation scores (BED).")

    apout = ap.add_argument_group("SQANTI3 output options")
    apout.add_argument('--saturation', action="store_true", default=False,
                       help='Include saturation curves in report.')
    apout.add_argument('--isoform_hits', action='store_true',
                       help='Report all FSM/ISM hits in a separate file.')
    apout.add_argument('--ratio_TSS_metric', default='max',
                       choices=['max', 'mean', 'median', '3quartile'],
                       help='Metric for ratio_TSS column. Default: max.')

    app = ap.add_argument_group("SQANTI3 performance options")
    app.add_argument('-t', '--mapping_cpus', default=10, type=int,
                     help='Threads for alignment. Default: 10.')
    app.add_argument('-n', '--chunks', default=10, type=int,
                     help='Number of chunks for analysis. Default: 10.')
    app.add_argument("-l", "--log_level", default="INFO",
                     choices=["ERROR", "WARNING", "INFO", "DEBUG"],
                     help="Set logging level. Default: INFO.")

    apm = ap.add_argument_group("SQANTI3 optional arguments")
    apm.add_argument('--is_fusion', action="store_true",
                     help="Input are fusion isoforms (GTF required).")
    apm.add_argument('-e', '--expression',
                     help='Expression matrix (e.g., Kallisto tsv).')
    apm.add_argument('-c', '--coverage',
                     help='Junction coverage files (single, comma-separated, or pattern).')
    apm.add_argument('-w', '--window', default=20, type=int,
                     help='Window size for Adenine content screen. Default: 20.')
    apm.add_argument('-fl', '--fl_count',
                     help='Full-length PacBio abundance file')
    apm.add_argument('-v', '--version', action='version',
                     version='sqanti-sc ' + str(__version__))
    apm.add_argument('--isoAnnotLite', action='store_true',
                     help='Run isoAnnot Lite for tappAS compatibility.')
    apm.add_argument('--gff3',
                     help='Precomputed tappAS GFF3 for functional transfer.')

    args = ap.parse_args()

    try:
        df = fill_design_table(args)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    get_files_runSQANTI3(args, df)

    if not args.SKIPHASH:
        make_UJC_hash(args, df)

    add_cell_data(args, df)

    generate_report(args, df)


if __name__ == "__main__":
    main()
