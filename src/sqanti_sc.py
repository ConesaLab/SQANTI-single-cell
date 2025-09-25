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
import numpy as np
import warnings
from pandas.errors import PerformanceWarning

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
__version__ = '0.2.1'  # Python 3.11

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
        update_logger(qc_logger, file_acc_dir, "qc", args.log_level)

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
        return f"{base}" + "_".join(junctions[2:])

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        print(f"**** Calculating UJCs for {file_acc}...", file=sys.stdout)

        # Ensure the corrected GTF contains gene_id attributes on every exon/CDS line so that
        # downstream `gtftools` does not fail.
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


def calculate_metrics_per_cell(args, df):
    """
    Python implementation of the R calculate_metrics_per_cell.
    Writes <out_dir>/<file_acc>/<sampleID>_SQANTI_cell_summary.txt.gz
    and returns nothing. Maintains feature parity for plotting.
    """
    structural_categories = [
        'full-splice_match', 'incomplete-splice_match', 'novel_in_catalog',
        'novel_not_in_catalog', 'genic', 'antisense', 'fusion', 'intergenic',
        'genic_intron'
    ]
    cat_to_root = {
        'full-splice_match': 'FULLSPLICEMATCH',
        'incomplete-splice_match': 'INCOMPLETESPLICEMATCH',
        'novel_in_catalog': 'NOVELINCATALOG',
        'novel_not_in_catalog': 'NOVELNOTINCATALOG',
        'genic': 'GENIC', 'antisense': 'ANTISENSE', 'fusion': 'FUSION',
        'intergenic': 'INTERGENIC', 'genic_intron': 'GENICINTRON'
    }
    # Tag names used for the final output column nomenclature
    cat_to_tag = {
        'full-splice_match': 'FSM',
        'incomplete-splice_match': 'ISM',
        'novel_in_catalog': 'NIC',
        'novel_not_in_catalog': 'NNC',
        'genic': 'Genic',
        'antisense': 'Antisense',
        'fusion': 'Fusion',
        'intergenic': 'Intergenic',
        'genic_intron': 'Genic_intron'
    }
    # Common abbr<->category pairs reused across multiple loops
    abbr_pairs = [
        ('FSM', 'full-splice_match'),
        ('ISM', 'incomplete-splice_match'),
        ('NIC', 'novel_in_catalog'),
        ('NNC', 'novel_not_in_catalog')
    ]

    def final_count_name(cat):
        tag = cat_to_tag[cat]
        return 'Genic_Genomic' if cat == 'genic' else tag

    def final_count_name(cat):
        tag = cat_to_tag[cat]
        return 'Genic_Genomic' if cat == 'genic' else tag

    def safe_prop(numer, denom):
        numer = numer.astype(float)
        denom = denom.astype(float)
        with np.errstate(divide='ignore', invalid='ignore'):
            arr = np.where(denom > 0, (numer / denom) * 100.0, 0.0)
        return pd.Series(arr, index=numer.index)

    # Mute pandas PerformanceWarning in this function; we intentionally build
    # the result by adding many columns for clarity and parity with R output.
    warnings.simplefilter('ignore', PerformanceWarning)

    for _, r in df.iterrows():
        file_acc = r['file_acc']
        sampleID = r['sampleID']
        prefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{prefix}_classification.txt"
        junc_file = f"{prefix}_junctions.txt"
        out_summary = f"{prefix}_SQANTI_cell_summary.txt.gz"

        try:
            cls = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False)
        except Exception:
            continue
        try:
            junc = pd.read_csv(junc_file, sep='\t', dtype=str, low_memory=False)
        except Exception:
            junc = pd.DataFrame()

        # Ensure required cols
        for col in ['CB','isoform','associated_gene','structural_category','exons','length','ref_length',
                    'all_canonical','subcategory','chrom','UMI','jxn_string','associated_transcript',
                    'RTS_stage','predicted_NMD','within_CAGE_peak','polyA_motif_found','perc_A_downstream_TTS','diff_to_gene_TSS','coding']:
            if col not in cls.columns:
                cls[col] = np.nan
        # Types
        for c in ['exons','length','ref_length','perc_A_downstream_TTS','diff_to_gene_TSS']:
            cls[c] = pd.to_numeric(cls[c], errors='coerce')

        if args.mode == 'isoforms':
            cls['CB'] = cls['CB'].fillna('')
            cls = cls.assign(CB=cls['CB'].str.split(',')).explode('CB')
            cls['CB'] = cls['CB'].fillna('')

        cls_valid = cls[(cls['CB'].notna()) & (cls['CB'] != '')].copy()

        # Base counts
        total_reads = cls_valid.groupby('CB').size().rename('total_reads')
        if args.mode == 'isoforms':
            total_umi = pd.Series(0, index=total_reads.index, name='total_UMI')
        else:
            total_umi = cls_valid.groupby('CB')['UMI'].nunique().rename('total_UMI')
        reads_no_mono = cls_valid[cls_valid['exons'] != 1].groupby('CB').size().rename('total_reads_no_monoexon')
        summary = pd.DataFrame(total_reads).join(total_umi, how='outer').join(reads_no_mono, how='outer').fillna(0)

        # Category counts/props (use final names directly)
        cat_counts = cls_valid.groupby(['CB','structural_category']).size().unstack(fill_value=0)
        for c in structural_categories:
            if c not in cat_counts.columns:
                cat_counts[c] = 0
        cat_counts = cat_counts[structural_categories]
        cat_counts.columns = [final_count_name(c) for c in structural_categories]
        summary = summary.join(cat_counts, how='left').fillna(0)
        for c in structural_categories:
            tag = final_count_name(c)
            summary[f"{tag}_prop"] = safe_prop(summary[tag], summary['total_reads'])

        # Genes/models
        summary['Genes_in_cell'] = cls_valid.groupby('CB')['associated_gene'].nunique().reindex(summary.index, fill_value=0)
        summary['UJCs_in_cell'] = cls_valid[cls_valid['exons'] > 1].groupby('CB')['jxn_string'].nunique().reindex(summary.index, fill_value=0)

        # MT
        mt = cls_valid[cls_valid['chrom'] == 'MT'].groupby('CB').size().rename('MT_reads_count')
        summary = summary.join(mt, how='left').fillna({'MT_reads_count': 0})
        summary['MT_perc'] = safe_prop(summary['MT_reads_count'], summary['total_reads'])

        # Annotated/novel genes
        anno = (~cls_valid['associated_gene'].fillna('').str.startswith('novel'))
        summary['Annotated_genes'] = cls_valid[anno].groupby('CB')['associated_gene'].nunique().reindex(summary.index, fill_value=0)
        summary['Novel_genes'] = cls_valid[~anno].groupby('CB')['associated_gene'].nunique().reindex(summary.index, fill_value=0)

        # Junctions per type and props
        if not junc.empty:
            if 'CB' not in junc.columns or (junc['CB'].fillna('') == '').all():
                iso_to_cb = cls_valid[['isoform','CB']].dropna().drop_duplicates()
                junc = pd.merge(junc, iso_to_cb, on='isoform', how='left')
            jv = junc[(junc['CB'].notna()) & (junc['CB'] != '')].copy()
            if not jv.empty:
                jv['junction_type'] = jv['junction_category'].astype(str) + '_' + jv['canonical'].astype(str)
                counts = jv.groupby(['CB','junction_type']).size().unstack(fill_value=0)
                for tp in ['known_canonical','known_non_canonical','novel_canonical','novel_non_canonical']:
                    if tp not in counts.columns:
                        counts[tp] = 0
                counts['total_junctions'] = counts.sum(axis=1)
                counts = counts.rename(columns={
                    'known_canonical':'Known_canonical_junctions',
                    'known_non_canonical':'Known_non_canonical_junctions',
                    'novel_canonical':'Novel_canonical_junctions',
                    'novel_non_canonical':'Novel_non_canonical_junctions'
                })
                for src, dst in [
                    ('Known_canonical_junctions','Known_canonical_junctions_prop'),
                    ('Known_non_canonical_junctions','Known_non_canonical_junctions_prop'),
                    ('Novel_canonical_junctions','Novel_canonical_junctions_prop'),
                    ('Novel_non_canonical_junctions','Novel_non_canonical_junctions_prop')]:
                    counts[dst] = safe_prop(counts[src].reindex(counts.index, fill_value=0), counts['total_junctions'])
                summary = summary.join(counts, how='left').fillna(0)
            else:
                summary[['Known_canonical_junctions','Known_non_canonical_junctions','Novel_canonical_junctions','Novel_non_canonical_junctions','total_junctions',
                         'Known_canonical_junctions_prop','Known_non_canonical_junctions_prop','Novel_canonical_junctions_prop','Novel_non_canonical_junctions_prop']] = 0
        else:
            summary[['Known_canonical_junctions','Known_non_canonical_junctions','Novel_canonical_junctions','Novel_non_canonical_junctions','total_junctions',
                     'Known_canonical_junctions_prop','Known_non_canonical_junctions_prop','Novel_canonical_junctions_prop','Novel_non_canonical_junctions_prop']] = 0

        # Subcategories (proportions within category) with final names directly
        sublevels = {
            'full-splice_match': ['alternative_3end','alternative_3end5end','alternative_5end','reference_match','mono-exon'],
            'incomplete-splice_match': ['3prime_fragment','internal_fragment','5prime_fragment','intron_retention','mono-exon'],
            'novel_in_catalog': ['combination_of_known_junctions','combination_of_known_splicesites','intron_retention','mono-exon_by_intron_retention','mono-exon'],
            'novel_not_in_catalog': ['at_least_one_novel_splicesite','intron_retention'],
            'genic': ['mono-exon','multi-exon'], 'antisense': ['mono-exon','multi-exon'],
            'fusion': ['intron_retention','multi-exon'], 'intergenic': ['mono-exon','multi-exon'],
            'genic_intron': ['mono-exon','multi-exon']
        }
        def _normalize_sub_lv(lv: str) -> str:
            # keep raw tokens; make column-safe: replace '-' with '_'
            return lv.replace('-', '_')
        def subkey(cat, lv):
            tag = cat_to_tag[cat]
            return f"{tag}_{_normalize_sub_lv(lv)}_prop"
        for cat in structural_categories:
            denom = summary[final_count_name(cat)]
            sub = cls_valid[cls_valid['structural_category'] == cat]
            tbl = sub.groupby(['CB','subcategory']).size().unstack(fill_value=0) if not sub.empty else pd.DataFrame()
            for lv in sublevels[cat]:
                numer = tbl.get(lv, pd.Series(0, index=summary.index)).reindex(summary.index, fill_value=0)
                summary[subkey(cat, lv.replace('-', '_'))] = safe_prop(numer, denom).fillna(0)

        # Gene read count bins per CB (annotated vs novel)
        gene_counts = cls_valid.groupby(['CB','associated_gene']).size().rename('read_count').reset_index()
        gene_counts['gene_type'] = np.where(gene_counts['associated_gene'].fillna('').str.startswith('novel'), 'novel', 'annotated')
        bins = gene_counts.groupby(['CB','gene_type']).agg(
            bin1_count=('read_count', lambda s: (s == 1).sum()),
            bin2_3_count=('read_count', lambda s: ((s >= 2) & (s <= 3)).sum()),
            bin4_5_count=('read_count', lambda s: ((s >= 4) & (s <= 5)).sum()),
            bin6plus_count=('read_count', lambda s: (s >= 6).sum()),
            total_genes_in_type=('associated_gene','nunique')
        ).reset_index()
        def bin_props(df, gene_kind, out_prefix):
            out = pd.DataFrame(index=summary.index)
            keyed = df[df['gene_type'] == gene_kind].set_index('CB') if not df.empty else pd.DataFrame(index=summary.index)
            for label, src in [(f"{out_prefix}_bin1_perc", 'bin1_count'), (f"{out_prefix}_bin2_3_perc", 'bin2_3_count'), (f"{out_prefix}_bin4_5_perc", 'bin4_5_count'), (f"{out_prefix}_bin6plus_perc", 'bin6plus_count')]:
                if not keyed.empty and src in keyed.columns:
                    out[label] = safe_prop(keyed[src].reindex(summary.index, fill_value=0), keyed['total_genes_in_type'].reindex(summary.index, fill_value=0)).fillna(0)
                else:
                    out[label] = 0
            return out
        summary = summary.join(bin_props(bins, 'annotated', 'anno')).join(bin_props(bins, 'novel', 'novel'))

        # UJC bins per CB for multiexonic genes
        gene_ujc = cls_valid[cls_valid['exons'] > 1].groupby(['CB','associated_gene'])['jxn_string'].nunique().rename('ujc_count').reset_index()
        gene_ujc['gene_type'] = np.where(gene_ujc['associated_gene'].fillna('').str.startswith('novel'), 'novel', 'annotated')
        ujc_bins = gene_ujc.groupby(['CB','gene_type']).agg(
            ujc_bin1_count=('ujc_count', lambda s: (s == 1).sum()),
            ujc_bin2_3_count=('ujc_count', lambda s: ((s >= 2) & (s <= 3)).sum()),
            ujc_bin4_5_count=('ujc_count', lambda s: ((s >= 4) & (s <= 5)).sum()),
            ujc_bin6plus_count=('ujc_count', lambda s: (s >= 6).sum()),
            total_genes_in_type_ujc=('associated_gene','nunique')
        ).reset_index()
        def ujc_props(df, gene_kind, out_prefix):
            out = pd.DataFrame(index=summary.index)
            keyed = df[df['gene_type'] == gene_kind].set_index('CB') if not df.empty else pd.DataFrame(index=summary.index)
            for label, src in [(f"{out_prefix}_ujc_bin1_perc", 'ujc_bin1_count'), (f"{out_prefix}_ujc_bin2_3_perc", 'ujc_bin2_3_count'), (f"{out_prefix}_ujc_bin4_5_perc", 'ujc_bin4_5_count'), (f"{out_prefix}_ujc_bin6plus_perc", 'ujc_bin6plus_count')]:
                if not keyed.empty and src in keyed.columns:
                    out[label] = safe_prop(keyed[src].reindex(summary.index, fill_value=0), keyed['total_genes_in_type_ujc'].reindex(summary.index, fill_value=0)).fillna(0)
                else:
                    out[label] = 0
            return out
        summary = summary.join(ujc_props(ujc_bins, 'annotated', 'anno')).join(ujc_props(ujc_bins, 'novel', 'novel'))

        # Length distributions (avoid groupby-apply deprecation) with final names directly
        def compute_lenbins_by_cb(df_group):
            gb = df_group.groupby('CB')['length']
            return pd.DataFrame({
                'two_fifty': gb.apply(lambda s: (s <= 250).sum()),
                'five_hund': gb.apply(lambda s: ((s > 250) & (s <= 500)).sum()),
                'short': gb.apply(lambda s: ((s > 500) & (s <= 1000)).sum()),
                'mid': gb.apply(lambda s: ((s > 1000) & (s <= 2000)).sum()),
                'long': gb.apply(lambda s: (s > 2000).sum()),
            })
        lo = compute_lenbins_by_cb(cls_valid)
        for dst, src in [('Total_250b_length_prop','two_fifty'),('Total_500b_length_prop','five_hund'),('Total_short_length_prop','short'),('Total_mid_length_prop','mid'),('Total_long_length_prop','long')]:
            summary[dst] = safe_prop(lo[src].reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        mono = cls_valid[cls_valid['exons'] == 1]
        if not mono.empty:
            lm = compute_lenbins_by_cb(mono)
            for dst, src in [('Total_250b_length_mono_prop','two_fifty'),('Total_500b_length_mono_prop','five_hund'),('Total_short_length_mono_prop','short'),('Total_mid_length_mono_prop','mid'),('Total_long_length_mono_prop','long')]:
                summary[dst] = safe_prop(lm[src].reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        else:
            for nm in ['Total_250b_length_mono_prop','Total_500b_length_mono_prop','Total_short_length_mono_prop','Total_mid_length_mono_prop','Total_long_length_mono_prop']:
                summary[nm] = 0
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            denom = summary[final_count_name(cat)]
            sub = cls_valid[cls_valid['structural_category'] == cat]
            if sub.empty:
                for nm in ['250b_length_prop','500b_length_prop','short_length_prop','mid_length_prop','long_length_prop']:
                    summary[f"{tag}_{nm}"] = 0
                for nm in ['250b_length_mono_prop','500b_length_mono_prop','short_length_mono_prop','mid_length_mono_prop','long_length_mono_prop']:
                    summary[f"{tag}_{nm}"] = 0
            else:
                lc = compute_lenbins_by_cb(sub)
                for dst, src in [('250b_length_prop','two_fifty'),('500b_length_prop','five_hund'),('short_length_prop','short'),('mid_length_prop','mid'),('long_length_prop','long')]:
                    summary[f"{tag}_{dst}"] = safe_prop(lc[src].reindex(summary.index, fill_value=0), denom).fillna(0)
                subm = sub[sub['exons'] == 1]
                if subm.empty:
                    for nm in ['250b_length_mono_prop','500b_length_mono_prop','short_length_mono_prop','mid_length_mono_prop','long_length_mono_prop']:
                        summary[f"{tag}_{nm}"] = 0
                else:
                    lcm = compute_lenbins_by_cb(subm)
                    for dst, src in [('250b_length_mono_prop','two_fifty'),('500b_length_mono_prop','five_hund'),('short_length_mono_prop','short'),('mid_length_mono_prop','mid'),('long_length_mono_prop','long')]:
                        summary[f"{tag}_{dst}"] = safe_prop(lcm[src].reindex(summary.index, fill_value=0), denom).fillna(0)

        # Reference body coverage >=45%
        cls_valid['ref_body_cov_flag'] = (cls_valid['length'] / cls_valid['ref_length'] * 100.0) >= 45.0
        for cat in structural_categories:
            tag = cat_to_tag[cat]
            sub = cls_valid[cls_valid['structural_category'] == cat]
            denom = summary[final_count_name(cat)]
            cov = sub[sub['ref_body_cov_flag']].groupby('CB').size() if not sub.empty else pd.Series(dtype=float)
            summary[f"{tag}_ref_coverage_prop"] = safe_prop(cov.reindex(summary.index, fill_value=0), denom).fillna(0)

        # RTS
        rts = cls_valid[cls_valid['RTS_stage'].astype(str) == 'TRUE'].groupby('CB').size()
        summary['RTS_prop_in_cell'] = safe_prop(rts.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        for abbr, cat in abbr_pairs:
            root = cat_to_root[cat]
            denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB').size().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['RTS_stage'].astype(str) == 'TRUE')].groupby('CB').size().reindex(summary.index, fill_value=0)
            summary[f"{abbr}_RTS_prop"] = safe_prop(numer, denom).fillna(0)

        # Non-canonical
        noncanon = cls_valid[cls_valid['all_canonical'] == 'non_canonical'].groupby('CB').size()
        summary['Non_canonical_prop_in_cell'] = safe_prop(noncanon.reindex(summary.index, fill_value=0), summary['total_reads_no_monoexon']).fillna(0)
        for abbr, cat in abbr_pairs:
            denom = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1)].groupby('CB').size().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1) & (cls_valid['all_canonical'] == 'non_canonical')].groupby('CB').size().reindex(summary.index, fill_value=0)
            summary[f"{abbr}_noncanon_prop"] = safe_prop(numer, denom).fillna(0)

        # Intrapriming
        intr = cls_valid[cls_valid['perc_A_downstream_TTS'] >= 60].groupby('CB').size()
        summary['Intrapriming_prop_in_cell'] = safe_prop(intr.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        for abbr, cat in abbr_pairs:
            denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB').size().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['perc_A_downstream_TTS'] >= 60)].groupby('CB').size().reindex(summary.index, fill_value=0)
            summary[f"{abbr}_intrapriming_prop"] = safe_prop(numer, denom).fillna(0)

        # TSS support
        tss_sup = (cls_valid['diff_to_gene_TSS'].abs() <= 50)
        sup_cnt = cls_valid[tss_sup].groupby('CB').size()
        summary['TSSAnnotationSupport_prop'] = safe_prop(sup_cnt.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
        for abbr, cat in abbr_pairs:
            denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB').size().reindex(summary.index, fill_value=0)
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['diff_to_gene_TSS'].abs() <= 50)].groupby('CB').size().reindex(summary.index, fill_value=0)
            summary[f"{abbr}_TSSAnnotationSupport"] = safe_prop(numer, denom).fillna(0)

        # Annotated genes in cell prop
        summary['Annotated_genes_prop_in_cell'] = safe_prop(summary['Annotated_genes'], summary['Genes_in_cell']).fillna(0)
        for abbr, cat in [('FSM','full-splice_match'),('ISM','incomplete-splice_match'),('NIC','novel_in_catalog'),('NNC','novel_not_in_catalog')]:
            sub = cls_valid[cls_valid['structural_category'] == cat]
            total_genes = sub.groupby('CB')['associated_gene'].nunique()
            anno_genes = sub[~sub['associated_gene'].fillna('').str.startswith('novel')].groupby('CB')['associated_gene'].nunique()
            summary[f"{abbr}_anno_genes_prop"] = safe_prop(anno_genes.reindex(summary.index, fill_value=0), total_genes.reindex(summary.index, fill_value=0)).fillna(0)

        # Annotated junction strings (multiexonic, non-novel transcript)
        anno_models = cls_valid[(cls_valid['exons'] > 1) & (~cls_valid['associated_transcript'].fillna('').str.startswith('novel'))] \
            .groupby('CB')['jxn_string'].nunique()
        summary['Annotated_juction_strings_prop_in_cell'] = safe_prop(anno_models.reindex(summary.index, fill_value=0), summary['UJCs_in_cell']).fillna(0)

        # Canonical overall and per category
        canon_over = cls_valid[cls_valid['all_canonical'] == 'canonical'].groupby('CB').size()
        summary['Canonical_prop_in_cell'] = safe_prop(canon_over.reindex(summary.index, fill_value=0), summary['total_reads_no_monoexon']).fillna(0)
        for abbr, cat in abbr_pairs:
            denom = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1)].groupby('CB').size()
            numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['exons'] > 1) & (cls_valid['all_canonical'] == 'canonical')].groupby('CB').size()
            summary[f"{abbr}_canon_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)

        # ORF-related
        if not args.skipORF:
            nmd = cls_valid[cls_valid['predicted_NMD'].astype(str) == 'TRUE'].groupby('CB').size()
            summary['NMD_prop_in_cell'] = safe_prop(nmd.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
            for abbr, cat in abbr_pairs:
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB').size()
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['predicted_NMD'].astype(str) == 'TRUE')].groupby('CB').size()
                summary[f"{abbr}_NMD_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                sub = cls_valid[cls_valid['structural_category'] == cat]
                denom = sub.groupby('CB').size()
                cod = sub[sub['coding'] == 'coding'].groupby('CB').size()
                ncod = sub[sub['coding'] == 'non_coding'].groupby('CB').size()
                coding_tag = tag if tag in ['FSM','ISM','NIC','NNC'] else tag.lower()
                summary[f"Coding_{coding_tag}_prop"] = safe_prop(cod.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
                summary[f"Non_coding_{coding_tag}_prop"] = safe_prop(ncod.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
        else:
            summary['NMD_prop_in_cell'] = 0
            for cat in structural_categories:
                tag = cat_to_tag[cat]
                coding_tag = tag if tag in ['FSM','ISM','NIC','NNC'] else tag.lower()
                summary[f"Coding_{coding_tag}_prop"] = 0
                summary[f"Non_coding_{coding_tag}_prop"] = 100

        # CAGE
        if getattr(args, 'CAGE_peak', None):
            cage = cls_valid[cls_valid['within_CAGE_peak'].astype(str) == 'TRUE'].groupby('CB').size()
            summary['CAGE_peak_support_prop'] = safe_prop(cage.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
            for abbr, cat in abbr_pairs:
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB').size()
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['within_CAGE_peak'].astype(str) == 'TRUE')].groupby('CB').size()
                summary[f"{abbr}_CAGE_peak_support_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
        else:
            summary['CAGE_peak_support_prop'] = 0
            for abbr in ['FSM','ISM','NIC','NNC']:
                summary[f"{abbr}_CAGE_peak_support_prop"] = 0

        # PolyA motif
        if getattr(args, 'polyA_motif_list', None):
            pa = cls_valid[cls_valid['polyA_motif_found'].astype(str) == 'TRUE'].groupby('CB').size()
            summary['PolyA_motif_support_prop'] = safe_prop(pa.reindex(summary.index, fill_value=0), summary['total_reads']).fillna(0)
            for abbr, cat in abbr_pairs:
                denom = cls_valid[cls_valid['structural_category'] == cat].groupby('CB').size()
                numer = cls_valid[(cls_valid['structural_category'] == cat) & (cls_valid['polyA_motif_found'].astype(str) == 'TRUE')].groupby('CB').size()
                summary[f"{abbr}_PolyA_motif_support_prop"] = safe_prop(numer.reindex(summary.index, fill_value=0), denom.reindex(summary.index, fill_value=0)).fillna(0)
        else:
            summary['PolyA_motif_support_prop'] = 0
            for abbr in ['FSM','ISM','NIC','NNC']:
                summary[f"{abbr}_PolyA_motif_support_prop"] = 0

        # Finalize table: reset index and minimally rename totals
        summary = summary.reset_index()
        summary = summary.rename(columns={'total_reads': 'Reads_in_cell', 'total_UMI': 'UMIs_in_cell'})

        # Coerce numeric
        for c in summary.columns[1:]:
            summary[c] = pd.to_numeric(summary[c], errors='coerce').fillna(0)

        try:
            summary.to_csv(out_summary, sep='\t', index=False, compression='gzip')
            print(f"**** Cell summary written: {out_summary}", file=sys.stdout)
        except Exception as e:
            print(f"[ERROR] Failed writing {out_summary}: {e}", file=sys.stderr)

def write_gene_counts_by_cell(args, df):
    """
    For each sample, read its classification file and write a gene_counts.csv
    that summarizes counts per associated_gene per CB, including structural
    category counts, total read count, unique junction-hash counts, and an
    annotated-gene flag.

    Output: <out_dir>/<file_acc>/<sampleID>_gene_counts.csv
    Columns: associated_gene, CB, <structural categories>, total_read_count,
             unique_jxnHash_counts, flag_annotated_gene
    """
    structural_categories = [
        'antisense',
        'full-splice_match',
        'fusion',
        'genic',
        'genic_intron',
        'incomplete-splice_match',
        'intergenic',
        'novel_in_catalog',
        'novel_not_in_catalog',
    ]

    for _, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{outputPathPrefix}_classification.txt"

        if not os.path.isfile(class_file):
            print(f"[INFO] Classification file not found for {file_acc}. Skipping gene counts.",
                  file=sys.stdout)
            continue

        try:
            needed_cols = ['associated_gene', 'structural_category', 'jxnHash', 'CB']
            class_df = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False,
                                   usecols=lambda c: c in needed_cols or c == 'isoform')
            # Reduce memory and speed up groupbys
            for col in ['associated_gene', 'structural_category', 'CB']:
                if col in class_df.columns:
                    try:
                        class_df[col] = class_df[col].astype('category')
                    except Exception:
                        pass
        except Exception as e:
            print(f"[ERROR] Failed reading {class_file}: {e}", file=sys.stderr)
            continue

        # Ensure required columns exist; if not, create with empty values
        for col in ['associated_gene', 'structural_category', 'CB']:
            if col not in class_df.columns:
                class_df[col] = ''
        if 'jxnHash' not in class_df.columns:
            # If UJC hashing was skipped, fallback to isoform as proxy uniqueness within (gene,CB)
            class_df['jxnHash'] = class_df.get('isoform', pd.Series(index=class_df.index, dtype=str))

        # Keep rows with valid CBs only
        filtered = class_df[(class_df['CB'].notna()) & (class_df['CB'] != '')].copy()
        if filtered.empty:
            # Still write an empty file with header for consistency
            out_file = f"{outputPathPrefix}_gene_counts.csv"
            empty_cols = (['associated_gene', 'CB'] + structural_categories +
                          ['total_read_count', 'unique_jxnHash_counts', 'flag_annotated_gene'])
            pd.DataFrame(columns=empty_cols).to_csv(out_file, index=False)
            print(f"[INFO] No valid CB rows for {file_acc}. Wrote empty gene counts.", file=sys.stdout)
            continue

        # Pivot counts per structural category within (associated_gene, CB)
        grouped_size = (
            filtered
            .groupby(['associated_gene', 'CB', 'structural_category'], observed=True)
            .size()
            .unstack(fill_value=0)
        )

        # Ensure all expected categories are present as columns
        for cat in structural_categories:
            if cat not in grouped_size.columns:
                grouped_size[cat] = 0
        grouped_size = grouped_size[structural_categories]

        # Total read count per (gene, CB)
        total_counts = (
            filtered
            .groupby(['associated_gene', 'CB'], observed=True)
            .size()
            .rename('total_read_count')
        )

        # Unique junction-hash counts per (gene, CB)
        unique_jxn_counts = (
            filtered
            .groupby(['associated_gene', 'CB'], observed=True)['jxnHash']
            .nunique()
            .rename('unique_jxnHash_counts')
        )

        # Merge all aggregates
        agg_df = grouped_size.merge(total_counts, left_index=True, right_index=True)
        agg_df = agg_df.merge(unique_jxn_counts, left_index=True, right_index=True)
        agg_df = agg_df.reset_index()

        # Annotated gene flag (1 if not novel*)
        agg_df['flag_annotated_gene'] = agg_df['associated_gene'].apply(
            lambda g: 0 if isinstance(g, str) and g.startswith('novel') else 1
        )

        # Final column order
        final_cols = (['associated_gene', 'CB'] + structural_categories +
                      ['total_read_count', 'unique_jxnHash_counts', 'flag_annotated_gene'])
        agg_df = agg_df[final_cols]

        out_file = f"{outputPathPrefix}_gene_counts.csv"
        try:
            agg_df.to_csv(out_file, index=False)
            print(f"**** Gene counts by cell written: {out_file}", file=sys.stdout)
        except Exception as e:
            print(f"[ERROR] Failed writing {out_file}: {e}", file=sys.stderr)


def write_ujc_counts_by_cell(args, df):
    """
    For each sample, write <sampleID>_ujc_counts.csv where each row
    corresponds to a unique junction chain (jxnHash) within a CB.

    Columns:
      - jxnHash, CB
      - read_count (reads with this jxnHash in this CB)
      - associated_gene (majority label among reads in this (jxnHash, CB))
      - known_canonical, known_non_canonical, novel_canonical, novel_non_canonical
        (counts of junctions in the chain, derived from junctions file via isoform->jxnHash map)
      - structural_category (majority among reads in this (jxnHash, CB))
      - flag_MEI (1 if any read in (jxnHash, CB) has exons > 1, else 0)
      - flag_annotated_gene (1 if associated_gene does not start with 'novel')
    """
    for _, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{outputPathPrefix}_classification.txt"
        junc_file = f"{outputPathPrefix}_junctions.txt"

        if not os.path.isfile(class_file):
            print(f"[INFO] Classification file not found for {file_acc}. Skipping UJC counts.",
                  file=sys.stdout)
            continue

        try:
            needed_cls_cols = ['isoform', 'associated_gene', 'structural_category', 'exons', 'jxnHash', 'CB']
            cls_df = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False,
                                 usecols=lambda c: c in needed_cls_cols)
            for col in ['associated_gene', 'structural_category', 'CB']:
                if col in cls_df.columns:
                    try:
                        cls_df[col] = cls_df[col].astype('category')
                    except Exception:
                        pass
        except Exception as e:
            print(f"[ERROR] Failed reading {class_file}: {e}", file=sys.stderr)
            continue

        # Normalize types and filter valid CBs
        for col in ['associated_gene', 'structural_category', 'CB', 'jxnHash']:
            if col not in cls_df.columns:
                cls_df[col] = ''
        if 'exons' not in cls_df.columns:
            cls_df['exons'] = ''
        cls_df['exons_num'] = pd.to_numeric(cls_df['exons'], errors='coerce')
        cls_df = cls_df[(cls_df['CB'].notna()) & (cls_df['CB'] != '') & (cls_df['jxnHash'].notna()) & (cls_df['jxnHash'] != '')].copy()

        if cls_df.empty:
            out_file = f"{outputPathPrefix}_ujc_counts.csv"
            empty_cols = ['jxnHash', 'CB', 'read_count', 'associated_gene',
                          'known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical',
                          'structural_category', 'flag_MEI', 'flag_annotated_gene']
            pd.DataFrame(columns=empty_cols).to_csv(out_file, index=False)
            print(f"[INFO] No valid CB/jxnHash rows for {file_acc}. Wrote empty UJC counts.", file=sys.stdout)
            continue

        # Majority functions
        def majority_value(series):
            counts = series.value_counts(dropna=True)
            return counts.idxmax() if len(counts) > 0 else ''

        def flag_mei(series):
            # Here, all rows are multiexonic by construction; keep robust fallback
            return 1 if (series.fillna(0).astype(float) > 1).any() else 0

        # Split into multiexonic and monoexonic for different grouping logic
        multi_df = cls_df[cls_df['exons_num'] > 1].copy()
        mono_df = cls_df[cls_df['exons_num'] <= 1].copy()

        # Multiexonic: aggregate at (jxnHash, CB), majority labels
        base_agg_multi = multi_df.groupby(['jxnHash', 'CB'], observed=True).agg(
            read_count=('isoform', 'size'),
            associated_gene=('associated_gene', majority_value),
            structural_category=('structural_category', majority_value),
            flag_MEI=('exons_num', flag_mei)
        ).reset_index()

        # Monoexonic: aggregate at (jxnHash, CB, associated_gene, structural_category)
        if not mono_df.empty:
            base_agg_mono = mono_df.groupby(['jxnHash', 'CB', 'associated_gene', 'structural_category'], observed=True).agg(
                read_count=('isoform', 'size')
            ).reset_index()
            base_agg_mono['flag_MEI'] = 0
        else:
            base_agg_mono = pd.DataFrame(columns=['jxnHash','CB','associated_gene','structural_category','read_count','flag_MEI'])

        # Compute junction-type counts per jxnHash from junctions file if present
        jxn_counts = None
        if os.path.isfile(junc_file):
            try:
                junc_needed = ['isoform', 'junction_category', 'canonical']
                jdf = pd.read_csv(junc_file, sep='\t', dtype=str, low_memory=False,
                                  usecols=lambda c: c in junc_needed)
                iso_to_hash = cls_df[['isoform', 'jxnHash']].drop_duplicates()
                jdf = pd.merge(jdf, iso_to_hash, on='isoform', how='left')
                jdf = jdf[jdf['jxnHash'].notna() & (jdf['jxnHash'] != '')]

                # Determine type
                def _type(row):
                    if row['junction_category'] == 'known' and row['canonical'] == 'canonical':
                        return 'known_canonical'
                    if row['junction_category'] == 'known' and row['canonical'] == 'non_canonical':
                        return 'known_non_canonical'
                    if row['junction_category'] == 'novel' and row['canonical'] == 'canonical':
                        return 'novel_canonical'
                    if row['junction_category'] == 'novel' and row['canonical'] == 'non_canonical':
                        return 'novel_non_canonical'
                    return None

                jdf['jxn_type'] = jdf.apply(_type, axis=1)
                jdf = jdf[jdf['jxn_type'].notna()]
                jxn_counts = jdf.groupby(['jxnHash', 'jxn_type']).size().unstack(fill_value=0)
                for col in ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']:
                    if col not in jxn_counts.columns:
                        jxn_counts[col] = 0
                jxn_counts = jxn_counts[['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']]
                jxn_counts = jxn_counts.reset_index()
            except Exception as e:
                print(f"[WARNING] Failed to compute junction counts for {file_acc}: {e}", file=sys.stderr)

        # Attach junction-type counts to multiexonic subset; monoexonic has zeros
        if jxn_counts is not None and not base_agg_multi.empty:
            out_multi = pd.merge(base_agg_multi, jxn_counts, on='jxnHash', how='left')
            for col in ['known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical']:
                if col not in out_multi.columns:
                    out_multi[col] = 0
                out_multi[col] = out_multi[col].fillna(0).astype(int)
        else:
            out_multi = base_agg_multi.copy()
            out_multi['known_canonical'] = 0
            out_multi['known_non_canonical'] = 0
            out_multi['novel_canonical'] = 0
            out_multi['novel_non_canonical'] = 0

        # Monoexonic: explicitly zero junction counts
        if not base_agg_mono.empty:
            out_mono = base_agg_mono.copy()
            out_mono['known_canonical'] = 0
            out_mono['known_non_canonical'] = 0
            out_mono['novel_canonical'] = 0
            out_mono['novel_non_canonical'] = 0
        else:
            out_mono = base_agg_mono

        # Normalize columns and concatenate
        # Ensure both have required columns
        required_cols = ['jxnHash','CB','read_count','associated_gene','structural_category','flag_MEI',
                         'known_canonical','known_non_canonical','novel_canonical','novel_non_canonical']
        for df_norm in (out_multi, out_mono):
            for col in required_cols:
                if col not in df_norm.columns:
                    df_norm[col] = 0 if col in ['read_count','flag_MEI','known_canonical','known_non_canonical','novel_canonical','novel_non_canonical'] else ''
        out_df = pd.concat([out_multi[required_cols], out_mono[required_cols]], ignore_index=True)

        # Recompute flag_MEI to match SQANTI-reads semantics: most-expressed UJC per gene (per CB)
        out_df['flag_MEI'] = out_df.groupby(['associated_gene','CB'])['read_count'] \
            .transform(lambda s: (s == s.max()).astype(int))

        # Annotated gene flag
        out_df['flag_annotated_gene'] = out_df['associated_gene'].apply(
            lambda g: 0 if isinstance(g, str) and g.startswith('novel') else 1
        )

        # Final order
        final_cols = ['jxnHash', 'CB', 'read_count', 'associated_gene',
                      'known_canonical', 'known_non_canonical', 'novel_canonical', 'novel_non_canonical',
                      'structural_category', 'flag_MEI', 'flag_annotated_gene']
        out_df = out_df[final_cols]

        out_file = f"{outputPathPrefix}_ujc_counts.csv"
        try:
            out_df.to_csv(out_file, index=False)
            print(f"**** UJC counts by cell written: {out_file}", file=sys.stdout)
        except Exception as e:
            print(f"[ERROR] Failed writing {out_file}: {e}", file=sys.stderr)


def write_cv_by_cell(args, df):
    """
    For each sample, write <sampleID>_cv.csv reporting donor/acceptor reference
    splice sites statistics per CB (cell barcode).

    Output columns:
      - chrom, strand, coord, mean_abs_diff, std_abs_diff, cv, count,
        associated_gene, gene_count, flag_multi_gene, flag_annotated_gene,
        flag_donor, flag_acceptor, flag_mean_0, flag_std_0, flag_ref_match,
        flag_cv_0, flag_cv_gt_0, CB
    """
    def safe_to_numeric(series):
        return pd.to_numeric(series, errors='coerce')

    def aggregate_cv(df_grouped, coord_col, abs_diff_col):
        agg = df_grouped.agg(
            mean_abs_diff=(abs_diff_col, 'mean'),
            std_abs_diff=(abs_diff_col, 'std'),
            count=(abs_diff_col, 'size'),
            associated_gene=('associated_gene', lambda s: '|'.join(sorted(set(s.dropna())))),
            gene_count=('associated_gene', lambda s: len(set([g for g in s.dropna()])))
        ).reset_index()

        group_keys = ['chrom', 'strand', coord_col, 'CB']
        cv_vals = df_grouped[abs_diff_col].apply(
            lambda s: np.nan if pd.isna(s.mean()) or s.mean() == 0 else float(np.std(s, ddof=0)) / float(s.mean())
        ).reset_index(name='cv')
        agg = pd.merge(agg, cv_vals, on=group_keys, how='left')

        # Flags
        agg['flag_mean_0'] = agg['mean_abs_diff'].apply(lambda x: 1 if pd.notna(x) and x == 0 else 0)
        agg['flag_std_0'] = agg['std_abs_diff'].apply(lambda x: 0 if pd.isna(x) else (1 if x == 0 else 0))
        agg['flag_ref_match'] = agg['flag_mean_0']
        agg['flag_cv_0'] = agg.apply(lambda r: 1 if pd.notna(r['cv']) and r['cv'] == 0 and r['flag_mean_0'] != 1 else 0, axis=1)
        agg['flag_cv_gt_0'] = agg['cv'].apply(lambda x: 1 if pd.notna(x) and x > 0 else 0)
        agg['flag_multi_gene'] = agg['gene_count'].apply(lambda x: 1 if x > 1 else 0)

        def is_annotated(genestr):
            genes = [g for g in str(genestr).split('|') if g]
            if len(genes) == 0:
                return 0
            return 1 if all(not str(g).startswith('novel') for g in genes) else 0

        agg['flag_annotated_gene'] = agg['associated_gene'].apply(is_annotated)

        agg = agg.rename(columns={coord_col: 'coord'})
        return agg

    for _, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        class_file = f"{outputPathPrefix}_classification.txt"
        junc_file = f"{outputPathPrefix}_junctions.txt"

        if not os.path.isfile(class_file):
            print(f"[INFO] Classification file not found for {file_acc}. Skipping CV.", file=sys.stdout)
            continue
        if not os.path.isfile(junc_file):
            print(f"[INFO] Junctions file not found for {file_acc}. Skipping CV.", file=sys.stdout)
            continue

        try:
            # Read classification for gene and CB mapping
            class_needed = ['isoform', 'associated_gene', 'CB']
            cls_df = pd.read_csv(class_file, sep='\t', dtype=str, low_memory=False,
                                 usecols=lambda c: c in class_needed)
        except Exception as e:
            print(f"[ERROR] Failed reading {class_file}: {e}", file=sys.stderr)
            continue

        try:
            # Read junctions with required columns
            j_needed = ['isoform', 'chrom', 'strand', 'genomic_start_coord', 'genomic_end_coord',
                        'diff_to_Ref_start_site', 'diff_to_Ref_end_site']
            jdf = pd.read_csv(junc_file, sep='\t', dtype=str, low_memory=False,
                              usecols=lambda c: c in j_needed)
        except Exception as e:
            print(f"[ERROR] Failed reading {junc_file}: {e}", file=sys.stderr)
            continue

        # Merge in associated_gene and CB
        jdf = pd.merge(jdf, cls_df, on='isoform', how='left')

        # Keep rows with valid CB
        jdf = jdf[(jdf['CB'].notna()) & (jdf['CB'] != '')].copy()
        if jdf.empty:
            out_file = f"{outputPathPrefix}_cv.csv"
            empty_cols = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff', 'cv', 'count',
                          'associated_gene', 'gene_count', 'flag_multi_gene', 'flag_annotated_gene',
                          'flag_donor', 'flag_acceptor', 'flag_mean_0', 'flag_std_0', 'flag_ref_match',
                          'flag_cv_0', 'flag_cv_gt_0', 'CB']
            pd.DataFrame(columns=empty_cols).to_csv(out_file, index=False)
            print(f"[INFO] No valid CB rows for {file_acc}. Wrote empty CV table.", file=sys.stdout)
            continue

        # Normalize numeric fields
        jdf['genomic_start_coord'] = safe_to_numeric(jdf['genomic_start_coord'])
        jdf['genomic_end_coord'] = safe_to_numeric(jdf['genomic_end_coord'])
        jdf['diff_to_Ref_start_site'] = safe_to_numeric(jdf['diff_to_Ref_start_site'])
        jdf['diff_to_Ref_end_site'] = safe_to_numeric(jdf['diff_to_Ref_end_site'])

        # Drop rows where diffs are NaN (cannot compute reference junction)
        jdf = jdf.dropna(subset=['diff_to_Ref_start_site', 'diff_to_Ref_end_site'])
        if jdf.empty:
            out_file = f"{outputPathPrefix}_cv.csv"
            empty_cols = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff', 'cv', 'count',
                          'associated_gene', 'gene_count', 'flag_multi_gene', 'flag_annotated_gene',
                          'flag_donor', 'flag_acceptor', 'flag_mean_0', 'flag_std_0', 'flag_ref_match',
                          'flag_cv_0', 'flag_cv_gt_0', 'CB']
            pd.DataFrame(columns=empty_cols).to_csv(out_file, index=False)
            print(f"[INFO] No diff-to-ref data for {file_acc}. Wrote empty CV table.", file=sys.stdout)
            continue

        # Compute reference junction coordinates and abs diffs (vectorized)
        jdf['ref_junction_start'] = jdf['genomic_start_coord'] + jdf['diff_to_Ref_start_site']
        jdf['ref_junction_end'] = jdf['genomic_end_coord'] + jdf['diff_to_Ref_end_site']
        jdf['abs_diff_start'] = jdf['diff_to_Ref_start_site'].abs()
        jdf['abs_diff_end'] = jdf['diff_to_Ref_end_site'].abs()

        # Group by CB and reference splice sites
        start_group = jdf.groupby(['chrom', 'strand', 'ref_junction_start', 'CB'], dropna=False, observed=True)
        end_group = jdf.groupby(['chrom', 'strand', 'ref_junction_end', 'CB'], dropna=False, observed=True)

        cv_start = aggregate_cv(start_group, 'ref_junction_start', 'abs_diff_start')
        cv_end = aggregate_cv(end_group, 'ref_junction_end', 'abs_diff_end')

        # Donor/Acceptor flags according to strand and coord type
        cv_start['flag_donor'] = cv_start['strand'].apply(lambda s: 1 if s == '+' else 0)
        cv_start['flag_acceptor'] = cv_start['strand'].apply(lambda s: 1 if s == '-' else 0)
        cv_end['flag_donor'] = cv_end['strand'].apply(lambda s: 1 if s == '-' else 0)
        cv_end['flag_acceptor'] = cv_end['strand'].apply(lambda s: 1 if s == '+' else 0)

        out_df = pd.concat([cv_start, cv_end], ignore_index=True)

        # Final column order
        final_cols = ['chrom', 'strand', 'coord', 'mean_abs_diff', 'std_abs_diff', 'cv', 'count',
                      'associated_gene', 'gene_count', 'flag_multi_gene', 'flag_annotated_gene',
                      'flag_donor', 'flag_acceptor', 'flag_mean_0', 'flag_std_0', 'flag_ref_match',
                      'flag_cv_0', 'flag_cv_gt_0', 'CB']
        # Ensure all columns exist
        for col in final_cols:
            if col not in out_df.columns:
                out_df[col] = np.nan if col in ['mean_abs_diff', 'std_abs_diff', 'cv'] else 0
        out_df = out_df[final_cols]

        # Ensure coordinate prints without .0 (integer formatting)
        out_df['coord'] = pd.to_numeric(out_df['coord'], errors='coerce').astype('Int64')

        out_file = f"{outputPathPrefix}_cv.csv"
        try:
            out_df.to_csv(out_file, index=False)
            print(f"**** CV by cell written: {out_file}", file=sys.stdout)
        except Exception as e:
            print(f"[ERROR] Failed writing {out_file}: {e}", file=sys.stderr)


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
                # Pass precomputed cell summary
                cell_summary_file = f"{outputPathPrefix}_SQANTI_cell_summary.txt.gz"
                if os.path.isfile(cell_summary_file):
                    flags.extend(["--cell_summary", cell_summary_file])

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


def generate_multisample_report(args, df):
    """
    Generate a multisample report by aggregating per-sample SQANTI single-cell summaries.
    It searches for <out_dir>/<file_acc>/<sampleID>_SQANTI_cell_summary.txt.gz
    for each row in the design file and invokes the cohort R script.
    """
    try:
        total_samples = df.shape[0]
    except Exception:
        total_samples = 0

    if total_samples < 2:
        print("[INFO] Design has fewer than 2 samples. Skipping multisample report.",
              file=sys.stdout)
        return

    # Collect existing cell summary files
    cell_summaries = []
    for _, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        outputPathPrefix = os.path.join(args.out_dir, file_acc, sampleID)
        cell_summary = f"{outputPathPrefix}_SQANTI_cell_summary.txt.gz"
        if os.path.isfile(cell_summary):
            cell_summaries.append(os.path.abspath(cell_summary))
        else:
            print(f"[INFO] Cell summary not found for {file_acc} ({sampleID}). Skipping this sample.",
                  file=sys.stdout)

    if len(cell_summaries) < 2:
        print("[INFO] Fewer than 2 cell summaries found. Skipping multisample report.",
              file=sys.stdout)
        return

    # Build command to run the multisample R script
    prefix = getattr(args, 'multisample_report_prefix', 'SQANTI_sc_multisample_report')
    files_arg = ",".join(cell_summaries)
    out_dir = os.path.abspath(args.out_dir)
    report_fmt = args.report
    mode = args.mode

    cmd = (
        f"Rscript {utilitiesPath}/SQANTI-sc_multisample.R "
        f"--files \"{files_arg}\" --out_dir \"{out_dir}\" "
        f"--mode {mode} --report {report_fmt} --prefix \"{prefix}\""
    )

    print("**** Generating multisample SQANTI-sc report...", file=sys.stdout)
    try:
        subprocess.run(cmd, shell=True, check=True)
        print("**** Multisample SQANTI-sc report generated.", file=sys.stdout)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Multisample report generation failed: {e}", file=sys.stderr)

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
    # Multisample report options
    apsc.add_argument('--multisample_report', action='store_true', default=False,
                     help='Generate a multisample cohort report from per-sample cell summaries.')
    apsc.add_argument('--multisample_report_prefix', default='SQANTI_sc_multisample_report',
                     help='Output prefix for the multisample report (default: SQANTI_sc_multisample_report).')

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
    apout.add_argument('--write_per_cell_outputs', action='store_true', default=False,
                       help='When set, writes per-cell gene_counts.csv, ujc_counts.csv, and cv.csv for each sample.')

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

    # Compute per-sample cell summary
    calculate_metrics_per_cell(args, df)

    # Additional per CB outputs (heavy, optional)
    if getattr(args, 'write_per_cell_outputs', False):
        write_gene_counts_by_cell(args, df)
        write_ujc_counts_by_cell(args, df)
        write_cv_by_cell(args, df)
    else:
        print("[INFO] Skipping per-cell outputs (enable with --write_per_cell_outputs).", file=sys.stdout)

    generate_report(args, df)

    # Generate multi-sample report
    if getattr(args, 'multisample_report', False):
        generate_multisample_report(args, df)


if __name__ == "__main__":
    main()
