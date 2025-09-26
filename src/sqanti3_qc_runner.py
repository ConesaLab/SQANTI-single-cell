import os
import sys
import glob
from paths import sqantiqcPath
from src.commands import run_command
from src.module_logging import qc_logger, update_logger


def run_sqanti3_qc(args, df):
    def check_files_exist(*files):
        missing_files = [f for f in files if not os.path.isfile(f)]
        if missing_files:
            print(f'[ERROR] Missing file(s): {", ".join(missing_files)}',
                  file=sys.stderr)
            return False
        return True

    def build_sqanti_command(input_file, file_acc, sampleID, is_fastq=False):
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

        if getattr(args, 'novel_gene_prefix', None):
            cmd_parts.extend(["--novel_gene_prefix", args.novel_gene_prefix])

        if is_fastq:
            cmd_parts.append("--fasta")

        if args.aligner_choice == "gmap" and getattr(args, 'gmap_index', None):
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
            if getattr(args, arg_name, None):
                cmd_parts.extend([flag, getattr(args, arg_name)])

        return " ".join(cmd_parts)

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
                            "[ERROR] --orf_input must be provided for fusion unless --skipORF is specified."
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


