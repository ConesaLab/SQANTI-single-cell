import os
import sys
import glob
import pandas as pd
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

    def build_sqanti_command(input_file, file_acc, sampleID, is_fastq=False, coverage=None, SR_bam=None):
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
            "--report", "skip",
            "--force_id_ignore"
        ]

        if getattr(args, 'novel_gene_prefix', None):
            cmd_parts.extend(["--novel_gene_prefix", args.novel_gene_prefix])

        if is_fastq:
            cmd_parts.append("--fasta")

        if args.aligner_choice == "gmap" and getattr(args, 'gmap_index', None):
            cmd_parts.extend(["-x", args.gmap_index])

        flag_args = [
            ('genename', '--genename'),
            ('isoAnnotLite', '--isoAnnotLite'),
            ('isoform_hits', '--isoform_hits')
        ]
        for arg, flag in flag_args:
            if getattr(args, arg):
                cmd_parts.append(flag)

        # Handle skipORF logic:
        # If mode is 'reads', ALWAYS skip ORF.
        # If mode is 'isoforms', use the user's --skipORF flag (if provided).
        if args.mode == 'reads':
            cmd_parts.append('--skipORF')
        elif getattr(args, 'skipORF', False):
            cmd_parts.append('--skipORF')

        optional_files = [
            ("CAGE_peak", "--CAGE_peak"),
            ("polyA_motif_list", "--polyA_motif_list"),
            ("polyA_peak", "--polyA_peak"),
            ("phyloP_bed", "--phyloP_bed"),
            ("orf_input", "--orf_input"),
            ("expression", "--expression"),
            ("gff3", "--gff3"),
        ]
        for arg_name, flag in optional_files:
            if getattr(args, arg_name, None):
                cmd_parts.extend([flag, getattr(args, arg_name)])

        # Handle per-sample orthogonal data
        if coverage:
            # Handle coverage FOFN (File of File Names) support
            if os.path.isfile(coverage):
                # Heuristic: Check if file looks like a STAR junction file (based on extension or content)
                is_star_file = False
                if coverage.endswith(".out") or coverage.endswith(".tab") or coverage.endswith("SJ.out.tab"):
                    is_star_file = True
                
                if not is_star_file:
                    try:
                        with open(coverage, 'r') as f:
                            lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
                        
                        if lines:
                            # Verify if first line looks like a STAR junction record (9 columns)
                            first_line_parts = lines[0].split()
                            if len(first_line_parts) == 9:
                                is_star_file = True
                            else:
                                # It is likely a FOFN
                                parsed_coverage = ",".join(lines)
                                if args.verbose:
                                    print(f"[INFO] Parsed coverage FOFN '{coverage}' into {len(lines)} files.", file=sys.stdout)
                                coverage = parsed_coverage
                    except Exception as e:
                        if args.verbose:
                            print(f"[WARNING] Failed to parse coverage file as FOFN: {e}. Treating as single file.", file=sys.stderr)

            cmd_parts.extend(["--coverage", coverage])
        
        if SR_bam:
            cmd_parts.extend(["--SR_bam", SR_bam])


        return " ".join(cmd_parts)

    if not check_files_exist(args.refFasta, args.refGTF):
        print("[ERROR] Genome or annotation file is missing.", file=sys.stderr)
        sys.exit(-1)

    for index, row in df.iterrows():
        file_acc = row['file_acc']
        sampleID = row['sampleID']
        
        # Extract per-sample short reads orthogonal data (ONLY for Isoforms mode)
        coverage_file = None
        sr_bam_file = None
        if args.mode == 'isoforms':
            coverage_file = row['coverage'] if 'coverage' in row and pd.notna(row['coverage']) and row['coverage'] else None
            sr_bam_file = row['SR_bam'] if 'SR_bam' in row and pd.notna(row['SR_bam']) and row['SR_bam'] else None

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
            cmd = build_sqanti_command(gtf_file, file_acc, sampleID, coverage=coverage_file, SR_bam=sr_bam_file)
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
            cmd = build_sqanti_command(fastq_file, file_acc, sampleID, is_fastq=True, coverage=coverage_file, SR_bam=sr_bam_file)
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
                    tmp_fastq, file_acc, sampleID, is_fastq=True, coverage=coverage_file, SR_bam=sr_bam_file
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


