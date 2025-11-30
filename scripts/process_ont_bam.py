#!/usr/bin/env python3
"""
process_ont_bam.py

This script processes ONT single-cell BAM files for SQANTI-sc.
It performs the following steps for each sample:
1.  Runs `umi_tools dedup` with specific flags for ONT single-cell data.
2.  Runs `spliced_bam2gff` to convert the processed BAM to GFF format.

Usage:
    # Single file
    python process_ont_bam.py --bam /path/to/sample.tagged.bam [-o ./results]

    # Batch processing (FOFN)
    python process_ont_bam.py --fofn bam_file_list.txt [-o ./results]

FOFN Format (Text file):
    /path/to/sample1.tagged.bam
    /path/to/sample2.tagged.bam
    ...
"""

import argparse
import os
import subprocess
import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Batch process ONT BAM files for SQANTI-sc.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--bam", help="Path to a single BAM file to process")
    group.add_argument("--fofn", help="Path to a File of File Names (FOFN) containing one BAM path per line")
    
    parser.add_argument("-o", "--output_dir", default=".", help="Directory to save outputs (default: current directory)")
    parser.add_argument("--umi_tools_path", default="umi_tools", help="Path to umi_tools executable (default: umi_tools)")
    parser.add_argument("--spliced_bam2gff_path", default="spliced_bam2gff", help="Path to spliced_bam2gff executable (default: spliced_bam2gff)")
    parser.add_argument("--dry_run", action="store_true", help="Print commands without executing them")
    parser.add_argument("--cpus", type=int, default=1, help="Number of threads per job (where applicable)")
    return parser.parse_args()

def check_tools(umi_tools_path, spliced_bam2gff_path):
    """Check if required tools are available."""
    if subprocess.call(f"which {umi_tools_path}", shell=True, stdout=subprocess.DEVNULL) != 0:
        logger.warning(f"{umi_tools_path} not found. Ensure it is installed and activated.")
    
    if subprocess.call(f"which {spliced_bam2gff_path}", shell=True, stdout=subprocess.DEVNULL) != 0:
        logger.warning(f"{spliced_bam2gff_path} not found. Ensure it is installed or provide full path.")

def get_sample_id(filepath):
    """Derive sample ID from filename.
    Example: /path/to/sample1.tagged.bam -> sample1
    """
    basename = os.path.basename(filepath)
    # Remove common extensions
    if basename.endswith(".tagged.bam"):
        return basename.replace(".tagged.bam", "")
    elif basename.endswith(".bam"):
        return basename.replace(".bam", "")
    else:
        return os.path.splitext(basename)[0]

def process_sample(bam_file, args):
    """Run the pipeline for a single sample."""
    bam_file = os.path.abspath(bam_file)
    if not os.path.exists(bam_file):
        logger.error(f"Input file not found: {bam_file}")
        return False

    sample_id = get_sample_id(bam_file)
    sample_dir = os.path.join(args.output_dir, sample_id)
    
    if not args.dry_run:
        os.makedirs(sample_dir, exist_ok=True)
    
    logger.info(f"Processing sample: {sample_id}")

    # Define file paths
    dedup_bam = os.path.join(sample_dir, f"{sample_id}.dedup.bam")
    output_gff = os.path.join(sample_dir, f"{sample_id}.gff")
    log_file = os.path.join(sample_dir, f"{sample_id}.log")

    # 1. Run umi_tools
    # Command provided by user:
    # umi_tools dedup --stdin=*.tagged.bam --stdout=*.dedup.bam --extract-umi-method=tag 
    # --umi-tag=UB --cell-tag=CB --method=unique --per-gene --gene-tag=GE --per-cell --keep-tag=CB,UB,GE,TX
    umi_cmd = (
        f"{args.umi_tools_path} dedup "
        f"--stdin={bam_file} "
        f"--stdout={dedup_bam} "
        f"--extract-umi-method=tag "
        f"--umi-tag=UB "
        f"--cell-tag=CB "
        f"--method=unique "
        f"--per-gene "
        f"--gene-tag=GE "
        f"--per-cell "
        f"--keep-tag=CB,UB,GE,TX"
    )
    
    # 2. Run spliced_bam2gff
    # Command provided by user:
    # spliced_bam2gff -t 1000000 -M *.dedup.bam > *.gff
    gff_cmd = (
        f"{args.spliced_bam2gff_path} "
        f"-t 1000000 -M {dedup_bam} > {output_gff}"
    )

    # Execute
    cmds = [umi_cmd, gff_cmd]
    
    for cmd in cmds:
        logger.info(f"Running: {cmd}")
        if not args.dry_run:
            try:
                # Run command and redirect output to log
                with open(log_file, "a") as log:
                    # Write command to log for record
                    log.write(f"COMMAND: {cmd}\n")
                    log.flush()
                    
                    subprocess.run(
                        cmd, 
                        shell=True, 
                        check=True, 
                        stdout=log if ">" not in cmd else None, # Don't double redirect if cmd has >
                        stderr=subprocess.STDOUT
                    )
            except subprocess.CalledProcessError as e:
                logger.error(f"Command failed: {cmd}")
                logger.error(f"Check log file: {log_file}")
                return False
    
    logger.info(f"Successfully processed {sample_id}")
    return True

def main():
    args = parse_args()
    
    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)
        check_tools(args.umi_tools_path, args.spliced_bam2gff_path)

    samples_to_process = []

    # Handle Single BAM
    if args.bam:
        samples_to_process.append(args.bam)

    # Handle FOFN
    if args.fofn:
        try:
            with open(args.fofn, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        samples_to_process.append(line)
        except Exception as e:
            logger.error(f"Error reading FOFN: {e}")
            sys.exit(1)

    if not samples_to_process:
        logger.error("No input files found.")
        sys.exit(1)

    logger.info(f"Found {len(samples_to_process)} samples to process.")

    # Process each sample
    success_count = 0
    for bam_file in samples_to_process:
        if process_sample(bam_file, args):
            success_count += 1

    logger.info(f"Finished. Successfully processed {success_count}/{len(samples_to_process)} samples.")

if __name__ == "__main__":
    main()
