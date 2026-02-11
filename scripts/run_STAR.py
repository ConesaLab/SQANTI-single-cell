#!/usr/bin/env python3
"""
run_STAR.py

This script processes short-read RNA-seq data by mapping them to a reference genome using STAR.

Features:
- Handles both Single-End (SE) and Paired-End (PE) data
- Supports multiple samples/replicates via a sample sheet
- Auto-detects and building of STAR genome index if missing (provided reference FASTA)
- Automatic handling of compressed (.gz) files

Usage:
    # Basic usage
    python run_STAR.py --refFasta /path/to/genome.fa \\
                       --samples samples.txt \\
                       --outDir ./results

    # Usage with existing index
    python run_STAR.py --genomeDir /path/to/STAR_index \\
                       --samples samples.txt \\
                       --outDir ./results

Sample Sheet Format (Tab-delimited):
    SampleID    Read1               [Read2]
    SampleA     sampleA_R1.fq.gz    sampleA_R2.fq.gz
    SampleB     sampleB.fq.gz

Output:
    A subdirectory is created for each SampleID containing:
    - Aligned.sortedByCoord.out.bam (Main output)
    - Log.final.out
    - SJ.out.tab
"""

import os
import sys
import argparse
import subprocess

def run_command(cmd, dry_run=False):
    """
    Executes a shell command.
    """
    if dry_run:
        print(f"[DRY-RUN] {cmd}")
        return

    print(f"[CMD] {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {cmd}", file=sys.stderr)
        sys.exit(1)

def build_star_cmd(genome_dir, files, output_prefix, threads, compressed=False):
    """
    Constructs the STAR command with SQANTI3 parameters.
    """
    # Base STAR command with hardcoded SQANTI3 parameters
    base_command = [
        'STAR',
        '--runThreadN', str(threads),
        '--genomeDir', genome_dir,
        '--readFilesIn'] + files + [
        '--outFileNamePrefix', output_prefix,
        '--alignSJoverhangMin', '8',
        '--alignSJDBoverhangMin', '1',
        '--outFilterType', 'BySJout',
        '--outSAMunmapped', 'Within',
        '--outFilterMultimapNmax', '20',
        '--outFilterMismatchNoverLmax', '0.04',
        '--outFilterMismatchNmax', '999',
        '--alignIntronMin', '20',
        '--alignIntronMax', '1000000',
        '--alignMatesGapMax', '1000000',
        '--sjdbScore', '1',
        '--genomeLoad', 'NoSharedMemory',
        '--outSAMtype', 'BAM', 'SortedByCoordinate',
        '--twopassMode', 'Basic'
    ]

    if compressed:
        base_command.extend(['--readFilesCommand', 'zcat'])

    return ' '.join(base_command)

def build_star_index(genome_dir, ref_fasta, threads):
    """
    Builds the STAR genome index.
    """
    cmd = [
        'STAR',
        '--runMode', 'genomeGenerate',
        '--genomeDir', genome_dir,
        '--genomeFastaFiles', ref_fasta,
        '--runThreadN', str(threads)
    ]
    return ' '.join(cmd)

def main():
    parser = argparse.ArgumentParser(description="Map short reads using STAR with SQANTI3 parameters.")
    parser.add_argument("--genomeDir", help="Path to STAR genome index directory (optional, defaults to outDir/STAR_index).")
    parser.add_argument("--samples", required=True, help="Path to tab-delimited sample sheet (SampleID R1 [R2]).")
    parser.add_argument("--outDir", required=True, help="Output directory.")
    parser.add_argument("--refFasta", help="Path to reference FASTA (optional, used to build index if genomeDir missing).")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads (default: 8).")
    parser.add_argument("--dry_run", action="store_true", help="Print commands without executing.")
    args = parser.parse_args()

    # Determine genome directory
    if not args.genomeDir:
        args.genomeDir = os.path.join(args.outDir, "STAR_index")
    
    # Check/Build Genome Index
    if not os.path.exists(args.genomeDir):
        if args.refFasta:
            print(f"[INFO] Genome directory {args.genomeDir} not found. Building index from {args.refFasta}...")
            if not args.dry_run:
                os.makedirs(args.genomeDir, exist_ok=True)
            
            cmd = build_star_index(args.genomeDir, args.refFasta, args.threads)
            run_command(cmd, args.dry_run)
        else:
            print(f"[ERROR] Genome directory {args.genomeDir} not found and --refFasta not provided.", file=sys.stderr)
            sys.exit(1)

    # Process sample sheet
    with open(args.samples, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 2:
                print(f"[WARNING] Skipping invalid line: {line}", file=sys.stderr)
                continue

            # Skip header if present
            if parts[0] == "SampleID":
                continue

            sample_id = parts[0]
            read1 = parts[1]
            read2 = parts[2] if len(parts) > 2 else None

            # Determine if inputs are compressed
            compressed = read1.endswith(".gz")
            if read2 and (read2.endswith(".gz") != compressed):
                print(f"[ERROR] Sample {sample_id}: Mixed compressed/uncompressed reads not supported.", file=sys.stderr)
                sys.exit(1)

            # Prepare file list
            files = [read1]
            if read2:
                files.append(read2)

            # Set up sample output directory
            sample_out_dir = os.path.join(args.outDir, sample_id)
            if not sample_out_dir.endswith(os.sep):
                sample_out_dir += os.sep  # STAR needs trailing slash for prefix being a directory
            
            if not args.dry_run:
                os.makedirs(sample_out_dir, exist_ok=True)

            # Build and run command
            cmd = build_star_cmd(args.genomeDir, files, sample_out_dir, args.threads, compressed)
            print(f"Processing sample: {sample_id}")
            run_command(cmd, args.dry_run)

    print("All samples processed.")

if __name__ == "__main__":
    main()
