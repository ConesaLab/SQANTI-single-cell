#!/usr/bin/env python3
"""
make_pacbio_matrix.py

This script generates single-cell count matrices (Market Exchange Format) from PacBio Iso-Seq data.
It requires two inputs per sample:
1.  `*.group.txt`: Output from Iso-Seq collapse (maps isoforms to molecules).
2.  `*.bam`: Unmapped BAM (uBAM) containing deduplicated molecules with Cell Barcode (CB) tags.

The script outputs three files per sample:
-   features.tsv: List of isoforms ids.
-   barcodes.tsv: List of cell barcodes.
-   matrix.mtx: Count matrix.

Usage:
    python make_pacbio_matrix.py --sample_sheet samples.csv [-o ./results]

Sample Sheet Format (CSV):
    sample_id,group_file,bam_file
    sample1,/path/to/sample1.group.txt,/path/to/sample1.dedup.bam
"""

import argparse
import csv
import os
import subprocess
import sys
import logging
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Generate single-cell matrices from PacBio data.")
    parser.add_argument("--sample_sheet", required=True, help="Path to CSV file with columns: sample_id, group_file, bam_file")
    parser.add_argument("-o", "--output_dir", default=".", help="Directory to save outputs (default: current directory)")
    parser.add_argument("--samtools_path", default="samtools", help="Path to samtools executable (default: samtools)")
    return parser.parse_args()

def check_tools(samtools_path):
    """Check if samtools is available."""
    if subprocess.call(f"which {samtools_path}", shell=True, stdout=subprocess.DEVNULL) != 0:
        logger.error(f"{samtools_path} not found. Ensure it is installed or provide full path.")
        sys.exit(1)

def parse_group_file(group_file):
    """
    Parse the group file to map molecule IDs to Isoform IDs.
    Format:
    PB.1.1  molecule/1,molecule/2
    """
    molecule_to_isoform = {}
    isoforms = set()
    
    logger.info(f"Parsing group file: {group_file}")
    try:
        with open(group_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                
                isoform_id = parts[0]
                molecules = parts[1].split(',')
                
                isoforms.add(isoform_id)
                for mol in molecules:
                    molecule_to_isoform[mol] = isoform_id
                    
    except Exception as e:
        logger.error(f"Error reading group file {group_file}: {e}")
        return None, None

    logger.info(f"Found {len(isoforms)} isoforms and {len(molecule_to_isoform)} molecules.")
    return molecule_to_isoform, sorted(list(isoforms))

def parse_bam_and_count(bam_file, molecule_to_isoform, samtools_path):
    """
    Stream BAM file to find CB tags for molecules and count (Isoform, CB).
    """
    logger.info(f"Parsing BAM file: {bam_file}")
    
    # Store counts: counts[(isoform_idx, barcode_idx)] = count
    # We need indices for matrix output, so we'll collect raw data first then map to indices.
    # To save memory, we can map barcodes on the fly if we want, but let's collect all valid (Isoform, CB) pairs first.
    
    # Dictionary: barcode -> {isoform_id -> count}
    barcode_counts = defaultdict(lambda: defaultdict(int))
    
    cmd = [samtools_path, "view", bam_file]
    
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        
        for line in process.stdout:
            parts = line.split('\t')
            read_name = parts[0]
            
            # Check if this read is in our group file
            if read_name not in molecule_to_isoform:
                continue
            
            isoform_id = molecule_to_isoform[read_name]
            
            # Find CB tag
            cb = None
            for field in parts[11:]:
                if field.startswith("CB:Z:"):
                    cb = field[5:]
                    break
            
            if cb:
                barcode_counts[cb][isoform_id] += 1
                
        process.wait()
        if process.returncode != 0:
            logger.error(f"samtools view failed for {bam_file}")
            return None

    except Exception as e:
        logger.error(f"Error processing BAM {bam_file}: {e}")
        return None

    return barcode_counts

def write_outputs(sample_dir, isoforms, barcode_counts):
    """Write features.tsv, barcodes.tsv, and matrix.mtx."""
    
    # 1. Prepare Lists
    # Isoforms are already sorted from parse_group_file
    # Barcodes need to be sorted
    barcodes = sorted(list(barcode_counts.keys()))
    
    # Create mappings for indices (1-based for Matrix Market)
    isoform_to_idx = {iso: i+1 for i, iso in enumerate(isoforms)}
    barcode_to_idx = {bc: i+1 for i, bc in enumerate(barcodes)}
    
    # 2. Write features.tsv
    # Format: id \t name \t type
    features_path = os.path.join(sample_dir, "features.tsv")
    with open(features_path, 'w') as f:
        for iso in isoforms:
            f.write(f"{iso}\t{iso}\tGene Expression\n")
            
    # 3. Write barcodes.tsv
    barcodes_path = os.path.join(sample_dir, "barcodes.tsv")
    with open(barcodes_path, 'w') as f:
        for bc in barcodes:
            f.write(f"{bc}\n")
            
    # 4. Write matrix.mtx
    matrix_path = os.path.join(sample_dir, "matrix.mtx")
    
    # Collect all non-zero entries
    entries = []
    for bc, iso_counts in barcode_counts.items():
        bc_idx = barcode_to_idx[bc]
        for iso, count in iso_counts.items():
            iso_idx = isoform_to_idx[iso]
            entries.append((iso_idx, bc_idx, count))
            
    # Sort entries (optional but good practice: by col (barcode) then row (isoform))
    entries.sort(key=lambda x: (x[1], x[0]))
    
    with open(matrix_path, 'w') as f:
        # Header
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write('%metadata_json: {"software_version": "ont-single-cell","format_version": 2}\n')
        # Dimensions: rows (features), cols (barcodes), entries
        f.write(f"{len(isoforms)} {len(barcodes)} {len(entries)}\n")
        
        for iso_idx, bc_idx, count in entries:
            f.write(f"{iso_idx} {bc_idx} {count}\n")
            
    logger.info(f"Generated outputs in {sample_dir}")

def process_sample(sample_id, group_file, bam_file, args):
    sample_dir = os.path.join(args.output_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    
    logger.info(f"Processing sample: {sample_id}")
    
    # 1. Parse Group File
    molecule_to_isoform, isoforms = parse_group_file(group_file)
    if not molecule_to_isoform:
        return False
        
    # 2. Parse BAM and Count
    barcode_counts = parse_bam_and_count(bam_file, molecule_to_isoform, args.samtools_path)
    if barcode_counts is None:
        return False
        
    if not barcode_counts:
        logger.warning(f"No valid molecule-barcode pairs found for {sample_id}")
        return False
        
    # 3. Write Outputs
    write_outputs(sample_dir, isoforms, barcode_counts)
    return True

def main():
    args = parse_args()
    check_tools(args.samtools_path)
    os.makedirs(args.output_dir, exist_ok=True)

    # Read sample sheet
    try:
        with open(args.sample_sheet, 'r') as f:
            reader = csv.DictReader(f)
            required_cols = {'sample_id', 'group_file', 'bam_file'}
            if not required_cols.issubset(reader.fieldnames):
                logger.error(f"Sample sheet must contain columns: {required_cols}")
                sys.exit(1)
            
            samples = list(reader)
    except Exception as e:
        logger.error(f"Error reading sample sheet: {e}")
        sys.exit(1)

    logger.info(f"Found {len(samples)} samples to process.")

    success_count = 0
    for sample in samples:
        if process_sample(sample['sample_id'], sample['group_file'], sample['bam_file'], args):
            success_count += 1

    logger.info(f"Finished. Successfully processed {success_count}/{len(samples)} samples.")

if __name__ == "__main__":
    main()
