#!/usr/bin/env python3
"""
make_isoquant_matrix.py

This script generates single-cell count matrices (Matrix Market Format) from spl-IsoQuant output.
It processes `*.allinfo` files, which contain read-to-isoform assignments and their classification (known/novel/etc.).

By default, the script mimics spl-IsoQuant's strict quantification logic, counting only "known" (consistent) reads.
Users can optionally include "novel" (inconsistent) or "ambiguous" reads.

Usage:
    # Single file
    python make_isoquant_matrix.py --allinfo /path/to/sample.allinfo [-o ./results]

    # Batch processing (FOFN)
    python make_isoquant_matrix.py --fofn allinfo_file_list.txt [-o ./results]

FOFN Format (Text file):
    /path/to/sample1.allinfo
    /path/to/sample2.allinfo
    ...
"""

import argparse
import os
import sys
import logging
import scipy.io
from collections import defaultdict
from scipy.sparse import coo_matrix

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Generate single-cell matrices from spl-IsoQuant allinfo files.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--allinfo", help="Path to a single allinfo file to process")
    group.add_argument("--fofn", help="Path to a File of File Names (FOFN) containing one allinfo path per line")
    
    parser.add_argument("-o", "--output_dir", default=".", help="Directory to save outputs (default: current directory)")
    
    parser.add_argument("--include_novel", action="store_true", default=False,
                        help="If set, includes 'novel' (inconsistent) assignments. Default: False (strict).")
    parser.add_argument("--include_ambiguous", action="store_true", default=False,
                        help="If set, includes 'ambiguous' assignments. Default: False (strict).")
    
    return parser.parse_args()

def get_sample_id(filepath):
    """Derive sample ID from filename.
    Removes common extensions: .allinfo, .UMI_filtered.ED4.allinfo, etc.
    """
    basename = os.path.basename(filepath)
    # Naive cleanup: remove known extensions
    # Order matters: longest match first
    extensions = [".UMI_filtered.ED4.allinfo", ".allinfo"]
    for ext in extensions:
        if basename.endswith(ext):
            return basename.replace(ext, "")
    return os.path.splitext(basename)[0]

def process_sample(allinfo_file, args):
    """Process a single allinfo file and generate outputs."""
    allinfo_file = os.path.abspath(allinfo_file)
    if not os.path.exists(allinfo_file):
        logger.error(f"Input file not found: {allinfo_file}")
        return False
        
    sample_id = get_sample_id(allinfo_file)
    sample_dir = os.path.join(args.output_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    
    logger.info(f"Processing sample: {sample_id}")
    
    # Data structures
    # counts[(transcript_id, barcode)] = count
    counts = defaultdict(int)
    # features[transcript_id] = (gene_id, transcript_type)
    features_metadata = {}
    all_barcodes = set()
    
    try:
        with open(allinfo_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                
                # Columns based on spl-IsoQuant src/barcode_calling/umi_filtering.py
                # 0: read_id, 1: gene_id, 3: barcode, 9: read_type, 11: transcript_id, 12: transcript_type
                gene_id = parts[1]
                barcode = parts[3]
                read_type = parts[9]
                transcript_id = parts[11]
                transcript_type = parts[12] if len(parts) > 12 else "unknown" # Safety check

                # FILTERING LOGIC
                valid_read = False
                if read_type == 'known':
                    valid_read = True
                elif read_type == 'novel' and args.include_novel:
                    valid_read = True
                elif 'ambiguous' in read_type and args.include_ambiguous:
                    valid_read = True
                
                if not valid_read:
                    continue
                
                counts[(transcript_id, barcode)] += 1
                all_barcodes.add(barcode)
                features_metadata[transcript_id] = (gene_id, transcript_type)
                
    except Exception as e:
        logger.error(f"Error reading file {allinfo_file}: {e}")
        return False

    if not counts:
        logger.warning(f"No valid counts found for {sample_id}. Check if filtering criteria are too strict.")
        return False
        
    logger.info(f"Found {len(counts)} non-zero entries, {len(all_barcodes)} barcodes, {len(features_metadata)} features.")
    
    # Write outputs
    write_matrix_market(sample_dir, counts, features_metadata, all_barcodes)
    return True

def write_matrix_market(output_dir, counts, features_metadata, all_barcodes):
    """Writes standard Matrix Market files to output_dir."""
    
    # Sort for deterministic output
    sorted_barcodes = sorted(list(all_barcodes))
    sorted_features = sorted(list(features_metadata.keys()))
    
    # Map to 0-based indices
    barcode_to_idx = {bc: i for i, bc in enumerate(sorted_barcodes)}
    feature_to_idx = {ft: i for i, ft in enumerate(sorted_features)}
    
    data = []
    row_ind = []
    col_ind = []
    
    for (tid, bc), count in counts.items():
        if tid in feature_to_idx and bc in barcode_to_idx:
            row_ind.append(feature_to_idx[tid])
            col_ind.append(barcode_to_idx[bc])
            data.append(count)
    
    # Create matrix
    mat = coo_matrix((data, (row_ind, col_ind)), shape=(len(sorted_features), len(sorted_barcodes)))
    
    # 1. matrix.mtx
    matrix_path = os.path.join(output_dir, "matrix.mtx")
    scipy.io.mmwrite(matrix_path, mat)
    
    # 2. barcodes.tsv
    barcodes_path = os.path.join(output_dir, "barcodes.tsv")
    with open(barcodes_path, 'w') as f:
        for bc in sorted_barcodes:
            f.write(f"{bc}\n")
            
    # 3. features.tsv
    features_path = os.path.join(output_dir, "features.tsv")
    with open(features_path, 'w') as f:
        for tid in sorted_features:
            gene_id, t_type = features_metadata[tid]
            f.write(f"{tid}\t{gene_id}\t{t_type}\n")
            
    logger.info(f"Generated outputs in {output_dir}")

def main():
    args = parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    samples_to_process = []
    
    # Handle Single File
    if args.allinfo:
        samples_to_process.append(args.allinfo)
        
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
    
    success_count = 0
    for f in samples_to_process:
        if process_sample(f, args):
            success_count += 1
            
    logger.info(f"Finished. Successfully processed {success_count}/{len(samples_to_process)} samples.")

if __name__ == "__main__":
    main()
