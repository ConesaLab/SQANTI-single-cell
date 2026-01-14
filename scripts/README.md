# SQANTI-sc Helper Scripts

This directory contains standalone scripts and utilities to help users prepare their data for SQANTI-sc.

## Available Scripts

*   `Generate_IsoQuant_Matrix.py`: Generates single-cell count matrices (Matrix Market format) from `spl-IsoQuant` output (`*.allinfo`). Supports filtering for known/consistent reads.
*   `make_pacbio_matrix.py`: Generates single-cell count matrices from PacBio Iso-Seq data (group file and BAM).
*   `process_ont_bam.py`: Batch processes ONT single-cell BAM files, performing UMI deduplication and GFF conversion.

## Usage

These scripts are designed to be run independently before the main SQANTI-sc pipeline.
