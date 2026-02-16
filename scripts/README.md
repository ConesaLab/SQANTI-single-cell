# SQANTI-sc Helper Scripts

This directory contains standalone scripts and utilities to help users prepare their data for SQANTI-sc.

## Available Scripts

*   `make_isoquant_matrix.py`: Generates single-cell count matrices (Matrix Market format) from `spl-IsoQuant` output (`*.allinfo`). Supports filtering for known/consistent reads.
*   `make_pacbio_matrix.py`: Generates single-cell count matrices from PacBio Iso-Seq data (group file and BAM).
*   `process_ont_bam.py`: Batch processes ONT single-cell BAM files, performing UMI deduplication and GFF conversion.
*   `run_STAR.py`: Automates STAR mapping of short-read RNA-seq data with SQANTI3-specific parameters. Supports single and paired-end reads, handling multiple samples/replicates.

## Usage
These scripts are designed to be run independently before the main SQANTI-sc pipeline.
