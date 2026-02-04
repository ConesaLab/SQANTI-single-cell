<p align="center">
  <img src="img/Conesalab_logos_V3_singlecell.svg" width="400" title="SQANTI-sc logo">
</p>

# SQANTI-single cell

**SQANTI-single cell (SQANTI-sc)** is a pipeline for the structural and quality control of long-read single-cell transcriptomics datasets. It extends the capabilities of [SQANTI3](https://github.com/ConesaLab/SQANTI3) and [SQANTI-reads](https://github.com/ConesaLab/SQANTI3) frameworks to provide cell-level structural and quality control metrics. 


Table of Contents:
- [Prerequisites & Installation](#prerequisites--installation)
    - [0. Install Anaconda](#0-install-anaconda)
    - [1. Install SQANTI-sc](#1-install-sqanti-sc)
    - [2. Install SQANTI3](#2-install-sqanti3)
    - [3. Set up the Environment](#3-set-up-the-environment)
    - [4. Docker Support (Coming Soon)](#4-docker-support-coming-soon)
- [Getting Ready](#getting-ready)
- [Arguments and parameters in SQANTI-sc](#arguments-and-parameters-in-sqanti-sc)
- [Running SQANTI-sc](#running-sqanti-sc)
    - [1. Reads Mode](#1-reads-mode---mode-reads)
    - [2. Isoforms Mode](#2-isoforms-mode---mode-isoforms)
- [Providing orthogonal data to SQANTI-sc](#providing-orthogonal-data-to-sqanti-sc)
- [Cell clustering](#clustering)
- [Understanding the output of SQANTI-sc](#understanding-the-output-of-sqanti-sc)
    - [1. SQANTI3-based Outputs](#1-sqanti3-based-outputs)
    - [2. SQANTI-reads-based Outputs](#2-sqanti-reads-based-outputs)
    - [3. SQANTI-sc Specific Outputs](#3-sqanti-sc-specific-outputs)

<a name="prerequisites--installation"></a>

## Prerequisites & Installation


<a name="0-install-anaconda"></a>

### 0. Install Anaconda
Make sure you have installed Anaconda. If not, you can download the generic installer for Linux [here](http://docs.continuum.io/anaconda/install/#linux-install).

<a name="1-install-sqanti-sc"></a>

### 1. Install SQANTI-sc
You can install SQANTI-sc by downloading the source code or cloning the repository.

**Option A: Download Source Code (Recommended for General Users)**
We recommend this option for general users who want to use the stable version of the tool.
Download the latest release from the [Releases page](https://github.com/ConesaLab/SQANTI-sc/releases) (if available) or download the repository as a ZIP file.
```bash
wget https://github.com/ConesaLab/SQANTI-sc/archive/refs/heads/main.zip
unzip main.zip
mv SQANTI-sc-main SQANTI-sc
cd SQANTI-sc
```

**Option B: Clone Repository (For Developers)**
If you intend to contribute to the development of SQANTI-sc, please clone the repository. This option sets up a git repository and is NOT recommended for general users unless you plan to submit pull requests or track the latest development changes.
```bash
git clone https://github.com/ConesaLab/SQANTI-sc.git
cd SQANTI-sc
```

<a name="2-install-sqanti3"></a>

### 2. Install SQANTI3
SQANTI-sc requires a functional installation of SQANTI3.
Please follow the [SQANTI3 Installation Instructions](https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation) to install SQANTI3.

**Important: SQANTI3 Location.**
SQANTI-sc needs to know where SQANTI3 is installed. You have two options:

*  **Option A (Simpler):** Place the `SQANTI3` folder inside the `SQANTI-sc` directory.
*  **Option B (Flexible):** Install SQANTI3 anywhere and set the `SQANTI3_DIR` environment variable:
```bash
export SQANTI3_DIR=/path/to/your/SQANTI3/directory
```

<a name="3-set-up-the-environment"></a>

### 3. Set up the Environment
We provide a unified Conda environment file `SQANTI-sc_env.yml` that includes all dependencies for both SQANTI3 and SQANTI-sc (including R packages for reporting).

```bash
conda env create -f SQANTI-sc_env.yml
conda activate SQANTI-sc_env
```

<a name="4-docker-support-coming-soon"></a>

### 4. Docker Support (Coming Soon)
Future releases of SQANTI-sc will be containerized and available on DockerHub. Currently, please use the Conda installation method.

<a name="getting-ready"></a>

## Getting Ready

Activate the SQANTI-sc conda environment:
```bash
conda activate SQANTI-sc_env
```

<a name="arguments-and-parameters-in-sqanti-sc"></a>

## Arguments and parameters in SQANTI-sc

The SQANTI-sc quality control script accepts the following arguments:

```bash
usage: sqanti_sc.py [-h] --mode {reads,isoforms} --design DESIGN --refGTF REFGTF --refFasta REFFASTA 
                    [--out_dir OUT_DIR] [--input_dir INPUT_DIR] [--report {pdf,html,both,skip}] 
                    [-@ SAMTOOLS_CPUS] [--verbose] [--skip_hash] [--ignore_cell_summary]
                    [--multisample_report] [--multisample_report_prefix PREFIX]
                    [--write_per_cell_outputs] [--run_clustering]
                    [--normalization {log1p,sqrt,pearson}] [--n_neighbors N] [--n_pc N]
                    [--resolution RES] [--n_top_genes N] [--clustering_method {leiden,louvain,kmeans}]
                    [--n_clusters N]
                    [--min_ref_len MIN_REF_LEN] [--force_id_ignore] [--genename] 
                    [--novel_gene_prefix NOVEL_GENE_PREFIX] [--ref_cov_min_pct PCT]
                    [--short_reads SHORT_READS] [--SR_bam SR_BAM] 
                    [--aligner_choice {minimap2,uLTRA,gmap,deSALT}] [-x GMAP_INDEX]
                    [--skipORF] [--orf_input ORF_INPUT]
                    [--CAGE_peak CAGE_PEAK] [--polyA_motif_list POLYA_MOTIF_LIST] 
                    [--polyA_peak POLYA_PEAK] [--phyloP_bed PHYLOP_BED] 
                    [-e EXPRESSION] [-c COVERAGE] 
                    [--isoAnnotLite] [--gff3 GFF3]
                    [--isoform_hits] [--ratio_TSS_metric {max,mean,median,3quartile}] 
                    [-t CPUS] [-n CHUNKS] [-l {ERROR,WARNING,INFO,DEBUG}] 
                    [--is_fusion] [-v]
```

<details>
  <summary>Arguments explanation</summary>

```text
Structural and Quality Analysis of Single-Cell Isoforms

options:
  -h, --help            show this help message and exit

Required arguments:
  --refFasta REFFASTA   Reference genome (Fasta format)
  --refGTF REFGTF       Reference annotation file (GTF format)
  --design DESIGN, -de DESIGN
                        Design file (CSV) containing sampleID, file_acc, and cell/abundance metadata
  --mode {reads,isoforms}, -m {reads,isoforms}
                        Input data mode: 'reads' for long reads, 'isoforms' for collapsed transcripts

SQANTI-sc specific options:
  --out_dir OUT_DIR, -d OUT_DIR
                        Output directory (Default: current dir)
  --input_dir INPUT_DIR, -i INPUT_DIR
                        Input directory for files referenced in design (Default: current dir)
  --report {pdf,html,both,skip}
                        Report format (Default: pdf)
  --samtools_cpus CPUS, -@ CPUS
                        Threads for samtools sorting/indexing (Default: 10)
  --verbose             Print detailed logs during execution
  --skip_hash           Skip UJC hashing calculation step
  --ignore_cell_summary Don't save cell summary table in report
  --multisample_report  Generate a multisample cohort report
  --multisample_report_prefix PREFIX
                        Prefix for multisample report (Default: SQANTI_sc_multisample_report)
  --write_per_cell_outputs
                        Write per-cell gene/UJC counts and CV matrices

SQANTI-sc Clustering and UMAP options:
  --run_clustering      Run cell clustering and UMAP analysis
  --normalization {log1p,sqrt,pearson}
                        Normalization method (Default: log1p)
  --n_neighbors N       Number of neighbors for UMAP (Default: 15)
  --n_pc N              Number of principal components (Default: 30)
  --resolution RES      Resolution for Leiden clustering (Default: 0.5)
  --n_top_genes N       Number of highly variable genes (Default: 2000)
  --clustering_method {leiden,louvain,kmeans}
                        Clustering method (Default: leiden)
  --n_clusters N        Number of clusters for K-means (Default: 10)

SQANTI3 Customization and filtering:
  --min_ref_len MIN_REF_LEN
                        Minimum reference transcript length (default: 0 bp)
  --force_id_ignore     Allow the usage of transcript IDs non related with PacBio's nomenclature
  --genename            Use gene_name tag from GTF to define genes. Default: gene_id
  --novel_gene_prefix PREFIX
                        Prefix for novel isoforms
  --ref_cov_min_pct PCT
                        Minimum % of reference transcript length a read must cover (Default: 45.0)

Aligner and mapping options:
  --aligner_choice {minimap2,uLTRA,gmap,deSALT}
                        Select aligner for FASTA input (default: minimap2)
  -x GMAP_INDEX, --gmap_index GMAP_INDEX
                        Path to gmap_build index. Mandatory if using GMAP.

ORF prediction:
  --skipORF             Skip ORF prediction
  --orf_input ORF_INPUT Input fasta to run ORF on.

SQANTI3 Orthogonal data inputs:
  --short_reads SHORT_READS
                        FOFN of short-read RNA-Seq files (FASTA/FASTQ) for validation
  --SR_bam SR_BAM       Directory or FOFN of short-read BAM files for validation
  --CAGE_peak CAGE_PEAK FANTOM5 Cage Peak (BED format)
  --polyA_motif_list LIST
                        Ranked list of polyA motifs (text)
  --polyA_peak POLYA_PEAK
                        PolyA Peak (BED format)
  --phyloP_bed PHYLOP_BED
                        PhyloP BED for conservation scores
  -e EXPRESSION, --expression EXPRESSION
                        Expression matrix (e.g. Kallisto tsv)
  -c COVERAGE, --coverage COVERAGE
                        Junction coverage files (comma-separated or pattern)

Functional annotation:
  --isoAnnotLite        Run isoAnnot Lite for tappAS compatibility
  --gff3 GFF3           Precomputed tappAS species specific GFF3 file

Output & Performance options:
  --isoform_hits        Report all FSM/ISM isoform hits in a separate file
  --ratio_TSS_metric {max,mean,median,3quartile}
                        Metric for ratio_TSS column (default: max)
  -t CPUS, --mapping_cpus CPUS
                        Number of threads used during alignment (default: 10)
  -n CHUNKS, --chunks CHUNKS
                        Number of chunks to split analysis (default: 10)
  -l {ERROR,WARNING,INFO,DEBUG}, --log_level LEVEL
                        Set the logging level (Default: INFO)

Optional arguments:
  --is_fusion           Input are fusion isoforms (GTF required)
  -v, --version         Display program version number
```
</details>

<a name="running-sqanti-sc"></a>

## Running SQANTI-sc

The tool operates in two main modes: **`reads`** and **`isoforms`**.

<a name="1-reads-mode---mode-reads"></a>

### 1. Reads Mode (`--mode reads`)
The **Reads Mode** is an adaptation of the **[SQANTI-reads](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI-reads)** tool specifically designed for single-cell data. It provides a comprehensive structural and quality control assessment of long-read single-cell RNA-seq data **at the read level**. This step can be useful for validating your data structure and quality before committing to the more complex steps of isoform identification and quantification.

**Key Applications:**
*   **Pre-analysis QC**: Provides actionable insights into the quality and structural properties of the library, such as splicing accuracy and the presence of technical artifacts. This information allows users to make informed decisions about downstream isoform identification strategies.
*   **Multisample Benchmarking**: Adapts core SQANTI-reads metrics for single-cell resolution. The full suite of quality metrics enables robust comparisons of library complexity and technical variability across different experimental conditions, sequencing platforms, preprocessing pipelines, etc. 

In this mode, SQANTI-sc treats each individual read as a "transcript" instance.

#### Compatible Pipelines
*   **PacBio Iso-Seq**: The [Iso-Seq single cell pipeline](https://isoseq.how/umi/cli-workflow.html) produces deduplicated unaligned reads in BAM format.
*   **Oxford Nanopore**: The [EPI2ME wf-single-cell pipeline](https://epi2me.nanoporetech.com/epi2me-docs/workflows/wf-single-cell/) producess undeduplicated aligned reads in BAM format. To use these reads as input for SQANTI-sc, first they must be deduplicated and converted to GTF/GFF. You can find a utility script named `process_ont_bam.py` that performs these steps for you automatically in the scripts/ folder.
*   **Any long-read pipeline** producing deduplicated unaligned reads in BAM or FASTA/FASTQ format; or deduplicated aligned reads in GTF/GFF format.

#### Reads Input Formats
*   **uBAM**: Unaligned reads. **Note**: SQANTI-sc will automatically convert BAM files to FASTQ for mapping and processing.
*   **FASTA/FASTQ**: Unaligned reads. If reads are provided in FASTA or FASTQ format, SQANTI-sc will automatically add the (`-fasta`) option to the SQANTI3 runner and perform a mapping step. Check the [SQANTI3 documentation](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control) for mapping options. 
*   **GTF/GFF**: Transcript annotations. We recommend providing reads in GTF/GFF format, as it enables better control over mapping parameters and faster runtimes. 

#### Minimal Input (Mandatory)
SQANTI-sc requires a specific set of input files to run. The main entry point is a **Design File** (CSV) that maps each sample to its corresponding reads and cell barcode information.

*   **Design File (`--design`)**
A comma-separated values (CSV) file containing the metadata for your samples.

| Column | Description |
| :--- | :--- |
| `sampleID` | **Required**. A unique identifier for each of the samples (e.g., `Sample1`) as we want them to be represented in the output. |
| `file_acc` | **Required**. The file prefix used to locate your input reads files. SQANTI-sc will search the `--input_dir` for files starting with this prefix (e.g., `PB_S1` matches `PB_S1.bam`, `PB_S1.fastq`, or `PB_S1.gtf`). These will also be the names of the output directories where the output files corresponding to each sample will be located.|
| `cell_association` | **Required**. Path to the file linking reads to cell barcodes for each sample. <br> There are 2 options for this file: <br> 1. A **TSV file** mapping Read IDs to Cell Barcodes. The minimum columns required are `id` (identifier of the read) and `cb` (cell barcode of the read).  <br> 2. A **uBAM/BAM file** containing `CB` (Cell Barcode) and `XM`/`UB` (UMI) tags. Easier to give if your input reads are already in uBAM format (the uBAM file can be the same as the one used as input reads).|

**Example `design_reads.csv`:**
```csv
sampleID,file_acc,cell_association
Sample1,PB_S1,/data/Sample1/barcodes.tsv
Sample2,PB_S2,/data/Sample2/sample2.bam
```

**Example `cell_association` file (barcodes.tsv):**
```tsv
id	cb	umi
m64012_250421_000242/120719489/ccs/10460_11196	AAACCCAAGGTTCCTA	TATGCCCGGTAT
m64012_250421_000242/17565024/ccs/13918_14203	GGGTTTAAGGTTCCTA	GCGCGCAATTCA
```
*Note: The `umi` column is optional. If included, the UMIs will be added to the reads classification under the column `UMI`.*

*   **Reference Genome (`--refFasta`)**: FASTA format (e.g., `hg38.fa`). The Chromosome/scaffold names must exactly match those in the reference annotation.
*   **Reference Annotation (`--refGTF`)**: GTF format (e.g., `gencode.v38.annotation.gtf`). Used to classify transcripts and assess novelty. Make sure that it matches the reference genome's coordinate system. You can find reference transcriptomes for different species in [GENCODE](https://www.gencodegenes.org) or [CHESS](https://ccb.jhu.edu/chess/).

---

<a name="2-isoforms-mode---mode-isoforms"></a>

### 2. Isoforms Mode (`--mode isoforms`)
The **Isoforms Mode** is designed for the in-depth characterization of unique transcript isoforms across single cells. Unlike Reads Mode, which focuses on individual reads, this mode takes collapsed or assembled consensus transcript models as input and performs quality control at the single-cell level.

**Key Applications:**
*   **Single-Cell Isoform Characterization**: Analyzes isoform expression across cell populations using user-provided abundance data, enabling the exploration of alternative splicing patterns and isoform usage heterogeneity.
*   **Structural Classification**: Classifies each isoform against the reference annotation following the standard SQANTI3 classification scheme, providing a detailed quality assessment of the transcriptome assembly.

#### Compatible Inputs (WIP)
*   **Iso-Seq Collapse**: Output from PacBio's [Iso-Seq Collapse](https://isoseq.how/classification/isoseq-collapse.html). To generate the isoform count matrices from the outputs provided by this tool you can use the script named `make_pacbio_matrix.py`, located in the scripts/ folder.
*   **Spl-IsoQuant**: Output from [spl-IsoQuant](https://github.com/algbio/spl-IsoQuant), IsoQuant version developed for upstream analysis of single-cell and spatial long-read data analysis. To generate the isoform count matrices from the outputs provided by this tool you can use the script named `make_isoquant_matrix.py`, located in the scripts/ folder.
*   **Any other tool** producing a GTF/GFF or FASTA representing specific transcript isoforms and the corresponding quantification of said isoforms in each cell barcode. 

#### Isoforms Input Formats
*   **GTF (`*.gtf`)**: **Recommended**. Transcript annotations. Like SQANTI3, this is the default and recommended format that SQANTI-sc expects. 
*   **FASTA/FASTQ (`*.fasta`, `*.fastq`)**: Transcript sequences. If provided, SQANTI-sc will map them to the reference genome (like in the Reads Mode).

#### Minimal Input (Mandatory)
Follows similar inputs to Reads mode. 

*   **Design File (`--design`)**
A comma-separated values (CSV) file containing the metadata for your samples.

| Column | Description |
| :--- | :--- |
| `sampleID` | **Required**. A unique identifier for each of the samples (e.g., `Sample1`) as we want them to be represented in the output. |
| `file_acc` | **Required**. The file prefix used to locate your input transcript models files (matches `{file_acc}.gtf` or `{file_acc}.fasta`). These will also be the names of the output directories where the output files corresponding to each sample will be located. |
| `cell_association` | **Conditional**. Path to the file linking Isoform IDs to Cell Barcodes (TSV) *if no abundance matrix is present*. |
| `abundance` | **Conditional**. Path to a folder containing quantification data in **Market Exchange (MEX) format**. The folder **MUST** contain three files: `matrix.mtx`, `features.tsv` (with the transcript model IDs used in the input), and `barcodes.tsv`.|
| `coverage` | **Optional**. Path to the STAR splice junction output file (`SJ.out.tab`). Used to validate splice junctions. |
| `SR_bam` | **Optional**. Path to a sorted and indexed BAM file of short reads. Used to validate TSS. |

*Note*: You can provide either a `cell_association` file or `abundance` directory as your cell barcode-isoform association file. If you want to perform quality control with the quantification of the expression of the isoforms (recommended) you will need the count matrix, but is not mandatory to run the isoforms mode.

**Example `design_isoforms.csv`:**
```csv
sampleID,file_acc,abundance,coverage,SR_bam
Sample1,Iso_S1,/data/S1_counts/,/path/to/S1_SJ.out.tab,/path/to/S1_sorted.bam
Sample2,Iso_S2,/data/S2_counts/,/path/to/S2_SJ.out.tab,/path/to/S2_sorted.bam
```

*   **Reference Files**: Same as Reads Mode (`--refFasta`, `--refGTF`).

#### SQANTI-sc usage Example
```bash
python sqanti_sc.py \
    --mode isoforms \
    --design design_isoforms.csv \
    --refFasta /genomes/human/hg38.fa \
    --refGTF /genomes/human/hg38.gtf \
    --input_dir ./isoforms \
    --out_dir ./results \
    --CAGE_peak ./SQANTI3/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed \
    --polyA_motif_list ./SQANTI3/data/polyA_motifs/mouse_and_human.polyA_motif.txt \
    --report both \
    --run_clustering \
    --multisample_report
```
<a name="providing-orthogonal-data-to-sqanti-sc"></a>

## Providing orthogonal data to SQANTI-sc


SQANTI-sc accepts orthogonal data to assist in the quality control and filtering of artifactual transcript models. 
* CAGE peak data (`--CAGE_peak`) for Transcription Start Site (TSS) validation.
* PolyA information (`--polyA_motif_list`, `--polyA_peak`) for Transcription Termination Site (TTS) validation. 
* Short Reads data (Splice Junctions and BAMs) provided via the **Design File**.

### Short Reads Validation

SQANTI-sc supports orthogonal validation of Splice Junctions and TSS using short reads. However, unlike SQANTI3, it does **not** accept raw FASTQ inputs (`--short_reads`). Instead, users must perform the short-read alignment externally (e.g., using [STAR](https://github.com/alexdobin/STAR)) and provide the resulting files via new columns in the **Design File**.

This strategy treats short-read data as a "bulk proxy" to validate the single-cell long-read isoforms. We recommend using deeper bulk RNA-seq data from the same tissue/condition to validate the single-cell library. To learn more about the metrics related to short-reads coverage, visit [SQANTI3 documentation](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control#SR).

**Design File Columns for Short Reads:**

| Column | Description |
| :--- | :--- |
| `coverage` | **Optional**. Path to the STAR splice junction output file (usually `SJ.out.tab`). Used to validate splice junctions. Equivalent to `--coverage` flag from SQANTI3. |
| `SR_bam` | **Optional**. Path to a sorted and indexed BAM file of short reads. Used to calculate the TSS ratio and validate 5' ends. Equivalent to `--SR_bam` flag from SQANTI3.|

> NOTE:
> **Why matching 10x Short Reads are not supported?**
>
> SQANTI-sc relies on "bulk-like" short read coverage to validate splice junctions across the full length of transcripts. Common single-cell short-read application (e.g., 10x Genomics) typically generate end-biased reads that do not cover the full transcript body, making them unsuitable for validating internal splice junctions.
> 
> Therefore, we recommend using **Bulk RNA-seq** samples (condition or tissue-comparable) for this validation step.


<a name="clustering"></a>

## Cell clustering (`--run_clustering`)

SQANTI-sc includes an optional step to perform cell clustering based on the expression of genes. This step is useful for exploring the heterogeneity of the cell population. SQANTI-sc uses the [Scanpy](https://scanpy.readthedocs.io/en/stable/) toolkit for cell clustering. 

To enable this step, you must use the flag `--run_clustering`.

### Customization options
The clustering process can be customized using the following arguments:

*   **`--normalization`**: Normalization method to apply to the count matrix before clustering. Options: `log1p` (default, log(CPM+1)), `sqrt` (square root), `pearson` (Pearson residuals).
*   **`--n_neighbors`**: Number of neighbors for UMAP construction (Default: 15).
*   **`--n_pc`**: Number of principal components to use for clustering (Default: 30).
*   **`--resolution`**: Resolution parameter for Leiden clustering (Default: 0.5). Higher values result in more clusters.
*   **`--n_top_genes`**: Number of highly variable genes to use for dimensionality reduction (Default: 2000).
*   **`--clustering_method`**: Algorithm to use for clustering. Options: `leiden` (default), `louvain`, `kmeans`.
*   **`--n_clusters`**: Number of clusters to force if using K-means clustering (Default: 10).


<a name="understanding-the-output-of-sqanti-sc"></a>

## Understanding the output of SQANTI-sc

The majority of the outputs follow the same logic and structure as [SQANTI3](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3) and [SQANTI-reads](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI-reads), with modifications to include single-cell information and new files specific to this pipeline.

<a name="1-sqanti3-based-outputs"></a>

### 1. SQANTI3-based Outputs
Standard SQANTI3 output files are generated for each sample, but they include additional columns to track cell barcode of origin.

*   **`*_classification.txt`**: The main output file containing structural classification and quality attributes.
    *   New columns:
        *   `CB`: Cell Barcode associated with the read/isoform. In isoforms mode, this column is a comma-separated list of the cell barcodes in which the isoform appears.
        *   `UMI`: Unique Molecular Identifier (**Reads Mode only**, if UMI information available (through cell association file)).
        *   `jxn_string` / `jxnHash` (**Reads Mode only**): Unique Junction Chain (UJC) representation of the transcript model structure. These columns, introduced by [SQANTI-reads](https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI-reads), allow for read grouping based on splice junctions structure. This step can be skipped using the `--skip_hash` option.
    * Changed columns:
        *   `FL`: Full-length count. In reads mode, this column will always diplay values of 1, as each row of the classification is suppossed to represent a unique UMI. In isoforms mode it will also display values of 1 except if quantification information is provided, in which case the column will be a comma-separated list of numeric values that represent the number of counts of that transcript model in each of the cells in the `CB` column, following the same order.
*   **`*_junctions.txt`**: File containing all splice junctions identified.
    *   New columns:
        *   `CB`: Cell Barcode associated with the junction. In isoforms mode, this column can contain multiple comma-separated cell barcodes if the junction pertains to a transcript model present in multiple cell barcodes. 

<a name="2-sqanti-reads-based-outputs"></a>

### 2. SQANTI-reads-based Outputs
The majority of SQANTI-reads-specific otuputs are not output by SQANTI-sc, with the exception of some tables that add the cell barcode dimension. These per-cell matrices follow the logic of **SQANTI-reads** and are generated optionally if the `--write_per_cell_outputs` flag is used. They provide detailed metrics at the single-cell level.

*   **`*_gene_counts.csv`**: Counts of reads/isoforms per gene per cell, broken down by structural category.
*   **`*_ujc_counts.csv`**: Counts of Unique Junction Chains (UJCs) per cell, including their structural classification and novelty status.
*   **`*_cv.csv`**: Coefficient of Variation (CV) metrics for splice junctions per cell, useful for identifying splicing variability.

<a name="3-sqanti-sc-specific-outputs"></a>

### 3. SQANTI-sc Specific Outputs

*   **`*_report.html` / `*.pdf`**: Comprehensive quality control reports summarizing the data at the single-cell level.
*   **`clustering/umap_results.csv`**: (If `--run_clustering` is active) Contains the UMAP coordinates and cluster assignments for each cell barcode.
*   **`*_SQANTI_cell_summary.txt.gz`**: A GZIP-compressed tab-delimited file containing a wide array of quality control metrics aggregated per cell. This is the core file for downstream analysis of cellular transcriptome quality.

#### Glossary of Cell Summary columns

The output `_SQANTI_cell_summary.txt.gz` has the following possible fields:

* **`CB`** : Cell Barcode identifier.  
* **`Reads_in_cell`** / **`Transcripts_in_cell`** : Total number of reads (Reads Mode) or transcripts (Isoforms Mode) associated with the cell.  
* **`UMIs_in_cell`** : Total number of unique Molecular Identifiers (UMIs) detected (Reads Mode only) in the cell.  
* **`total_reads_no_monoexon`** / **`total_transcripts_no_monoexon`** : Count of reads/transcripts excluding mono-exons.  
* **`FSM`**, **`ISM`**, **`NIC`**, **`NNC`**, **`Genic_Genomic`**, **`Antisense`**, **`Fusion`**, **`Intergenic`**, **`Genic_intron`** : Absolute counts of reads/transcripts belonging to each SQANTI3 structural category.
* **`[Category]_prop`** (e.g., `FSM_prop`, `NIC_prop`) : The proportion of reads/transcripts in the cell belonging to each structural category.
* **`Genes_in_cell`** : Number of unique genes detected in the cell.
* **`UJCs_in_cell`** : Number of Unique Junction Chains (UJCs) detected (multi-exon only).
* **`MT_reads_count`** : Number of reads/transcripts mapping to mitochondrial genes.  
* **`MT_perc`** : Percentage of reads/transcripts mapping to mitochondrial genes.  
* **`Annotated_genes`** : Number of known (annotated) genes detected. 
* **`Novel_genes`** : Number of novel genes detected.
* **`*_junctions`** (e.g., `Known_canonical_junctions`, `Novel_canonical_junctions`, `Novel_non_canonical_junctions`, `Known_non_canonical_junctions`) : Counts of splice junctions by type.  
* **`total_junctions`** : Total number of splice junctions identified in the cell.  
* **`*_junctions_prop`** : Proportions of each junction type relative to total junctions.  
* **`[Category]_[Subcategory]_prop`** : Proportion of reads/transcripts within a category that belong to a specific subcategory (e.g., `FSM_alternative_3end_prop`, `ISM_intron_retention_prop`, `Genic_mono_exon_prop`).  
* **`anno_bin*_perc`** : Proportion of annotated genes with expression levels (reads or transcript model counts) in specific bins (`bin1`: 1 count, `bin2_4`: 2-4 counts, `bin5_9`: 5-9 counts, `bin10plus`: >=10 counts).  
* **`novel_bin*_perc`** : Proportion of novel genes with expression levels in specific bins (same as above).  
* **`*_ujc_bin*_perc`** : Similar bins calculated for Unique Junction Chains (UJCs) (Reads Mode only): `bin1` (1 count), `bin2_3` (2-3 counts), `bin4_5` (4-5 counts), `bin6plus` (>=6 counts).  
* **`[Total|Category]_*_length_[mono]_prop`** : Proportion of reads falling into specific length bins (Total and per-category, including mono-exonic specific). Length bins: `<250bp`, `250-500bp`, `500-1000bp` (`short`), `1000-2000bp` (`mid`), `>2000bp` (`long`).  
* **`*_ref_coverage_prop`** : Proportion of reads/transcripts in each category covering at least `ref_cov_min_pct` of the reference transcript length.  
* **`ref_cov_min_pct`** : The minimum reference coverage percentage used for the above metric (parameter value).  
* **`RTS_prop_in_cell`** : Proportion of reads flagged as Reverse Transcriptase Switching (RTS) artifacts.  
* **`[Category]_RTS_prop`** : Proportion of RTS artifacts within each structural category.  
* **`Non_canonical_prop_in_cell`** : Proportion of reads with non-canonical splicing.  
* **`[Category]_noncanon_prop`** : Proportion of non-canonical splicing within each structural category.  
* **`Intrapriming_prop_in_cell`** : Proportion of reads flagged as intra-priming artifacts.  
* **`[Category]_intrapriming_prop`** : Proportion of intra-priming artifacts within each structural category.  
* **`TSSAnnotationSupport_prop`** : Proportion of reads where the TSS is within 50bp of an annotated TSS.  
* **`[Category]_TSSAnnotationSupport`** : Proportion of TSS support within each structural category.  
* **`Annotated_genes_prop_in_cell`** : Proportion of genes detected that are annotated.  
* **`Annotated_juction_strings_prop_in_cell`** : Proportion of UJCs that match annotated junction combinations.  
* **`Canonical_prop_in_cell`** : Proportion of reads/transcripts with canonical splicing.  
* **`[Category]_canon_prop`** : Proportion of canonical splicing within each structural category.  
* **`NMD_prop_in_cell`** : Proportion of reads predicted to be NMD candidates (if ORF prediction is on).  
* **`[Category]_NMD_prop`** : Proportion of NMD candidates within each structural category.  
* **`[Category]_[non]_coding_prop`** : Proportion of coding vs non-coding transcripts within each category.  
* **`CAGE_peak_support_prop`** : Proportion of reads supported by CAGE peaks (if provided).  
* **`[Category]_CAGE_peak_support_prop`** : Proportion of CAGE support within each structural category.  
* **`PolyA_motif_support_prop`** : Proportion of reads with identified PolyA motifs (if provided).  
* **`[Category]_PolyA_motif_support_prop`** : Proportion of PolyA support within each structural category.  
