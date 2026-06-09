#!/usr/bin/env Rscript
# make_recount3_sj.R
#
# Download junction data from recount3 and convert to STAR SJ.out.tab format
# for use as --coverage input to SQANTI-sc.
#
# Junction coverage (column 7 of SJ.out.tab) is the TOTAL number of reads
# supporting that junction, summed across all samples -- standard STAR semantics,
# so SQANTI3's --min_cov keeps its usual meaning ("minimum read coverage of the
# weakest junction in an isoform"). The orthogonal SAMPLE support (in how many
# samples the junction appears, and that as a % of the study) is written to a
# companion "<output>.sample_support.tsv" sidecar, which SQANTI-sc folds into the
# classification as the min_sample_cov / min_sample_pct columns.
#
# Human data uses GTEx (tissue names). Mouse data uses SRA project IDs.
# For multi-tissue SRA projects, use --filter_col / --filter_val to restrict to
# the samples of interest before summing; use --list_metadata to discover which
# colData columns and values are available for a given project.
#
# Usage:
#   Human GTEx:      Rscript make_recount3_sj.R --organism human --tissue BLOOD \
#                        --output blood_gtex.SJ.out.tab
#   Mouse SRA:       Rscript make_recount3_sj.R --organism mouse --project SRP123456 \
#                        --output mouse_brain.SJ.out.tab
#   Filter samples:  Rscript make_recount3_sj.R --organism mouse --project SRP199494 \
#                        --filter_col sra.sample_attributes --filter_val Heart \
#                        --output heart_tms.SJ.out.tab
#   List projects:   Rscript make_recount3_sj.R --list_tissues [--organism mouse]
#   List metadata:   Rscript make_recount3_sj.R --organism mouse --project SRP199494 \
#                        --list_metadata
#
# Requirements:
#   install.packages("BiocManager")
#   BiocManager::install(c("recount3", "optparse"))

suppressPackageStartupMessages({
    library(optparse)
    library(recount3)
    library(Matrix)
})

option_list <- list(
    make_option("--organism", type = "character", default = "human",
                help = "Organism: 'human' or 'mouse' [default: human]"),
    make_option("--tissue", type = "character", default = NULL,
                help = "GTEx tissue name for human (e.g. BLOOD). Use --list_tissues to see all options."),
    make_option("--project", type = "character", default = NULL,
                help = "SRA project ID for mouse or human SRA data (e.g. SRP123456)."),
    make_option("--output", type = "character", default = NULL,
                help = paste("Output file [default: <tissue/project>_recount3.SJ.out.tab].",
                             "Written uncompressed (STAR SJ.out.tab format) so it can be",
                             "passed straight to SQANTI3 --coverage; a .gz extension",
                             "compresses it but must be decompressed before use.")),
    make_option("--min_samples", type = "integer", default = 1,
                help = "Min samples a junction must appear in to be kept [default: 1]"),
    make_option("--list_tissues", action = "store_true", default = FALSE,
                help = "List available GTEx tissues (human) or SRA projects (mouse) and exit"),
    make_option("--n_projects", type = "integer", default = 50L,
                help = "Max number of mouse SRA projects to show with --list_tissues [default: 50]. Use 0 for all."),
    make_option("--list_metadata", action = "store_true", default = FALSE,
                help = paste("Print the biological sample metadata columns available for the",
                             "given --project or --tissue and exit. Use this to discover which",
                             "--filter_col and --filter_val values to use.")),
    make_option("--filter_col", type = "character", default = NULL,
                help = paste("colData column to filter samples on before summing",
                             "(e.g. sra.sample_title). Must be paired with --filter_val.",
                             "Run --list_metadata first to see available columns.")),
    make_option("--filter_val", type = "character", default = NULL,
                help = paste("Value to match in --filter_col (case-insensitive substring).",
                             "e.g. 'Heart' matches 'Heart' and 'Heart - Left Ventricle'.")),
    make_option("--strip_chr", action = "store_true", default = FALSE,
                help = paste("Strip 'chr' prefix from chromosome names.",
                             "Only needed if your reference uses Ensembl-style names (1, 2, X).",
                             "Default: keep 'chr' to match GENCODE/UCSC references (chr1, chr2, chrX)."))
)

opt <- parse_args(OptionParser(
    option_list = option_list,
    description  = paste(
        "Download recount3 junction data and write a STAR SJ.out.tab file.",
        "Column 7 encodes total read coverage summed across (filtered) samples."
    )
))

# ── List available tissues/projects and exit ─────────────────────────────────
if (opt$list_tissues) {
    message("Fetching available projects from recount3 (", opt$organism, ")...")
    projects <- available_projects(organism = opt$organism)
    if (opt$organism == "human") {
        gtex <- projects[projects$file_source == "gtex", ]
        cat("Available GTEx tissues (human):\n")
        cat(paste0("  ", sort(unique(gtex$project)), "\n"), sep = "")
    } else {
        sra    <- projects[projects$file_source == "sra", ]
        sra    <- sra[order(-sra$n_samples), ]
        # Pin Tabula Muris Senis (SRP199494) to the top as a recommended reference.
        tms_id  <- "SRP199494"
        tms_idx <- which(sra$project == tms_id)
        if (length(tms_idx) > 0)
            sra <- rbind(sra[tms_idx, ], sra[-tms_idx, ])
        n_show <- if (opt$n_projects == 0L) nrow(sra) else min(opt$n_projects, nrow(sra))
        cat(sprintf("Available mouse SRA projects (%d total, showing top %d by sample count):\n",
                    nrow(sra), n_show))
        cat(sprintf("  %-42s  %s\n", "Project", "Samples"), sep = "")
        cat(sprintf("  %-42s  %s\n", "-------", "-------"), sep = "")
        for (i in seq_len(n_show)) {
            label <- sra$project[i]
            if (label == tms_id) label <- paste0(label, " (Tabula Muris Senis bulk)")
            cat(sprintf("  %-42s  %d\n", label, sra$n_samples[i]), sep = "")
        }
        cat(sprintf("\nTo see all projects:  Rscript make_recount3_sj.R --list_tissues --organism mouse --n_projects 0 > projects.txt\n"))
        cat("To find your tissue: https://www.ncbi.nlm.nih.gov/sra or https://rna.recount.bio/\n")
    }
    quit(status = 0)
}

# ── Validate arguments ────────────────────────────────────────────────────────
if (opt$organism == "human" && is.null(opt$tissue) && is.null(opt$project)) {
    stop("For human, specify --tissue (GTEx) or --project (SRA). ",
         "Use --list_tissues to see GTEx tissue names.")
}
if (opt$organism == "mouse" && is.null(opt$project)) {
    stop("For mouse, specify --project (SRA project ID). ",
         "Use --list_tissues to see available projects.")
}
if (xor(is.null(opt$filter_col), is.null(opt$filter_val))) {
    stop("--filter_col and --filter_val must be used together.")
}

label      <- if (!is.null(opt$tissue)) opt$tissue else opt$project
opt$output <- if (is.null(opt$output)) paste0(label, "_recount3.SJ.out.tab") else opt$output

# ── Load project ──────────────────────────────────────────────────────────────
message("Fetching available projects...")
projects <- available_projects(organism = opt$organism)

if (!is.null(opt$tissue)) {
    proj <- projects[projects$file_source == "gtex" &
                     projects$project     == opt$tissue, ]
    if (nrow(proj) == 0)
        stop("Tissue '", opt$tissue, "' not found. Use --list_tissues to see options.")
} else {
    proj <- projects[projects$project == opt$project, ]
    if (nrow(proj) == 0)
        stop("Project '", opt$project, "' not found. Use --list_tissues to see options.")
    proj <- proj[1, ]
}

# ── Helper: build the URL for the per-sample metadata TSV ─────────────────────
# Pattern observed from recount3 caching logs:
#   {base}/{organism}/data_sources/{source}/metadata/{last2}/{project}/
#   {source}.{source}.{project}.MD.gz
# Returns NULL for non-SRA sources (GTEx uses a different layout).
sra_metadata_url <- function(organism, proj_row) {
    if (proj_row$file_source != "sra") return(NULL)
    src  <- proj_row$file_source
    pid  <- proj_row$project
    last2 <- substr(pid, nchar(pid) - 1, nchar(pid))
    sprintf("http://duffel.rail.bio/recount3/%s/data_sources/%s/metadata/%s/%s/%s.%s.%s.MD.gz",
            organism, src, last2, pid, src, src, pid)
}

# Download metadata TSV to a tempfile; returns data.frame or NULL on failure.
read_sra_metadata <- function(url) {
    tmp <- tempfile(fileext = ".gz")
    ok  <- tryCatch({
        download.file(url, tmp, mode = "wb", quiet = TRUE); TRUE
    }, error = function(e) FALSE)
    if (!ok) { unlink(tmp); return(NULL) }
    df <- tryCatch(
        read.table(gzfile(tmp), header = TRUE, sep = "\t",
                   comment.char = "", quote = "", fill = TRUE),
        error = function(e) NULL
    )
    unlink(tmp)
    df
}

# ── List sample metadata and exit ─────────────────────────────────────────────
# Downloads only the sample-level metadata TSV (a few MB at most) — no count
# matrix needed. Bypasses create_rse() so this works on login nodes too.
# GTEx projects are already tissue-specific so --list_metadata just confirms that.
if (opt$list_metadata) {
    if (!is.null(opt$tissue)) {
        cat(sprintf("\nGTEx project '%s' (%d samples) is already tissue-specific.\n",
                    proj$project, proj$n_samples))
        cat("--filter_col / --filter_val are not needed for GTEx projects.\n")
        cat("To inspect GTEx demographic metadata (sex, age, etc.) the full\n")
        cat("count matrix must be loaded via create_rse() — not supported with\n")
        cat("--list_metadata because it requires a compute node.\n")
        quit(status = 0)
    }

    source_id  <- proj$file_source
    project_id <- proj$project
    meta_url   <- sra_metadata_url(opt$organism, proj)
    if (is.null(meta_url))
        stop("--list_metadata is only supported for SRA projects.")

    message("Downloading sample metadata for ", project_id,
            " (", proj$n_samples, " samples)...")
    meta <- read_sra_metadata(meta_url)
    if (is.null(meta))
        stop("Could not download metadata from:\n  ", meta_url)

    cat(sprintf("\nSample metadata for %s  (%d samples)\n", proj$project, nrow(meta)))
    cat("─────────────────────────────────────────────────────────────────────\n")
    cat("Use --filter_col <COL> --filter_val <VALUE> to subset samples.\n")
    cat("Matching is case-insensitive substring, so partial values work.\n\n")

    # In create_rse() colData, the raw metadata columns are prefixed with the
    # source name (e.g. "sra.sample_attributes"). Prepend it here so the column
    # names shown to the user match exactly what --filter_col expects.
    # Also skip obvious non-biological columns (IDs, sizes, dates, counts).
    skip_pattern <- paste0("^(",
        "rail_id|external_id|sample_acc|experiment_acc|submission_acc|",
        "sample_bases|sample_spots|run_published|size|run_total_bases|",
        "run_total_spots|num_spots|read_info|run_alias|sample_name|",
        "sample_title|library_name",
        ")$")

    shown <- 0L
    for (col in colnames(meta)) {
        if (grepl(skip_pattern, col)) next
        uvals <- sort(unique(as.character(meta[[col]])))
        uvals <- uvals[!is.na(uvals) & nchar(uvals) > 0]
        if (length(uvals) <= 1L) next    # same value for every sample — uninformative
        shown <- shown + 1L
        display_col <- paste0(source_id, ".", col)   # e.g. sra.sample_attributes
        if (length(uvals) <= 12L) {
            cat(sprintf("  %-45s  [%3d unique]  %s\n",
                        display_col, length(uvals), paste(uvals, collapse = " | ")))
        } else {
            cat(sprintf("  %-45s  [%3d unique]  %s ... (truncated)\n",
                        display_col, length(uvals), paste(head(uvals, 5), collapse = " | ")))
        }
    }

    if (shown == 0L)
        cat("  No informative metadata columns found.\n")

    cat(sprintf("\n%d informative column(s) shown (IDs, sizes, dates and recount_qc.* omitted).\n",
                shown))
    cat(sprintf("\nExample usage:\n"))
    cat(sprintf("  Rscript make_recount3_sj.R --%s %s --filter_col <COL> --filter_val <VALUE> --output out.SJ.out.tab\n",
                if (!is.null(opt$tissue)) "tissue" else "project",
                if (!is.null(opt$tissue)) opt$tissue  else opt$project))
    quit(status = 0)
}

# ── Pre-validate --filter_col before the heavy junction download ───────────────
# For SRA projects, fetch the small metadata TSV (a few MB) and check the
# column name now — avoids waiting for the 100MB+ junction matrix download only
# to fail on a typo.  Silently skipped if the URL is unavailable (e.g. GTEx).
if (!is.null(opt$filter_col)) {
    meta_url <- sra_metadata_url(opt$organism, proj)
    if (!is.null(meta_url)) {
        message("Pre-validating --filter_col against sample metadata...")
        meta_preview <- read_sra_metadata(meta_url)
        if (!is.null(meta_preview)) {
            source_id <- proj$file_source
            col_raw   <- sub(paste0("^", source_id, "\\."), "", opt$filter_col)
            if (!col_raw %in% colnames(meta_preview)) {
                avail <- paste(paste0(source_id, ".", colnames(meta_preview)),
                               collapse = "\n    ")
                stop(sprintf(
                    "Column '%s' not found in sample metadata for %s.\nAvailable columns:\n    %s\nRun --list_metadata to see values.",
                    opt$filter_col, proj$project, avail))
            }
        }
    }
}

# ── Patch create_rse_manual to handle missing-sample data gaps ────────────────
# Some recount3 projects have samples in the metadata whose junction counts
# were never produced (pipeline failure). create_rse_manual matches samples via
# rail_id and crashes when some are absent from the count file (NA subscripts
# in Matrix 1.6+, dimension mismatch otherwise). The fix: drop NA matches from
# BOTH m AND metadata before the subscript so dimensions stay consistent.
local({
    ns  <- asNamespace("recount3")
    old <- tryCatch(get("create_rse_manual", envir = ns, inherits = FALSE),
                    error = function(e) NULL)
    if (is.null(old)) return(invisible(NULL))

    b <- body(old)
    # Find and patch the line:  counts <- counts[, m, drop = FALSE]
    # It sits right after the match() and stopifnot() for rail_id.
    patch_fn <- function(node) {
        if (is.call(node)) {
            txt <- deparse(node, width.cutoff = 500L)
            target <- "counts <- counts[, m, drop = FALSE]"
            if (identical(txt, target)) {
                return(quote({
                    if (anyNA(m)) {
                        n_skip <- sum(is.na(m))
                        warning(sprintf(
                            "%d sample(s) absent from the junction count file and skipped: %s",
                            n_skip,
                            paste(metadata$external_id[is.na(m)], collapse = ", ")),
                            call. = FALSE)
                        metadata <- metadata[!is.na(m), ]
                        m <- m[!is.na(m)]
                    }
                    counts <- counts[, m, drop = FALSE]
                }))
            }
            node[] <- lapply(node, patch_fn)
        }
        node
    }
    body(old) <- patch_fn(b)
    tryCatch(utils::assignInNamespace("create_rse_manual", old, ns = "recount3"),
             error = function(e) NULL)
})

message("Downloading junction data: ", proj$project,
        " (", proj$n_samples, " samples)...")
rse_jxn <- create_rse(proj, type = "jxn")

# ── Optional: filter to a subset of samples ───────────────────────────────────
# Restricts the count matrix to samples matching --filter_col / --filter_val
# before summing. pct_samples in the sidecar is then relative to the filtered
# cohort (not the full study), which is the meaningful denominator for
# tissue-matched validation.
if (!is.null(opt$filter_col)) {
    meta <- as.data.frame(colData(rse_jxn))
    if (!opt$filter_col %in% colnames(meta))
        stop(sprintf(
            "Column '%s' not found in sample metadata. ",
            "Run with --list_metadata to see available columns.",
            opt$filter_col))
    keep_samples <- grepl(opt$filter_val, as.character(meta[[opt$filter_col]]),
                          ignore.case = TRUE)
    n_before <- ncol(rse_jxn)
    rse_jxn  <- rse_jxn[, keep_samples]
    message(sprintf("Sample filter  --filter_col '%s'  --filter_val '%s':  %d / %d samples kept.",
                    opt$filter_col, opt$filter_val, sum(keep_samples), n_before))
    if (ncol(rse_jxn) == 0L)
        stop("No samples matched the filter. Check --filter_col / --filter_val, ",
             "or run --list_metadata to inspect the available metadata.")
}

# ── Aggregate per junction: total read coverage AND sample support ────────────
# Two complementary metrics come out of the junction x sample count matrix:
#   * total_reads : reads summed across all (filtered) samples -> SJ.out.tab
#                   column 7, so SQANTI3's min_cov is the standard "min read
#                   coverage" semantics.
#   * n_samples   : number of (filtered) samples with count > 0 -> written to
#                   the sidecar and folded into classification
#                   (min_sample_cov / _pct) by SQANTI-sc.
message("Aggregating across ", ncol(rse_jxn), " samples...")
counts_mat          <- assay(rse_jxn, "counts")
n_total_samples     <- ncol(counts_mat)
total_reads_per_jxn <- rowSums(counts_mat)
n_samples_with_jxn  <- rowSums(counts_mat > 0)

keep <- n_samples_with_jxn >= opt$min_samples
message(sprintf("Junctions total: %d  |  passing --min_samples %d: %d",
                length(keep), opt$min_samples, sum(keep)))

rse_jxn             <- rse_jxn[keep, ]
total_reads_per_jxn <- total_reads_per_jxn[keep]
n_samples_with_jxn  <- n_samples_with_jxn[keep]

# ── Build SJ.out.tab ──────────────────────────────────────────────────────────
rr <- rowRanges(rse_jxn)
rd <- as.data.frame(rowData(rse_jxn))

chrom       <- as.character(seqnames(rr))
jxn_start   <- start(rr)
jxn_end     <- end(rr)
strand_char <- as.character(strand(rr))

# STAR strand codes: 0=undefined, 1=+, 2=-
strand_code <- ifelse(strand_char == "+", 1L,
               ifelse(strand_char == "-", 2L, 0L))

# STAR motif codes: 0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG,
#                   4=CT/GC, 5=AT/AC, 6=GT/AT
motif_map  <- c(GTAG = 1L, CTAC = 2L, GCAG = 3L,
                CTGC = 4L, ATAC = 5L, GTAT = 6L)
motif_str  <- paste0(rd$left_motif, rd$right_motif)
motif_code <- ifelse(motif_str %in% names(motif_map),
                     motif_map[motif_str], 0L)

annotated  <- as.integer(rd$annotated)

# Keep 'chr' prefix by default (matches GENCODE/UCSC references).
# Pass --strip_chr only if your reference uses Ensembl-style names (1, 2, X).
if (opt$strip_chr) {
    chrom <- sub("^chr", "", chrom)
}

sj <- data.frame(
    chrom     = chrom,
    start     = jxn_start,
    end       = jxn_end,
    strand    = strand_code,
    motif     = motif_code,
    annotated = annotated,
    coverage  = as.integer(total_reads_per_jxn),  # col 7: total reads across (filtered) samples
    col8      = 0L,  # multi-mapping reads (not available from recount3)
    col9      = 0L,  # max overhang (not available from recount3)
    stringsAsFactors = FALSE
)

sj <- sj[order(sj$chrom, sj$start, sj$end), ]

# ── Write output ──────────────────────────────────────────────────────────────
# Real STAR SJ.out.tab files are uncompressed, and SQANTI3 reads coverage with a
# plain open() (no gzip). Write uncompressed by default so the file works as a
# direct --coverage input; honour an explicit .gz extension for archival only
# (such a file must be decompressed before use with SQANTI3).
message("Writing ", nrow(sj), " junctions to ", opt$output)
con <- if (grepl("\\.gz$", opt$output)) gzfile(opt$output, "w") else file(opt$output, "w")
write.table(sj, con, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
close(con)

# ── Sidecar: per-junction sample support ──────────────────────────────────────
# SJ.out.tab column 7 can only carry one value per junction (read coverage), so
# sample support (count + percentage) goes to a companion TSV that SQANTI-sc
# joins onto the SQANTI3 junctions file by genomic coordinate and reduces per
# isoform into the min_sample_cov / min_sample_pct classification columns.
# start/end are the SAME 1-based intron coordinates written to SJ.out.tab
# columns 2-3, which SQANTI3 stores verbatim as genomic_start_coord /
# genomic_end_coord, so the downstream join needs no offset.
# pct_samples is relative to the filtered sample count (n_total_samples), so
# when --filter_col is used the percentage reflects the chosen cohort.
sidecar <- data.frame(
    chrom       = chrom,  # already chr-stripped/kept above
    start       = jxn_start,
    end         = jxn_end,
    strand      = strand_char,                                  # '+', '-' or '*'
    n_samples   = as.integer(n_samples_with_jxn),
    pct_samples = round(100 * n_samples_with_jxn / n_total_samples, 4),
    stringsAsFactors = FALSE
)
# STAR (and SQANTI3) place undefined-strand junctions on BOTH strands; mirror
# that so the coordinate join finds them whatever strand the isoform is on.
und <- sidecar$strand == "*"
if (any(und)) {
    sidecar <- rbind(sidecar[!und, ],
                     transform(sidecar[und, ], strand = "+"),
                     transform(sidecar[und, ], strand = "-"))
}
sidecar <- sidecar[order(sidecar$chrom, sidecar$start, sidecar$end), ]

sidecar_path <- paste0(opt$output, ".sample_support.tsv")
filter_note  <- if (!is.null(opt$filter_col))
    sprintf(" (filtered to %s='%s')", opt$filter_col, opt$filter_val) else ""
message(sprintf("Writing sample-support sidecar (%d rows, %d%s samples) to %s",
                nrow(sidecar), n_total_samples, filter_note, sidecar_path))
write.table(sidecar, sidecar_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("Done.")
message("Pass to SQANTI-sc with:  --coverage ", opt$output)
message("Note: SJ.out.tab column 7 is total read coverage (so SQANTI3 --min_cov ",
        "= min read coverage). Per-junction sample support lives in ",
        sidecar_path, "; SQANTI-sc folds it into min_sample_cov / min_sample_pct.")
