#!/usr/bin/env Rscript
# make_recount3_sj.R
#
# Download junction data from recount3 and convert to STAR SJ.out.tab format
# for use as --coverage input to SQANTI-sc.
#
# Junction coverage (column 7 of SJ.out.tab) is encoded as the NUMBER OF SAMPLES
# in which that junction was observed (count > 0), NOT total read counts. This
# makes the SQANTI-sc --min_cov threshold directly interpretable: --min_cov 10
# means "junction must appear in at least 10 independent samples."
#
# Human data uses GTEx (tissue names). Mouse data uses SRA project IDs.
#
# Usage:
#   Human GTEx:   Rscript make_recount3_sj.R --organism human --tissue BLOOD \
#                     --output blood_gtex.SJ.out.tab.gz
#   Mouse SRA:    Rscript make_recount3_sj.R --organism mouse --project SRP123456 \
#                     --output mouse_brain.SJ.out.tab.gz
#   List tissues: Rscript make_recount3_sj.R --list_tissues [--organism mouse]
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
                help = "Output file [default: <tissue/project>_recount3.SJ.out.tab.gz]"),
    make_option("--min_samples", type = "integer", default = 1,
                help = "Min samples a junction must appear in to be kept [default: 1]"),
    make_option("--list_tissues", action = "store_true", default = FALSE,
                help = "List available GTEx tissues (human) or SRA projects (mouse) and exit"),
    make_option("--n_projects", type = "integer", default = 50L,
                help = "Max number of mouse SRA projects to show with --list_tissues [default: 50]. Use 0 for all.")
)

opt <- parse_args(OptionParser(
    option_list = option_list,
    description  = paste(
        "Download recount3 junction data and write a STAR SJ.out.tab file.",
        "Column 7 encodes number of samples supporting each junction."
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
        n_show <- if (opt$n_projects == 0L) nrow(sra) else min(opt$n_projects, nrow(sra))
        cat(sprintf("Available mouse SRA projects (%d total, showing top %d by sample count):\n",
                    nrow(sra), n_show))
        cat(sprintf("  %-15s  %s\n", "Project", "Samples"), sep = "")
        cat(sprintf("  %-15s  %s\n", "-------", "-------"), sep = "")
        for (i in seq_len(n_show)) {
            cat(sprintf("  %-15s  %d\n", sra$project[i], sra$n_samples[i]), sep = "")
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

label      <- if (!is.null(opt$tissue)) opt$tissue else opt$project
opt$output <- if (is.null(opt$output)) paste0(label, "_recount3.SJ.out.tab.gz") else opt$output

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

message("Downloading junction data: ", proj$project,
        " (", proj$n_samples, " samples)...")
rse_jxn <- create_rse(proj, type = "jxn")

# ── Aggregate: number of samples with count > 0 ───────────────────────────────
message("Aggregating across ", ncol(rse_jxn), " samples...")
counts_mat          <- assay(rse_jxn, "counts")
n_samples_with_jxn  <- rowSums(counts_mat > 0)

keep <- n_samples_with_jxn >= opt$min_samples
message(sprintf("Junctions total: %d  |  passing --min_samples %d: %d",
                length(keep), opt$min_samples, sum(keep)))

rse_jxn            <- rse_jxn[keep, ]
n_samples_with_jxn <- n_samples_with_jxn[keep]

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

sj <- data.frame(
    chrom     = chrom,
    start     = jxn_start,
    end       = jxn_end,
    strand    = strand_code,
    motif     = motif_code,
    annotated = annotated,
    n_samples = as.integer(n_samples_with_jxn),  # col 7: sample support
    col8      = 0L,  # multi-mapping reads (not available from recount3)
    col9      = 0L,  # max overhang (not available from recount3)
    stringsAsFactors = FALSE
)

sj <- sj[order(sj$chrom, sj$start, sj$end), ]

# ── Write output ──────────────────────────────────────────────────────────────
message("Writing ", nrow(sj), " junctions to ", opt$output)
gz <- gzfile(opt$output, "w")
write.table(sj, gz, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
close(gz)

message("Done.")
message("Pass to SQANTI-sc with:  --coverage ", opt$output)
message("Note: --min_cov threshold now means 'min number of GTEx/SRA samples",
        " supporting the junction'.")
