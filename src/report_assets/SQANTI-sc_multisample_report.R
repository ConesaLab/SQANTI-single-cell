#!/usr/env/bin Rscript

  # ##### SQANTI single-cell multisample report generation #####

### Author: Carlos Blanco

#********************** Packages

# !/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(stringr)
  library(rmarkdown)
})

# Ensure RColorConesa is available; DO NOT install from GitHub to prevent HPC timeouts.
if (!requireNamespace("RColorConesa", quietly = TRUE)) {
  message("[WARNING] RColorConesa not installed. Using fallback color palettes.")
  get_conesa_fallback <- function(n) {
    base_cols <- c("#15918A", "#F58A53", "#FDC659", "#74CDF0", "#FDA3D1", "#9F7BB8", "#EE446F")
    if (n <= length(base_cols)) return(base_cols[1:n])
    grDevices::colorRampPalette(base_cols)(n)
  }
  scale_fill_conesa <- function(palette = "complete", drop = FALSE, ...) {
    ggplot2::discrete_scale("fill", "conesa_fallback", get_conesa_fallback, drop = drop, ...)
  }
  scale_color_conesa <- function(palette = "complete", drop = FALSE, ...) {
    ggplot2::discrete_scale("colour", "conesa_fallback", get_conesa_fallback, drop = drop, ...)
  }
} else {
  library(RColorConesa)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  # Simple flag parser: expects --key value pairs, and a single string for --files (comma-separated)
  res <- list(
    files = NULL, class_files = NULL, out_dir = ".", mode = "reads",
    report = "pdf", prefix = "SQANTI_sc_multi_report",
    color_group = NULL, shape_group = NULL, shade_group = NULL,
    pca_features = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (startsWith(key, "--")) {
      k <- substring(key, 3)
      if (k %in% c("files", "class_files", "out_dir", "mode", "report", "prefix",
                    "color_group", "shape_group", "shade_group", "pca_features")) {
        if (i + 1 <= length(args)) {
          res[[k]] <- args[i + 1]
          i <- i + 2
          next
        } else {
          stop(sprintf("Missing value for flag %s", key))
        }
      } else {
        stop(sprintf("Unknown flag: %s", key))
      }
    } else {
      stop(sprintf("Unexpected argument: %s", key))
    }
  }
  if (is.null(res$files) || !nzchar(res$files)) {
    stop("--files must be provided (comma-separated list of cell summary files)")
  }
  if (!(res$report %in% c("pdf", "html", "both"))) {
    stop("--report must be one of: pdf, html, both")
  }
  res
}

# Helper: lighten a hex color toward white by a given amount (0 = unchanged, 1 = white)
lighten_hex <- function(hex, amount) {
  rgb_vals <- grDevices::col2rgb(hex) / 255
  new_rgb  <- rgb_vals + amount * (1 - rgb_vals)
  grDevices::rgb(new_rgb[1, ], new_rgb[2, ], new_rgb[3, ])
}

# Helper: parse a comma-separated group vector from a CLI arg; returns NULL if absent/mismatched
parse_group_vec <- function(param_val, n_files) {
  if (is.null(param_val) || !nzchar(param_val)) return(NULL)
  vec <- trimws(unlist(strsplit(param_val, ",", fixed = TRUE)))
  if (length(vec) != n_files) {
    message(sprintf(
      "[WARNING] Group vector length (%d) does not match number of files (%d). Ignoring.",
      length(vec), n_files
    ))
    return(NULL)
  }
  vec
}

safe_read_summary <- function(fpath) {
  # read.table supports gz automatically
  df <- tryCatch(
    {
      if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::fread(fpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, data.table = FALSE)
      } else {
        read.table(fpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
      }
    },
    error = function(e) {
      message(sprintf("[ERROR] Failed to read summary %s: %s", fpath, e$message))
      return(NULL)
    }
  )
  if (is.null(df)) {
    return(NULL)
  }
  if (ncol(df) < 2 || !("CB" %in% colnames(df))) {
    message(sprintf("[WARNING] Summary %s does not have expected structure; skipping", fpath))
    return(NULL)
  }
  # Coerce numeric columns (col 2..n) to numeric
  if (ncol(df) >= 2) {
    for (j in 2:ncol(df)) {
      df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
    }
  }
  # Derive sampleID from filename: <sampleID>_SQANTI_cell_summary.txt.gz
  base <- basename(fpath)
  sample <- sub("_SQANTI_cell_summary\\.txt(\\.gz)?$", "", base)
  df$sampleID <- sample
  df
}

# Helper: tidy feature names for figure titles
format_feature_display_name <- function(feature) {
  cleaned <- feature
  cleaned <- gsub(
    "(_prop_in_cell|_perc_in_cell|_prop|_perc|_percentage|_pct|_ratio|_counts?|_in_cell|_per_cell|_value)$",
    "", cleaned, ignore.case = TRUE
  )
  cleaned <- gsub("([0-9]+)b", "\\1 bp", cleaned, ignore.case = TRUE)
  cleaned <- stringr::str_replace_all(cleaned, "_+", " ")
  cleaned <- stringr::str_squish(cleaned)
  if (cleaned == "") cleaned <- feature
  tokens <- unlist(strsplit(cleaned, " ", fixed = FALSE))
  if (length(tokens) == 0) {
    return(feature)
  }
  replacements <- c(
    "Fsm" = "FSM",
    "Ism" = "ISM",
    "Nic" = "NIC",
    "Nnc" = "NNC",
    "Rts" = "RTS",
    "Tss" = "TSS",
    "Nmd" = "NMD",
    "Cage" = "CAGE",
    "Ujcs" = "UJCs",
    "Ujc" = "UJC",
    "Umis" = "UMIs",
    "Umi" = "UMI",
    "Mt" = "MT",
    "Cb" = "CB",
    "Pc" = "PC",
    "Qc" = "QC",
    "Orf" = "ORF",
    "Tpm" = "TPM"
  )
  normalize_token <- function(tok) {
    if (tok == "") {
      return(tok)
    }
    if (grepl("-", tok, fixed = TRUE)) {
      parts <- strsplit(tok, "-", fixed = TRUE)[[1]]
      parts <- vapply(parts, normalize_token, character(1), USE.NAMES = FALSE)
      return(paste(parts, collapse = "-"))
    }
    if (tok == toupper(tok)) {
      return(tok)
    }
    if (grepl("^[0-9]+(\\.[0-9]+)?$", tok)) {
      return(tok)
    }
    if (grepl("^[0-9]+bp$", tok, ignore.case = TRUE)) {
      return(gsub("bp$", "bp", tok, ignore.case = TRUE))
    }
    key <- stringr::str_to_title(tok)
    if (key %in% names(replacements)) {
      return(replacements[[key]])
    }
    stringr::str_to_title(tok)
  }
  tokens <- vapply(tokens, normalize_token, character(1), USE.NAMES = FALSE)
  stringr::str_squish(paste(tokens, collapse = " "))
}

# Helper: derive axis labels and scaling behaviour from feature metadata
infer_feature_metadata <- function(feature, values) {
  name_lower <- tolower(feature)
  # A per-read/per-transcript rate: the numerator names the domain, so the
  # suffix is stripped before the rules below run, and the values (<= 1) must
  # not trip the "looks like a fraction, call it a percentage" branch.
  rate_unit <- stringr::str_match(name_lower, "_per_(read|transcript)$")[, 2]
  if (!is.na(rate_unit)) name_lower <- sub("_per_(read|transcript)$", "", name_lower)
  domain <- "Value"
  if (stringr::str_detect(name_lower, "length")) domain <- "Length"
  if (stringr::str_detect(name_lower, "length") && stringr::str_detect(name_lower, "prop|perc|pct|ratio|fraction")) {
    if (exists("params") && params$mode == "isoforms") domain <- "Transcripts" else domain <- "Reads"
  }
  if (stringr::str_detect(name_lower, "read")) {
    if (exists("params") && params$mode == "isoforms") domain <- "Transcripts" else domain <- "Reads"
  }
  if (stringr::str_detect(name_lower, "\\bmt\\b") || stringr::str_detect(name_lower, "^mt_")) domain <- "Reads"
  if (stringr::str_detect(name_lower, "exon") && !stringr::str_detect(name_lower, "monoexon")) domain <- "Exons"
  if (stringr::str_detect(name_lower, "intron")) domain <- "Introns"
  if (stringr::str_detect(name_lower, "coverage")) domain <- "Coverage"
  if (stringr::str_detect(name_lower, "isoform")) domain <- "Isoforms"
  if (stringr::str_detect(name_lower, "transcript")) domain <- "Transcripts"
  if (stringr::str_detect(name_lower, "gene")) domain <- "Genes"
  if (stringr::str_detect(name_lower, "^anno_.*_bin") || stringr::str_detect(name_lower, "^novel_.*_bin")) domain <- "Genes"
  if (stringr::str_detect(name_lower, "umi")) domain <- "UMIs"
  if (stringr::str_detect(name_lower, "junction")) domain <- "Junctions"
  if (stringr::str_detect(name_lower, "^ujcs?(_|$)")) domain <- "UJCs"
  structural_keywords <- c(
    "(?<![a-z])fsm(?![a-z])", "(?<![a-z])ism(?![a-z])",
    "(?<![a-z])nic(?![a-z])", "(?<![a-z])nnc(?![a-z])",
    "genic_genomic", "genic genomic", "(?<![a-z])genic(?![a-z])",
    "antisense", "fusion", "intergenic",
    "genic_intron", "genic intron"
  )
  if (any(stringr::str_detect(name_lower, structural_keywords))) {
    if (domain != "Genes") {
      if (exists("params") && params$mode == "isoforms") domain <- "Transcripts" else domain <- "Reads"
    }
  }
  is_prop_keyword <- stringr::str_detect(name_lower, "prop|perc|pct|ratio|fraction")
  finite_vals <- values[is.finite(values)]
  unit <- if (is_prop_keyword) "%" else "count"
  scale_to_percent <- FALSE
  if (length(finite_vals) > 0) {
    maxv <- max(finite_vals)
    minv <- min(finite_vals)
    if (unit == "%" && maxv <= 1.5) {
      scale_to_percent <- TRUE
    }
    if (unit == "count" && !is_prop_keyword && is.na(rate_unit) && maxv <= 1.5 && minv >= 0) {
      unit <- "%"
      scale_to_percent <- TRUE
    }
    if (domain == "Length" && unit == "count") {
      unit <- "bp"
    }
  } else {
    if (unit == "%") scale_to_percent <- TRUE
  }
  if (unit == "%" && domain == "Value") {
    domain <- "Percentage"
  }
  value_label <- switch(unit,
    "%" = if (domain %in% c("Percentage", "Value")) "Value, %" else sprintf("%s, %%", domain),
    "bp" = sprintf("%s, bp", domain),
    sprintf("%s, count", domain)
  )
  if (domain == "Value" && unit == "count") {
    value_label <- "Value"
  }
  value_label <- infer_junction_display_label(feature, value_label)
  # Last, so the junction rule cannot append ", junctions" to UJCs per read --
  # that rate counts junction CHAINS per read, not junctions. Its own unit, not
  # "count": counts span orders of magnitude and get a log axis, while these
  # ratios sit inside (0, 1].
  if (!is.na(rate_unit)) {
    unit <- "rate"
    value_label <- sprintf("%s per %s", domain, rate_unit)
  }
  list(
    display_name = format_feature_display_name(feature),
    value_label = value_label,
    scale_to_percent = scale_to_percent,
    unit = unit,
    domain = domain
  )
}

infer_junction_display_label <- function(feature_name, current_label) {
  lower_name <- tolower(feature_name)
  
  # Gene bin features should not get junction annotations even if they contain "ujc"
  if (grepl("^anno_.*_bin", lower_name) || grepl("^novel_.*_bin", lower_name)) {
    return(current_label)
  }

  junction_keywords <- c("junction", "junctions", "splice", "sj", "canonical", "noncanonical", "ujc", "ujcs")
  contains_junction <- any(vapply(junction_keywords, function(kw) grepl(kw, lower_name, fixed = TRUE), logical(1)))
  if (!contains_junction) {
    return(current_label)
  }
  if (grepl("junct", current_label, ignore.case = TRUE)) {
    return(current_label)
  }
  if (grepl("%", current_label, fixed = TRUE)) {
    return(sub("%", " (junctions, %)", current_label, fixed = TRUE))
  }
  if (grepl("count", current_label, ignore.case = TRUE) || grepl("junction", current_label, ignore.case = TRUE)) {
    return(current_label)
  }
  paste0(current_label, ", junctions")
}

get_conesa_palette_colors <- function(n, palette = "complete") {
  if (n <= 0) {
    return(character(0))
  }

  col_vec <- NULL
  if (requireNamespace("RColorConesa", quietly = TRUE)) {
    col_vec <- tryCatch(
      {
        scale_obj <- RColorConesa::scale_fill_conesa(palette = palette)
        if (!is.null(scale_obj$palette) && is.function(scale_obj$palette)) {
          scale_obj$palette(n)
        } else {
          NULL
        }
      },
      error = function(e) NULL
    )

    if (is.null(col_vec) || length(col_vec) == 0) {
      col_vec <- tryCatch(
        {
          ns <- asNamespace("RColorConesa")
          if (exists("conesa_palettes", envir = ns, inherits = FALSE)) {
            pal_list <- get("conesa_palettes", envir = ns, inherits = FALSE)
            pal_entry <- pal_list[[palette]]
            if (is.function(pal_entry)) {
              pal_entry(n)
            } else if (is.vector(pal_entry)) {
              unname(pal_entry)
            } else {
              NULL
            }
          } else {
            NULL
          }
        },
        error = function(e) NULL
      )
    }

    if ((is.null(col_vec) || length(col_vec) == 0) && exists("palette_conesa", envir = asNamespace("RColorConesa"), inherits = FALSE)) {
      pal_fun <- get("palette_conesa", envir = asNamespace("RColorConesa"), inherits = FALSE)
      col_vec <- tryCatch(pal_fun(palette, n), error = function(e) NULL)
    }
  }

  if (is.null(col_vec) || length(col_vec) == 0) {
    fallback <- c(
      "#15918A", "#F58A53", "#FDC659", "#74CDF0",
      "#FDA3D1", "#9F7BB8", "#EE446F"
    )
    if (n <= length(fallback)) {
      return(fallback[1:n])
    } else {
      return(grDevices::colorRampPalette(fallback)(n))
    }
  }

  if (n <= length(col_vec)) {
    return(col_vec[1:n])
  } else {
    return(grDevices::colorRampPalette(col_vec)(n))
  }
}

# Helper: embed a small zoomed inset into the top-right whitespace of a fixed 0-100% plot.
# Conditions: max_val > 0 and max_val < 20 (data stays below 20%, leaving whitespace to fill).
# Shown in both HTML and PDF. Uses annotation_custom — no extra packages beyond ggplot2 + grid.
attach_zoom_inset <- function(main_plot, plot_df, x_var, y_var,
                               max_val, min_val, is_html,
                               use_violin = TRUE, fill_col = NULL) {
  if (max_val == 0 || max_val >= 20) return(main_plot)
  if (!is.factor(plot_df[[x_var]])) plot_df[[x_var]] <- factor(plot_df[[x_var]])
  n_grps <- nlevels(plot_df[[x_var]])

  gp_inset <- ggplot(plot_df, aes(x = .data[[x_var]], y = .data[[y_var]]))

  if (!is.null(fill_col)) {
    if (use_violin) {
      gp_inset <- gp_inset +
        geom_violin(fill = fill_col, colour = fill_col, alpha = 0.7,
                    scale = "width", trim = TRUE, linewidth = 0.2)
    }
    gp_inset <- gp_inset +
      geom_boxplot(fill = fill_col, colour = "grey30", width = 0.08,
                   outlier.shape = NA, alpha = 0.3, lwd = 0.2)
  } else {
    if (use_violin) {
      gp_inset <- gp_inset +
        geom_violin(aes(fill = .data[[x_var]], colour = .data[[x_var]]),
                    alpha = 0.7, scale = "width", trim = TRUE, linewidth = 0.2)
    }
    gp_inset <- gp_inset +
      geom_boxplot(aes(fill = .data[[x_var]]), colour = "grey30",
                   width = 0.08, outlier.shape = NA, alpha = 0.3, lwd = 0.2) +
      scale_fill_conesa(palette = "complete", drop = FALSE) +
      scale_color_conesa(palette = "complete", guide = "none", drop = FALSE)
  }

  gp_inset <- gp_inset +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1,
                 colour = "red", stroke = 0.45) +
    scale_y_continuous(
      limits = c(0, 100),
      labels = function(x) as.integer(x),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position  = "none",
      axis.title       = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.ticks.length = unit(2, "pt"),
      plot.margin      = margin(2, 3, 2, 2, "pt"),
      plot.background  = element_rect(fill = "white", colour = "grey60", linewidth = 0.4),
      panel.background = element_rect(fill = "white")
    )

  main_plot + annotation_custom(
    ggplotGrob(gp_inset),
    xmin = n_grps + 0.6 - n_grps * 0.5, xmax = n_grps + 0.6,
    ymin = 45,  ymax = 99
  )
}

# ── Curated feature registry ─────────────────────────────────────────────────
#
# The per-cell summary carries >300 numeric columns, and most of them are
# decompositions of one another: abundance bins sum to the gene count, the
# per-category junction breakdown sums to total_junctions, per-subcategory
# proportions sum to their parent category. Feeding every numeric column to
# prcomp(scale. = TRUE) weights each feature *family* by how many columns it
# happens to have rather than by how much it says, so a 36-column junction
# breakdown outvotes transcript length 36:1 and sets the direction of PC1.
#
# The registry keeps one feature per concept, prefers rates over raw counts
# (counts track sequencing depth, which otherwise becomes PC1) and tags each
# feature with an interpretable block. Decomposed families are not dropped from
# the report -- they keep their own dedicated figures; they just stop voting.
#
# conditional = TRUE marks features that only carry information when the
# matching SQANTI3 run flag was used. Those columns are written unconditionally
# by cell_metrics.py, filled with a constant sentinel when the flag is absent
# (e.g. NMD_prop_in_cell = 0, PolyA_motif_support_prop = 0), so a plain
# "column exists in every sample" test does not catch them -- see
# curated_feature_table().
# Entries are named vectors: name = cell-summary column, value = display label.
# Curating the feature set includes curating how it reads, so labels are
# explicit here rather than derived from the column name.
#
# Block names are constants because the report layout keys off them (see
# QC_BLOCKS_PLOTTED_ELSEWHERE): a literal string in both places would let a
# rename silently duplicate or drop a whole section of figures.
QC_BLOCK_YIELD      <- "Yield & detection"
QC_BLOCK_LENGTH     <- "Transcript length"
QC_BLOCK_STRUCTURAL <- "Structural categories"
QC_BLOCK_JUNCTIONS  <- "Splice junction composition"
QC_BLOCK_GOOD       <- "Good-quality features"
QC_BLOCK_BAD        <- "Bad-quality features"

# Blocks whose per-feature violins already exist earlier in the report: the
# entity count and MT share under library size, the two per-unit rates under
# gene characterisation, the per-cell median length under length, and all nine
# structural categories under structural categories. Those figures are richer
# than a generic panel (per-category colour, zoom insets, the combined
# comparison figure), so the sections built from this registry skip them rather
# than drawing a second, plainer copy.
QC_BLOCKS_PLOTTED_ELSEWHERE <- c(QC_BLOCK_YIELD, QC_BLOCK_LENGTH, QC_BLOCK_STRUCTURAL)

curated_feature_registry <- function(mode) {
  is_iso <- identical(mode, "isoforms")
  count_col <- if (is_iso) "Transcripts_in_cell" else "Reads_in_cell"
  count_lab <- if (is_iso) "Transcripts per cell" else "Reads per cell"
  unit <- if (is_iso) "transcript" else "read"
  unit_pl <- if (is_iso) "Transcripts" else "Reads"
  # `conditional` lists the features in a block that are only computed when a
  # SQANTI3 run flag was passed; the rest are always present. It is per feature,
  # not per block, because the good/bad quality blocks mix the two.
  #
  # Each mode contributes its own per-cell DIVERSITY feature to "Yield & detection",
  # so that block has the same number of rows either way. Reads mode divides
  # UJCs_in_cell (distinct junction chains) by the read count; isoforms mode divides
  # Isoforms_in_cell (distinct collapsed models) by the FL sum. The denominator is
  # depth in both cases -- Transcripts_in_cell is the FL sum, i.e. transcript copies,
  # which is why it is not itself a diversity feature. Only one of the two source
  # columns exists per mode, so naming just the applicable one here keeps
  # curated_feature_table() from reporting the other as an exclusion every run.
  div_feat <- if (is_iso) {
    stats::setNames(paste("Isoforms per", unit), paste0("Isoforms_per_", unit))
  } else {
    stats::setNames(paste("UJCs per", unit), paste0("UJCs_per_", unit))
  }
  reg <- list(
    # The two per-unit features are derived (see derived_feature_defs). Raw
    # counts are not comparable across samples because they track sequencing
    # depth; dividing by the cell's own read count puts them on a detection-
    # efficiency scale. The depth itself is kept as an explicit feature so an
    # odd yield row can be read against it.
    #
    # Genes_in_cell and Novel_genes are deliberately absent. SQANTI3 mints a
    # fresh novelGene_* id for almost every unassignable row (739 distinct ids
    # across 741 rows in the reads test set), so Novel_genes is a count of
    # unassignable reads rather than of genes, and Genes_in_cell -- which
    # counts distinct ids over ALL rows -- inherits that inflation in
    # proportion to the artifact rate. Annotated_genes is the clean one.
    #
    # MT_perc sits here rather than under a quality heading: it is
    # MT_reads / total_reads, a description of what was captured, and cell
    # types differ in mitochondrial expression for ordinary biological
    # reasons. A high value is not by itself a quality verdict.
    #
    # div_feat (defined above) is this mode's diversity feature.
    list(block = QC_BLOCK_YIELD, features = c(
      stats::setNames(count_lab, count_col),
      stats::setNames(paste("Annotated genes per", unit), paste0("Annotated_genes_per_", unit)),
      div_feat,
      MT_perc = "Mitochondrial (%)"
    )),
    list(block = QC_BLOCK_LENGTH, features = c(
      Median_length_per_cell = "Median length (bp)"
    )),
    list(block = QC_BLOCK_STRUCTURAL, features = c(
      FSM_prop           = "FSM (%)",
      ISM_prop           = "ISM (%)",
      NIC_prop           = "NIC (%)",
      NNC_prop           = "NNC (%)",
      Genic_Genomic_prop = "Genic genomic (%)",
      Antisense_prop     = "Antisense (%)",
      Fusion_prop        = "Fusion (%)",
      Intergenic_prop    = "Intergenic (%)",
      Genic_intron_prop  = "Genic intron (%)"
    )),
    # The only per-JUNCTION features in the summary: these divide by
    # total_junctions, everything else in the registry divides by the cell's
    # read (or FL) count. Left as composition rather than split into good and
    # bad, because novel canonical is genuinely ambiguous -- a real novel
    # junction and a mapping artifact both land there.
    list(block = QC_BLOCK_JUNCTIONS, features = c(
      Known_canonical_junctions_prop     = "Known canonical SJ (%)",
      Known_non_canonical_junctions_prop = "Known non-canonical SJ (%)",
      Novel_canonical_junctions_prop     = "Novel canonical SJ (%)",
      Novel_non_canonical_junctions_prop = "Novel non-canonical SJ (%)"
    )),
    # No subcategory shares (*_reference_match_prop, *_3prime_fragment_prop and
    # the rest of cell_metrics.py's sublevels table). There are ~30 of them and
    # no principled way to pick a subset, while taking all of them would make
    # subcategories the largest family in the registry and reintroduce the
    # column-count weighting this curation exists to remove. They also divide
    # by their parent CATEGORY rather than the cell's reads, so over a handful
    # of ISMs the ratio swings on nothing. --pca_features is the escape hatch
    # for anyone who wants them.
    #
    # Good/bad hold exactly what the per-sample report already declares as such
    # (SQANTI-sc_report.R, "Good quality features" / "Bad quality features"),
    # so one taxonomy covers both reports and this file invents no verdicts of
    # its own. Every feature in both blocks is a proportion of the cell's reads.
    list(block = QC_BLOCK_GOOD,
         conditional = c("PolyA_motif_support_prop", "CAGE_peak_support_prop",
                         "TSS_ratio_validated_prop", "srjunctions_support_prop"),
         features = c(
      TSSAnnotationSupport_prop = "TSS annotation support (%)",
      PolyA_motif_support_prop  = "PolyA motif support (%)",
      CAGE_peak_support_prop    = "CAGE peak support (%)",
      TSS_ratio_validated_prop  = "TSS ratio validated (%)",
      srjunctions_support_prop  = "Short-read SJ support (%)"
    )),
    list(block = QC_BLOCK_BAD,
         conditional = "NMD_prop_in_cell",
         features = c(
      RTS_prop_in_cell           = "RT-switching (%)",
      Intrapriming_prop_in_cell  = "Intrapriming (%)",
      # NOT a junction proportion: the fraction of multi-exon reads whose
      # all_canonical field is non_canonical, i.e. reads carrying at least one
      # non-canonical junction. Named accordingly so it is not read as a share
      # of junctions. Canonical_prop_in_cell is its exact complement over the
      # same denominator, so including both would double-count.
      Non_canonical_prop_in_cell = paste(unit_pl, "with non-canonical SJ (%)"),
      NMD_prop_in_cell           = "Predicted NMD (%)"
    ))
  )
  do.call(rbind, lapply(reg, function(b) {
    data.frame(
      feature = names(b$features), label = unname(b$features),
      block = b$block,
      conditional = names(b$features) %in%
        (if (is.null(b$conditional)) character(0) else b$conditional),
      stringsAsFactors = FALSE
    )
  }))
}

# Features computed from existing columns rather than read from the summary.
#
# These MUST be derived per cell and only then aggregated: the median of the
# per-cell ratios is not the ratio of the per-cell medians. A cell with 10
# genes from 10 reads and one with 1 gene from 100 reads have per-cell ratios
# of 1.00 and 0.01 (median 0.505), while dividing the medians gives 0.10.
derived_feature_defs <- function(mode) {
  is_iso <- identical(mode, "isoforms")
  count_col <- if (is_iso) "Transcripts_in_cell" else "Reads_in_cell"
  unit <- if (is_iso) "transcript" else "read"
  defs <- list()
  defs[[paste0("Annotated_genes_per_", unit)]] <- list(
    inputs = c("Annotated_genes", count_col),
    fun = function(d) d[["Annotated_genes"]] / d[[count_col]]
  )
  defs[[paste0("UJCs_per_", unit)]] <- list(
    inputs = c("UJCs_in_cell", count_col),
    fun = function(d) d[["UJCs_in_cell"]] / d[[count_col]]
  )
  # Each mode has exactly one per-cell diversity column, and it is the other mode's
  # that is missing: UJCs_in_cell is reads-only (it needs the jxn_string that
  # annotate_with_ujc_hash writes, which isoforms mode skips), Isoforms_in_cell is
  # isoforms-only. Defining both unconditionally is safe -- add_derived_features()
  # drops whichever def has missing inputs -- and keeps the two symmetric here.
  defs[[paste0("Isoforms_per_", unit)]] <- list(
    inputs = c("Isoforms_in_cell", count_col),
    fun = function(d) d[["Isoforms_in_cell"]] / d[[count_col]]
  )
  defs
}

# Append the derived columns to the merged per-cell table. Returns the
# augmented table plus the names that were actually added, so the availability
# check can treat them as present-in-every-sample (their INPUTS were).
add_derived_features <- function(multi, lst, mode) {
  defs <- derived_feature_defs(mode)
  common_cols <- if (length(lst) >= 1) Reduce(intersect, lapply(lst, colnames)) else character(0)
  added <- character(0)
  for (nm in names(defs)) {
    d <- defs[[nm]]
    if (!all(d$inputs %in% common_cols) || !all(d$inputs %in% colnames(multi))) next
    v <- suppressWarnings(as.numeric(d$fun(multi)))
    v[!is.finite(v)] <- NA_real_   # a zero denominator must not become Inf
    multi[[nm]] <- v
    added <- c(added, nm)
  }
  list(multi = multi, added = added)
}

# Read a user-supplied feature whitelist (one feature name per line; blank
# lines and lines starting with '#' ignored).
read_feature_whitelist <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    message(sprintf("[WARNING] --pca_features file not found, using the curated default: %s", path))
    return(NULL)
  }
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) == 0) {
    message("[WARNING] --pca_features file is empty, using the curated default.")
    return(NULL)
  }
  data.frame(
    feature = lines,
    label = vapply(lines, format_feature_display_name, character(1), USE.NAMES = FALSE),
    block = "User-specified", conditional = FALSE,
    stringsAsFactors = FALSE
  )
}

# Resolve the registry against the data actually loaded. Returns the retained
# feature/block table (possibly zero rows) and messages every exclusion.
curated_feature_table <- function(lst, multi, mode, features_file = NULL,
                                  derived = character(0)) {
  reg <- read_feature_whitelist(features_file)
  if (is.null(reg)) reg <- curated_feature_registry(mode)

  # (1) present in every original summary. A column missing from one sample can
  # reflect a different run flag rather than a true zero, so it is not usable
  # for a cross-sample comparison. Derived columns qualify when their inputs
  # passed the same test (add_derived_features only emits those).
  common_cols <- c(
    if (length(lst) >= 1) Reduce(intersect, lapply(lst, colnames)) else character(0),
    derived
  )
  usable <- reg$feature %in% common_cols &
    reg$feature %in% colnames(multi) &
    vapply(reg$feature, function(f) {
      f %in% colnames(multi) && is.numeric(multi[[f]]) && !all(is.na(multi[[f]]))
    }, logical(1))
  if (any(!usable)) {
    message(sprintf(
      "[INFO] Curated feature(s) not available in all samples, excluded: %s",
      paste(sort(reg$feature[!usable]), collapse = ", ")
    ))
  }
  reg <- reg[usable, , drop = FALSE]
  if (nrow(reg) == 0) return(reg)

  # (2) conditional features: the column exists everywhere, but is a constant
  # sentinel in the samples whose run did not compute it. If a conditional
  # feature does not vary within EVERY sample, the between-sample differences
  # are a run-flag artifact (they would otherwise dominate PC1), so drop it.
  varies_within_each_sample <- function(feat) {
    parts <- split(multi[[feat]], multi$sampleID)
    all(vapply(parts, function(v) {
      v <- v[is.finite(v)]
      length(v) > 1 && stats::sd(v) > 0
    }, logical(1)))
  }
  cond_idx <- which(reg$conditional)
  if (length(cond_idx) > 0) {
    keep_cond <- vapply(reg$feature[cond_idx], varies_within_each_sample, logical(1))
    drop_idx <- cond_idx[!keep_cond]
    if (length(drop_idx) > 0) {
      message(sprintf(
        paste0("[INFO] Conditional feature(s) constant within at least one sample ",
               "(the corresponding SQANTI3 run flag was not used everywhere), excluded: %s"),
        paste(sort(reg$feature[drop_idx]), collapse = ", ")
      ))
      # NB: reg[-integer(0), ] drops every row, so only negate a non-empty index.
      reg <- reg[-drop_idx, , drop = FALSE]
    }
  }

  reg$block <- factor(reg$block, levels = unique(reg$block))
  reg
}

# Helper: per-sample QC overview -- curated features x samples, coloured by how
# far each sample sits from the cohort. This is the between-sample view in
# directly readable form: you see which BLOCK a sample departs in, not just
# that it sits somewhere on a PC.
#
# Colour is (sample median - cohort median) / pooled within-sample CELL spread,
# computed row by row. Standardising per row is what lets bp, counts and
# percentages share one colour scale; dividing by cell spread rather than by the
# spread of the medians themselves is what lets a row say that a metric
# separates nothing (see qc_cell_deviation).
#
# The centre has to be the cohort median: no sample can be designated the
# baseline, and one bad library would otherwise drag a mean toward itself.
#
# robust_zscore() below is the pre-#1a colour and is still the fallback whenever
# per-cell values are unavailable; it also supplies the side-by-side value in
# the HTML tooltip.
#
# Deliberately NOT log-transformed first. A z-score is already invariant to
# scale, so features differing only by a constant factor get identical scores;
# log2(x + 1) breaks that, because the +1 pseudocount is large relative to a
# 0.5% value and compresses exactly the rare features we want treated fairly.
# Zeros are not a problem here either -- they are ordinary values to a z-score,
# unlike a fold-change, which would need x / 0.
#
# Flat-row guard: rows varying by less than this fraction of their centre are
# reported as flat. Expressed as a coefficient of variation so it stays
# relative rather than tied to any feature's units.
QC_OVERVIEW_FLAT_CV <- 0.02

# Above this many samples the per-tile medians stop being legible, so the
# heatmap drops them and colour carries the plot on its own. The exact values
# remain available in <prefix>_pca_feature_medians.tsv.
QC_OVERVIEW_MAX_LABELLED_SAMPLES <- 8

# Colour scale limit, in cell MADs. Two samples on opposite sides of the centre
# are 2x this far apart, so at 4 their cell distributions no longer overlap at
# all: past that point a larger number does not change the reading, and letting
# it stretch the gradient would compress every other row into near-white -- the
# failure this figure exists to avoid. Measured on the 4-sample cohort, the
# choice is flat between 3 and 6 (one row of 24 moves), so it is not finely
# tuned. Deliberately NOT calibrated to that cohort's own break at 5-6: two
# tools differ far more than the biological replicates this must also serve.
QC_OVERVIEW_DEV_CAP <- 4

# Restates a normal's 10-90 range on the scale stats::mad()'s default
# constant (1.4826) puts MAD on, so the fallback below is commensurable with
# the MAD it replaces rather than a second, silently different unit.
QC_OVERVIEW_P10_P90_TO_MAD <- 2.563

# Robust z with a mean/sd fallback. Returns all-zero for a row that does not
# meaningfully vary: z divides by the spread, so a cohort agreeing to within a
# fraction of a percent would otherwise have that noise stretched across the
# full colour range and read as a real difference.
robust_zscore <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  flat <- rep(0, length(v))
  centre <- stats::median(v, na.rm = TRUE)
  spread <- stats::mad(v, center = centre, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) {
    centre <- mean(v, na.rm = TRUE)
    spread <- stats::sd(v, na.rm = TRUE)
  }
  if (!is.finite(spread) || spread <= 0 || !is.finite(centre)) return(flat)
  if (centre > 0 && spread / centre < QC_OVERVIEW_FLAT_CV) return(flat)
  z <- (v - centre) / spread
  z[!is.finite(z)] <- 0
  z
}

# Spread of ONE sample's cell values for one feature, in MAD units. The
# quantile branch catches zero-inflation: a feature that is 0 in more than half
# a sample's cells has a cell MAD of exactly 0. NA rather than 0 when even the
# 10-90 range is empty, so such a sample drops out of the pool below instead of
# flattening a row the other samples can still separate on.
qc_cell_spread <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  v <- v[is.finite(v)]
  if (length(v) < 2) return(NA_real_)
  s <- stats::mad(v)
  if (!is.finite(s) || s <= 0) {
    q <- stats::quantile(v, c(0.10, 0.90), names = FALSE)
    s <- (q[2] - q[1]) / QC_OVERVIEW_P10_P90_TO_MAD
  }
  if (!is.finite(s) || s <= 0) return(NA_real_)
  s
}

qc_cell_spread_table <- function(multi, feats) {
  feats <- feats[feats %in% colnames(multi)]
  if (length(feats) == 0 || !"sampleID" %in% colnames(multi)) return(NULL)
  parts <- split(seq_len(nrow(multi)), as.character(multi$sampleID))
  do.call(rbind, lapply(feats, function(f) {
    data.frame(
      feature = f, sampleID = names(parts),
      cell_spread = vapply(parts, function(i) qc_cell_spread(multi[[f]][i]), numeric(1)),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }))
}

# How far each sample sits from the cohort centre, in cells' worth of spread.
#
# Pooled with the median rather than used per sample: dividing each sample by
# its own spread lets a noisy sample shrink its own tile and leaves one row's
# tiles in different units. Cohen's d pools for the same reason.
#
# No flat-CV guard here. That guard exists because the old denominator collapsed
# together with the numerator, which cannot happen once the denominator comes
# from within samples; applying it would veto exactly the small-but-consistent
# differences this score exists to detect.
qc_cell_deviation <- function(v, spreads) {
  v <- suppressWarnings(as.numeric(v))
  flat <- list(z = rep(0, length(v)), pooled = NA_real_)
  centre <- stats::median(v, na.rm = TRUE)
  pooled <- stats::median(suppressWarnings(as.numeric(spreads)), na.rm = TRUE)
  if (!is.finite(centre) || !is.finite(pooled) || pooled <= 0) return(flat)
  z <- (v - centre) / pooled
  z[!is.finite(z)] <- 0
  list(z = z, pooled = pooled)
}

# Order samples for the columns. Grouping beats clustering when the design
# declares one: the user brought a hypothesis ("do my conditions separate?")
# and reordering by similarity would answer a different question. With no
# groups, cluster on the z-scored features so similar samples sit together.
#
# Returns the column order plus the group of each sample (NULL when there is no
# grouping), because ordering alone leaves the groups invisible -- the caller
# facets on them so the boundaries are drawn and named.
order_overview_samples <- function(zmat, sample_levels, sample_group_map = NULL) {
  grp <- NULL
  if (!is.null(sample_group_map) && "color_group" %in% colnames(sample_group_map)) {
    g <- sample_group_map$color_group[match(sample_levels, sample_group_map$sampleID)]
    if (any(nzchar(g) & !is.na(g))) grp <- ifelse(is.na(g) | !nzchar(g), "ungrouped", g)
  }
  if (!is.null(grp)) {
    ord <- order(grp, seq_along(sample_levels))
    return(list(order = sample_levels[ord], groups = grp[ord]))
  }
  if (length(sample_levels) < 3) return(list(order = sample_levels, groups = NULL))
  ord <- tryCatch({
    # Columns are samples in zmat, so cluster the transpose.
    d <- stats::dist(t(zmat[, sample_levels, drop = FALSE]))
    if (!all(is.finite(d))) stop("non-finite distances")
    sample_levels[stats::hclust(d, method = "average")$order]
  }, error = function(e) sample_levels)
  list(order = ord, groups = NULL)
}

# Is the interactive HTML renderer usable? ggiraph draws the real ggplot to
# SVG, so the free-space facets, the block strips and the theme survive as-is;
# plotly would re-implement the plot and drop space = "free", collapsing the
# 9-row and 1-row blocks to equal heights. Absence is not fatal -- the static
# figure is still rendered.
qc_overview_interactive_ok <- function() {
  requireNamespace("ggiraph", quietly = TRUE)
}

build_qc_overview_plot <- function(agg_median, feature_map, sample_levels,
                                   is_html = FALSE, sample_group_map = NULL,
                                   interactive = FALSE, cell_spread = NULL) {
  feats <- as.character(feature_map$feature)
  feats <- feats[feats %in% colnames(agg_median)]
  n_samples <- nrow(agg_median)
  if (length(feats) < 2 || n_samples < 2) return(NULL)

  use_cell <- is.data.frame(cell_spread) && nrow(cell_spread) > 0 &&
    all(c("feature", "sampleID", "cell_spread") %in% colnames(cell_spread))

  long <- do.call(rbind, lapply(feats, function(f) {
    v <- suppressWarnings(as.numeric(agg_median[[f]]))
    own <- rep(NA_real_, length(v))
    if (use_cell) {
      i <- match(paste(f, agg_median$sampleID),
                 paste(cell_spread$feature, cell_spread$sampleID))
      own <- suppressWarnings(as.numeric(cell_spread$cell_spread))[i]
    }
    z_mad <- robust_zscore(v)
    dev <- if (use_cell) qc_cell_deviation(v, own) else NULL
    data.frame(
      sampleID = agg_median$sampleID, feature = f, value = v,
      z = if (is.null(dev)) z_mad else dev$z,
      z_mad = z_mad,
      cell_spread = own,
      pooled_spread = if (is.null(dev)) NA_real_ else dev$pooled,
      stringsAsFactors = FALSE
    )
  }))
  # The row headline. Not the maximum: (0, 0, 0, 2.5) and (-2.5, -0.8, 0.8, 2.5)
  # share a maximum but not a range, and only the second is the cohort splitting.
  long$row_range <- ave(long$z, long$feature, FUN = function(z) diff(range(z)))
  long <- dplyr::left_join(long, feature_map, by = "feature")
  # NB: keep this distinct from the joined `label` (the curated feature name),
  # which supplies the y-axis levels below.
  long$tile_text <- vapply(long$value, function(v) {
    if (!is.finite(v)) return("")
    # trimws: formatC pads to the requested width, which would push the
    # centred tile label off-centre.
    trimws(if (abs(v) >= 1000) formatC(v, format = "d", big.mark = ",")
           else formatC(signif(v, 3), format = "fg", digits = 3))
  }, character(1))

  # Rows top-to-bottom in registry order, blocks stacked in registry order.
  lab_levels <- feature_map$label[match(feats, feature_map$feature)]
  long$feature_lab <- factor(long$label, levels = rev(unique(lab_levels)))

  sample_levels <- sample_levels[sample_levels %in% long$sampleID]
  zmat <- stats::reshape(long[, c("feature", "sampleID", "z")], idvar = "feature",
                         timevar = "sampleID", direction = "wide")
  colnames(zmat) <- sub("^z\\.", "", colnames(zmat))
  col_order <- order_overview_samples(zmat, sample_levels, sample_group_map)
  long$sampleID <- factor(long$sampleID, levels = col_order$order)
  if (!is.null(col_order$groups)) {
    long$col_group <- factor(col_order$groups[match(long$sampleID, col_order$order)],
                             levels = unique(col_order$groups))
  }

  # Everything below scales with the sample count: the row count is fixed by
  # the registry, but the columns are however many samples the user brought.
  show_values <- n_samples <= QC_OVERVIEW_MAX_LABELLED_SAMPLES
  base_size <- if (is_html) 13 else 12
  if (n_samples > 12) base_size <- base_size - 2
  if (n_samples > 24) base_size <- base_size - 2
  x_angle <- if (n_samples > 6) 45 else 30
  tile_border <- if (n_samples > 20) 0.2 else 0.5

  interactive <- isTRUE(interactive) && qc_overview_interactive_ok()

  p <- ggplot(long, aes(x = sampleID, y = feature_lab, fill = z))
  if (interactive) {
    # The hover text carries what the tile cannot: the exact median AND the
    # z-score behind the colour, which is otherwise only readable off the
    # legend. This is also what lets the numbers be dropped from wide cohorts
    # without losing them.
    #
    # Single line on purpose. ggiraph 0.8.12 rewrites "\n" to <br/> and then
    # XML-escapes it, so the tooltip renders a literal "<br/>"; an explicit
    # <br/> is mangled the same way. Separators are the only form that survives.
    num2 <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 2), "n/a")
    num3 <- function(x) ifelse(is.finite(x), formatC(signif(x, 3), format = "fg", digits = 3), "n/a")
    long$tooltip <- if (use_cell) {
      # No apostrophes: ggiraph escapes them to &#39; and the surrounding markup
      # escapes the & again, so the browser shows a literal "&#39;".
      sprintf(
        paste0("%s | %s | median: %s | deviation: %s cell MADs | sample spread: %s",
               " | pooled spread: %s | row range: %s | robust z: %s"),
        as.character(long$sampleID), long$label, long$tile_text,
        num2(long$z), num3(long$cell_spread), num3(long$pooled_spread),
        num2(long$row_range), num2(long$z_mad)
      )
    } else {
      sprintf(
        "%s | %s | median: %s | robust z: %s",
        as.character(long$sampleID), long$label, long$tile_text, num2(long$z)
      )
    }
    p <- ggplot(long, aes(x = sampleID, y = feature_lab, fill = z)) +
      ggiraph::geom_tile_interactive(
        aes(tooltip = tooltip, data_id = interaction(sampleID, feature)),
        colour = "grey80", linewidth = tile_border
      )
  } else {
    # Grey borders, not white: when every feature is flat across samples the
    # fill is white everywhere and white borders make the grid disappear.
    p <- p + geom_tile(colour = "grey80", linewidth = tile_border)
  }
  if (show_values) {
    p <- p + geom_text(aes(label = tile_text),
                       size = if (is_html) 3.4 else 3.0, colour = "grey15")
  }
  # Only what cannot be inferred from the figure itself: what the printed
  # number is. The colour scale is named in the legend and explained in the
  # documentation; repeating it here just crowds the page.
  caption <- if (show_values) {
    "Tile label: per-sample median."
  } else if (interactive) {
    if (use_cell) "Hover a tile for its median and deviation."
    else "Hover a tile for its median and z-score."
  } else NULL

  # Colour limit adapts to the cohort, up to the cap. A fixed +/-cap would fix
  # the outlier case but break the ordinary one: a cohort whose largest
  # deviation is 0.7 would be drawn entirely in the palest tenth of the
  # gradient. Taking the observed maximum uses the full gradient when nothing
  # is extreme, and only clamps once a sample exceeds the cap.
  zmax <- suppressWarnings(max(abs(long$z), na.rm = TRUE))
  if (!is.finite(zmax) || zmax <= 0) zmax <- 1
  lim <- min(QC_OVERVIEW_DEV_CAP, zmax)
  is_capped <- zmax > QC_OVERVIEW_DEV_CAP

  # Columns faceted by group when there is one, so the boundaries are drawn and
  # the conditions named; row blocks are faceted either way. Block names are
  # wrapped to ~12 chars so each word stacks onto its own line (e.g. "Splice /
  # junction / composition") -- keeps the horizontal strip readable while cutting
  # its width ~45% vs the unwrapped label. Vertical (rotated) strip text was
  # rejected: block heights are data-dependent and a 1-row block (e.g. Transcript
  # length) cannot fit a rotated word without shrinking the font to ~4pt.
  block_labeller <- labeller(block = label_wrap_gen(width = 12))
  p <- p + if (is.null(long$col_group)) {
    facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y",
               labeller = block_labeller)
  } else {
    facet_grid(block ~ col_group, scales = "free", space = "free", switch = "y",
               labeller = block_labeller)
  }

  p +
    scale_fill_gradient2(
      low = "#0072B2", mid = "white", high = "#E69F00", midpoint = 0,
      name = if (use_cell) "Deviation\n(cell MADs)" else "Robust\nz-score",
      limits = c(-lim, lim), oob = scales::squish,
      # Only advertise clamping when something is actually being clamped.
      breaks = if (is_capped) c(-lim, -lim / 2, 0, lim / 2, lim) else waiver(),
      labels = if (is_capped) {
        c(sprintf("<= -%g", lim), -lim / 2, 0, lim / 2, sprintf(">= %g", lim))
      } else waiver()
    ) +
    theme_classic(base_size = base_size) +
    labs(title = "Per-sample QC overview", caption = caption, x = NULL, y = NULL) +
    theme(
      # "plot" rather than the default "panel": centred on the whole figure,
      # not on the tile grid, which the block strips and legend push off-centre.
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = base_size - 2, hjust = 0.5, colour = "grey30"),
      axis.text.x = element_text(size = base_size, angle = x_angle, hjust = 1),
      axis.text.y = element_text(size = base_size - 1),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey92", colour = NA),
      # The SVG renderer draws strip text ~12% wider than the layout measured
      # it, so block names first got clipped mid-word ("Splice junction
      # compositi") and then, with the clip off, spilled out of the grey box.
      # ggplot sizes the strip as measured_width + margins and right-aligns the
      # label at rect_right - r, so the overshoot is (actual - measured) - r:
      # a generous right margin absorbs it outright, and clip = "off" keeps a
      # larger-than-expected metric from ever truncating a word.
      strip.clip = "off",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = base_size - 1,
                                       hjust = 1, margin = margin(r = 30, l = 6)),
      strip.text.x = element_text(face = "bold", size = base_size - 1,
                                  margin = margin(t = 4, b = 4)),
      panel.spacing.y = unit(4, "pt"),
      panel.spacing.x = unit(6, "pt"),
      legend.position = "right"
    )
}

# Figure canvas for the QC overview. The row count is fixed by the registry, so
# height is stable; width has to grow with the sample count or the columns are
# squeezed into slivers. Capped so a 60-sample cohort does not produce a PDF
# page metres wide.
qc_overview_fig_dims <- function(n_samples, n_features) {
  list(
    width  = max(10, min(26, 6 + 0.55 * n_samples)),
    height = max(7, min(20, 3 + 0.42 * n_features))
  )
}

# Wrap the plot as a self-sizing SVG widget for the HTML report.
#
# The static PNG is rendered at a fixed inch size and then shown at that size,
# which is why a 24-row figure needed scrolling. An SVG can be drawn at a
# comfortable aspect ratio and then scaled by the browser without going soft,
# so opts_sizing(rescale = TRUE) is what makes it fit at 100% zoom.
#
# The drawing aspect is deliberately wider and shorter than the PDF canvas.
# rescale = TRUE fits the SVG to the container WIDTH, so the rendered height is
# container_width * (height_svg / width_svg). The content column is NOT the
# full 1300px of .main-container -- the floating TOC takes a chunk, leaving
# roughly 800px -- so the aspect can be close to square and still fit: ~800px
# wide at this ratio is ~720px tall. The stylesheet caps the SVG at 82vh, which
# is what actually protects short windows, so this only needs to avoid being
# gratuitously tall rather than guarantee the fit on its own.
QC_OVERVIEW_SVG_MAX_ASPECT <- 0.9

qc_overview_widget <- function(plot, n_samples, n_features) {
  if (is.null(plot) || !qc_overview_interactive_ok()) return(NULL)
  dims <- qc_overview_fig_dims(n_samples, n_features)
  tryCatch(
    ggiraph::girafe(
      ggobj = plot,
      width_svg = dims$width,
      height_svg = min(dims$height, QC_OVERVIEW_SVG_MAX_ASPECT * dims$width),
      options = list(
        ggiraph::opts_sizing(rescale = TRUE, width = 1),
        ggiraph::opts_tooltip(
          css = paste0("background-color:#333; color:#fff; padding:6px 8px;",
                       "border-radius:4px; font-size:12px; max-width:340px;"),
          opacity = 0.95
        ),
        ggiraph::opts_hover(css = "stroke:#333; stroke-width:1.5px;"),
        ggiraph::opts_toolbar(saveaspng = TRUE)
      )
    ),
    error = function(e) {
      message(sprintf("[WARNING] Could not build the interactive QC overview: %s", e$message))
      NULL
    }
  )
}

# Helper: build the violin + boxplot for one curated feature. Takes a row of the
# curated feature table, so the panel is titled with the registry's own label --
# the same string the QC overview heatmap puts on that feature's row. The block
# only routes the panel to a section (see main()); it is not drawn, because no
# other figure in the report subtitles itself and the section heading already
# names the block.
build_curated_feature_plot <- function(multi, feature_row, sample_levels, is_html = FALSE) {
  feature_name <- as.character(feature_row$feature)[1]
  if (!feature_name %in% colnames(multi)) {
    return(NULL)
  }
  values <- multi[[feature_name]]
  if (!is.numeric(values)) {
    return(NULL)
  }
  plot_df <- multi %>%
    select(sampleID, value = all_of(feature_name)) %>%
    mutate(value = as.numeric(value)) %>%
    filter(is.finite(value))
  if (nrow(plot_df) == 0) {
    return(NULL)
  }
  info <- infer_feature_metadata(feature_name, plot_df$value)
  use_log <- (info$unit == "count")
  if (info$scale_to_percent) {
    plot_df <- plot_df %>% mutate(value = value * 100)
  }
  if (use_log) {
    plot_df <- plot_df %>% filter(value > 0)
  }
  if (length(sample_levels) == 0) {
    sample_levels <- unique(plot_df$sampleID)
  }
  plot_df <- plot_df %>% mutate(sampleID = factor(sampleID, levels = sample_levels))
  uniqueness <- plot_df %>%
    group_by(sampleID) %>%
    summarise(unique_vals = n_distinct(value), .groups = "drop")
  use_violin <- any(uniqueness$unique_vals > 1)
  # The label is why this reads "Reads with non-canonical SJ" rather than the
  # column name's "Non canonical", which #44 established is a different claim.
  # Its trailing unit is dropped -- the y axis already carries it.
  feature_label <- as.character(feature_row$label)[1]
  feature_title <- if (length(feature_label) == 1 && !is.na(feature_label) && nzchar(feature_label)) {
    trimws(sub("\\s*\\((%|bp)\\)$", "", feature_label))
  } else {
    info$display_name
  }
  gp <- ggplot(plot_df, aes(x = sampleID, y = value, fill = sampleID, colour = sampleID))
  if (use_violin) {
    gp <- gp + geom_violin(trim = TRUE, scale = "width", alpha = 0.7, linewidth = 0.3)
  }
  gp <- gp +
    geom_boxplot(
      width = 0.05, outlier.shape = NA,
      alpha = 0.3,
      colour = "grey20"
    ) +
    stat_summary(
      fun = mean, geom = "point", shape = 4, size = 1,
      colour = "red", stroke = 0.45
    ) +
    scale_fill_conesa(palette = "complete", drop = FALSE) +
    scale_color_conesa(palette = "complete", guide = "none", drop = FALSE) +
    labs(
      title = feature_title,
      x = "Sample",
      y = info$value_label
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text.x = element_text(size = 16, angle = 35, hjust = 1),
      axis.text.y = element_text(size = 16)
    )
  if (use_log) {
    gp <- gp + scale_y_log10(labels = scales::comma)
  } else if (info$unit == "%") {
    gp <- gp + coord_cartesian(ylim = c(0, 100))
    finite_vals <- plot_df$value[is.finite(plot_df$value)]
    max_val <- if (length(finite_vals) > 0) max(finite_vals) else 0
    min_val <- if (length(finite_vals) > 0) min(finite_vals) else 0
    gp <- attach_zoom_inset(gp, plot_df, "sampleID", "value",
                             max_val, min_val, is_html, use_violin)
  }
  return(gp)
}

# Reusable helper: violin + boxplot comparing a single metric across samples
build_sample_comparison_plot <- function(data, col_name, title, y_label,
                                         sample_levels, scale_percent = FALSE,
                                         log_scale = FALSE, is_html = FALSE,
                                         y_breaks = NULL) {
  if (!col_name %in% colnames(data)) return(NULL)
  plot_df <- data %>%
    select(sampleID, value = all_of(col_name)) %>%
    mutate(value = as.numeric(value)) %>%
    filter(is.finite(value))
  if (nrow(plot_df) == 0) return(NULL)
  if (scale_percent) plot_df$value <- plot_df$value * 100
  if (isTRUE(log_scale)) plot_df <- plot_df %>% filter(value > 0)
  plot_df$sampleID <- factor(plot_df$sampleID, levels = sample_levels)

  uniqueness <- plot_df %>%
    group_by(sampleID) %>%
    summarise(uv = n_distinct(value), .groups = "drop")
  use_violin <- any(uniqueness$uv > 1)

  gp <- ggplot(plot_df, aes(x = sampleID, y = value, fill = sampleID, colour = sampleID))
  if (use_violin) {
    gp <- gp + geom_violin(trim = TRUE, scale = "width", alpha = 0.7, linewidth = 0.3)
  }
  gp <- gp +
    geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.3, colour = "grey20") +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1,
                 colour = "red", stroke = 0.45) +
    scale_fill_conesa(palette = "complete", drop = FALSE) +
    scale_color_conesa(palette = "complete", guide = "none", drop = FALSE) +
    labs(title = title, x = "Sample", y = y_label) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text.x = element_text(size = 16, angle = 35, hjust = 1),
      axis.text.y = element_text(size = 16)
    )
  is_pct <- grepl("%", y_label) && !isTRUE(log_scale)
  if (isTRUE(log_scale)) {
    gp <- gp + scale_y_log10(
      labels = scales::comma,
      breaks = if (!is.null(y_breaks)) y_breaks else waiver()
    )
  } else if (is_pct) {
    gp <- gp + coord_cartesian(ylim = c(0, 100))
    finite_vals <- plot_df$value[is.finite(plot_df$value)]
    max_val <- if (length(finite_vals) > 0) max(finite_vals) else 0
    min_val <- if (length(finite_vals) > 0) min(finite_vals) else 0
    gp <- attach_zoom_inset(gp, plot_df, "sampleID", "value",
                             max_val, min_val, is_html, use_violin)
  }
  gp
}

main <- function() {
  params <- parse_args()
  is_html_output <- params$report %in% c("html", "both")

  files <- unlist(strsplit(params$files, ",", fixed = TRUE))
  files <- trimws(files)
  files <- files[nchar(files) > 0]
  if (length(files) < 2) {
    message("[INFO] Fewer than 2 files provided. Nothing to do.")
    quit(status = 0)
  }

  # Derive sample IDs from filenames before loading (needed to align group vectors)
  all_sample_ids <- sub("_SQANTI_cell_summary\\.txt(\\.gz)?$", "", basename(files))

  # Parse optional group vectors (parallel to files)
  color_groups_vec <- parse_group_vec(params$color_group, length(files))
  shape_groups_vec <- parse_group_vec(params$shape_group, length(files))
  shade_groups_vec  <- parse_group_vec(params$shade_group,  length(files))

  # Read all summaries, tracking which files loaded successfully
  raw_lst <- lapply(files, safe_read_summary)
  valid_mask <- !sapply(raw_lst, is.null)
  lst <- raw_lst[valid_mask]

  # Apply the same mask to sample IDs and group vectors
  all_sample_ids   <- all_sample_ids[valid_mask]
  if (!is.null(color_groups_vec)) color_groups_vec <- color_groups_vec[valid_mask]
  if (!is.null(shape_groups_vec)) shape_groups_vec <- shape_groups_vec[valid_mask]
  if (!is.null(shade_groups_vec))  shade_groups_vec  <- shade_groups_vec[valid_mask]

  if (length(lst) < 2) {
    message("[INFO] Fewer than 2 valid summaries after reading. Skipping.")
    quit(status = 0)
  }

  # Build sample group map for PCA annotation
  sample_group_map <- data.frame(sampleID = all_sample_ids, stringsAsFactors = FALSE)
  if (!is.null(color_groups_vec)) sample_group_map$color_group <- color_groups_vec
  if (!is.null(shape_groups_vec)) sample_group_map$shape_group <- shape_groups_vec
  if (!is.null(shade_groups_vec))  sample_group_map$shade_group  <- shade_groups_vec

  # Harmonize columns: union of all names; fill missing with 0 for numeric, "" otherwise
  all_cols <- Reduce(union, lapply(lst, colnames))
  # Ensure CB and sampleID exist in final order front
  all_cols <- unique(c("CB", setdiff(all_cols, "CB")))
  all_cols <- unique(c(all_cols, "sampleID"))

  norm_list <- lapply(lst, function(df) {
    missing <- setdiff(all_cols, colnames(df))
    for (m in missing) {
      df[[m]] <- if (m == "CB" || m == "sampleID") "" else 0
    }
    # Reorder
    df <- df[, all_cols]
    df
  })
  multi <- bind_rows(norm_list)

  # Sample order for tables and plots: same as design / --files order (cell summary path order),
  # not alphabetical. dplyr::summarise() group order can differ from row order.
  ids_in_data <- unique(multi$sampleID[!is.na(multi$sampleID)])
  sample_levels_global <- all_sample_ids[all_sample_ids %in% ids_in_data]

  render_pdf <- params$report %in% c("pdf", "both")
  render_html <- params$report %in% c("html", "both")

  # Basic cohort-level aggregates per sample
  count_col <- if (params$mode == "isoforms") "Transcripts_in_cell" else "Reads_in_cell"

  per_sample_stats <- multi %>%
    group_by(sampleID) %>%
    summarise(
      cells = n_distinct(CB),
      median_reads = median(.data[[count_col]], na.rm = TRUE),
      median_umis = if ("UMIs_in_cell" %in% names(.)) median(UMIs_in_cell, na.rm = TRUE) else NA,
      median_genes = median(Genes_in_cell, na.rm = TRUE),
      median_annotated = median(Annotated_genes, na.rm = TRUE),
      median_ujc = if ("UJCs_in_cell" %in% names(.)) median(UJCs_in_cell, na.rm = TRUE) else NA
    )

  # Compute per-sample length statistics from the cell summary in both modes.
  # Median_length_per_cell is one value per real cell (the same source as the
  # per-cell median-length violin), so the summary table agrees with the violin.
  # We deliberately do NOT re-derive this from the classification file: that
  # would give a pooled median over all reads/transcripts (and in isoforms mode
  # the comma-separated FL breaks as.integer() outright), which disagrees with
  # the violin. The raw length *distribution* still reads the classification
  # file below, since per-read lengths cannot be recovered from a per-cell median.
  len_stats <- NULL
  if ("Median_length_per_cell" %in% colnames(multi)) {
    len_stats <- multi %>%
      filter(!is.na(Median_length_per_cell)) %>%
      group_by(sampleID) %>%
      summarise(
        median_length = median(Median_length_per_cell, na.rm = TRUE),
        iqr_length    = IQR(Median_length_per_cell, na.rm = TRUE),
        .groups = "drop"
      )
  }
  if (!is.null(len_stats)) {
    per_sample_stats <- left_join(per_sample_stats, len_stats, by = "sampleID")
  } else {
    per_sample_stats$median_length <- NA_real_
    per_sample_stats$iqr_length <- NA_real_
  }

  # Match Samples Summary row order to sample_levels_global (figures use the same levels)
  per_sample_stats <- per_sample_stats %>%
    mutate(sampleID = factor(.data$sampleID, levels = sample_levels_global)) %>%
    dplyr::arrange(.data$sampleID) %>%
    mutate(sampleID = as.character(.data$sampleID))

  entity_label_plural <- if (params$mode == "isoforms") "Transcripts" else "Reads"

  summary_tbl <- per_sample_stats %>%
    mutate(across(where(is.numeric), ~ round(., 1))) %>%
    transmute(
      Sample = sampleID,
      `Cell\nBarcodes` = cells,
      `Median\nReads/Cell` = median_reads,
      `Median\nUMIs/Cell` = median_umis,
      `Median Annotated\nGenes/Cell` = median_annotated,
      `Median\nUJCs/Cell` = median_ujc,
      `Median\nLength\n(bp)` = median_length,
      `Length\nIQR\n(bp)` = iqr_length
    )

  # Rename reads/transcripts column dynamically to match the mode
  colnames(summary_tbl)[colnames(summary_tbl) == "Median\nReads/Cell"] <- paste0("Median\n", entity_label_plural, "/Cell")

  if (params$mode == "isoforms") {
    summary_tbl <- summary_tbl %>% select(-`Median\nUMIs/Cell`, -`Median\nUJCs/Cell`)
  }

  # Drop length columns if no classification files were provided
  if (all(is.na(summary_tbl$`Median\nLength\n(bp)`))) {
    summary_tbl <- summary_tbl %>% select(-`Median\nLength\n(bp)`, -`Length\nIQR\n(bp)`)
  }

  summary_tbl_html <- summary_tbl
  colnames(summary_tbl_html) <- gsub("\\n", "<br>", colnames(summary_tbl_html), fixed = TRUE)

  assign("multi_per_sample_stats", per_sample_stats, envir = .GlobalEnv)
  assign("multi_summary_tbl_pdf", summary_tbl, envir = .GlobalEnv)
  assign("multi_summary_tbl_html", summary_tbl_html, envir = .GlobalEnv)
  assign("entity_label", if (params$mode == "isoforms") "Transcript" else "Read", envir = .GlobalEnv)
  assign("entity_label_plural", entity_label_plural, envir = .GlobalEnv)
  assign("mode", params$mode, envir = .GlobalEnv)

  # Derived before any plotting: the per-unit rates are shown as figures in the
  # sections below as well as feeding the curated feature set.
  derived_res <- add_derived_features(multi, lst, params$mode)
  multi <- derived_res$multi

  # Each count is followed by its per-unit rate. A raw count tracks how deeply
  # the cell was sequenced, so it cannot be compared across samples of different
  # depth; the rate is what survives that. These are the same derived columns the
  # curated feature set and the QC overview use, and they are linear, not log:
  # they sit inside (0, 1] rather than spanning orders of magnitude.
  add_rate_spec <- function(specs, nm, col, y) {
    if (col %in% colnames(multi)) specs[[nm]] <- list(col = col, y = y)
    specs
  }
  rate_ent <- if (params$mode == "isoforms") "Transcript" else "Read"
  rate_unit <- tolower(rate_ent)

  # -------- Per-Cell Library Size plots --------
  # Section membership follows the per-sample report's taxonomy: the entity
  # count, UMIs and the mode's distinct-structure count are all measures of what
  # the cell's library holds. Anything gene-related lives in gene
  # characterisation below, including the mitochondrial share.
  library_size_specs <- list()
  library_size_specs[[paste0(entity_label_plural, " per Cell")]] <- list(
    col = count_col, y = paste0(entity_label_plural, ", count"), log = TRUE
  )
  if (params$mode == "reads" && "UMIs_in_cell" %in% colnames(multi)) {
    library_size_specs[["UMIs per Cell"]] <- list(col = "UMIs_in_cell", y = "UMIs, count", log = TRUE)
  }
  if (params$mode == "reads" && "UJCs_in_cell" %in% colnames(multi)) {
    library_size_specs[["UJCs per Cell"]] <- list(col = "UJCs_in_cell", y = "UJCs, count", log = TRUE)
  }
  if (params$mode == "isoforms" && "Isoforms_in_cell" %in% colnames(multi)) {
    library_size_specs[["Isoforms per Cell"]] <- list(
      col = "Isoforms_in_cell", y = "Isoforms, count", log = TRUE
    )
  }
  rate_div <- if (params$mode == "isoforms") "Isoforms" else "UJCs"
  library_size_specs <- add_rate_spec(
    library_size_specs, paste0(rate_div, " per ", rate_ent),
    paste0(rate_div, "_per_", rate_unit), paste0(rate_div, " per ", rate_unit)
  )


  multi_library_size_plots <- list()
  for (nm in names(library_size_specs)) {
    spec <- library_size_specs[[nm]]
    gp <- build_sample_comparison_plot(
      multi, spec$col,
      title = paste0("Per Sample ", nm, " Distribution"),
      y_label = spec$y,
      sample_levels = sample_levels_global,
      scale_percent = isTRUE(spec$scale_pct),
      log_scale = isTRUE(spec$log),
      is_html = is_html_output
    )
    if (!is.null(gp)) multi_library_size_plots[[nm]] <- gp
  }
  if (length(multi_library_size_plots) > 0) {
    assign("multi_library_size_plots", multi_library_size_plots, envir = .GlobalEnv)
  }

  # -------- Gene Characterization plots --------
  # Everything gene-related, mitochondrial content included: MT_perc is
  # MT-gene reads over total reads, so it belongs with the genes rather than
  # with library size.
  gene_char_specs <- list()
  if ("Genes_in_cell" %in% colnames(multi)) {
    gene_char_specs[["Genes per Cell"]] <- list(col = "Genes_in_cell", y = "Genes, count", log = TRUE)
  }
  if ("Annotated_genes" %in% colnames(multi)) {
    gene_char_specs[["Annotated Genes per Cell"]] <- list(col = "Annotated_genes", y = "Genes, count", log = TRUE)
  }
  gene_char_specs <- add_rate_spec(
    gene_char_specs, paste0("Annotated Genes per ", rate_ent),
    paste0("Annotated_genes_per_", rate_unit), paste0("Genes per ", rate_unit)
  )
  if ("MT_perc" %in% colnames(multi)) {
    gene_char_specs[["Mitochondrial % per Cell"]] <- list(
      col = "MT_perc", y = paste0(entity_label_plural, ", %")
    )
  }

  multi_gene_char_plots <- list()
  for (nm in names(gene_char_specs)) {
    spec <- gene_char_specs[[nm]]
    gp <- build_sample_comparison_plot(
      multi, spec$col,
      title = paste0("Per Sample ", nm, " Distribution"),
      y_label = spec$y,
      sample_levels = sample_levels_global,
      log_scale = isTRUE(spec$log),
      is_html = is_html_output
    )
    if (!is.null(gp)) multi_gene_char_plots[[nm]] <- gp
  }
  if (length(multi_gene_char_plots) > 0) {
    assign("multi_gene_char_plots", multi_gene_char_plots, envir = .GlobalEnv)
  }

  # Weighted quantile via linear interpolation on the cumulative weight.
  # Lets the bulk length distribution reflect FL (expression) without exploding
  # rows: a model of length L with total FL n counts as n observations of L.
  weighted_quantile <- function(x, w, probs) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    x <- x[ok]; w <- w[ok]
    if (length(x) == 0) return(rep(NA_real_, length(probs)))
    o <- order(x); x <- x[o]; w <- w[o]
    cw <- (cumsum(w) - 0.5 * w) / sum(w)
    stats::approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
  }

  # -------- Length Distribution plots --------
  multi_length_plots_local <- list()
  if (!is.null(params$class_files) && nzchar(params$class_files)) {
    class_file_paths <- trimws(unlist(strsplit(params$class_files, ",", fixed = TRUE)))
    class_file_paths <- class_file_paths[nchar(class_file_paths) > 0 & file.exists(class_file_paths)]

    if (length(class_file_paths) >= 2) {
      # Read length, and in isoforms mode the FL abundance, so the distribution
      # can be FL-weighted (expression-level), consistent with the per-cell
      # median length and with reads mode where each read is one observation.
      # FL is comma-separated per cell in isoforms mode, so sum it to a per-model
      # total; as.integer() on the raw field would NA out (the length version of
      # the median-length bug) and silently drop to unweighted.
      len_dfs <- lapply(class_file_paths, function(f) {
        sel <- if (params$mode == "isoforms") c("length", "FL") else c("length")
        df <- tryCatch(
          data.table::fread(f, select = sel, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, data.table = FALSE),
          error = function(e) { message("[WARNING] Could not read ", f, ": ", e$message); NULL }
        )
        if (is.null(df) || nrow(df) == 0) return(NULL)
        df$length <- suppressWarnings(as.numeric(df$length))
        if ("FL" %in% colnames(df)) {
          df$w <- vapply(strsplit(as.character(df$FL), ",", fixed = TRUE),
                         function(x) sum(as.numeric(x), na.rm = TRUE), numeric(1))
          df$FL <- NULL
        } else {
          df$w <- 1  # reads mode: one row per read is already expression-level
        }
        df <- df[is.finite(df$length) & df$length > 0 & is.finite(df$w) & df$w > 0, , drop = FALSE]
        sample_id <- sub("_classification\\.txt(\\.gz)?$", "", basename(f))
        df$sampleID <- sample_id
        df
      })
      len_combined <- bind_rows(Filter(Negate(is.null), len_dfs))

      if (nrow(len_combined) > 0) {
        len_combined$sampleID <- factor(len_combined$sampleID, levels = sample_levels_global)

        uniqueness <- len_combined %>%
          group_by(sampleID) %>%
          summarise(uv = n_distinct(length), .groups = "drop")
        use_violin <- any(uniqueness$uv > 1)

        # FL-weighted box statistics + mean, drawn with stat = "identity" so the
        # box and mean point match the (also FL-weighted) violin density.
        box_stats <- len_combined %>%
          group_by(sampleID) %>%
          summarise(
            lower  = weighted_quantile(length, w, 0.25),
            middle = weighted_quantile(length, w, 0.50),
            upper  = weighted_quantile(length, w, 0.75),
            wmean  = sum(length * w) / sum(w),
            .groups = "drop"
          )
        box_stats <- len_combined %>%
          left_join(box_stats, by = "sampleID") %>%
          group_by(sampleID) %>%
          summarise(
            lower  = lower[1],
            middle = middle[1],
            upper  = upper[1],
            wmean  = wmean[1],
            # whiskers: most extreme lengths within 1.5*IQR of the box
            ymin = min(length[length >= lower[1] - 1.5 * (upper[1] - lower[1])]),
            ymax = max(length[length <= upper[1] + 1.5 * (upper[1] - lower[1])]),
            .groups = "drop"
          )
        box_stats$sampleID <- factor(box_stats$sampleID, levels = sample_levels_global)

        gp_len <- ggplot(len_combined, aes(x = sampleID, y = length, fill = sampleID, colour = sampleID))
        if (use_violin) {
          gp_len <- gp_len + geom_violin(aes(weight = w), trim = TRUE, scale = "width",
                                         alpha = 0.7, linewidth = 0.3)
        }
        gp_len <- gp_len +
          geom_boxplot(
            data = box_stats,
            aes(x = sampleID, ymin = ymin, lower = lower, middle = middle,
                upper = upper, ymax = ymax, fill = sampleID),
            stat = "identity", width = 0.05, alpha = 0.3, colour = "grey20",
            inherit.aes = FALSE
          ) +
          geom_point(
            data = box_stats, aes(x = sampleID, y = wmean),
            shape = 4, size = 1, colour = "red", stroke = 0.45, inherit.aes = FALSE
          ) +
          scale_fill_conesa(palette = "complete", drop = FALSE) +
          scale_color_conesa(palette = "complete", guide = "none", drop = FALSE) +
          scale_y_log10(
            breaks = c(50, 100, 250, 500, 1000, 2000, 5000, 10000, 20000, 50000),
            labels = scales::comma
          ) +
          annotation_logticks(sides = "l") +
          coord_cartesian(clip = "off") +
          labs(
            title = paste0("Per Sample ", entity_label_plural, " Length Distribution"),
            x = "Sample", y = "Length, bp"
          ) +
          theme_classic(base_size = 14) +
          theme(
            legend.position = "none",
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 18),
            axis.text.x = element_text(size = 16, angle = 35, hjust = 1),
            axis.text.y = element_text(size = 16),
            plot.margin = margin(t = 5, r = 5, b = 5, l = 10, unit = "pt")
          )
        multi_length_plots_local[["Bulk Length Distribution"]] <- gp_len
      }
    }
  }

  if ("Median_length_per_cell" %in% colnames(multi)) {
    gp <- build_sample_comparison_plot(
      multi, "Median_length_per_cell",
      title = paste0("Per Sample Median Length per Cell Distribution"),
      y_label = "Length, bp",
      sample_levels = sample_levels_global,
      log_scale = TRUE,
      is_html = is_html_output,
      # Denser log-spaced ticks adapted to the (narrow) per-cell median range.
      # A fixed break vector like the bulk plot's would fall mostly outside this
      # range and leave only 2-3 ticks; breaks_log() fills the actual range.
      y_breaks = scales::breaks_log(n = 8)
    )
    if (!is.null(gp)) multi_length_plots_local[["Median Length per Cell"]] <- gp
  }

  if (length(multi_length_plots_local) > 0) {
    assign("multi_length_plots", multi_length_plots_local, envir = .GlobalEnv)
  }


  # Output path
  out_dir <- params$out_dir
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  pdf_out <- file.path(out_dir, paste0(params$prefix, ".pdf"))

  if (render_html) {
    message("[INFO] HTML report requested; plots will be rendered via R Markdown template.")
  }

  # Structural category proportions, if available
  cat_cols <- c("FSM_prop", "ISM_prop", "NIC_prop", "NNC_prop", "Genic_Genomic_prop", "Antisense_prop", "Fusion_prop", "Intergenic_prop", "Genic_intron_prop")
  have_cats <- all(cat_cols %in% colnames(multi))
  if (have_cats) {
    cats_long <- multi %>%
      select(all_of(c("sampleID", cat_cols))) %>%
      pivot_longer(cols = all_of(cat_cols), names_to = "category", values_to = "prop") %>%
      mutate(
        sampleID = factor(sampleID, levels = sample_levels_global),
        category = factor(category,
          levels = c(
            "FSM_prop", "ISM_prop", "NIC_prop", "NNC_prop",
            "Genic_Genomic_prop", "Antisense_prop", "Fusion_prop", "Intergenic_prop", "Genic_intron_prop"
          ),
          labels = c(
            "FSM", "ISM", "NIC", "NNC",
            "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic Intron"
          )
        )
      ) %>%
      mutate(prop = suppressWarnings(as.numeric(prop)))

    category_levels <- levels(cats_long$category)
    sample_levels_all <- levels(cats_long$sampleID)
    cat_to_col <- c(
      "FSM" = "#6BAED6",
      "ISM" = "#FC8D59",
      "NIC" = "#78C679",
      "NNC" = "#EE6A50",
      "Genic Genomic" = "#969696",
      "Antisense" = "#66C2A4",
      "Fusion" = "goldenrod1",
      "Intergenic" = "darksalmon",
      "Genic Intron" = "#41B6C4"
    )

    category_plots <- lapply(category_levels, function(cat_lab) {
      dfp <- cats_long %>% filter(category == cat_lab)
      cat_col <- unname(cat_to_col[[as.character(cat_lab)]])
      if (is.null(cat_col) || is.na(cat_col)) cat_col <- "grey60"
      box_outline_col <- if (as.character(cat_lab) == "Genic Genomic") "grey90" else "grey20"
      violin_fill <- grDevices::adjustcolor(cat_col, alpha.f = 0.7)
      gg <- ggplot(dfp, aes(x = sampleID, y = prop)) +
        geom_violin(fill = violin_fill, color = cat_col, linewidth = 0.3, width = 0.8, trim = TRUE, scale = "width") +
        geom_boxplot(width = 0.05, outlier.shape = NA, fill = cat_col, color = box_outline_col, alpha = 0.3) +
        stat_summary(fun = mean, geom = "point", shape = 4, size = 1, colour = "red", stroke = 0.9) +
        scale_y_continuous(limits = c(0, 100), expand = expansion(add = c(1, 0))) +
        theme_classic(base_size = 13) +
        labs(title = paste0("Per Sample ", cat_lab, " Reads Distribution Across Cells"), x = "Sample", y = "Reads, %") +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 18),
          axis.text.y = element_text(size = 16)
        )
      cat_finite <- dfp$prop[is.finite(dfp$prop)]
      cat_max <- if (length(cat_finite) > 0) max(cat_finite) else 0
      cat_min <- if (length(cat_finite) > 0) min(cat_finite) else 0
      attach_zoom_inset(gg, dfp, "sampleID", "prop",
                        cat_max, cat_min, is_html_output,
                        use_violin = TRUE, fill_col = violin_fill)
    })
    names(category_plots) <- as.character(category_levels)

    assign("multi_structural_category_violin_plots", category_plots, envir = .GlobalEnv)

  }

  multi_pca_scores_plot_local <- NULL

  multi_pca_group_plot_local <- NULL

  multi_pca_scree_plot_local <- NULL

  multi_pca_top_loadings_plots_local <- NULL

  multi_sj_composition_plots_local <- NULL

  multi_good_quality_plots_local <- NULL

  multi_bad_quality_plots_local <- NULL

  multi_extra_feature_plots_local <- NULL

  multi_qc_overview_plot_local <- NULL

  # -------- Curated feature set (drives the QC overview and the PCA) --------
  # Selecting "all numeric columns" weights each feature family by its column
  # count and lets run-flag sentinels in as real signal; curated_feature_table()
  # resolves the registry against the data and reports every exclusion.
  feature_map <- curated_feature_table(lst, multi, params$mode, params$pca_features,
                                       derived = derived_res$added)
  pca_cols <- as.character(feature_map$feature)
  if (length(pca_cols) > 0) {
    message(sprintf(
      "[INFO] Using %d curated feature(s) across %d block(s): %s",
      length(pca_cols), nlevels(feature_map$block),
      paste(levels(feature_map$block), collapse = ", ")
    ))
  } else {
    message("[WARNING] No curated features available; skipping the QC overview and PCA.")
  }

  # -------- Per-feature distributions for the blocks that lack a section -----
  # Every curated feature gets a panel, in registry order, EXCEPT the blocks
  # already drawn earlier in the report (QC_BLOCKS_PLOTTED_ELSEWHERE). The panels
  # used to be the union of the top-10 |loading| on PC1 and PC2, which let the
  # PCA decide which distributions the user is shown -- and loading rank answers
  # "which features define the largest axis of variance", not "which features
  # differ between samples", which is what these figures are read for.
  block_plots <- list()
  for (idx in seq_len(nrow(feature_map))) {
    blk <- as.character(feature_map$block[idx])
    if (blk %in% QC_BLOCKS_PLOTTED_ELSEWHERE) next
    gp_feat <- build_curated_feature_plot(multi, feature_map[idx, ], sample_levels_global,
                                          is_html = is_html_output)
    if (is.null(gp_feat)) {
      message(sprintf("[INFO] Skipping curated feature %s due to missing or constant data.",
                      feature_map$feature[idx]))
      next
    }
    block_plots[[blk]] <- c(block_plots[[blk]],
                            stats::setNames(list(gp_feat), as.character(feature_map$label[idx])))
  }
  if (length(block_plots) > 0) {
    message(sprintf("[INFO] Feature distribution section(s) built for: %s",
                    paste(names(block_plots), collapse = ", ")))
  }
  multi_sj_composition_plots_local <- block_plots[[QC_BLOCK_JUNCTIONS]]
  multi_good_quality_plots_local <- block_plots[[QC_BLOCK_GOOD]]
  multi_bad_quality_plots_local <- block_plots[[QC_BLOCK_BAD]]
  # Anything left over is a block the layout does not know about -- in practice
  # the "User-specified" block a --pca_features whitelist produces. It still gets
  # figures, in one catch-all section, rather than silently having none.
  extra <- block_plots[setdiff(names(block_plots),
                               c(QC_BLOCK_JUNCTIONS, QC_BLOCK_GOOD, QC_BLOCK_BAD))]
  multi_extra_feature_plots_local <- do.call(c, unname(extra))
  for (nm in c("multi_sj_composition_plots", "multi_good_quality_plots",
               "multi_bad_quality_plots", "multi_extra_feature_plots")) {
    v <- get(paste0(nm, "_local"))
    if (length(v) > 0) assign(nm, v, envir = .GlobalEnv)
  }

  # Aggregate per-sample medians across retained features
  agg_median <- multi %>%
    group_by(sampleID) %>%
    summarise(across(all_of(pca_cols), ~ median(., na.rm = TRUE)), .groups = "drop")
  agg_median$sampleID <- factor(agg_median$sampleID, levels = sample_levels_global)
  agg_median <- agg_median[order(agg_median$sampleID), , drop = FALSE]
  agg_median$sampleID <- as.character(agg_median$sampleID)

  # Write the per-sample feature medians table so users can inspect / reuse PCA input
  medians_out <- file.path(out_dir, paste0(params$prefix, "_pca_feature_medians.tsv"))
  tryCatch(
    {
      write.table(agg_median, file = medians_out, sep = "\t", row.names = FALSE, quote = FALSE)
      message(sprintf("[INFO] PCA feature medians written to: %s", medians_out))
    },
    error = function(e) message(sprintf("[WARNING] Could not write PCA feature medians: %s", e$message))
  )

  # -------- Per-sample QC overview heatmap --------
  # Two renders from one builder: the PDF gets plain geom_tile, the HTML gets
  # the interactive tiles wrapped in a self-sizing SVG widget. The PDF must
  # stay static, so it cannot simply reuse the interactive object.
  if (length(pca_cols) >= 2 && nrow(agg_median) >= 2) {
    cell_spread_tbl <- tryCatch(
      qc_cell_spread_table(multi, pca_cols),
      error = function(e) {
        message(sprintf("[WARNING] Could not compute per-cell spreads: %s", e$message))
        NULL
      }
    )
    build_overview <- function(interactive) {
      tryCatch(
        build_qc_overview_plot(agg_median, feature_map, sample_levels_global,
                               is_html = is_html_output,
                               sample_group_map = sample_group_map,
                               interactive = interactive,
                               cell_spread = cell_spread_tbl),
        error = function(e) {
          message(sprintf("[WARNING] Could not build the QC overview heatmap: %s", e$message))
          NULL
        }
      )
    }
    multi_qc_overview_plot_local <- build_overview(interactive = FALSE)
    if (!is.null(multi_qc_overview_plot_local)) {
      assign("multi_qc_overview_plot", multi_qc_overview_plot_local, envir = .GlobalEnv)
      assign("multi_qc_overview_dims",
             qc_overview_fig_dims(nrow(agg_median), length(pca_cols)), envir = .GlobalEnv)
    }
    if (is_html_output) {
      if (!qc_overview_interactive_ok()) {
        message("[INFO] ggiraph not available; the QC overview will be a static figure.")
      }
      widget <- qc_overview_widget(build_overview(interactive = TRUE),
                                   nrow(agg_median), length(pca_cols))
      if (!is.null(widget)) assign("multi_qc_overview_widget", widget, envir = .GlobalEnv)
    }
  }

  if (nrow(agg_median) >= 2 && ncol(agg_median) >= 2) {
    # Drop features with zero variance across samples
    feat_sds <- sapply(agg_median %>% select(-sampleID), function(x) stats::sd(x, na.rm = TRUE))
    feat_keep <- names(feat_sds)[is.finite(feat_sds) & !is.na(feat_sds) & feat_sds > 0]

    if (length(feat_keep) >= 2) {
      mat <- as.matrix(agg_median[, feat_keep, drop = FALSE])
      rownames(mat) <- agg_median$sampleID
      pca_fit <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
      var_expl <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)

      # A) PC1–PC2 scatter (first among PCA plots)
      if (ncol(pca_fit$x) >= 2) {
        scores <- as.data.frame(pca_fit$x)
        scores$sampleID <- rownames(scores)
        gp_scores <- ggplot(scores, aes(x = PC1, y = PC2, colour = sampleID, label = sampleID)) +
          geom_point(size = 4.5, alpha = 0.95, shape = 19, stroke = 0) +
          scale_color_conesa(palette = "complete") +
          theme_classic(base_size = 14) +
          labs(
            title = "PCA Plot Based on sampleID",
            x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
            y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2])
          ) +
          scale_x_continuous(labels = function(x) sprintf("%.2f", x)) +
          scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
          theme(
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 13),
            legend.key = element_blank(),
            legend.margin = margin(t = 16),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16)
          ) +
          guides(colour = guide_legend(override.aes = list(size = 5, alpha = 0.95, stroke = 0)))
        multi_pca_scores_plot_local <- gp_scores
        assign("multi_pca_scores_plot", gp_scores, envir = .GlobalEnv)

        # A2) Group-annotated PCA (only when group columns are present)
        has_color_grp <- "color_group" %in% colnames(sample_group_map) &&
                         any(nzchar(sample_group_map$color_group))
        has_shape_grp <- "shape_group" %in% colnames(sample_group_map) &&
                         any(nzchar(sample_group_map$shape_group))
        has_shade_grp <- "shade_group"  %in% colnames(sample_group_map) &&
                         any(nzchar(sample_group_map$shade_group))

        if (has_color_grp || has_shape_grp || has_shade_grp) {
          scores_grp <- dplyr::left_join(scores, sample_group_map, by = "sampleID")

          # If color_group was not supplied, use a single placeholder so that
          # color_is_trivial = TRUE fires and suppresses the uninformative legend.
          if (!has_color_grp) scores_grp$color_group <- "group"

          unique_color_vals <- sort(unique(scores_grp$color_group))
          n_color_vals      <- length(unique_color_vals)
          base_hues         <- get_conesa_palette_colors(n_color_vals)
          names(base_hues)  <- unique_color_vals

          shape_palette <- c(19L, 17L, 15L, 18L, 8L, 10L, 13L)

          pca_x_label <- sprintf("PC1 (%.1f%%)", 100 * var_expl[1])
          pca_y_label <- sprintf("PC2 (%.1f%%)", 100 * var_expl[2])

          pca_group_theme <- list(
            theme_classic(base_size = 14),
            scale_x_continuous(labels = function(x) sprintf("%.2f", x)),
            scale_y_continuous(labels = function(x) sprintf("%.2f", x)),
            labs(title = "PCA Plot - Group Annotation",
                 x = pca_x_label, y = pca_y_label),
            theme(
              legend.position  = "bottom",
              legend.text      = element_text(size = 12),
              legend.key       = element_blank(),
              legend.margin    = margin(t = 12),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 18),
              axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 16)
            )
          )

          # When color_group has only 1 unique value the colour channel carries no
          # information: suppress its legend so only the informative channels show.
          color_is_trivial <- (n_color_vals == 1)
          color_guide <- guide_legend(override.aes = list(size = 5, shape = 15, stroke = 0))

          if (has_shade_grp) {
            shade_levels    <- sort(unique(scores_grp$shade_group))
            n_shades        <- length(shade_levels)
            lighten_amounts <- seq(0.45, 0, length.out = n_shades)

            if (color_is_trivial) {
              # Shade alone distinguishes — drop the colour prefix from legend keys
              color_map <- setNames(
                vapply(lighten_amounts, function(a) lighten_hex(base_hues[[1]], a), character(1)),
                shade_levels
              )
              scores_grp$color_shade_key <- factor(scores_grp$shade_group, levels = shade_levels)
              color_scale <- scale_colour_manual(values = color_map, name = "Shade group", drop = FALSE)
            } else {
              color_map <- character(0)
              for (cg in unique_color_vals) {
                for (si in seq_along(shade_levels)) {
                  key <- paste(cg, shade_levels[si], sep = " | ")
                  color_map[key] <- lighten_hex(base_hues[[cg]], lighten_amounts[si])
                }
              }
              legend_key_order <- unlist(lapply(unique_color_vals, function(cg) {
                paste(cg, shade_levels, sep = " | ")
              }))
              scores_grp$color_shade_key <- factor(
                paste(scores_grp$color_group, scores_grp$shade_group, sep = " | "),
                levels = legend_key_order
              )
              color_scale <- scale_colour_manual(values = color_map, name = "Group | Shade",
                                                 drop = FALSE)
            }

            if (has_shape_grp) {
              unique_shape_vals <- sort(unique(scores_grp$shape_group))
              shape_map <- setNames(shape_palette[seq_len(min(length(unique_shape_vals), length(shape_palette)))],
                                    unique_shape_vals)
              scores_grp$shape_group <- factor(scores_grp$shape_group, levels = unique_shape_vals)
              gp_grp <- ggplot(scores_grp,
                aes(x = PC1, y = PC2, colour = color_shade_key, shape = shape_group)) +
                geom_point(size = 4.5, alpha = 0.95, stroke = 0.5) +
                color_scale +
                scale_shape_manual(values = shape_map, name = "Shape group") +
                guides(colour = color_guide,
                       shape  = guide_legend(override.aes = list(size = 5, colour = "grey20")))
            } else {
              gp_grp <- ggplot(scores_grp,
                aes(x = PC1, y = PC2, colour = color_shade_key)) +
                geom_point(size = 4.5, alpha = 0.95, shape = 19L, stroke = 0) +
                color_scale +
                guides(colour = color_guide)
            }

          } else {
            # No shade — plain hues only
            scores_grp$color_group <- factor(scores_grp$color_group, levels = unique_color_vals)
            color_scale <- scale_colour_manual(values = base_hues, name = "Color group", drop = FALSE)

            if (has_shape_grp) {
              unique_shape_vals <- sort(unique(scores_grp$shape_group))
              shape_map <- setNames(shape_palette[seq_len(min(length(unique_shape_vals), length(shape_palette)))],
                                    unique_shape_vals)
              scores_grp$shape_group <- factor(scores_grp$shape_group, levels = unique_shape_vals)
              gp_grp <- ggplot(scores_grp,
                aes(x = PC1, y = PC2, colour = color_group, shape = shape_group)) +
                geom_point(size = 4.5, alpha = 0.95, stroke = 0.5) +
                color_scale +
                scale_shape_manual(values = shape_map, name = "Shape group") +
                guides(
                  colour = if (color_is_trivial) "none" else color_guide,
                  shape  = guide_legend(override.aes = list(size = 5, colour = "grey20"))
                )
            } else {
              gp_grp <- ggplot(scores_grp,
                aes(x = PC1, y = PC2, colour = color_group)) +
                geom_point(size = 4.5, alpha = 0.95, shape = 19L, stroke = 0) +
                color_scale +
                guides(colour = if (color_is_trivial) "none" else color_guide)
            }
          }

          gp_grp <- Reduce(`+`, c(list(gp_grp), pca_group_theme))
          multi_pca_group_plot_local <- gp_grp
          assign("multi_pca_group_plot", gp_grp, envir = .GlobalEnv)
        }

      }

      # B) Scree plot (second)
      k <- min(length(var_expl), 10)
      scree_df <- data.frame(
        PC = factor(paste0("PC", seq_len(k)), levels = paste0("PC", seq_len(k))),
        Proportion = var_expl[seq_len(k)],
        Cumulative = cumsum(var_expl)[seq_len(k)]
      )
      gp_scree <- ggplot(scree_df, aes(x = PC)) +
        geom_col(aes(y = Proportion, fill = "Proportion"), width = 0.8, colour = NA) +
        geom_point(aes(y = Cumulative, colour = "Cumulative"), size = 2.2) +
        geom_line(aes(y = Cumulative, colour = "Cumulative", group = 1), linewidth = 0.6) +
        scale_fill_manual(values = c("Proportion" = "#6BAED6"), name = "") +
        scale_color_manual(values = c("Cumulative" = "#4D4D4D"), name = "") +
        theme_classic(base_size = 14) +
        labs(title = "PCA scree plot", y = "Variance explained", x = "Principal component") +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.position = "bottom",
          legend.margin = margin(t = 20)
        )
      multi_pca_scree_plot_local <- gp_scree
      assign("multi_pca_scree_plot", gp_scree, envir = .GlobalEnv)

      # C) Top loadings for PC1 and PC2 (third)
      if (ncol(pca_fit$rotation) >= 2) {
        rot <- as.data.frame(pca_fit$rotation)
        rot$variable <- rownames(rot)
        top_n <- 10L
        pick_top <- function(colname, n = top_n) {
          ord <- order(abs(rot[[colname]]), decreasing = TRUE)
          df <- head(rot[ord, c("variable", colname)], n)
          colnames(df) <- c("variable", "loading")
          df %>%
            mutate(
              variable = as.character(variable),
              PC = colname,
              rank = dplyr::row_number(),
              sign = if_else(loading >= 0, "Positive", "Negative"),
              abs_loading = abs(loading)
            )
        }
        top_pc1 <- pick_top("PC1")
        top_pc2 <- pick_top("PC2")

        top_pc1_plot <- top_pc1
        top_pc2_plot <- top_pc2
        top_pc1_plot$variable <- factor(top_pc1_plot$variable, levels = rev(top_pc1_plot$variable))
        top_pc2_plot$variable <- factor(top_pc2_plot$variable, levels = rev(top_pc2_plot$variable))
        top_pc1_plot$sign <- factor(top_pc1_plot$sign, levels = c("Positive", "Negative"))
        top_pc2_plot$sign <- factor(top_pc2_plot$sign, levels = c("Positive", "Negative"))

        gp_load1 <- ggplot(top_pc1_plot, aes(x = variable, y = abs_loading, fill = sign)) +
          geom_col(width = 0.7) +
          coord_flip() +
          scale_fill_manual(values = c("Positive" = "#E69F00", "Negative" = "#0072B2"), name = "Sign", limits = c("Positive", "Negative"), drop = FALSE) +
          theme_classic(base_size = 14) +
          labs(title = "Top 10 loadings: PC1", x = "Feature", y = "Absolute loading") +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            legend.position = "bottom"
          )
        gp_load2 <- ggplot(top_pc2_plot, aes(x = variable, y = abs_loading, fill = sign)) +
          geom_col(width = 0.7) +
          coord_flip() +
          scale_fill_manual(values = c("Positive" = "#E69F00", "Negative" = "#0072B2"), name = "Sign", limits = c("Positive", "Negative"), drop = FALSE) +
          theme_classic(base_size = 14) +
          labs(title = "Top 10 loadings: PC2", x = "Feature", y = "Absolute loading") +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            legend.position = "bottom"
          )
        loadings_plots <- list(PC1 = gp_load1, PC2 = gp_load2)
        multi_pca_top_loadings_plots_local <- loadings_plots
        assign("multi_pca_top_loadings_plots", loadings_plots, envir = .GlobalEnv)
      }
    }
  }

  if (render_pdf) {
    pdf(pdf_out, paper = "a4r", width = 14, height = 11)
    grid.newpage()
    title_text <- if (params$mode == "isoforms") "SQANTI-single cell\nmulti-sample isoforms report" else "SQANTI-single cell\nmulti-sample reads report"
    cover <- textGrob(title_text,
      gp = gpar(fontface = "italic", fontsize = 40, col = "orangered")
    )
    grid.draw(cover)

    theme_pdf_paper <- theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 20, face = "bold")
    )

    tbl_theme <- ttheme_default(
      core = list(fg_params = list(cex = 1.4, hjust = 0.5, x = 0.5)),
      colhead = list(fg_params = list(cex = 1.4, fontface = "bold", hjust = 0.5, x = 0.5))
    )
    tbl_grob <- tableGrob(summary_tbl, rows = NULL, theme = tbl_theme)
    title_grob <- textGrob("Per cell summary of samples", gp = gpar(fontface = "italic", fontsize = 28))
    grid.newpage()
    pushViewport(viewport(x = 0.5, y = 0.95))
    grid.draw(title_grob)
    popViewport()
    pushViewport(viewport(x = 0.5, y = 0.5))
    grid.draw(tbl_grob)
    popViewport()

    # Straight after the summary table, same as the HTML: the overview is what
    # tells the reader which of the distributions below are worth opening.
    # No theme_pdf_paper here: its larger base sizes blow the tile text out of
    # the tiles. The PDF page is a fixed A4 landscape for the whole device, so
    # a wide cohort is absorbed by the font scaling inside the builder rather
    # than by a wider canvas.
    if (!is.null(multi_qc_overview_plot_local)) {
      print(multi_qc_overview_plot_local)
    }

    if (length(multi_library_size_plots) > 0) {
      for (plt in multi_library_size_plots) print(plt + theme_pdf_paper)
    }
    if (length(multi_gene_char_plots) > 0) {
      for (plt in multi_gene_char_plots) print(plt + theme_pdf_paper)
    }
    if (length(multi_length_plots_local) > 0) {
      for (plt in multi_length_plots_local) print(plt + theme_pdf_paper)
    }

    if (have_cats) {
      for (gp in category_plots) {
        print(gp + theme_pdf_paper)
      }
    }

    for (plts in list(multi_sj_composition_plots_local, multi_good_quality_plots_local,
                      multi_bad_quality_plots_local, multi_extra_feature_plots_local)) {
      for (gp in plts) print(gp + theme_pdf_paper)
    }

    if (!is.null(multi_pca_scores_plot_local)) {
      print(multi_pca_scores_plot_local + theme_pdf_paper)
    }
    if (!is.null(multi_pca_group_plot_local)) {
      print(multi_pca_group_plot_local + theme_pdf_paper)
    }
    if (!is.null(multi_pca_scree_plot_local)) {
      print(multi_pca_scree_plot_local + theme_pdf_paper)
    }
    if (!is.null(multi_pca_top_loadings_plots_local)) {
      if (!is.null(multi_pca_top_loadings_plots_local[["PC1"]])) {
        print(multi_pca_top_loadings_plots_local[["PC1"]] + theme_pdf_paper)
      }
      if (!is.null(multi_pca_top_loadings_plots_local[["PC2"]])) {
        print(multi_pca_top_loadings_plots_local[["PC2"]] + theme_pdf_paper)
      }
    }

    dev.off()
    message(sprintf("**** Multisample report written: %s", pdf_out))
  }

  if (render_html) {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    script_arg <- cmd_args[grep("--file=", cmd_args)]
    if (length(script_arg) > 0) {
      script_path <- substring(script_arg, 8L)
      script_dir <- dirname(normalizePath(script_path))
    } else {
      script_dir <- getwd()
    }

    rmd_file <- file.path(script_dir, "SQANTI-sc_multisample_report.Rmd")
    css_file <- file.path(script_dir, "style-multisample.css")
    html_output_file <- file.path(out_dir, paste0(params$prefix, ".html"))

    if (!file.exists(rmd_file)) {
      stop("HTML report template not found: ", rmd_file)
    }

    if (file.exists(css_file)) {
      file.copy(css_file, dirname(html_output_file), overwrite = TRUE)
    }

    rmarkdown::render(
      rmd_file,
      output_file = html_output_file,
      intermediates_dir = dirname(html_output_file),
      knit_root_dir = dirname(html_output_file),
      envir = globalenv(),
      quiet = TRUE
    )

    # Cleanup: remove the copied CSS file
    css_output <- file.path(dirname(html_output_file), basename(css_file))
    if (file.exists(css_output)) {
      file.remove(css_output)
    }

    message("HTML report generated: ", html_output_file)
  }
}

main()
