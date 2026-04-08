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
  library(ggdist)
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
  res <- list(files = NULL, class_files = NULL, out_dir = ".", mode = "reads", report = "pdf", prefix = "SQANTI_sc_multi_report")
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (startsWith(key, "--")) {
      k <- substring(key, 3)
      if (k %in% c("files", "class_files", "out_dir", "mode", "report", "prefix")) {
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

# Helper: tidy feature names for titles and subtitles
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
  domain <- "Value"
  if (stringr::str_detect(name_lower, "length")) domain <- "Length"
  if (stringr::str_detect(name_lower, "length") && stringr::str_detect(name_lower, "prop|perc|pct|ratio|fraction")) {
    if (exists("params") && params$mode == "isoforms") domain <- "Transcripts" else domain <- "Reads"
  }
  if (stringr::str_detect(name_lower, "read")) {
    if (exists("params") && params$mode == "isoforms") domain <- "Transcripts" else domain <- "Reads"
  }
  if (stringr::str_detect(name_lower, "\\bmt\\b") || stringr::str_detect(name_lower, "^mt_")) domain <- "Reads"
  if (stringr::str_detect(name_lower, "exon")) domain <- "Exons"
  if (stringr::str_detect(name_lower, "intron")) domain <- "Introns"
  if (stringr::str_detect(name_lower, "coverage")) domain <- "Coverage"
  if (stringr::str_detect(name_lower, "isoform")) domain <- "Isoforms"
  if (stringr::str_detect(name_lower, "transcript")) domain <- "Transcripts"
  if (stringr::str_detect(name_lower, "gene")) domain <- "Genes"
  if (stringr::str_detect(name_lower, "umi")) domain <- "UMIs"
  if (stringr::str_detect(name_lower, "junction")) domain <- "Junctions"
  structural_keywords <- c(
    "fsm", "ism", "nic", "nnc",
    "genic_genomic", "genic genomic", "genic",
    "antisense", "fusion", "intergenic",
    "genic_intron", "genic intron"
  )
  if (any(stringr::str_detect(name_lower, structural_keywords))) {
    if (exists("params") && params$mode == "isoforms") domain <- "Transcripts" else domain <- "Reads"
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
    if (unit == "count" && !is_prop_keyword && maxv <= 1.5 && minv >= 0) {
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

# Helper: build per-feature violin + boxplot for a PCA loading
build_loading_feature_plot <- function(multi, feature_info, sample_levels) {
  feature_name <- as.character(feature_info$variable)[1]
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
  if (info$scale_to_percent) {
    plot_df <- plot_df %>% mutate(value = value * 100)
  }
  if (length(sample_levels) == 0) {
    sample_levels <- unique(plot_df$sampleID)
  }
  plot_df <- plot_df %>% mutate(sampleID = factor(sampleID, levels = sample_levels))
  uniqueness <- plot_df %>%
    group_by(sampleID) %>%
    summarise(unique_vals = n_distinct(value), .groups = "drop")
  use_violin <- any(uniqueness$unique_vals > 1)
  loading_value <- as.numeric(feature_info$loading)[1]
  loading_rank <- as.integer(feature_info$rank)[1]
  pc_label <- as.character(feature_info$PC)[1]
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
      title = sprintf("Per Sample %s Distribution Across Cells", info$display_name),
      subtitle = sprintf("%s loading rank #%d (loading = %.3f)", pc_label, loading_rank, loading_value),
      x = "Sample",
      y = info$value_label
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      axis.title = element_text(size = 15),
      axis.text.x = element_text(size = 13, angle = 35, hjust = 1),
      axis.text.y = element_text(size = 13)
    )
  return(gp)
}

# Reusable helper: violin + boxplot comparing a single metric across samples
build_sample_comparison_plot <- function(data, col_name, title, y_label,
                                         sample_levels, scale_percent = FALSE) {
  if (!col_name %in% colnames(data)) return(NULL)
  plot_df <- data %>%
    select(sampleID, value = all_of(col_name)) %>%
    mutate(value = as.numeric(value)) %>%
    filter(is.finite(value))
  if (nrow(plot_df) == 0) return(NULL)
  if (scale_percent) plot_df$value <- plot_df$value * 100
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
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 15),
      axis.text.x = element_text(size = 13, angle = 35, hjust = 1),
      axis.text.y = element_text(size = 13)
    )
  gp
}

main <- function() {
  params <- parse_args()

  files <- unlist(strsplit(params$files, ",", fixed = TRUE))
  files <- trimws(files)
  files <- files[nchar(files) > 0]
  if (length(files) < 2) {
    message("[INFO] Fewer than 2 files provided. Nothing to do.")
    quit(status = 0)
  }

  # Read all summaries
  lst <- lapply(files, safe_read_summary)
  lst <- Filter(Negate(is.null), lst)
  if (length(lst) < 2) {
    message("[INFO] Fewer than 2 valid summaries after reading. Skipping.")
    quit(status = 0)
  }

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

  sample_levels_global <- unique(multi$sampleID[!is.na(multi$sampleID)])

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
      median_novel = median(Novel_genes, na.rm = TRUE),
      median_ujc = if ("UJCs_in_cell" %in% names(.)) median(UJCs_in_cell, na.rm = TRUE) else NA,
      median_mt = median(MT_perc, na.rm = TRUE)
    )

  entity_label_plural <- if (params$mode == "isoforms") "Transcripts" else "Reads"

  summary_tbl <- per_sample_stats %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    transmute(
      Sample = sampleID,
      `Cell\nBarcodes` = cells,
      `Median\nReads` = median_reads,
      `Median\nUMIs` = median_umis,
      `Median\nAnnotated\nGenes` = median_annotated,
      `Median\nNovel\nGenes` = median_novel,
      `Median\nUJCs` = median_ujc,
      `Median\nMitochondrial\nReads` = median_mt
    )

  # Rename columns dynamically
  colnames(summary_tbl)[colnames(summary_tbl) == "Median\nReads"] <- paste0("Median\n", entity_label_plural)
  colnames(summary_tbl)[colnames(summary_tbl) == "Median\nMitochondrial\nReads"] <- paste0("Median\nMitochondrial\n", entity_label_plural)

  if (params$mode == "isoforms") {
    summary_tbl <- summary_tbl %>% select(-`Median\nUMIs`, -`Median\nUJCs`)
  }

  summary_tbl_html <- summary_tbl
  colnames(summary_tbl_html) <- gsub("\\n", "<br>", colnames(summary_tbl_html), fixed = TRUE)

  assign("multi_per_sample_stats", per_sample_stats, envir = .GlobalEnv)
  assign("multi_summary_tbl_pdf", summary_tbl, envir = .GlobalEnv)
  assign("multi_summary_tbl_html", summary_tbl_html, envir = .GlobalEnv)
  assign("entity_label", if (params$mode == "isoforms") "Transcript" else "Read", envir = .GlobalEnv)
  assign("entity_label_plural", entity_label_plural, envir = .GlobalEnv)
  assign("mode", params$mode, envir = .GlobalEnv)

  # -------- Per-Cell Library Size plots --------
  library_size_specs <- list()
  library_size_specs[[paste0(entity_label_plural, " per Cell")]] <- list(
    col = count_col, y = paste0(entity_label_plural, ", count")
  )
  if (params$mode == "reads" && "UMIs_in_cell" %in% colnames(multi)) {
    library_size_specs[["UMIs per Cell"]] <- list(col = "UMIs_in_cell", y = "UMIs, count")
  }
  if ("Genes_in_cell" %in% colnames(multi)) {
    library_size_specs[["Genes per Cell"]] <- list(col = "Genes_in_cell", y = "Genes, count")
  }
  if ("MT_perc" %in% colnames(multi)) {
    library_size_specs[["Mitochondrial % per Cell"]] <- list(
      col = "MT_perc", y = paste0(entity_label_plural, ", %")
    )
  }

  multi_library_size_plots <- list()
  for (nm in names(library_size_specs)) {
    spec <- library_size_specs[[nm]]
    gp <- build_sample_comparison_plot(
      multi, spec$col,
      title = paste0("Per Sample ", nm, " Distribution"),
      y_label = spec$y,
      sample_levels = sample_levels_global,
      scale_percent = isTRUE(spec$scale_pct)
    )
    if (!is.null(gp)) multi_library_size_plots[[nm]] <- gp
  }
  if (length(multi_library_size_plots) > 0) {
    assign("multi_library_size_plots", multi_library_size_plots, envir = .GlobalEnv)
  }

  # -------- Gene Characterization plots --------
  gene_char_specs <- list()
  if ("Annotated_genes" %in% colnames(multi)) {
    gene_char_specs[["Annotated Genes per Cell"]] <- list(col = "Annotated_genes", y = "Genes, count")
  }
  if (params$mode == "reads" && "UJCs_in_cell" %in% colnames(multi)) {
    gene_char_specs[["UJCs per Cell"]] <- list(col = "UJCs_in_cell", y = "UJCs, count")
  }

  multi_gene_char_plots <- list()
  for (nm in names(gene_char_specs)) {
    spec <- gene_char_specs[[nm]]
    gp <- build_sample_comparison_plot(
      multi, spec$col,
      title = paste0("Per Sample ", nm, " Distribution"),
      y_label = spec$y,
      sample_levels = sample_levels_global
    )
    if (!is.null(gp)) multi_gene_char_plots[[nm]] <- gp
  }
  if (length(multi_gene_char_plots) > 0) {
    assign("multi_gene_char_plots", multi_gene_char_plots, envir = .GlobalEnv)
  }

  # -------- Bulk-level Length Distribution --------
  multi_length_distribution_plot_local <- NULL
  if (!is.null(params$class_files) && nzchar(params$class_files)) {
    class_file_paths <- trimws(unlist(strsplit(params$class_files, ",", fixed = TRUE)))
    class_file_paths <- class_file_paths[nchar(class_file_paths) > 0 & file.exists(class_file_paths)]

    if (length(class_file_paths) >= 2) {
      select_cols <- if (params$mode == "isoforms") c("length", "FL") else c("length")
      len_dfs <- lapply(class_file_paths, function(f) {
        df <- tryCatch(
          data.table::fread(f, select = select_cols, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, data.table = FALSE),
          error = function(e) { message("[WARNING] Could not read ", f, ": ", e$message); NULL }
        )
        if (is.null(df) || nrow(df) == 0) return(NULL)
        df$length <- suppressWarnings(as.numeric(df$length))
        df <- df[is.finite(df$length) & df$length > 0, , drop = FALSE]
        sample_id <- sub("_classification\\.txt(\\.gz)?$", "", basename(f))
        df$sampleID <- sample_id
        df
      })
      len_combined <- bind_rows(Filter(Negate(is.null), len_dfs))

      if (nrow(len_combined) > 0) {
        if (params$mode == "isoforms" && "FL" %in% colnames(len_combined)) {
          len_combined$FL <- suppressWarnings(as.integer(len_combined$FL))
          len_combined$FL[is.na(len_combined$FL) | len_combined$FL < 1] <- 1L
          len_combined <- len_combined[rep(seq_len(nrow(len_combined)), len_combined$FL), c("sampleID", "length")]
        }
        len_combined$sampleID <- factor(len_combined$sampleID, levels = sample_levels_global)

        uniqueness <- len_combined %>%
          group_by(sampleID) %>%
          summarise(uv = n_distinct(length), .groups = "drop")
        use_violin <- any(uniqueness$uv > 1)

        gp_len <- ggplot(len_combined, aes(x = sampleID, y = length, fill = sampleID, colour = sampleID))
        if (use_violin) {
          gp_len <- gp_len + geom_violin(trim = TRUE, scale = "width", alpha = 0.7, linewidth = 0.3)
        }
        gp_len <- gp_len +
          geom_boxplot(width = 0.05, outlier.shape = NA, alpha = 0.3, colour = "grey20") +
          stat_summary(fun = mean, geom = "point", shape = 4, size = 1,
                       colour = "red", stroke = 0.45) +
          scale_fill_conesa(palette = "complete", drop = FALSE) +
          scale_color_conesa(palette = "complete", guide = "none", drop = FALSE) +
          labs(
            title = paste0("Per Sample ", entity_label_plural, " Length Distribution"),
            x = "Sample", y = "Length, bp"
          ) +
          theme_classic(base_size = 14) +
          theme(
            legend.position = "none",
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 15),
            axis.text.x = element_text(size = 13, angle = 35, hjust = 1),
            axis.text.y = element_text(size = 13)
          )
        multi_length_distribution_plot_local <- gp_len
        assign("multi_length_distribution_plot", gp_len, envir = .GlobalEnv)
      }
    }
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
        sampleID = factor(sampleID, levels = unique(multi$sampleID)),
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

    p_cats <- ggplot(cats_long, aes(x = category, y = prop, fill = sampleID, colour = sampleID, group = sampleID)) +
      ggdist::stat_slabinterval(
        side = "left",
        position = position_dodge(width = 0.6),
        density = "unbounded",
        bw = "nrd0",
        normalize = "groups", scale = 0.85, adjust = 2.5, trim = TRUE,
        show_point = FALSE, show_interval = FALSE, # slab only (fill only)
        slab_colour = NA, alpha = 0.7
      ) +
      ggdist::stat_slabinterval(
        side = "left",
        position = position_dodge(width = 0.6),
        density = "unbounded",
        bw = "nrd0",
        normalize = "groups", scale = 0.85, adjust = 2.5, trim = TRUE,
        show_point = FALSE, show_interval = FALSE, # outline-only
        mapping = aes(slab_colour = after_scale(colour)),
        fill = NA, slab_linewidth = 0.05, show.legend = FALSE
      ) +
      # Draw medians as short horizontal lines centered within each dodged slab
      stat_summary(
        data = cats_long, aes(group = sampleID),
        fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", width = 0.1,
        position = position_dodge(width = 0.6),
        color = "black", linewidth = 0.1, alpha = 1
      ) +
      # Add a small lower expansion so slabs don't touch the x-axis
      scale_y_continuous(limits = c(0, 100), expand = expansion(add = c(1, 0))) +
      scale_fill_conesa(palette = "complete") +
      scale_color_conesa(palette = "complete", guide = "none") +
      guides(
        fill = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 0.95, colour = NA, stroke = 0)),
        linetype = "none", alpha = "none", size = "none", colour = "none"
      ) +
      theme_classic(base_size = 13) +
      labs(
        title = "Structural Category Proportions by Sample",
        x = "Structural category", y = "Reads, %"
      ) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      )

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
      ggplot(dfp, aes(x = sampleID, y = prop)) +
        geom_violin(fill = violin_fill, color = cat_col, linewidth = 0.3, width = 0.8, trim = TRUE) +
        geom_boxplot(width = 0.05, outlier.shape = NA, fill = cat_col, color = box_outline_col, alpha = 0.3) +
        stat_summary(fun = mean, geom = "point", shape = 4, size = 1, colour = "red", stroke = 0.9) +
        scale_y_continuous(limits = c(0, 100), expand = expansion(add = c(1, 0))) +
        theme_classic(base_size = 13) +
        labs(title = paste0("Per Sample ", cat_lab, " Reads Distribution Across Cells"), x = "Sample", y = "Reads, %") +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1, size = 13),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 13)
        )
    })
    names(category_plots) <- as.character(category_levels)

    assign("multi_structural_category_combined_plot", p_cats, envir = .GlobalEnv)
    assign("multi_structural_category_violin_plots", category_plots, envir = .GlobalEnv)

  }

  multi_pca_scores_plot_local <- NULL

  multi_pca_scree_plot_local <- NULL

  multi_pca_top_loadings_plots_local <- NULL

  multi_pca_loading_distribution_plots_local <- list()

  # -------- PCA (all numeric features, per-sample medians) --------
  # 1) Select all numeric columns from the cell summary
  num_cols <- names(multi)[sapply(multi, function(x) is.numeric(x) && !all(is.na(x)))]
  # 2) Aggregate per-sample medians across all numeric features
  agg_median <- multi %>%
    group_by(sampleID) %>%
    summarise(across(all_of(num_cols), ~ median(., na.rm = TRUE)), .groups = "drop")

  if (nrow(agg_median) >= 2 && ncol(agg_median) >= 2) {
    # 3) Drop features with zero variance across samples
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
          geom_point(size = 3.8, alpha = 0.95, shape = 19, stroke = 0) +
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
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 15),
            axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 13)
          ) +
          guides(colour = guide_legend(override.aes = list(size = 5, alpha = 0.95, stroke = 0)))
        multi_pca_scores_plot_local <- gp_scores
        assign("multi_pca_scores_plot", gp_scores, envir = .GlobalEnv)

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
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 15),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),
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
          scale_fill_manual(values = c("Positive" = "#78C679", "Negative" = "#EE6A50"), name = "Sign", limits = c("Positive", "Negative"), drop = FALSE) +
          theme_classic(base_size = 14) +
          labs(title = "Top 10 loadings: PC1", x = "Feature", y = "Absolute loading") +
          theme(
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 15),
            axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 13),
            legend.position = "bottom"
          )
        gp_load2 <- ggplot(top_pc2_plot, aes(x = variable, y = abs_loading, fill = sign)) +
          geom_col(width = 0.7) +
          coord_flip() +
          scale_fill_manual(values = c("Positive" = "#78C679", "Negative" = "#EE6A50"), name = "Sign", limits = c("Positive", "Negative"), drop = FALSE) +
          theme_classic(base_size = 14) +
          labs(title = "Top 10 loadings: PC2", x = "Feature", y = "Absolute loading") +
          theme(
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 15),
            axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 13),
            legend.position = "bottom"
          )
        loadings_plots <- list(PC1 = gp_load1, PC2 = gp_load2)
        multi_pca_top_loadings_plots_local <- loadings_plots
        assign("multi_pca_top_loadings_plots", loadings_plots, envir = .GlobalEnv)

        # D) Distribution plots for top-loading features on PC1/PC2
        sample_levels <- sample_levels_global
        loading_plot_info <- bind_rows(top_pc1, top_pc2) %>%
          distinct(variable, .keep_all = TRUE)
        loading_distribution_plots <- list()
        if (nrow(loading_plot_info) > 0) {
          for (idx in seq_len(nrow(loading_plot_info))) {
            gp_loading <- build_loading_feature_plot(multi, loading_plot_info[idx, ], sample_levels)
            feat_name <- loading_plot_info$variable[idx]
            if (is.null(gp_loading)) {
              message(sprintf("[INFO] Skipping PCA loading feature %s due to missing or constant data.", feat_name))
            } else {
              loading_distribution_plots[[feat_name]] <- gp_loading
            }
          }
        }
        multi_pca_loading_distribution_plots_local <- loading_distribution_plots
        assign("multi_pca_loading_distribution_plots", loading_distribution_plots, envir = .GlobalEnv)

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

    if (length(multi_library_size_plots) > 0) {
      for (plt in multi_library_size_plots) print(plt)
    }
    if (length(multi_gene_char_plots) > 0) {
      for (plt in multi_gene_char_plots) print(plt)
    }
    if (!is.null(multi_length_distribution_plot_local)) {
      print(multi_length_distribution_plot_local)
    }

    if (have_cats) {
      for (gp in category_plots) {
        print(gp)
      }
      print(p_cats)
    }

    if (!is.null(multi_pca_scores_plot_local)) {
      print(multi_pca_scores_plot_local)
    }
    if (!is.null(multi_pca_scree_plot_local)) {
      print(multi_pca_scree_plot_local)
    }
    if (!is.null(multi_pca_top_loadings_plots_local)) {
      gp_load1 <- multi_pca_top_loadings_plots_local[["PC1"]]
      gp_load2 <- multi_pca_top_loadings_plots_local[["PC2"]]
      if (!is.null(gp_load1) && !is.null(gp_load2)) {
        legend_df <- data.frame(
          variable = c("pos", "neg"),
          abs_loading = c(1, 1),
          sign = factor(c("Positive", "Negative"), levels = c("Positive", "Negative"))
        )
        legend_plot <- ggplot(legend_df, aes(x = variable, y = abs_loading, fill = sign)) +
          geom_col() +
          scale_fill_manual(values = c("Positive" = "#78C679", "Negative" = "#EE6A50"), name = "Sign", limits = c("Positive", "Negative"), drop = FALSE) +
          theme_void(base_size = 14) +
          theme(legend.position = "bottom")
        legend_grob <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")
        row_plots <- arrangeGrob(gp_load1 + theme(legend.position = "none"), gp_load2 + theme(legend.position = "none"), ncol = 2)
        grid.arrange(row_plots, legend_grob, ncol = 1, heights = c(0.86, 0.14))
      } else {
        for (plt in multi_pca_top_loadings_plots_local) {
          print(plt)
        }
      }
    }
    if (length(multi_pca_loading_distribution_plots_local) > 0) {
      for (nm in names(multi_pca_loading_distribution_plots_local)) {
        print(multi_pca_loading_distribution_plots_local[[nm]])
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
