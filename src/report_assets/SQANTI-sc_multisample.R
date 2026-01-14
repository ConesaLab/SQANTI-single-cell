#!/usr/env/bin Rscript

############################################################
##### SQANTI single-cell multisample report generation #####
############################################################



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
  library(plotly)
  library(stringr)
  library(rmarkdown)
})

# Ensure RColorConesa is available; if not, install from GitHub, then load
if (!requireNamespace("RColorConesa", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("ConesaLab/RColorConesa")
}
library(RColorConesa)


parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  # Simple flag parser: expects --key value pairs, and a single string for --files (comma-separated)
  res <- list(files = NULL, out_dir = ".", mode = "reads", report = "pdf", prefix = "SQANTI_sc_multi_report")
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (startsWith(key, "--")) {
      k <- substring(key, 3)
      if (k %in% c("files", "out_dir", "mode", "report", "prefix")) {
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
      read.table(fpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
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
  suffixes <- c(
    "_prop_in_cell$", "_perc_in_cell$", "_prop$", "_perc$",
    "_percentage$", "_pct$", "_ratio$", "_count$", "_counts$",
    "_in_cell$", "_per_cell$", "_value$"
  )
  for (pattern in suffixes) {
    cleaned <- gsub(pattern, "", cleaned, ignore.case = TRUE)
  }
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

to_rgba <- function(col, alpha = 1) {
  if (is.null(col) || is.na(col) || !nzchar(col)) {
    return(sprintf("rgba(0,0,0,%.3f)", alpha))
  }
  rgb <- grDevices::col2rgb(col)
  sprintf("rgba(%d,%d,%d,%.3f)", rgb[1], rgb[2], rgb[3], alpha)
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
      "#6BAED6", "#FC8D59", "#78C679", "#EE6A50", "#969696",
      "#66C2A4", "#FFD92F", "#E78AC3", "#A6D854", "#8DA0CB",
      "#E5C494", "#B3B3B3"
    )
    col_vec <- fallback
  }

  rep_len(col_vec, n)
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
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text.x = element_text(size = 14, angle = 35, hjust = 1),
      axis.text.y = element_text(size = 14)
    )
  plot_df_html <- plot_df
  plot_df_html$sampleID <- as.character(plot_df_html$sampleID)
  hover_tmpl <- sprintf("Sample: %%{x}<br>%s: %%{y:.3f}<extra></extra>", info$value_label)

  sample_levels_html <- if (!is.null(sample_levels) && length(sample_levels) > 0) sample_levels else unique(plot_df_html$sampleID)
  palette_cols <- get_conesa_palette_colors(length(sample_levels_html), palette = "complete")
  sample_color_map <- setNames(palette_cols, sample_levels_html)

  plt_html <- plotly::plot_ly()
  for (idx in seq_along(sample_levels_html)) {
    sample_nm <- sample_levels_html[[idx]]
    sample_df <- plot_df_html %>% filter(sampleID == sample_nm)
    if (nrow(sample_df) == 0) next
    col_val <- sample_color_map[[sample_nm]]
    violin_fill <- to_rgba(col_val, 0.7)
    box_fill <- to_rgba(col_val, 0.3)
    line_col <- to_rgba(col_val, 1)
    plt_html <- plt_html %>%
      plotly::add_trace(
        data = sample_df,
        x = ~sampleID,
        y = ~value,
        type = "violin",
        name = sample_nm,
        legendgroup = sample_nm,
        showlegend = TRUE,
        hovertemplate = hover_tmpl,
        fillcolor = violin_fill,
        line = list(color = line_col, width = 1.1),
        spanmode = "hard",
        scalemode = "width",
        width = 0.85,
        points = "none",
        box = list(visible = FALSE),
        meanline = list(visible = FALSE)
      ) %>%
      plotly::add_trace(
        data = sample_df,
        x = ~sampleID,
        y = ~value,
        type = "box",
        name = paste0(sample_nm, " (IQR)"),
        legendgroup = sample_nm,
        showlegend = FALSE,
        hoverinfo = "skip",
        fillcolor = box_fill,
        line = list(color = to_rgba("#333333", 1), width = 1),
        boxpoints = FALSE,
        width = 0.05
      )
  }

  mean_df <- plot_df_html %>%
    group_by(sampleID) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
  if (nrow(mean_df) > 0) {
    plt_html <- plt_html %>% plotly::add_trace(
      data = mean_df,
      x = ~sampleID,
      y = ~mean_value,
      type = "scatter",
      mode = "markers",
      name = "Mean",
      legendgroup = "Mean",
      hovertemplate = hover_tmpl,
      marker = list(symbol = "x-thin", size = 8, color = "red", line = list(width = 0)),
      showlegend = FALSE
    )
  }

  plt_html <- plt_html %>% plotly::layout(
    title = list(
      text = sprintf("<b>Per Sample %s Distribution Across Cells</b>", info$display_name),
      x = 0.5,
      xanchor = "center",
      font = list(size = 18, family = "Arial")
    ),
    xaxis = list(
      title = list(text = "Sample", font = list(size = 16, family = "Arial")),
      tickfont = list(size = 14, family = "Arial"),
      tickangle = 45,
      showline = TRUE,
      linecolor = "#000000",
      linewidth = 1.1,
      mirror = FALSE,
      zeroline = FALSE,
      standoff = 26,
      automargin = TRUE
    ),
    yaxis = list(
      title = list(text = info$value_label, font = list(size = 16, family = "Arial")),
      tickfont = list(size = 14, family = "Arial"),
      showline = TRUE,
      linecolor = "#000000",
      linewidth = 1.1,
      zeroline = FALSE,
      standoff = 8,
      automargin = TRUE
    ),
    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = -0.25,
      yanchor = "top",
      font = list(size = 14, family = "Arial"),
      title = list(text = "")
    ),
    margin = list(t = 60, b = 250, l = 130, r = 80),
    hovermode = "closest",
    paper_bgcolor = "rgba(0,0,0,0)",
    plot_bgcolor = "rgba(0,0,0,0)",
    font = list(family = "Arial", size = 14),
    height = 700
  )

  list(ggplot = gp, plotly = plt_html)
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
      mean_reads = mean(.data[[count_col]], na.rm = TRUE),
      mean_umis = if ("UMIs_in_cell" %in% names(.)) mean(UMIs_in_cell, na.rm = TRUE) else NA,
      median_reads = median(.data[[count_col]], na.rm = TRUE),
      mean_genes = mean(Genes_in_cell, na.rm = TRUE),
      median_genes = median(Genes_in_cell, na.rm = TRUE),
      mean_annotated = mean(Annotated_genes, na.rm = TRUE),
      mean_novel = mean(Novel_genes, na.rm = TRUE),
      mean_ujc = if ("UJCs_in_cell" %in% names(.)) mean(UJCs_in_cell, na.rm = TRUE) else NA,
      median_ujc = if ("UJCs_in_cell" %in% names(.)) median(UJCs_in_cell, na.rm = TRUE) else NA,
      mean_mt = mean(MT_perc, na.rm = TRUE)
    )

  entity_label_plural <- if (params$mode == "isoforms") "Transcripts" else "Reads"

  summary_tbl <- per_sample_stats %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    transmute(
      Sample = sampleID,
      `Cell\nBarcodes` = cells,
      `Average\nReads` = mean_reads,
      `Average\nUMIs` = mean_umis,
      `Average\nAnnotated\nGenes` = mean_annotated,
      `Average\nNovel\nGenes` = mean_novel,
      `Average\nUJCs` = mean_ujc,
      `Average\nMitochondrial\nReads` = mean_mt
    )

  # Rename columns dynamically
  colnames(summary_tbl)[colnames(summary_tbl) == "Average\nReads"] <- paste0("Average\n", entity_label_plural)
  colnames(summary_tbl)[colnames(summary_tbl) == "Average\nMitochondrial\nReads"] <- paste0("Average\nMitochondrial\n", entity_label_plural)

  if (params$mode == "isoforms") {
    summary_tbl <- summary_tbl %>% select(-`Average\nUMIs`, -`Average\nUJCs`)
  }

  summary_tbl_html <- summary_tbl
  colnames(summary_tbl_html) <- gsub("\\n", "<br>", colnames(summary_tbl_html), fixed = TRUE)

  assign("multi_per_sample_stats", per_sample_stats, envir = .GlobalEnv)
  assign("multi_summary_tbl_pdf", summary_tbl, envir = .GlobalEnv)
  assign("multi_summary_tbl_html", summary_tbl_html, envir = .GlobalEnv)
  assign("entity_label", if (params$mode == "isoforms") "Transcript" else "Read", envir = .GlobalEnv)
  assign("entity_label_plural", entity_label_plural, envir = .GlobalEnv)
  assign("mode", params$mode, envir = .GlobalEnv)

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
            "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic intron"
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
      theme_classic(base_size = 12) +
      labs(
        title = "Structural Category Proportions by Sample",
        x = "Structural category", y = "Reads, %"
      ) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
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
      "Genic intron" = "#41B6C4"
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
        theme_classic(base_size = 14) +
        labs(title = paste0("Per Sample ", cat_lab, " Reads Distribution Across Cells"), x = "Sample", y = "Reads, %") +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14)
        )
    })
    names(category_plots) <- as.character(category_levels)

    category_plots_html <- lapply(category_levels, function(cat_lab) {
      dfp <- cats_long %>%
        filter(category == cat_lab) %>%
        mutate(sampleID = factor(sampleID, levels = sample_levels_all)) %>%
        filter(is.finite(prop))

      cat_col <- unname(cat_to_col[[as.character(cat_lab)]])
      if (is.null(cat_col) || is.na(cat_col)) cat_col <- "#6C757D"
      violin_fill <- to_rgba(cat_col, 0.7)
      box_fill <- to_rgba(cat_col, 0.3)
      line_col <- to_rgba(cat_col, 1)

      plt <- plotly::plot_ly()

      for (sample_nm in levels(dfp$sampleID)) {
        sample_df <- dfp %>% filter(sampleID == sample_nm)
        if (nrow(sample_df) == 0) next

        plt <- plt %>%
          plotly::add_trace(
            data = sample_df,
            x = ~sampleID,
            y = ~prop,
            type = "violin",
            name = NULL,
            showlegend = FALSE,
            hovertemplate = "Sample: %{x}<br>Reads, %: %{y:.3f}<extra></extra>",
            fillcolor = violin_fill,
            line = list(color = line_col, width = 1.1),
            spanmode = "hard",
            scalemode = "width",
            width = 0.85,
            points = "none",
            box = list(visible = FALSE),
            meanline = list(visible = FALSE)
          ) %>%
          plotly::add_trace(
            data = sample_df,
            x = ~sampleID,
            y = ~prop,
            type = "box",
            name = NULL,
            showlegend = FALSE,
            hoverinfo = "skip",
            fillcolor = box_fill,
            line = list(color = to_rgba("#333333", 1), width = 1),
            boxpoints = FALSE,
            width = 0.05
          )
      }

      mean_df <- dfp %>%
        group_by(sampleID) %>%
        summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")
      if (nrow(mean_df) > 0) {
        plt <- plt %>% plotly::add_trace(
          data = mean_df,
          x = ~sampleID,
          y = ~mean_prop,
          type = "scatter",
          mode = "markers",
          name = NULL,
          showlegend = FALSE,
          hovertemplate = "Sample: %{x}<br>Reads, %: %{y:.3f}<extra></extra>",
          marker = list(symbol = "x-thin", size = 8, line = list(width = 0), color = "red")
        )
      }

      plt %>%
        plotly::layout(
          showlegend = FALSE,
          title = list(
            text = sprintf("<b>Per Sample %s Reads Distribution Across Cells</b>", cat_lab),
            x = 0.5,
            xanchor = "center",
            font = list(size = 18, family = "Arial")
          ),
          xaxis = list(
            title = list(text = "Sample", font = list(size = 16, family = "Arial")),
            tickfont = list(size = 14, family = "Arial"),
            tickangle = 45,
            showline = TRUE,
            linecolor = "#000000",
            linewidth = 1.1,
            mirror = FALSE,
            zeroline = FALSE
          ),
          yaxis = list(
            title = list(text = "Reads, %", font = list(size = 16, family = "Arial")),
            tickfont = list(size = 14, family = "Arial"),
            range = c(0, 100),
            showline = TRUE,
            linecolor = "#000000",
            linewidth = 1.1,
            zeroline = FALSE
          ),
          margin = list(t = 60, b = 160, l = 110, r = 60),
          hovermode = "closest",
          paper_bgcolor = "rgba(0,0,0,0)",
          plot_bgcolor = "rgba(0,0,0,0)",
          font = list(family = "Arial", size = 14),
          height = 560
        )
    })
    names(category_plots_html) <- as.character(category_levels)

    assign("multi_structural_category_combined_plot", p_cats, envir = .GlobalEnv)
    assign("multi_structural_category_violin_plots", category_plots, envir = .GlobalEnv)
    assign("multi_structural_category_violin_plots_html", category_plots_html, envir = .GlobalEnv)
  }

  multi_pca_scores_plot_local <- NULL
  multi_pca_scores_plot_html_local <- NULL
  multi_pca_scree_plot_local <- NULL
  multi_pca_scree_plot_html_local <- NULL
  multi_pca_top_loadings_plots_local <- NULL
  multi_pca_top_loadings_plots_html_local <- NULL
  multi_pca_loading_distribution_plots_local <- list()
  multi_pca_loading_distribution_plots_html_local <- list()
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
          theme_classic(base_size = 16) +
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
            legend.text = element_text(size = 12),
            legend.key = element_blank(),
            legend.margin = margin(t = 16),
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14)
          ) +
          guides(colour = guide_legend(override.aes = list(size = 5, alpha = 0.95, stroke = 0)))
        multi_pca_scores_plot_local <- gp_scores
        assign("multi_pca_scores_plot", gp_scores, envir = .GlobalEnv)

        sample_levels_pca <- sample_levels_global[sample_levels_global %in% scores$sampleID]
        if (length(sample_levels_pca) == 0) {
          sample_levels_pca <- unique(scores$sampleID)
        }
        palette_cols_pca <- get_conesa_palette_colors(length(sample_levels_pca), palette = "complete")
        sample_color_map_pca <- setNames(palette_cols_pca, sample_levels_pca)
        scores_plotly <- plotly::plot_ly()
        for (sample_nm in sample_levels_pca) {
          sample_df <- scores %>% filter(sampleID == sample_nm)
          if (nrow(sample_df) == 0) next
          scores_plotly <- scores_plotly %>% plotly::add_trace(
            data = sample_df,
            x = ~PC1,
            y = ~PC2,
            type = "scatter",
            mode = "markers",
            name = sample_nm,
            text = ~sampleID,
            hovertemplate = "Sample: %{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>",
            marker = list(size = 12, color = sample_color_map_pca[[sample_nm]], line = list(width = 0))
          )
        }
        scores_plotly <- scores_plotly %>%
          plotly::layout(
            title = list(text = "<b>PCA Plot Based on sampleID</b>", font = list(size = 20, family = "Arial")),
            xaxis = list(
              title = list(text = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]), font = list(size = 16, family = "Arial")),
              tickfont = list(size = 14, family = "Arial"), standoff = 12, automargin = TRUE
            ),
            yaxis = list(
              title = list(text = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]), font = list(size = 16, family = "Arial")),
              tickfont = list(size = 14, family = "Arial"), standoff = 10, automargin = TRUE
            ),
            legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25, yanchor = "top", title = list(text = ""), font = list(size = 14, family = "Arial")),
            paper_bgcolor = "rgba(0,0,0,0)",
            plot_bgcolor = "rgba(0,0,0,0)",
            font = list(family = "Arial", size = 16),
            margin = list(t = 60, b = 120, l = 90, r = 60),
            height = 520
          )
        multi_pca_scores_plot_html_local <- scores_plotly
        assign("multi_pca_scores_plot_html", scores_plotly, envir = .GlobalEnv)
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
        theme_classic(base_size = 16) +
        labs(title = "PCA scree plot", y = "Variance explained", x = "Principal component") +
        theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "bottom",
          legend.margin = margin(t = 20)
        )
      multi_pca_scree_plot_local <- gp_scree
      assign("multi_pca_scree_plot", gp_scree, envir = .GlobalEnv)

      scree_plotly <- plotly::plot_ly(
        scree_df,
        x = ~PC,
        y = ~Proportion,
        type = "bar",
        name = "Proportion",
        marker = list(color = "#6BAED6")
      ) %>%
        plotly::add_trace(
          y = ~Cumulative,
          type = "scatter",
          mode = "lines+markers",
          name = "Cumulative",
          hovertemplate = "PC: %{x}<br>Cumulative: %{y:.3f}<extra></extra>",
          line = list(color = "#4D4D4D", width = 2),
          marker = list(color = "#4D4D4D", size = 9, line = list(width = 0))
        ) %>%
        plotly::layout(
          title = list(text = "<b>PCA scree plot</b>", font = list(size = 18, family = "Arial")),
          yaxis = list(
            title = list(text = "Variance explained", font = list(size = 16, family = "Arial")), tickfont = list(size = 14),
            standoff = 6, automargin = TRUE
          ),
          xaxis = list(
            title = list(text = "Principal component", font = list(size = 16, family = "Arial")), tickfont = list(size = 14),
            standoff = 14, automargin = TRUE
          ),
          legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25, yanchor = "top", title = list(text = "")),
          paper_bgcolor = "rgba(0,0,0,0)",
          plot_bgcolor = "rgba(0,0,0,0)",
          margin = list(t = 60, b = 120, l = 80, r = 40)
        )
      multi_pca_scree_plot_html_local <- scree_plotly
      assign("multi_pca_scree_plot_html", scree_plotly, envir = .GlobalEnv)

      # C) Top loadings for PC1 and PC2 (third)
      if (ncol(pca_fit$rotation) >= 2) {
        rot <- as.data.frame(pca_fit$rotation)
        rot$variable <- rownames(rot)
        top_n <- 10L
        pick_top <- function(colname) {
          ord <- order(abs(rot[[colname]]), decreasing = TRUE)
          head(rot[ord, c("variable", colname)], top_n)
        }
        top_pc1 <- pick_top("PC1")
        colnames(top_pc1) <- c("variable", "loading")
        top_pc1 <- top_pc1 %>%
          mutate(
            variable = as.character(variable),
            PC = "PC1",
            rank = dplyr::row_number(),
            sign = if_else(loading >= 0, "Positive", "Negative"),
            abs_loading = abs(loading)
          )
        top_pc2 <- pick_top("PC2")
        colnames(top_pc2) <- c("variable", "loading")
        top_pc2 <- top_pc2 %>%
          mutate(
            variable = as.character(variable),
            PC = "PC2",
            rank = dplyr::row_number(),
            sign = if_else(loading >= 0, "Positive", "Negative"),
            abs_loading = abs(loading)
          )

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
          theme_classic(base_size = 16) +
          labs(title = "Top 10 loadings: PC1", x = "Feature", y = "Absolute loading") +
          theme(
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            legend.position = "bottom"
          )
        gp_load2 <- ggplot(top_pc2_plot, aes(x = variable, y = abs_loading, fill = sign)) +
          geom_col(width = 0.7) +
          coord_flip() +
          scale_fill_manual(values = c("Positive" = "#78C679", "Negative" = "#EE6A50"), name = "Sign", limits = c("Positive", "Negative"), drop = FALSE) +
          theme_classic(base_size = 16) +
          labs(title = "Top 10 loadings: PC2", x = "Feature", y = "Absolute loading") +
          theme(
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            legend.position = "bottom"
          )
        loadings_plots <- list(PC1 = gp_load1, PC2 = gp_load2)
        multi_pca_top_loadings_plots_local <- loadings_plots
        assign("multi_pca_top_loadings_plots", loadings_plots, envir = .GlobalEnv)

        pc1_plot_html <- plotly::plot_ly(
          top_pc1_plot,
          x = ~abs_loading,
          y = ~variable,
          color = ~sign,
          colors = c("Positive" = "#78C679", "Negative" = "#EE6A50"),
          type = "bar",
          orientation = "h",
          customdata = ~sign,
          hovertemplate = "Feature: %{y}<br>|loading|: %{x:.3f}<br>Sign: %{customdata}<extra></extra>"
        ) %>%
          plotly::layout(
            title = list(text = "<b>Top 10 loadings: PC1</b>", font = list(size = 18, family = "Arial")),
            xaxis = list(
              title = list(text = "Absolute loading", font = list(size = 16, family = "Arial")),
              tickfont = list(size = 14),
              automargin = TRUE,
              standoff = 10
            ),
            yaxis = list(
              title = list(text = "Feature", font = list(size = 16, family = "Arial")),
              tickfont = list(size = 14),
              automargin = TRUE,
              standoff = 6
            ),
            legend = list(
              orientation = "h",
              x = 0.5,
              xanchor = "center",
              y = -0.25,
              yanchor = "top",
              font = list(size = 14, family = "Arial"),
              title = list(text = "Sign")
            ),
            paper_bgcolor = "rgba(0,0,0,0)",
            plot_bgcolor = "rgba(0,0,0,0)",
            margin = list(t = 60, b = 110, l = 190, r = 70),
            font = list(family = "Arial", size = 14),
            height = 460
          )

        pc2_plot_html <- plotly::plot_ly(
          top_pc2_plot,
          x = ~abs_loading,
          y = ~variable,
          color = ~sign,
          colors = c("Positive" = "#78C679", "Negative" = "#EE6A50"),
          type = "bar",
          orientation = "h",
          customdata = ~sign,
          hovertemplate = "Feature: %{y}<br>|loading|: %{x:.3f}<br>Sign: %{customdata}<extra></extra>"
        ) %>%
          plotly::layout(
            title = list(text = "<b>Top 10 loadings: PC2</b>", font = list(size = 18, family = "Arial")),
            xaxis = list(
              title = list(text = "Absolute loading", font = list(size = 16, family = "Arial")),
              tickfont = list(size = 14),
              automargin = TRUE,
              standoff = 10
            ),
            yaxis = list(
              title = list(text = "", font = list(size = 16, family = "Arial")),
              tickfont = list(size = 14),
              automargin = TRUE,
              standoff = 6
            ),
            legend = list(
              orientation = "h",
              x = 0.5,
              xanchor = "center",
              y = -0.25,
              yanchor = "top",
              font = list(size = 14, family = "Arial"),
              title = list(text = "Sign")
            ),
            paper_bgcolor = "rgba(0,0,0,0)",
            plot_bgcolor = "rgba(0,0,0,0)",
            margin = list(t = 60, b = 110, l = 190, r = 70),
            font = list(family = "Arial", size = 14),
            height = 460
          )

        loadings_plots_html <- list(PC1 = pc1_plot_html, PC2 = pc2_plot_html)
        multi_pca_top_loadings_plots_html_local <- loadings_plots_html
        assign("multi_pca_top_loadings_plots_html", loadings_plots_html, envir = .GlobalEnv)

        positive_col <- "#78C679"
        negative_col <- "#EE6A50"
        pc1_levels <- rev(as.character(top_pc1_plot$variable))
        pc2_levels <- rev(as.character(top_pc2_plot$variable))

        combined_loadings_html <- plotly::plot_ly()
        for (panel in c("PC1", "PC2")) {
          axis_suffix <- if (panel == "PC1") "" else "2"
          panel_df <- if (panel == "PC1") top_pc1_plot else top_pc2_plot
          panel_df <- panel_df %>% mutate(variable = as.character(variable))
          for (sgn in c("Positive", "Negative")) {
            sgn_df <- panel_df %>% filter(sign == sgn)
            if (nrow(sgn_df) == 0) next
            combined_loadings_html <- combined_loadings_html %>%
              plotly::add_trace(
                data = sgn_df,
                x = ~abs_loading,
                y = ~variable,
                type = "bar",
                orientation = "h",
                name = sgn,
                legendgroup = sgn,
                showlegend = (panel == "PC1"),
                marker = list(color = if (sgn == "Positive") positive_col else negative_col),
                hovertemplate = paste0("Feature: %{y}<br>|loading|: %{x:.3f}<br>Sign: ", sgn, "<extra></extra>"),
                xaxis = paste0("x", axis_suffix),
                yaxis = paste0("y", axis_suffix)
              )
          }
        }

        combined_loadings_html <- combined_loadings_html %>%
          plotly::layout(
            barmode = "stack",
            xaxis = list(
              domain = c(0, 0.35),
              title = list(text = "Absolute loading", font = list(size = 14, family = "Arial")),
              tickfont = list(size = 14),
              automargin = TRUE,
              standoff = 12
            ),
            yaxis = list(
              domain = c(0, 1),
              title = list(text = "Feature", standoff = 30, automargin = TRUE, font = list(size = 14, family = "Arial")),
              tickfont = list(size = 14),
              automargin = TRUE,
              categoryorder = "array",
              categoryarray = pc1_levels
            ),
            xaxis2 = list(
              domain = c(0.65, 1),
              title = list(text = "Absolute loading", font = list(size = 14, family = "Arial")),
              tickfont = list(size = 14),
              automargin = TRUE,
              standoff = 12,
              anchor = "y2"
            ),
            yaxis2 = list(
              domain = c(0, 1),
              title = list(text = "", font = list(size = 14, family = "Arial")),
              tickfont = list(size = 14),
              automargin = TRUE,
              categoryorder = "array",
              categoryarray = pc2_levels,
              anchor = "x2"
            ),
            legend = list(
              orientation = "h",
              x = 0.5,
              xanchor = "center",
              y = -0.25,
              yanchor = "top",
              font = list(size = 14, family = "Arial"),
              title = list(text = "Sign")
            ),
            annotations = list(
              list(text = "<b>Top 10 loadings: PC1</b>", x = 0.13, y = 1.08, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size = 16, family = "Arial")),
              list(text = "<b>Top 10 loadings: PC2</b>", x = 0.9, y = 1.08, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size = 16, family = "Arial"))
            ),
            paper_bgcolor = "rgba(0,0,0,0)",
            plot_bgcolor = "rgba(0,0,0,0)",
            font = list(family = "Arial", size = 14),
            margin = list(t = 60, b = 140, l = 220, r = 160),
            height = 640
          )
        assign("multi_pca_top_loadings_combined_html", combined_loadings_html, envir = .GlobalEnv)

        # D) Distribution plots for top-loading features on PC1/PC2
        sample_levels <- sample_levels_global
        loading_plot_info <- bind_rows(top_pc1, top_pc2) %>%
          distinct(variable, .keep_all = TRUE)
        loading_distribution_plots <- list()
        loading_distribution_plots_html <- list()
        if (nrow(loading_plot_info) > 0) {
          for (idx in seq_len(nrow(loading_plot_info))) {
            gp_loading <- build_loading_feature_plot(multi, loading_plot_info[idx, ], sample_levels)
            feat_name <- loading_plot_info$variable[idx]
            if (is.null(gp_loading)) {
              message(sprintf("[INFO] Skipping PCA loading feature %s due to missing or constant data.", feat_name))
            } else {
              loading_distribution_plots[[feat_name]] <- gp_loading$ggplot
              loading_distribution_plots_html[[feat_name]] <- gp_loading$plotly
            }
          }
        }
        multi_pca_loading_distribution_plots_local <- loading_distribution_plots
        multi_pca_loading_distribution_plots_html_local <- loading_distribution_plots_html
        assign("multi_pca_loading_distribution_plots", loading_distribution_plots, envir = .GlobalEnv)
        assign("multi_pca_loading_distribution_plots_html", loading_distribution_plots_html, envir = .GlobalEnv)
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

    message("Generating HTML report...")
    message("Rmd file: ", rmd_file)
    message("Output file: ", html_output_file)

    rmarkdown::render(
      rmd_file,
      output_file = html_output_file,
      envir = globalenv(),
      quiet = TRUE
    )

    message("HTML report generated: ", html_output_file)
  }
}

main()
