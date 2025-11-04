#!/usr/env/bin Rscript

############################################################
##### SQANTI single-cell multisample report generation #####
############################################################



### Author: Carlos Blanco

#********************** Packages

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(ggdist)
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
        stop(sprintf("Unknown flag %s", key))
      }
    } else {
      stop(sprintf("Unexpected argument: %s", key))
    }
  }
  if (is.null(res$files) || nchar(res$files) == 0) {
    stop("--files must be provided (comma-separated list of cell summary files)")
  }
  if (!(res$report %in% c("pdf", "html", "both"))) {
    stop("--report must be one of: pdf, html, both")
  }
  res
}

safe_read_summary <- function(fpath) {
  # read.table supports gz automatically
  df <- tryCatch({
    read.table(fpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    message(sprintf("[ERROR] Failed to read summary %s: %s", fpath, e$message))
    return(NULL)
  })
  if (is.null(df)) return(NULL)
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

  # Basic cohort-level aggregates per sample
  per_sample_stats <- multi %>%
    group_by(sampleID) %>%
    summarise(
      cells = n_distinct(CB),
      mean_reads = mean(Reads_in_cell, na.rm = TRUE),
      median_reads = median(Reads_in_cell, na.rm = TRUE),
      mean_genes = mean(Genes_in_cell, na.rm = TRUE),
      median_genes = median(Genes_in_cell, na.rm = TRUE),
      mean_annotated = mean(Annotated_genes, na.rm = TRUE),
      mean_novel = mean(Novel_genes, na.rm = TRUE),
      mean_ujc = mean(UJCs_in_cell, na.rm = TRUE),
      median_ujc = median(UJCs_in_cell, na.rm = TRUE),
      mean_mt = mean(MT_perc, na.rm = TRUE)
    )

  # Output path
  out_dir <- params$out_dir
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  pdf_out <- file.path(out_dir, paste0(params$prefix, ".pdf"))

  if (params$report %in% c("html", "both")) {
    message("[INFO] This R script produces a PDF. HTML implementation is WIP.")
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
          levels = c("FSM_prop", "ISM_prop", "NIC_prop", "NNC_prop",
                     "Genic_Genomic_prop", "Antisense_prop", "Fusion_prop", "Intergenic_prop", "Genic_intron_prop"),
          labels = c("FSM", "ISM", "NIC", "NNC",
                     "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic intron")
        )
      ) %>%
      mutate(prop = suppressWarnings(as.numeric(prop)))

    # Prepare a tiny dataset to generate clean legend keys (off-plot)
    legend_df <- cats_long %>%
      distinct(sampleID) %>%
      mutate(
        category = factor(levels(cats_long$category)[1], levels = levels(cats_long$category)),
        prop = 0
      )

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
      labs(title = "Structural Category Proportions by Sample",
           x = "Structural category", y = "Reads, %") +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            axis.title = element_text(size = 16),
            axis.text.y = element_text(size = 14),
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  }

  # Write PDF
  pdf(pdf_out, paper = "a4r", width = 14, height = 11)
  # Cover page
  grid.newpage()
  cover <- textGrob("SQANTI-single cell\nmulti-sample reads report",
                    gp = gpar(fontface = "italic", fontsize = 40, col = "orangered"))
  grid.draw(cover)
  # Per cell summary table (centered)
  summary_tbl <- per_sample_stats %>%
    mutate(across(where(is.numeric), ~round(., 3))) %>%
    transmute(
      Sample = sampleID,
      `Cell\nBarcodes` = cells,
      `Average\nReads` = mean_reads,
      `Average\nAnnotated\nGenes` = mean_annotated,
      `Average\nNovel\nGenes` = mean_novel,
      `Average\nUJCs` = mean_ujc,
      `Average\nMitochondrial\nReads` = mean_mt
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
  if (have_cats) {
    # --- Per-category violins: one page per category, all samples shown together (FIRST) ---
    category_levels <- levels(cats_long$category)
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
    for (cat_lab in category_levels) {
      dfp <- cats_long %>% filter(category == cat_lab)
      cat_col <- unname(cat_to_col[[as.character(cat_lab)]])
      if (is.null(cat_col) || is.na(cat_col)) cat_col <- "grey60"
      box_outline_col <- if (as.character(cat_lab) == "Genic Genomic") "grey90" else "grey20"
      gp <- ggplot(dfp, aes(x = sampleID, y = prop)) +
        geom_violin(fill = cat_col, color = cat_col, alpha = 0.7, width = 0.8, trim = TRUE) +
        geom_boxplot(width = 0.08, outlier.shape = NA, fill = cat_col, color = box_outline_col, alpha = 0.6) +
        scale_y_continuous(limits = c(0, 100), expand = expansion(add = c(1, 0))) +
        theme_classic(base_size = 14) +
        labs(title = paste0("Per Sample ", cat_lab, " Reads Distribution Across Cells"), x = "Sample", y = "Reads, %") +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14))
      print(gp)
    }
    
    # --- Structural category slabs across all samples (SECOND) ---
    print(p_cats)
  }
  
  # -------- PCA (all numeric features, per-sample medians) --------
  # 1) Select all numeric columns from the cell summary
  num_cols <- names(multi)[sapply(multi, function(x) is.numeric(x) && !all(is.na(x)))]
  # 2) Aggregate per-sample medians across all numeric features
  agg_median <- multi %>%
    group_by(sampleID) %>%
    summarise(across(all_of(num_cols), ~median(., na.rm = TRUE)), .groups = "drop")
  
  if (nrow(agg_median) >= 2 && ncol(agg_median) >= 2) {
    # 3) Drop features with zero variance across samples
    feat_sds <- sapply(agg_median %>% select(-sampleID), function(x) stats::sd(x, na.rm = TRUE))
    feat_keep <- names(feat_sds)[is.finite(feat_sds) & !is.na(feat_sds) & feat_sds > 0]
    
    if (length(feat_keep) >= 2) {
      mat <- as.matrix(agg_median[, feat_keep, drop = FALSE])
      rownames(mat) <- agg_median$sampleID
      pca_fit <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
      var_expl <- (pca_fit$sdev ^ 2) / sum(pca_fit$sdev ^ 2)
      
      # A) PC1–PC2 scatter (first among PCA plots)
      if (ncol(pca_fit$x) >= 2) {
        scores <- as.data.frame(pca_fit$x)
        scores$sampleID <- rownames(scores)
        gp_scores <- ggplot(scores, aes(x = PC1, y = PC2, colour = sampleID, label = sampleID)) +
          geom_point(size = 3, alpha = 0.95) +
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
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14)
          )
        print(gp_scores)
      }
      
      # B) Scree plot (second)
      k <- min(length(var_expl), 10)
      scree_df <- data.frame(
        PC = factor(paste0("PC", seq_len(k)), levels = paste0("PC", seq_len(k))),
        Proportion = var_expl[seq_len(k)],
        Cumulative = cumsum(var_expl)[seq_len(k)]
      )
      gp_scree <- ggplot(scree_df, aes(x = PC, y = Proportion)) +
        geom_col(fill = "#6BAED6", width = 0.8) +
        geom_point(aes(y = Cumulative), color = "grey20", size = 1.6) +
        geom_line(aes(y = Cumulative, group = 1), color = "grey20", linewidth = 0.3) +
        theme_classic(base_size = 14) +
        labs(title = "PCA scree plot", y = "Variance explained", x = "Principal component") +
        theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14)
        )
      print(gp_scree)
      
      # C) Top loadings for PC1 and PC2 (third)
      if (ncol(pca_fit$rotation) >= 2) {
        rot <- as.data.frame(pca_fit$rotation)
        rot$variable <- rownames(rot)
        top_n <- 10L
        pick_top <- function(colname) {
          ord <- order(abs(rot[[colname]]), decreasing = TRUE)
          head(rot[ord, c("variable", colname)], top_n)
        }
        top_pc1 <- pick_top("PC1"); colnames(top_pc1) <- c("variable", "loading"); top_pc1$PC <- "PC1"
        top_pc2 <- pick_top("PC2"); colnames(top_pc2) <- c("variable", "loading"); top_pc2$PC <- "PC2"

        # Encode sign via color and plot absolute magnitudes to save space
        top_pc1$sign <- ifelse(top_pc1$loading >= 0, "Positive", "Negative")
        top_pc1$abs_loading <- abs(top_pc1$loading)
        top_pc2$sign <- ifelse(top_pc2$loading >= 0, "Positive", "Negative")
        top_pc2$abs_loading <- abs(top_pc2$loading)
        top_pc1$sign <- factor(top_pc1$sign, levels = c("Positive", "Negative"))
        top_pc2$sign <- factor(top_pc2$sign, levels = c("Positive", "Negative"))

        top_pc1$variable <- factor(top_pc1$variable, levels = rev(top_pc1$variable))
        top_pc2$variable <- factor(top_pc2$variable, levels = rev(top_pc2$variable))

        gp_load1 <- ggplot(top_pc1, aes(x = variable, y = abs_loading, fill = sign)) +
          geom_col() + coord_flip() +
          scale_fill_manual(values = c("Positive" = "#78C679", "Negative" = "#EE6A50"), name = "Sign", limits = c("Positive", "Negative"), drop = FALSE) +
          theme_classic(base_size = 14) +
          labs(title = "Top 10 loadings: PC1", x = "Feature", y = "Absolute loading") +
          theme(
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14)
          )
        gp_load2 <- ggplot(top_pc2, aes(x = variable, y = abs_loading, fill = sign)) +
          geom_col() + coord_flip() +
          scale_fill_manual(values = c("Positive" = "#78C679", "Negative" = "#EE6A50"), name = "Sign", limits = c("Positive", "Negative"), drop = FALSE) +
          theme_classic(base_size = 14) +
          labs(title = "Top 10 loadings: PC2", x = "Feature", y = "Absolute loading") +
          theme(
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14)
          )
        # Build shared legend with both levels present, then arrange below plots
        legend_df <- data.frame(
          variable = c("pos", "neg"),
          abs_loading = c(1, 1),
          sign = factor(c("Positive", "Negative"), levels = c("Positive", "Negative"))
        )
        legend_plot <- ggplot(legend_df, aes(x = variable, y = abs_loading, fill = sign)) +
          geom_col() +
          scale_fill_manual(values = c("Positive" = "#78C679", "Negative" = "#EE6A50"), name = "Sign", limits = c("Positive", "Negative"), drop = FALSE) +
          theme_void(base_size = 14) + theme(legend.position = "bottom")
        legend_grob <- gtable::gtable_filter(ggplotGrob(legend_plot), "guide-box")
        row_plots <- arrangeGrob(gp_load1 + theme(legend.position = "none"), gp_load2 + theme(legend.position = "none"), ncol = 2)
        grid.arrange(row_plots, legend_grob, ncol = 1, heights = c(0.86, 0.14))
      }
    }
  }

  
  dev.off()

  message(sprintf("**** Multisample report written: %s", pdf_out))
}

main()
