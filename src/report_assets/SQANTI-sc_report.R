#!/usr/env/bin Rscript

######################################################
##### SQANTI single-cell reads report generation #####
######################################################

### Author: Juan Francisco Cervilla & Carlos Blanco

#********************** Packages

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(grid)
  library(gridExtra)
  library(rmarkdown)
  library(scales)
  library(data.table)
}))

# Prevent Rplots.pdf generation
pdf(NULL)

#********************** Taking arguments from python script

args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1]
junc.file <- args[2]
report.format <- args[3]
outputPathPrefix <- args[4]
mode <- args[5]

# Initialize ignore_cell_summary flag
ignore_cell_summary <- FALSE
include_ORF <- FALSE
CAGE_peak <- FALSE
polyA_motif_list <- FALSE
cell_summary_path <- NULL
ref_gtf_path <- NULL

# Check for optional arguments
if (length(args) > 5) {
  i <- 6
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--ignore_cell_summary") {
      ignore_cell_summary <- TRUE
      i <- i + 1
      next
    }
    if (arg == "--include_ORF") {
      include_ORF <- TRUE
      i <- i + 1
      next
    }
    if (arg == "--CAGE_peak") {
      CAGE_peak <- TRUE
      i <- i + 1
      next
    }
    if (arg == "--polyA_motif_list") {
      polyA_motif_list <- TRUE
      i <- i + 1
      next
    }
    if (arg == "--cell_summary") {
      if ((i + 1) <= length(args)) {
        cell_summary_path <- args[i + 1]
        i <- i + 2
        next
      } else {
        stop("--cell_summary requires a path argument")
      }
    }
    if (arg == "--clustering") {
      if ((i + 1) <= length(args)) {
        clustering_path <- args[i + 1]
        i <- i + 2
        next
      } else {
        stop("--clustering requires a path argument")
      }
    }
    if (arg == "--refGTF") {
      if ((i + 1) <= length(args)) {
        ref_gtf_path <- args[i + 1]
        i <- i + 2
        next
      } else {
        stop("--refGTF requires a path argument")
      }
    }
    i <- i + 1
  }
}

# Validate arguments
if (length(args) < 5) {
  stop("Incorrect number of arguments! Required: [classification file] [junc file] [report format] [outputPathPrefix] [mode]. Abort!")
}

if (!(report.format %in% c("pdf", "html", "both"))) {
  stop("Report format needs to be: pdf, html, or both. Abort!")
}

# Validate mode argument
if (!(mode %in% c("reads", "isoforms"))) {
  stop("Mode needs to be: reads or isoforms. Abort!")
}

# Set labels based on mode
if (mode == "isoforms") {
  entity_label <- "Transcript"
  entity_label_plural <- "Transcripts"
} else {
  entity_label <- "Read"
  entity_label_plural <- "Reads"
}

# Lowercase versions for inline text
entity_label_lower <- tolower(entity_label)
entity_label_plural_lower <- tolower(entity_label_plural)

# Print cell summary saving status
if (ignore_cell_summary) {
  print("Cell summary table will not be saved (--ignore_cell_summary flag is active).")
} else {
  print("Cell summary table will be saved.")
}

# Call the function with the appropriate Save parameter
save_option <- ifelse(ignore_cell_summary, "N", "Y")

# Define column names based on mode
if (mode == "isoforms") {
  count_col <- "Transcripts_in_cell"
  no_mono_col <- "total_transcripts_no_monoexon"
} else {
  count_col <- "Reads_in_cell"
  no_mono_col <- "total_reads_no_monoexon"
}

# Generate output file names with full paths
cell_summary_output <- file.path(paste0(outputPathPrefix, "_SQANTI_cell_summary"))
report_output <- file.path(paste0(outputPathPrefix, "_SQANTI_sc_report_", mode))
clustering_output <- file.path(dirname(outputPathPrefix), "clustering", "umap_results.csv")

# Define standard colors
fill_color_orange <- "#CC6633"

# Check for clustering results
gg_umap <- NULL
if (file.exists(clustering_output)) {
  print(paste("Found clustering results at:", clustering_output))
  tryCatch(
    {
      umap_df <- read.csv(clustering_output)
      umap_df$Cluster <- as.factor(umap_df$Cluster)

      gg_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
        geom_point(alpha = 0.6, size = 0.5) +
        theme_classic() +
        labs(title = "UMAP Projection", x = "UMAP 1", y = "UMAP 2") +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 16, face = "bold"),
          legend.position = "right"
        ) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    },
    error = function(e) {
      print(paste("Error reading clustering results:", e$message))
    }
  )
} else {
  print("No clustering results found.")
}

# ----------------------------------------------------------------
# Helper Functions (Global Scope)
# ----------------------------------------------------------------

# Helper: pivot selected columns to long and return factor-ordered long df
pivot_long <- function(df, cols) {
  out <- pivot_longer(df, cols = all_of(cols), names_to = "Variable", values_to = "Value") %>%
    select(Variable, Value)
  out$Variable <- factor(out$Variable, levels = cols)
  out
}

# Helper: generic violin + box + mean-cross plot with shared theme
build_violin_plot <- function(df_long,
                              title,
                              x_labels,
                              fill_map,
                              color_map = fill_map,
                              x_title = "",
                              y_label = paste(entity_label_plural, ", %", sep = ""),
                              legend = FALSE,
                              ylim = NULL,
                              override_outline_vars = character(0),
                              violin_alpha = 0.7,
                              box_alpha = 0.6,
                              box_width = 0.05,
                              x_tickangle = 45,
                              violin_outline_fill = FALSE,
                              box_outline_default = "grey20",
                              adjust = 1,
                              log_scale = FALSE,
                              ...) {
  # Clamp values for percentage / count plots
  if (grepl("%", y_label)) {
    df_long$Value <- pmin(pmax(df_long$Value, 0), 100)
  } else if (grepl("count", y_label, ignore.case = TRUE)) {
    df_long$Value <- pmax(df_long$Value, 0)
  }

  # Percentage violins: fixed 0–100 y-axis so plots are comparable, especially between categories
  if (is.null(ylim) && !isTRUE(log_scale) && grepl("%", y_label)) {
    ylim <- c(0, 100)
  }

  if (isTRUE(log_scale)) {
    df_long <- df_long[is.finite(df_long$Value) & df_long$Value > 0, , drop = FALSE]
  }

  p <- ggplot(df_long, aes(x = Variable, y = Value)) +
    {
      if (isTRUE(violin_outline_fill)) {
        geom_violin(aes(fill = Variable, color = Variable), alpha = violin_alpha, scale = "width", show.legend = legend, adjust = adjust, trim = TRUE)
      } else {
        geom_violin(aes(fill = Variable), color = "black", alpha = violin_alpha, scale = "width", show.legend = legend, adjust = adjust, trim = TRUE)
      }
    } +
    scale_fill_manual(values = fill_map, labels = x_labels) +
    {
      if (isTRUE(violin_outline_fill)) scale_color_manual(values = fill_map, guide = "none") else NULL
    } +
    scale_x_discrete(labels = x_labels) +
    labs(title = title, x = x_title, y = y_label) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16, angle = x_tickangle, hjust = ifelse(x_tickangle == 0, 0.5, 1)),
      legend.position = if (legend) "bottom" else "none"
    )

  # Add boxplots per variable with correct outline color (grey90 overrides)
  for (var in levels(df_long$Variable)) {
    var_df <- df_long[df_long$Variable == var, , drop = FALSE]
    box_col <- if (var %in% override_outline_vars) "grey90" else box_outline_default
    p <- p + geom_boxplot(
      data = var_df,
      aes(x = Variable, y = Value, fill = Variable),
      width = box_width, outlier.shape = NA, alpha = box_alpha, show.legend = FALSE, color = box_col, lwd = 0.3
    )
  }

  # Add mean markers on top (moved here to separate from violin layer and ensure it is on top of boxplots)
  p <- p + stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1, show.legend = FALSE)

  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }

  if (isTRUE(log_scale)) {
    p <- p + scale_y_log10(labels = scales::comma)
  }

  return(p)
}

generate_sqantisc_plots <- function(SQANTI_cell_summary, Classification_file, Junctions, report_output, generate_pdf = TRUE) {
  # Helper function to mix colors
  mix_color <- function(col, target, amount) {
    c_rgb <- col2rgb(col)
    t_rgb <- col2rgb(target)
    mix <- c_rgb * (1 - amount) + t_rgb * amount
    rgb(mix[1], mix[2], mix[3], maxColorValue = 255)
  }

  # Generate UMAP plots by structural category if UMAP exists
  if (exists("gg_umap") && !is.null(gg_umap)) {
    tryCatch(
      {
        umap_data <- gg_umap$data

        # Merge with SQANTI_cell_summary
        # umap_data has 'Barcode', SQANTI_cell_summary has 'CB'
        merged_umap <- inner_join(umap_data, SQANTI_cell_summary, by = c("Barcode" = "CB"))

        if (nrow(merged_umap) > 0) {
          gg_umap_by_category <<- list()

          # Define categories and their colors
          cat_colors <- c(
            "FSM_prop" = "#6BAED6",
            "ISM_prop" = "#FC8D59",
            "NIC_prop" = "#78C679",
            "NNC_prop" = "#EE6A50",
            "Genic_Genomic_prop" = "#969696",
            "Antisense_prop" = "#66C2A4",
            "Fusion_prop" = "goldenrod1",
            "Intergenic_prop" = "darksalmon",
            "Genic_intron_prop" = "#41B6C4"
          )

          cat_labels <- c(
            "FSM_prop" = "FSM",
            "ISM_prop" = "ISM",
            "NIC_prop" = "NIC",
            "NNC_prop" = "NNC",
            "Genic_Genomic_prop" = "Genic Genomic",
            "Antisense_prop" = "Antisense",
            "Fusion_prop" = "Fusion",
            "Intergenic_prop" = "Intergenic",
            "Genic_intron_prop" = "Genic Intron"
          )

          for (cat_col in names(cat_colors)) {
            if (cat_col %in% colnames(merged_umap)) {
              cat_color <- cat_colors[[cat_col]]
              cat_label <- cat_labels[[cat_col]]

              dark_color <- mix_color(cat_color, "black", 0.6)
              light_color <- mix_color(cat_color, "white", 0.8)

              p <- ggplot(merged_umap, aes(x = UMAP_1, y = UMAP_2, color = .data[[cat_col]])) +
                geom_point(alpha = 0.6, size = 0.5) +
                theme_classic() +
                labs(title = paste("UMAP - %", cat_label), x = "UMAP 1", y = "UMAP 2", color = paste0(entity_label_plural, ", %")) +
                theme(
                  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  axis.title = element_text(size = 18),
                  axis.text.x = element_text(size = 16),
                  axis.text.y = element_text(size = 16),
                  legend.title = element_text(size = 16, face = "bold"),
                  legend.position = "right",
                  legend.key.height = unit(3, "cm"),
                  legend.key.width = unit(1, "cm") # Thicker legend bar
                ) +
                scale_color_gradientn(colors = c(light_color, cat_color, dark_color)) + # Custom gradient
                guides(color = guide_colorbar(barwidth = 2.5, barheight = 15)) # Make legend bar thicker and taller

              gg_umap_by_category[[cat_label]] <<- p
            }
          }

          # ----------------------------------------------------------------
          # Common Helpers for Cluster Violin Plots
          # ----------------------------------------------------------------
          prepare_violin_data <- function(data, y_col) {
            df <- data[!is.na(data[[y_col]]), c("Cluster", y_col)]
            colnames(df) <- c("Variable", "Value")
            df$Variable <- as.factor(df$Variable)
            return(df)
          }

          unique_clusters <- levels(merged_umap$Cluster)
          cluster_colors <- scales::hue_pal()(length(unique_clusters))
          names(cluster_colors) <- unique_clusters

          # ----------------------------------------------------------------
          # Structural Categories Support by Cluster (Violin Plots)
          # ----------------------------------------------------------------
          gg_cat_cluster_plots <<- list()
          for (cat_col in names(cat_colors)) {
            if (cat_col %in% colnames(merged_umap)) {
              cat_color <- cat_colors[[cat_col]]
              cat_label <- cat_labels[[cat_col]]

              fixed_color_map <- rep(cat_color, length(unique_clusters))
              names(fixed_color_map) <- unique_clusters

              cat_data <- prepare_violin_data(merged_umap, cat_col)

              gg_cat_cluster_plots[[cat_label]] <<- build_violin_plot(
                df_long = cat_data,
                title = paste(cat_label, "Distribution"),
                x_labels = levels(cat_data$Variable),
                fill_map = fixed_color_map,
                y_label = paste(entity_label_plural, ", %", sep = ""),
                legend = FALSE,
                x_title = "Cluster",
                x_tickangle = 0,
                ylim = c(0, 100),
                violin_outline_fill = TRUE,
                violin_alpha = 0.7,
                box_alpha = 0.3
              )
            }
          }

          # ----------------------------------------------------------------
          # Length Distribution by Cluster (Violin/Box Plots)
          # ----------------------------------------------------------------
          tryCatch({
            cb_cluster_map <- umap_data[, c("Barcode", "Cluster"), drop = FALSE]
            cb_cluster_map <- cb_cluster_map[!is.na(cb_cluster_map$Barcode) & cb_cluster_map$Barcode != "", , drop = FALSE]

            cls_for_len <- Classification_file
            cls_for_len$length_num <- suppressWarnings(as.numeric(cls_for_len$length))

            if (mode == "isoforms" && "FL" %in% colnames(cls_for_len) && "CB" %in% colnames(cls_for_len)) {
              cls_for_len$CB_raw <- as.character(cls_for_len$CB)
              cls_for_len$FL_raw <- as.character(cls_for_len$FL)
              cls_for_len <- tidyr::separate_rows(cls_for_len, CB_raw, FL_raw, sep = ",")
              cls_for_len$FL_num   <- suppressWarnings(as.numeric(trimws(cls_for_len$FL_raw)))
              cls_for_len$FL_num[is.na(cls_for_len$FL_num) | cls_for_len$FL_num < 1] <- 1
              cls_for_len$CB_clean <- trimws(cls_for_len$CB_raw)
            } else if ("CB" %in% colnames(cls_for_len)) {
              cls_for_len$CB_clean <- as.character(cls_for_len$CB)
              cls_for_len$FL_num   <- 1
            } else {
              cls_for_len <- data.frame()
            }

            if (nrow(cls_for_len) > 0) {
              cls_for_len <- merge(cls_for_len, cb_cluster_map,
                                   by.x = "CB_clean", by.y = "Barcode", all.x = FALSE)
              cls_for_len <- cls_for_len[
                !is.na(cls_for_len$length_num) & cls_for_len$length_num > 0 &
                !is.na(cls_for_len$Cluster), , drop = FALSE]

              if (nrow(cls_for_len) > 0 && any(cls_for_len$FL_num > 1)) {
                rep_idx <- rep(seq_len(nrow(cls_for_len)), times = as.integer(cls_for_len$FL_num))
                cls_for_len <- cls_for_len[rep_idx, , drop = FALSE]
              }

              build_len_cluster_plot <- function(df_sub, title_str, colors = cluster_colors) {
                if (nrow(df_sub) == 0) return(NULL)
                df_sub$Cluster <- factor(df_sub$Cluster, levels = unique_clusters)
                ggplot(df_sub, aes(x = Cluster, y = length_num, fill = Cluster)) +
                  geom_violin(aes(color = Cluster), alpha = 0.7, scale = "width",
                              adjust = 1, trim = TRUE, show.legend = FALSE) +
                  scale_color_manual(values = colors, guide = "none") +
                  geom_boxplot(width = 0.05, alpha = 0.5, outlier.shape = NA,
                               color = "grey20", show.legend = FALSE) +
                  stat_summary(fun = mean, geom = "point", shape = 4,
                               size = 1, color = "red", stroke = 1, show.legend = FALSE) +
                  scale_fill_manual(values = colors) +
                  scale_y_log10(labels = scales::comma) +
                  labs(
                    title = title_str,
                    x     = "Cluster",
                    y     = "Length, bp"
                  ) +
                  theme_classic(base_size = 11) +
                  theme(
                    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                    axis.title = element_text(size = 18),
                    axis.text.y = element_text(size = 16),
                    axis.text.x = element_text(size = 16),
                    legend.position = "none"
                  )
              }

              all_label <- if (mode == "isoforms") "All Transcripts" else "All Reads"
              gg_len_cluster_plots <<- list()
              gg_len_cluster_plots[[all_label]] <<- build_len_cluster_plot(
                cls_for_len, paste(all_label, "Length Distribution by Cluster")
              )

              cat_sc_map <- c(
                "FSM"           = "full-splice_match",
                "ISM"           = "incomplete-splice_match",
                "NIC"           = "novel_in_catalog",
                "NNC"           = "novel_not_in_catalog",
                "Genic Genomic" = "genic",
                "Antisense"     = "antisense",
                "Fusion"        = "fusion",
                "Intergenic"    = "intergenic",
                "Genic Intron"  = "genic_intron"
              )
              for (lbl in names(cat_sc_map)) {
                # Find the matching cat_colors key (e.g. "FSM_prop" for lbl "FSM")
                cat_col_key <- names(cat_colors)[match(lbl, cat_labels)]
                cat_hex <- if (!is.na(cat_col_key) && cat_col_key %in% names(cat_colors))
                  cat_colors[[cat_col_key]] else "#888888"
                fixed_cat_colors <- rep(cat_hex, length(unique_clusters))
                names(fixed_cat_colors) <- unique_clusters

                sub_df <- cls_for_len[
                  !is.na(cls_for_len$structural_category) &
                    cls_for_len$structural_category == cat_sc_map[[lbl]], , drop = FALSE]
                if (nrow(sub_df) > 0) {
                  gg_len_cluster_plots[[lbl]] <<- build_len_cluster_plot(
                    sub_df, paste(lbl, "Length Distribution by Cluster"),
                    colors = fixed_cat_colors
                  )
                }
              }
            }
          }, error = function(e) {
            message("Could not build length-by-cluster plots: ", e$message)
          })

          # ----------------------------------------------------------------
          # Short Read Support by Cluster (Violin Plots)
          # ----------------------------------------------------------------
          # Use the new column name: srjunctions_support_prop
          if ("srjunctions_support_prop" %in% colnames(merged_umap) && sum(merged_umap$srjunctions_support_prop, na.rm = TRUE) > 0) {
            gg_sr_cluster_plots <<- list()

            # ----------------------------------------------------------------
            # TSS Ratio Validated Support by Cluster (Violin Plots)
            # ----------------------------------------------------------------
            if ("TSS_ratio_validated_prop" %in% colnames(merged_umap) && sum(merged_umap$TSS_ratio_validated_prop, na.rm = TRUE) > 0) {
              gg_tss_cluster_plots <<- list()

              # 1. All Transcripts Plot - TSS
              # Use Cluster Colors (same as Junctions Coverage)
              p_all_tss <- build_violin_plot(
                df_long = prepare_violin_data(merged_umap, "TSS_ratio_validated_prop"),
                title = "All Transcripts TSS Validation by Short Reads",
                x_labels = levels(merged_umap$Cluster),
                fill_map = cluster_colors,
                x_title = "Cluster",
                y_label = "TSS Ratio Validated, %",
                x_tickangle = 0,
                violin_outline_fill = TRUE,
                violin_alpha = 0.7,
                box_alpha = 0.3
              )
              gg_tss_cluster_plots[["All Transcripts"]] <<- p_all_tss

              # 2. Per-Category Plots - TSS
              # Use Category Color for ALL clusters
              for (cat_col in names(cat_colors)) {
                tag <- cat_labels[[cat_col]]
                prop_col <- paste0(tag, "_TSS_ratio_validated_prop")
                if (tag == "Genic Genomic") prop_col <- "Genic_TSS_ratio_validated_prop" # Handle Genic weirdness if needed
                if (tag == "Genic Intron") prop_col <- "Genic_intron_TSS_ratio_validated_prop"

                # Clean up tag to match column naming convention if straightforward
                simple_tag <- names(cat_labels)[which(cat_labels == tag)]
                simple_tag <- gsub("_prop", "", simple_tag)
                prop_col <- paste0(simple_tag, "_TSS_ratio_validated_prop")

                if (prop_col %in% colnames(merged_umap)) {
                  # Define single color map
                  current_cat_color <- cat_colors[[cat_col]]
                  fixed_color_map <- rep(current_cat_color, length(unique_clusters))
                  names(fixed_color_map) <- unique_clusters

                  p_cat_tss <- build_violin_plot(
                    df_long = prepare_violin_data(merged_umap, prop_col),
                    title = paste(tag, "TSS Validation by Short Reads"),
                    x_labels = levels(merged_umap$Cluster),
                    fill_map = fixed_color_map,
                    x_title = "Cluster",
                    y_label = "TSS Ratio Validated, %",
                    x_tickangle = 0,
                    violin_outline_fill = TRUE,
                    violin_alpha = 0.7,
                    box_alpha = 0.3
                  )
                  gg_tss_cluster_plots[[tag]] <<- p_cat_tss
                }
              }
            }

            # ----------------------------------------------------------------
            # Short Read Support by Cluster (Violin Plots)
            # ----------------------------------------------------------------
            # 1. Global Plot
            global_data <- prepare_violin_data(merged_umap, "srjunctions_support_prop")

            # Use build_violin_plot (which returns a ggplot object)
            gg_sr_cluster_plots[["All Transcripts"]] <<- build_violin_plot(
              df_long = global_data,
              title = "All Transcripts Junction Coverage by Short Reads",
              x_labels = levels(global_data$Variable),
              fill_map = cluster_colors,
              y_label = "Transcripts Supported, %",
              legend = FALSE,
              x_title = "Cluster",
              x_tickangle = 0,
              ylim = c(0, 100),
              violin_outline_fill = TRUE,
              violin_alpha = 0.7,
              box_alpha = 0.3
            )

            # 2. Per-Category Plots
            # For these, we want to maintain the specific Structural Category color Scheme?
            # The user said: "the colors used for the violins and boxes should be the same in each structural category and the color should be the corresponding to the structural category."
            # This implies that for the FSM plot, ALL clusters should be colored with the FSM color.

            for (cat_col in names(cat_labels)) {
              # cat_labels[[cat_col]] is e.g. "FSM", "Genic Genomic"
              tag <- cat_labels[[cat_col]]
              tag_clean <- gsub(" ", "_", tag)
              # New column name format: {TAG}_srjunctions_support_prop
              sr_col <- paste0(tag_clean, "_srjunctions_support_prop")

              if (sr_col %in% colnames(merged_umap)) {
                cat_data <- prepare_violin_data(merged_umap, sr_col)

                current_cat_color <- cat_colors[[cat_col]]
                fixed_color_map <- rep(current_cat_color, length(unique_clusters))
                names(fixed_color_map) <- unique_clusters

                title <- paste(tag, "Junction Coverage by Short Reads")

                gg_sr_cluster_plots[[tag]] <<- build_violin_plot(
                  df_long = cat_data,
                  title = title,
                  x_labels = levels(cat_data$Variable),
                  fill_map = fixed_color_map, # Per-category uses the category color for all clusters
                  y_label = "Transcripts Supported, %",
                  legend = FALSE,
                  x_title = "Cluster",
                  x_tickangle = 0,
                  ylim = c(0, 100),
                  violin_outline_fill = TRUE,
                  violin_alpha = 0.7,
                  box_alpha = 0.3
                )
              }
            }
          }
        }
      },
      error = function(e) {
        print(paste("Error generating UMAP by category plots:", e$message))
      }
    )
  }

  # ----------------------------------------------------------------
  # Helper: Build Continuous UMAP (for coloring by %)
  # ----------------------------------------------------------------
  build_continuous_umap <- function(data, color_col, title, color_base = "blue") {
    mix_color <- function(col, target, amount) {
      c_rgb <- col2rgb(col)
      t_rgb <- col2rgb(target)
      mix <- c_rgb * (1 - amount) + t_rgb * amount
      rgb(mix[1], mix[2], mix[3], maxColorValue = 255)
    }

    if (!all(c("UMAP_1", "UMAP_2") %in% colnames(data))) {
      return(NULL)
    }

    plot_data <- data[!is.na(data[[color_col]]), ]
    if (nrow(plot_data) == 0) {
      return(NULL)
    }

    dark_color <- mix_color(color_base, "black", 0.6)
    light_color <- mix_color(color_base, "white", 0.8)

    p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]])) +
      geom_point(alpha = 0.6, size = 0.5) +
      scale_color_gradientn(
        colors = c(light_color, color_base, dark_color),
        limits = c(0, max(plot_data[[color_col]], na.rm = TRUE))
      ) +
      guides(color = guide_colorbar(barwidth = 2.5, barheight = 15)) +
      labs(title = title, x = "UMAP 1", y = "UMAP 2", color = paste0(entity_label_plural, ", %")) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 16, face = "bold"),
        legend.key.height = unit(3, "cm"),
        legend.key.width = unit(1, "cm")
      )
    return(p)
  }

  # ----------------------------------------------------------------
  # Short Read (SJ) Validation UMAPs
  # ----------------------------------------------------------------
  gg_sr_umap_plots <<- list()
  if (exists("merged_umap") && "srjunctions_support_prop" %in% colnames(merged_umap)) {
    # Global
    gg_sr_umap_plots[["All Transcripts"]] <<- build_continuous_umap(
      merged_umap,
      "srjunctions_support_prop",
      "All Transcripts Junction Coverage by Short Reads",
      color_base = "#cd4f39"
    )

    # Per-Category
    for (cat_col in names(cat_labels)) {
      tag <- cat_labels[[cat_col]]
      tag_clean <- gsub(" ", "_", tag)
      sr_col <- paste0(tag_clean, "_srjunctions_support_prop")

      if (sr_col %in% colnames(merged_umap)) {
        gg_sr_umap_plots[[tag]] <<- build_continuous_umap(
          merged_umap,
          sr_col,
          paste(tag, "Junction Coverage by Short Reads"),
          color_base = cat_colors[[cat_col]]
        )
      }
    }
  }

  # ----------------------------------------------------------------
  # TSS Validation UMAPs
  # ----------------------------------------------------------------
  gg_tss_umap_plots <<- list()
  if (exists("merged_umap") && "TSS_ratio_validated_prop" %in% colnames(merged_umap)) {
    # Global
    gg_tss_umap_plots[["All Transcripts"]] <<- build_continuous_umap(
      merged_umap,
      "TSS_ratio_validated_prop",
      "All Transcripts TSS Validation by Short Reads",
      color_base = "#ffc125"
    )

    # Per-Category
    for (cat_col in names(cat_labels)) {
      tag <- cat_labels[[cat_col]]
      tag_clean <- gsub(" ", "_", tag)
      # Handle special cases if any (e.g. Genic Genomic)
      simple_tag <- names(cat_labels)[which(cat_labels == tag)]
      simple_tag <- gsub("_prop", "", simple_tag)
      tss_col <- paste0(simple_tag, "_TSS_ratio_validated_prop")

      if (tss_col %in% colnames(merged_umap)) {
        gg_tss_umap_plots[[tag]] <<- build_continuous_umap(
          merged_umap,
          tss_col,
          paste(tag, "TSS Validation by Short Reads"),
          color_base = cat_colors[[cat_col]]
        )
      }
    }
  }

  # Helper: grouped violins by bin with legend (Annotated/Novel) using ggplot, rebuildable for PDF
  # df must contain columns: bin, group, value
  build_grouped_violin_plot <- function(df,
                                        bin_levels,
                                        group_levels,
                                        title,
                                        fill_map,
                                        legend_labels,
                                        y_label = "Genes, %",
                                        ylim = c(0, 100),
                                        violin_alpha = 0.5,
                                        box_alpha = 0.3,
                                        box_width = 0.05,
                                        x_tickangle = 45,
                                        violin_width = 0.45,
                                        dodge_width = 0.8,
                                        violangap = 0.05,
                                        violingroupgap = 0.15,
                                        legend_title = NULL) {
    # Ensure factors
    df$bin <- factor(df$bin, levels = bin_levels)
    df$group <- factor(df$group, levels = group_levels)

    # Clamp values to ylim range
    df$value <- pmin(pmax(df$value, ylim[1]), ylim[2])

    p <- ggplot(df, aes(x = bin, y = value, fill = group)) +
      geom_violin(aes(color = group),
        alpha = violin_alpha,
        position = position_dodge(width = dodge_width),
        scale = "width",
        trim = TRUE,
        adjust = 1,
        linewidth = 0.3,
        show.legend = TRUE
      ) +
      geom_boxplot(aes(group = interaction(bin, group)),
        width = box_width,
        outlier.shape = NA,
        fill = NA,
        color = "grey20",
        alpha = box_alpha,
        position = position_dodge(width = dodge_width),
        show.legend = FALSE
      ) +
      stat_summary(aes(group = group),
        fun = mean, geom = "point", shape = 4, size = 1,
        colour = "red", stroke = 0.9,
        position = position_dodge(width = dodge_width),
        show.legend = FALSE
      ) +
      scale_fill_manual(values = fill_map, labels = legend_labels) +
      scale_color_manual(values = fill_map, guide = "none") +
      labs(title = title, x = "", y = y_label) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, 
          angle = if (!is.null(x_tickangle) && x_tickangle != 0) x_tickangle else 0,
          hjust = ifelse(!is.null(x_tickangle) && x_tickangle == 0, 0.5, 1)
        ),
        legend.position = "bottom",
        legend.title = element_blank()
      )

    if (!is.null(ylim)) {
      p <- p + coord_cartesian(ylim = ylim)
    }
    return(p)
  }
  `%||%` <- function(x, y) if (is.null(x)) y else x
  assign_plot <- function(name, plot) assign(name, plot, envir = .GlobalEnv)
  build_violin_from_long <- function(df_long, args) {
    do.call(build_violin_plot, c(list(df_long = df_long), args))
  }
  single_violin <- function(df, cfg) {
    var <- cfg$column
    df_long <- data.frame(Variable = factor(var, levels = var), Value = df[[var]])
    fill_map <- setNames(cfg$fill, var)
    base_args <- list(
      title = cfg$title,
      x_labels = cfg$x_labels %||% cfg$x_label,
      fill_map = fill_map,
      legend = cfg$legend %||% FALSE
    )
    if (!is.null(cfg$y_label)) base_args$y_label <- cfg$y_label
    plot_args <- c(base_args, cfg$plot_args %||% list())
    assign_plot(cfg$name, build_violin_from_long(df_long, plot_args))
  }
  pivot_violin <- function(df, cfg) {
    df_long <- pivot_long(df, cfg$columns)
    fill_map <- cfg$fill_map %||% setNames(rep(cfg$fill, length(cfg$columns)), cfg$columns)
    base_args <- list(
      title = cfg$title,
      x_labels = cfg$x_labels,
      fill_map = fill_map,
      legend = cfg$legend %||% FALSE
    )
    if (!is.null(cfg$y_label)) base_args$y_label <- cfg$y_label
    plot_args <- c(base_args, cfg$plot_args %||% list())
    assign_plot(cfg$name, build_violin_from_long(df_long, plot_args))
  }
  render_pdf_plot <- function(name) {
    if (exists(name)) print(get(name))
  }

  render_pdf_plot_centered <- function(name, width_frac = 0.45) {
    if (!exists(name)) {
      return(invisible(NULL))
    }
    p <- get(name)
    g <- if (inherits(p, "grob")) p else ggplotGrob(p)
    left_right <- (1 - width_frac) / 2
    grid.arrange(nullGrob(), g, nullGrob(), widths = c(left_right, width_frac, left_right), newpage = TRUE)
  }

  # Helper: build length-distribution violins for given column prefix
  # If mono=TRUE, uses *_length_mono_prop columns; otherwise *_length_prop
  build_len_violin_for_prefix <- function(df, prefix, title, fill_color, box_fill = NULL, mono = FALSE, box_outline_color = "grey20", violin_alpha = 0.5, box_alpha = 0.3, violin_outline_fill = FALSE, format = "ggplot") {
    suffix <- if (mono) "_length_mono_prop" else "_length_prop"
    cols <- paste0(prefix, c("_250b", "_500b", "_short", "_mid", "_long"), suffix)
    cols <- cols[cols %in% colnames(df)]
    if (length(cols) == 0) {
      return(NULL)
    }

    df_long <- pivot_long(df, cols)
    fill_map <- setNames(rep(fill_color, length(cols)), cols)
    color_map <- setNames(rep(fill_color, length(cols)), cols)
    x_labels <- c("0-250bp", "250-500bp", "500-1000bp", "1000-2000bp", ">2000bp")
    names(x_labels) <- cols

    build_violin_plot(
      df_long = df_long,
      title = title,
      x_labels = x_labels,
      fill_map = fill_map,
      color_map = color_map,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      legend = FALSE,
      violin_alpha = violin_alpha,
      box_alpha = box_alpha,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = violin_outline_fill,
      box_outline_default = box_outline_color
    )
  }

  # Meta-transcript body coverage profile plot
  # Shows where reads cover their reference transcripts (0%=5' to 100%=3'),
  # one line per reference transcript length bin.
  build_meta_coverage_plot <- function(cls_df, n_bins = 100) {
    required_cols <- c("diff_to_TSS", "diff_to_TTS", "ref_length", "length")
    if (!all(required_cols %in% colnames(cls_df))) {
      return(NULL)
    }

    keep_cols <- c(
      required_cols,
      if ("count" %in% colnames(cls_df)) "count"
    )
    df <- cls_df[, keep_cols, drop = FALSE]
    for (col in required_cols) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    if (!("count" %in% colnames(df))) df$count <- 1L
    df$count <- suppressWarnings(as.numeric(df$count))
    df$count[is.na(df$count) | df$count <= 0] <- 1
    df <- df[complete.cases(df[, required_cols]) & df$ref_length > 0, ]
    if (nrow(df) == 0) {
      return(NULL)
    }

    # Compute coverage fractions (fraction of ref body covered by each read).
    # diff_to_TSS/diff_to_TTS are now transcript-level (spliced) distances from
    # SQANTI3, so the same formula applies to both FSM and ISM:
    #   start_frac = -diff_to_TSS / ref_length   (positive diff = upstream overshoot)
    #   end_frac   = 1 + diff_to_TTS / ref_length (negative diff = early termination)
    df$start_frac <- pmax(0, pmin(1, -df$diff_to_TSS / df$ref_length))
    df$end_frac   <- pmax(0, pmin(1, 1 + df$diff_to_TTS / df$ref_length))

    # Safety: keep only rows where end > start
    df <- df[df$end_frac > df$start_frac, ]
    if (nrow(df) == 0) {
      return(NULL)
    }

    # Bin by reference transcript length
    df$len_bin <- cut(df$ref_length,
      breaks = c(0, 250, 500, 1000, 2000, 5000, Inf),
      labels = c("0-250bp", "250-500bp", "500bp-1kb", "1-2kb", "2-5kb", ">5kb"),
      right = FALSE
    )

    # For each bin, compute weighted % of reads covering each position
    bin_positions <- seq(0, 1, length.out = n_bins + 1)
    bin_mids <- (bin_positions[-length(bin_positions)] + bin_positions[-1]) / 2

    profile_list <- lapply(levels(df$len_bin), function(lb) {
      sub <- df[df$len_bin == lb, ]
      if (nrow(sub) == 0) {
        return(NULL)
      }
      total_weight <- sum(sub$count)
      coverage <- vapply(seq_len(n_bins), function(i) {
        pos <- bin_mids[i]
        sum(sub$count[sub$start_frac <= pos & sub$end_frac >= pos]) / total_weight * 100
      }, numeric(1))
      data.frame(
        position = bin_mids * 100,
        coverage = coverage,
        len_bin = lb,
        n_reads = total_weight,
        stringsAsFactors = FALSE
      )
    })
    profile_df <- do.call(rbind, Filter(Negate(is.null), profile_list))
    if (is.null(profile_df) || nrow(profile_df) == 0) {
      return(NULL)
    }

    profile_df$len_bin <- factor(profile_df$len_bin, levels = levels(df$len_bin))

    # RColorConesa main palette (7 discrete base colors)
    n_levels <- nlevels(df$len_bin)
    bin_colors <- c("#00B0A5", "#E1744E", "#FAC24A", "#6DC8E5", "#E7A5CB", "#9C8AB4", "#E44067")[seq_len(n_levels)]

    p <- ggplot(profile_df, aes(x = position, y = coverage, color = len_bin)) +
      geom_line(linewidth = 0.9, alpha = 0.85) +
      scale_color_manual(values = setNames(bin_colors, levels(profile_df$len_bin))) +
      scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(
        title = paste(entity_label, "Body Coverage Along Reference Transcript"),
        x = "Position along reference transcript (%)\n5' -> 3'",
        y = paste(entity_label_plural, "covering position, %"),
        color = "Reference length"
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 14)
      ) +
      guides(color = guide_legend(nrow = 2, byrow = TRUE))
    return(p)
  }

  # Build GTF reference length violin/boxplot comparison (works for both 'isoforms' and 'reads' modes)
  build_ref_vs_sample_lengths <- function(cls_df, ref_gtf, mode = "isoforms") {
    if (is.null(ref_gtf)) {
      return(NULL)
    }

    # Strip quotes if passed by the shell
    ref_gtf <- gsub('^"|"$', "", ref_gtf)

    if (!file.exists(ref_gtf)) {
      return(NULL)
    }

    # 1. Sample Transcripts
    if (!"length" %in% colnames(cls_df)) {
      return(NULL)
    }
    sample_lengths <- suppressWarnings(as.numeric(cls_df$length))
    sample_lengths <- sample_lengths[!is.na(sample_lengths) & sample_lengths > 0]

    if (length(sample_lengths) == 0) {
      return(NULL)
    }

    sample_label <- if (mode == "isoforms") "Sample Transcriptome" else "Sample Reads"
    sample_df <- data.frame(
      length = sample_lengths,
      Dataset = sample_label,
      stringsAsFactors = FALSE
    )

    # 2. Reference Transcripts via data.table rapid GTF parsing
    # Skip comment lines explicitly to prevent fread from guessing incorrect column number
    skip_lines <- 0
    con <- try(file(ref_gtf, "r"), silent = TRUE)
    if (!inherits(con, "try-error")) {
      while (TRUE) {
        line <- suppressWarnings(readLines(con, n = 1))
        if (length(line) == 0) break
        if (startsWith(line, "#")) {
          skip_lines <- skip_lines + 1
        } else {
          break
        }
      }
      close(con)
    }

    gtf <- tryCatch(
      data.table::fread(ref_gtf, sep = "\t", skip = skip_lines, header = FALSE, fill = TRUE, quote = ""),
      error = function(e) NULL
    )
    if (is.null(gtf) || nrow(gtf) == 0 || ncol(gtf) < 9) {
      return(NULL)
    }

    # V3 = feature, V4 = start, V5 = end, V9 = attributes
    exons <- gtf[V3 == "exon"]
    if (nrow(exons) == 0) {
      return(NULL)
    }

    exons[, start_pos := suppressWarnings(as.numeric(V4))]
    exons[, end_pos := suppressWarnings(as.numeric(V5))]
    exons <- exons[!is.na(start_pos) & !is.na(end_pos)]
    exons[, exon_length := end_pos - start_pos + 1]

    # Extract transcript_id via regex
    exons[, transcript_id := sub(".*transcript_id \"([^\"]+)\".*", "\\1", V9)]

    # Aggregate by transcript
    ref_lengths_dt <- exons[, .(length = sum(exon_length, na.rm = TRUE)), by = transcript_id]
    ref_lengths <- ref_lengths_dt$length[ref_lengths_dt$length > 0]

    if (length(ref_lengths) == 0) {
      return(NULL)
    }

    ref_df <- data.frame(
      length = ref_lengths,
      Dataset = "Reference Transcriptome",
      stringsAsFactors = FALSE
    )

    # 3. Combine and Plot
    plot_df <- rbind(sample_df, ref_df)
    dataset_levels <- c("Reference Transcriptome", sample_label)
    plot_df$Dataset <- factor(plot_df$Dataset, levels = dataset_levels)
    pal <- setNames(c("#1fa291", "#f5c05d"), dataset_levels)
    y_axis_label <- "Length, bp"
    plot_subtitle <- if (mode == "isoforms") "Reference transcriptome vs. Sample transcriptome" else "Reference transcriptome vs. Sample reads"

    p <- ggplot(plot_df, aes(x = Dataset, y = length, fill = Dataset)) +
      geom_violin(aes(color = Dataset), alpha = 0.7, scale = "width", adjust = 1, trim = TRUE, show.legend = FALSE) +
      scale_color_manual(values = pal, guide = "none") +
      geom_boxplot(width = 0.05, alpha = 0.6, outlier.shape = NA, color = "grey20", show.legend = FALSE) +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1, show.legend = FALSE) +
      scale_y_log10(labels = scales::comma) +
      scale_fill_manual(values = pal) +
      labs(
        title = paste0("Transcript Length Distribution:\n", plot_subtitle),
        x = "",
        y = y_axis_label
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1)
      )

    return(p)
  }

  # Helper: build a per-category exon count profile (median + IQR across cells)
  build_exon_profile_plot <- function(df_prof, title, line_color, k_max = 20, y_label = paste(entity_label_plural, ", %", sep = ""), n_cells = NULL) {
    # Sanitize title
    title <- tryCatch(
      {
        t <- as.character(title)
        trimws(gsub("\n.*", "", t))
      },
      error = function(e) "Exon Profile"
    )

    # Canonical Fusion color override
    FUSION_COLOR <- "#F1C40F"
    detect_fusion <- function(df) {
      tryCatch(
        {
          ttl_has <- grepl("Fusion", title, ignore.case = TRUE)
          in_df <- "category" %in% colnames(df) && any(grepl("Fusion", df$category, ignore.case = TRUE))
          ttl_has || in_df
        },
        error = function(e) FALSE
      )
    }
    if (detect_fusion(df_prof)) {
      line_color <- FUSION_COLOR
    }

    # Helpers: color utilities
    lighten_hex <- function(hex, amount = 0.35) {
      rgb <- grDevices::col2rgb(hex)
      r <- as.integer(round(rgb[1] + (255 - rgb[1]) * amount))
      g <- as.integer(round(rgb[2] + (255 - rgb[2]) * amount))
      b <- as.integer(round(rgb[3] + (255 - rgb[3]) * amount))
      grDevices::rgb(r, g, b, maxColorValue = 255)
    }

    if (is.null(df_prof) || nrow(df_prof) == 0 || all(!is.finite(df_prof$median))) {
      p_empty <- ggplot() +
        labs(title = title, subtitle = "No data available for this category") +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray")
        )
      return(p_empty)
    }

    tick_breaks <- seq_len(k_max)
    label_last <- paste0("\u2265", k_max)
    ticktexts <- c(as.character(seq_len(k_max - 1)), label_last)

    # Choose center line: use mean if present, else median
    y_center <- if (!is.null(df_prof$mean)) df_prof$mean else df_prof$median

    # Build stat columns for the line/ribbon
    stat_cols <- intersect(colnames(df_prof), c("mean", "median"))
    line_stats <- if (length(stat_cols)) {
      df_prof %>%
        dplyr::select(k, dplyr::all_of(stat_cols)) %>%
        tidyr::pivot_longer(cols = dplyr::all_of(stat_cols), names_to = "stat", values_to = "value") %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::mutate(stat = factor(stat, levels = stat_cols))
    } else {
      NULL
    }

    iqr_fill <- lighten_hex(line_color, 0.55)

    p <- ggplot(df_prof, aes(x = k)) +
      geom_ribbon(aes(ymin = q1, ymax = q3, fill = "IQR"), alpha = 0.25, show.legend = TRUE, key_glyph = "rect") +
      scale_y_continuous(limits = c(0, 100)) +
      scale_x_continuous(breaks = tick_breaks, labels = ticktexts) +
      scale_fill_manual(values = c("IQR" = iqr_fill), name = "") +
      labs(title = title, x = "Number of exons", y = y_label) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "bottom"
      )

    if (!is.null(line_stats) && nrow(line_stats) > 0) {
      line_palette <- setNames(rep(line_color, length(stat_cols)), stat_cols)
      linetype_values <- c("mean" = "solid", "median" = "dashed")
      legend_linewidths <- c("mean" = 1.2, "median" = 1.0)
      if (!"mean" %in% stat_cols) {
        linetype_values <- linetype_values[names(linetype_values) != "mean"]
        legend_linewidths <- legend_linewidths[names(legend_linewidths) != "mean"]
      }
      if (!"median" %in% stat_cols) {
        linetype_values <- linetype_values[names(linetype_values) != "median"]
        legend_linewidths <- legend_linewidths[names(legend_linewidths) != "median"]
      }
      p <- p +
        geom_line(data = line_stats, aes(x = k, y = value, linetype = stat, color = stat, linewidth = stat)) +
        scale_color_manual(values = line_palette, name = "") +
        scale_linetype_manual(values = linetype_values, name = "") +
        scale_linewidth_manual(values = legend_linewidths, name = "")
    }

    return(p)
  }

  # Basic cell information

  common_plot_args <- list(
    violin_alpha = 0.5,
    box_alpha = 0.3,
    box_width = 0.05,
    violin_outline_fill = FALSE,
    box_outline_default = "black",
    log_scale = TRUE,
    adjust = 1
  )

  # 1. Number of Reads Across Cells
  cfg_reads <- list(
    column = count_col,
    name = "gg_reads_in_cells",
    title = paste0("Number of ", entity_label_plural, "\nAcross Cells"),
    fill = "#CC6633",
    y_label = paste(entity_label_plural, ", count", sep = ""),
    x_label = "Cells",
    plot_args = common_plot_args
  )
  single_violin(SQANTI_cell_summary, cfg_reads)

  # 2. Number of UMIs Across Cells (only if not isoforms mode)
  if (mode != "isoforms") {
    cfg_umis <- list(
      column = "UMIs_in_cell",
      name = "gg_umis_in_cells",
      title = "Number of UMIs\nAcross Cells",
      fill = "#CC6633",
      y_label = "UMIs, count",
      x_label = "Cells",
      plot_args = common_plot_args
    )
    single_violin(SQANTI_cell_summary, cfg_umis)
  }

  # 3. Number of Genes Across Cells
  cfg_genes <- list(
    column = "Genes_in_cell",
    name = "gg_genes_in_cells",
    title = "Number of Genes\nAcross Cells",
    fill = "#CC6633",
    y_label = "Genes, count",
    x_label = "Cells",
    plot_args = common_plot_args
  )
  single_violin(SQANTI_cell_summary, cfg_genes)

  # 4. Number of Unique Junction Chains Across Cells
  if (mode != "isoforms" && "UJCs_in_cell" %in% names(SQANTI_cell_summary) && !all(is.na(SQANTI_cell_summary$UJCs_in_cell)) && max(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE) > 0) {
    cfg_ujcs <- list(
      column = "UJCs_in_cell",
      name = "gg_JCs_in_cell",
      title = "Number of Unique Junction\nChains Across Cells",
      fill = "#CC6633",
      y_label = "UJCs, count",
      x_label = "Cells",
      plot_args = common_plot_args
    )
    single_violin(SQANTI_cell_summary, cfg_ujcs)
  }

  pivot_defaults <- list(
    violin_alpha = 0.5,
    box_alpha = 0.3,
    box_width = 0.05,
    violin_outline_fill = FALSE,
    box_outline_default = "black"
  )
  pivot_violin(SQANTI_cell_summary, list(
    name = "gg_annotation_of_genes_in_cell",
    columns = c("Annotated_genes", "Novel_genes"),
    title = "Number of Known/Novel Genes\nAcross Cells",
    x_labels = c("Annotated Genes", "Novel Genes"),
    y_label = paste(entity_label_plural, ", count", sep = ""),
    fill_map = c("Annotated_genes" = fill_color_orange, "Novel_genes" = fill_color_orange),
    plot_args = c(pivot_defaults, list(log_scale = TRUE))
  ))

  if ("Genes_in_cell" %in% colnames(SQANTI_cell_summary)) {
    SQANTI_cell_summary$Annotated_genes_perc <- ifelse(
      SQANTI_cell_summary$Genes_in_cell > 0,
      100 * SQANTI_cell_summary$Annotated_genes / SQANTI_cell_summary$Genes_in_cell,
      0
    )
    SQANTI_cell_summary$Novel_genes_perc <- ifelse(
      SQANTI_cell_summary$Genes_in_cell > 0,
      100 * SQANTI_cell_summary$Novel_genes / SQANTI_cell_summary$Genes_in_cell,
      0
    )

    pivot_violin(SQANTI_cell_summary, list(
      name = "gg_annotation_of_genes_percent_in_cell",
      columns = c("Annotated_genes_perc", "Novel_genes_perc"),
      title = "Percentage of Known/Novel Genes Across Cells",
      x_labels = c("Annotated Genes", "Novel Genes"),
      y_label = "Genes, %",
      fill_map = c("Annotated_genes_perc" = fill_color_orange, "Novel_genes_perc" = fill_color_orange),
      plot_args = pivot_defaults
    ))
  }

  # 5. Percentage of Reads/Transcripts from Known/Novel Genes Across Cells
  # (Enabled for both reads and isoforms modes)
  {
    classification_valid <- Classification_file[Classification_file$CB != "unassigned" & !is.na(Classification_file$CB), ]

    if (nrow(classification_valid) > 0) {
      # Function to expand FL and CB columns into a long format for correct counting per cell
      expand_isoform_counts <- function(df, mode) {
        if (mode == "reads") {
          return(df %>% group_by(CB) %>% summarise(count = n(), .groups = "drop"))
        } else {
          # Isoforms mode: Each row has comma-separated FL (counts) and CB (barcodes)
          # We need to split them and sum counts per barcode

          # Initialize lists to store expanded data
          all_cbs <- character()
          all_counts <- numeric()

          # Iterate through rows (this might be slow for huge files, but safe)
          # A vectorised approach would be better if possible, but strsplit returns list
          fl_list <- strsplit(as.character(df$FL), ",")
          cb_list <- strsplit(as.character(df$CB), ",")

          # Check if lengths match (they should)
          if (length(fl_list) != length(cb_list)) {
            stop("Mismatch in row counts between FL and CB columns")
          }

          # Use mapply to create a data frame of all counts
          # This creates a list of data frames, one per isoform
          expanded_list <- mapply(function(fl, cb) {
            if (length(fl) != length(cb)) {
              # Warning or skip? For now, we assume they match as per SQANTI specs
              return(NULL)
            }
            data.frame(CB = cb, count = as.numeric(fl), stringsAsFactors = FALSE)
          }, fl_list, cb_list, SIMPLIFY = FALSE)

          # Bind all tiny data frames
          long_df <- do.call(rbind, expanded_list)

          # Now group by CB and sum
          return(long_df %>% group_by(CB) %>% summarise(count = sum(count, na.rm = TRUE), .groups = "drop"))
        }
      }

      annotated_reads_per_cell <- classification_valid %>%
        filter(!grepl("^novel", associated_gene))

      annotated_reads_per_cell <- expand_isoform_counts(annotated_reads_per_cell, mode) %>%
        rename(Annotated_genes_reads = count)

      novel_reads_per_cell <- classification_valid %>%
        filter(grepl("^novel", associated_gene))

      novel_reads_per_cell <- expand_isoform_counts(novel_reads_per_cell, mode) %>%
        rename(Novel_genes_reads = count)

      SQANTI_cell_summary <- SQANTI_cell_summary %>%
        left_join(annotated_reads_per_cell, by = "CB") %>%
        left_join(novel_reads_per_cell, by = "CB")

      SQANTI_cell_summary$Annotated_genes_reads[is.na(SQANTI_cell_summary$Annotated_genes_reads)] <- 0
      SQANTI_cell_summary$Novel_genes_reads[is.na(SQANTI_cell_summary$Novel_genes_reads)] <- 0

      # Revert to original denominator (Total Transcripts in Cell) now that numerators are correct
      SQANTI_cell_summary$Annotated_reads_perc <- 100 * SQANTI_cell_summary$Annotated_genes_reads / SQANTI_cell_summary[[count_col]]
      SQANTI_cell_summary$Novel_reads_perc <- 100 * SQANTI_cell_summary$Novel_genes_reads / SQANTI_cell_summary[[count_col]]

      SQANTI_cell_summary$Annotated_reads_perc <- ifelse(is.na(SQANTI_cell_summary$Annotated_reads_perc) | is.infinite(SQANTI_cell_summary$Annotated_reads_perc), 0, SQANTI_cell_summary$Annotated_reads_perc)
      SQANTI_cell_summary$Novel_reads_perc <- ifelse(is.na(SQANTI_cell_summary$Novel_reads_perc) | is.infinite(SQANTI_cell_summary$Novel_reads_perc), 0, SQANTI_cell_summary$Novel_reads_perc)

      pivot_violin(SQANTI_cell_summary, list(
        name = "gg_annotation_of_reads_in_cell",
        columns = c("Annotated_reads_perc", "Novel_reads_perc"),
        title = paste("Percentage of", entity_label_plural, "from Known/Novel Genes Across Cells"),
        x_labels = c("Annotated Genes", "Novel Genes"),
        y_label = paste(entity_label_plural, ", %", sep = ""),
        fill_map = c("Annotated_reads_perc" = fill_color_orange, "Novel_reads_perc" = fill_color_orange),
        plot_args = pivot_defaults
      ))
    } else {
      message("Warning: No valid classification data found. Skipping read expression by gene annotation plot.")
      gg_annotation_of_reads_in_cell <<- ggplot() +
        labs(title = "Plot not available") +
        theme_minimal()
      layout(
        title = paste("Percentage of", entity_label_plural, "from Known/Novel Genes Across Cells"),
        annotations = list(
          text = paste(entity_label, "expression by gene annotation\nnot available"),
          showarrow = FALSE,
          font = list(size = 16, color = "gray")
        )
      )
    }
  }

  single_violin(SQANTI_cell_summary, list(
    name = "gg_MT_perc",
    column = "MT_perc",
    title = paste("Mitochondrial", entity_label_plural, "Across Cells"),
    x_labels = c("Cell"),
    y_label = paste(entity_label_plural, ", %", sep = ""),
    fill = "#CC6633",
    plot_args = common_plot_args
  ))

  ### Gene Distribution by Read Count Bins (configurable gene bins) ###
  ####################################################################

  # Define gene read-count bins and labels
  gene_bin_label <- function(n) {
    if (is.na(n)) {
      return(NA_character_)
    }
    if (n == 1) {
      return("1")
    }
    if (n >= 2 && n <= 5) {
      return("2-5")
    }
    if (n >= 6 && n <= 9) {
      return("6-9")
    }
    return(">=10")
  }
  gene_bin_levels <- c("1", "2-5", "6-9", ">=10")

  # Build per-cell per-gene read counts from classification.
  # In isoforms mode, must explode the comma-separated CB/FL columns so each
  # (cell, isoform) pair is weighted by its FL count, then sum per (CB, gene).
  # In reads mode, each row is one read so n() is correct.
  if (mode == "isoforms" && "FL" %in% colnames(Classification_file) && "CB" %in% colnames(Classification_file)) {
    genes_by_cb_base <- Classification_file %>%
      filter(!is.na(CB), CB != "unassigned", !is.na(associated_gene)) %>%
      select(CB, FL, associated_gene)

    genes_by_cb_base$CB_raw <- as.character(genes_by_cb_base$CB)
    genes_by_cb_base$FL_raw <- as.character(genes_by_cb_base$FL)
    genes_by_cb_base <- tidyr::separate_rows(genes_by_cb_base, CB_raw, FL_raw, sep = ",")
    genes_by_cb_base$FL_num <- suppressWarnings(as.numeric(trimws(genes_by_cb_base$FL_raw)))
    genes_by_cb_base$FL_num[is.na(genes_by_cb_base$FL_num) | genes_by_cb_base$FL_num < 0] <- 0
    genes_by_cb_base$CB_clean <- trimws(genes_by_cb_base$CB_raw)

    genes_by_cb <- genes_by_cb_base %>%
      filter(CB_clean != "" & CB_clean != "unassigned" & FL_num > 0) %>%
      group_by(CB = CB_clean, associated_gene) %>%
      summarise(reads_per_gene = sum(FL_num), .groups = "drop") %>%
      mutate(
        gene_type = ifelse(grepl("^novel", associated_gene), "Novel", "Annotated"),
        bin = vapply(reads_per_gene, gene_bin_label, character(1))
      ) %>%
      filter(!is.na(bin))
  } else {
    genes_by_cb <- Classification_file %>%
      filter(!is.na(CB), CB != "unassigned", !is.na(associated_gene)) %>%
      group_by(CB, associated_gene) %>%
      summarise(reads_per_gene = n(), .groups = "drop") %>%
      mutate(
        gene_type = ifelse(grepl("^novel", associated_gene), "Novel", "Annotated"),
        bin = vapply(reads_per_gene, gene_bin_label, character(1))
      ) %>%
      filter(!is.na(bin))
  }

  # Percent of genes per bin within each CB and gene type
  read_bins_data <- genes_by_cb %>%
    group_by(CB, gene_type, bin) %>%
    summarise(num_genes = n(), .groups = "drop") %>%
    group_by(CB, gene_type) %>%
    mutate(percentage = 100 * num_genes / sum(num_genes)) %>%
    ungroup() %>%
    tidyr::complete(CB, gene_type, bin = gene_bin_levels, fill = list(num_genes = 0, percentage = 0))

  read_bins_data$bin <- factor(read_bins_data$bin, levels = gene_bin_levels)
  read_bins_data$gene_type <- factor(read_bins_data$gene_type, levels = c("Annotated", "Novel"))

  if (mode == "isoforms") {
    gg_read_bins <<- build_grouped_violin_plot(
      df = read_bins_data %>% transmute(bin = as.character(bin), group = as.character(gene_type), value = percentage),
      bin_levels = gene_bin_levels,
      group_levels = c("Annotated", "Novel"),
      title = paste("Distribution of Known/Novel Genes by", entity_label, "Count Bins Across Cells"),
      fill_map = c("Annotated" = "#e37744", "Novel" = "#78C679"),
      legend_labels = c("Annotated" = "Annotated", "Novel" = "Novel"),
      y_label = "Genes, %",
      ylim = c(0, 100),
      violin_alpha = 0.5,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_width = 0.28,
      dodge_width = 1.0
    )
  }

  # Combined (all genes together): one violin per bin
  if (mode == "reads") {
    # Filter for Annotated genes only
    read_bins_all <- genes_by_cb %>%
      filter(gene_type == "Annotated") %>%
      group_by(CB, bin) %>%
      summarise(num_genes = n(), .groups = "drop") %>%
      group_by(CB) %>%
      mutate(percentage = 100 * num_genes / sum(num_genes)) %>%
      ungroup() %>%
      tidyr::complete(CB, bin = gene_bin_levels, fill = list(num_genes = 0, percentage = 0))

    plot_title_all <- paste("Distribution of Annotated Genes by", entity_label, "Count Bins Across Cells")
  } else {
    # All genes (Annotated + Novel)
    read_bins_all <- genes_by_cb %>%
      group_by(CB, bin) %>%
      summarise(num_genes = n(), .groups = "drop") %>%
      group_by(CB) %>%
      mutate(percentage = 100 * num_genes / sum(num_genes)) %>%
      ungroup() %>%
      tidyr::complete(CB, bin = gene_bin_levels, fill = list(num_genes = 0, percentage = 0))

    plot_title_all <- paste("Distribution of Genes by", entity_label, "Count Bins Across Cells")
  }

  read_bins_all$bin <- factor(read_bins_all$bin, levels = gene_bin_levels)

  {
    df_long <- data.frame(
      Variable = factor(read_bins_all$bin, levels = gene_bin_levels),
      Value = read_bins_all$percentage
    )
    fill_map <- setNames(rep("#CC6633", length(gene_bin_levels)), gene_bin_levels)
    gg_read_bins_all <<- build_violin_plot(
      df_long,
      title = plot_title_all,
      x_labels = as.character(gene_bin_levels),
      fill_map = fill_map,
      y_label = "Genes, %",
      legend = FALSE,
      ylim = c(0, 100),
      violin_alpha = 0.5,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      box_outline_default = "black",
      violin_outline_fill = FALSE
    )
  }

  # New plot: Distribution of Known Genes by Unique Isoform Count Bins Across Cells (Isoforms mode)
  if (mode == "isoforms") {
    iso_bins_annot <- genes_by_cb %>%
      filter(gene_type == "Annotated") %>%
      group_by(CB, bin) %>%
      summarise(num_genes = n(), .groups = "drop") %>%
      group_by(CB) %>%
      mutate(percentage = 100 * num_genes / sum(num_genes)) %>%
      ungroup() %>%
      tidyr::complete(CB, bin = gene_bin_levels, fill = list(num_genes = 0, percentage = 0))

    iso_bins_annot$bin <- factor(iso_bins_annot$bin, levels = gene_bin_levels)

    df_long_iso <- data.frame(
      Variable = factor(iso_bins_annot$bin, levels = gene_bin_levels),
      Value = iso_bins_annot$percentage
    )
    # Use Annotated color #e37744
    fill_map_iso <- setNames(rep("#e37744", length(gene_bin_levels)), gene_bin_levels)

    gg_isoform_bins <<- build_violin_plot(
      df_long_iso,
      title = "Distribution of Known Genes by Unique Isoform Count Bins Across Cells",
      x_labels = as.character(gene_bin_levels),
      fill_map = fill_map_iso,
      y_label = "Genes, %",
      legend = FALSE,
      ylim = c(0, 100),
      violin_alpha = 0.5,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      box_outline_default = "black",
      violin_outline_fill = FALSE
    )
  }

  # UJC bins (combined) using jxn strings per gene per CB
  if (mode != "isoforms") {
    ujc_bin_label <- function(n) {
      if (is.na(n)) {
        return(NA_character_)
      }
      if (n == 1) {
        return("1")
      }
      if (n >= 2 && n <= 3) {
        return("2-3")
      }
      if (n >= 4 && n <= 5) {
        return("4-5")
      }
      return(">=6")
    }
    ujc_bin_levels <- c("1", "2-3", "4-5", ">=6")

    # Check if jxn_string exists (it won't if --skip_hash was used)
    if ("jxn_string" %in% colnames(Classification_file)) {
      ujc_by_cb <- Classification_file %>%
        filter(!is.na(CB), CB != "unassigned", !is.na(associated_gene), exons > 1) %>%
        group_by(CB, associated_gene) %>%
        summarise(ujc_per_gene = dplyr::n_distinct(jxn_string), .groups = "drop") %>%
        mutate(bin = vapply(ujc_per_gene, ujc_bin_label, character(1))) %>%
        filter(!is.na(bin))
    } else {
      ujc_by_cb <- data.frame()
    }

    if (nrow(ujc_by_cb) > 0) {
      # For reads mode, filter for Annotated genes only
      if (mode == "reads") {
        ujc_by_cb <- ujc_by_cb %>%
          mutate(gene_type = ifelse(grepl("^novel", associated_gene), "Novel", "Annotated")) %>%
          filter(gene_type == "Annotated")

        plot_title_ujc_all <- "Distribution of Annotated Genes by UJC Count Bins Across Cells"
      } else {
        plot_title_ujc_all <- "Distribution of Genes by UJC Count Bins Across Cells"
      }

      ujc_bins_all <- ujc_by_cb %>%
        group_by(CB, bin) %>%
        summarise(num_genes = n(), .groups = "drop") %>%
        group_by(CB) %>%
        mutate(percentage = 100 * num_genes / sum(num_genes)) %>%
        ungroup() %>%
        tidyr::complete(CB, bin = ujc_bin_levels, fill = list(num_genes = 0, percentage = 0))

      ujc_bins_all$bin <- factor(ujc_bins_all$bin, levels = ujc_bin_levels)

      {
        df_long <- data.frame(
          Variable = factor(ujc_bins_all$bin, levels = ujc_bin_levels),
          Value = ujc_bins_all$percentage
        )
      }
      fill_map <- setNames(rep("#CC6633", length(ujc_bin_levels)), ujc_bin_levels)
      gg_ujc_bins_all <<- build_violin_plot(
        df_long,
        title = plot_title_ujc_all,
        x_labels = as.character(ujc_bin_levels),
        fill_map = fill_map,
        y_label = "Genes, %",
        legend = FALSE,
        ylim = c(0, 100),
        violin_alpha = 0.5,
        box_alpha = 0.3,
        box_width = 0.05,
        x_tickangle = 45,
        box_outline_default = "black",
        violin_outline_fill = FALSE
      )
    }

    # Create UJC bins data
    ujc_bins_data <- data.frame(
      CB = rep(SQANTI_cell_summary$CB, 8),
      bin = rep(c("1", "2-3", "4-5", ">=6", "1", "2-3", "4-5", ">=6"), each = nrow(SQANTI_cell_summary)),
      gene_type = rep(c("Annotated", "Annotated", "Annotated", "Annotated", "Novel", "Novel", "Novel", "Novel"), each = nrow(SQANTI_cell_summary)),
      percentage = c(
        SQANTI_cell_summary$anno_ujc_bin1_perc,
        SQANTI_cell_summary$anno_ujc_bin2_3_perc,
        SQANTI_cell_summary$anno_ujc_bin4_5_perc,
        SQANTI_cell_summary$anno_ujc_bin6plus_perc,
        SQANTI_cell_summary$novel_ujc_bin1_perc,
        SQANTI_cell_summary$novel_ujc_bin2_3_perc,
        SQANTI_cell_summary$novel_ujc_bin4_5_perc,
        SQANTI_cell_summary$novel_ujc_bin6plus_perc
      )
    )

    # Handle NA and invalid values
    ujc_bins_data <- ujc_bins_data %>%
      mutate(percentage = ifelse(is.na(percentage) | is.infinite(percentage) | percentage < 0, 0, percentage))

    ujc_bins_data$bin <- factor(ujc_bins_data$bin, levels = c("1", "2-3", "4-5", ">=6"))
    ujc_bins_data$gene_type <- factor(ujc_bins_data$gene_type, levels = c("Annotated", "Novel"))

    # Only generate split plot if NOT in reads mode
    if (mode != "reads") {
      gg_ujc_bins <<- build_grouped_violin_plot(
        df = ujc_bins_data %>% transmute(bin = as.character(bin), group = as.character(gene_type), value = percentage),
        bin_levels = ujc_bin_levels,
        group_levels = c("Annotated", "Novel"),
        title = "Distribution of Known/Novel Genes by UJC Count Bins Across Cells",
        fill_map = c("Annotated" = "#e37744", "Novel" = "#78C679"),
        legend_labels = c("Annotated" = "Annotated", "Novel" = "Novel"),
        y_label = "Genes, %",
        ylim = c(0, 100),
        violin_alpha = 0.5,
        box_alpha = 0.3,
        box_width = 0.05,
        x_tickangle = 45,
        violin_width = 0.28,
        dodge_width = 1.0
      )
    }
  }

  # Mitochondrial percentage in cell
  {
    df_long <- data.frame(Variable = "MT_perc", Value = SQANTI_cell_summary$MT_perc)
    df_long$Variable <- factor(df_long$Variable, levels = "MT_perc")
    fill_map <- c("MT_perc" = fill_color_orange)
    x_labels <- c("Cell")
    gg_MT_perc <<- build_violin_plot(
      df_long,
      title = paste("Mitochondrial", entity_label_plural, "\nAcross Cells"),
      x_labels = x_labels,
      fill_map = fill_map,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      legend = FALSE,
      violin_alpha = 0.5,
      box_alpha = 0.3,
      box_width = 0.05,
      box_outline_default = "black",
      violin_outline_fill = FALSE,
      x_tickangle = 45
    )
  }

  #  Mono/multi-exon prop novel vs annotated genes

  ### Length distribution ###
  ###########################

  # Compact helpers for repeated per-category and length plots
  cat_tags <- c("FSM", "ISM", "NIC", "NNC", "Genic", "Antisense", "Fusion", "Intergenic", "Genic_intron")
  cat_labels_pretty <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic", "Antisense", "Fusion", "Intergenic", "Genic\nIntron")
  cat_fill_map <- c(FSM = "#6BAED6", ISM = "#FC8D59", NIC = "#78C679", NNC = "#EE6A50", Genic = "#969696", Antisense = "#66C2A4", Fusion = "goldenrod1", Intergenic = "darksalmon", Genic_intron = "#41B6C4")
  structural_category_map <- c(
    "full-splice_match" = "FSM",
    "incomplete-splice_match" = "ISM",
    "novel_in_catalog" = "NIC",
    "novel_not_in_catalog" = "NNC",
    "genic" = "Genic",
    "antisense" = "Antisense",
    "fusion" = "Fusion",
    "intergenic" = "Intergenic",
    "genic_intron" = "Genic_intron"
  )
  structural_category_levels <- unname(structural_category_map)
  # Build violin across categories and assign to a global name
  # Helper: build 9 tag column names from suffix (e.g. "_intrapriming_prop")
  cat_cols <- function(suffix) paste0(cat_tags, suffix)
  # Length plot generator and variable name mapping
  cat_var_base <- c(FSM = "FSM", ISM = "ISM", NIC = "NIC", NNC = "NNC", Genic = "genic", Antisense = "antisense", Fusion = "fusion", Intergenic = "intergenic", Genic_intron = "genic_intron")
  make_len_plot <- function(prefix, pretty, color, mono = FALSE) {
    var_nm <- if (mono) paste0("gg_", cat_var_base[[prefix]], "_mono_read_distr") else paste0("gg_", cat_var_base[[prefix]], "_read_distr")
    title_txt <- if (mono) paste0(pretty, " Mono-exonic Read Lengths Distribution Across Cells") else paste0(pretty, " Reads Length Distribution Across Cells")
    assign(var_nm, build_len_violin_for_prefix(
      SQANTI_cell_summary,
      prefix = prefix,
      title = title_txt,
      fill_color = color,
      box_fill = color,
      mono = mono,
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_outline_color = if (prefix %in% c("Genic")) "grey90" else "grey20",
      violin_outline_fill = TRUE
    ), envir = .GlobalEnv)
  }

  # Bulk distributions
  bulk_len_breaks <- c(50, 100, 250, 500, 1000, 2000, 5000, 10000, 20000, 50000)

  gg_bulk_all_reads <<- ggplot(Classification_file[Classification_file$length > 0, ], aes(x = length)) +
    geom_histogram(bins = 80, fill = "#CC6633", color = "black", alpha = 0.5) +
    scale_x_log10(breaks = bulk_len_breaks, labels = scales::comma) +
    annotation_logticks(sides = "b") +
    coord_cartesian(clip = "off") +
    labs(
      title = paste("All", entity_label, "Lengths Distribution"),
      x = "Length, bp",
      y = paste(entity_label_plural, ", count", sep = "")
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 40, r = 20, b = 5, l = 5, unit = "pt"),
      axis.title = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16, angle = 35, hjust = 1)
    )

  # Bulk read length distribution by structural category
  Classification_file$structural_category <- factor(
    Classification_file$structural_category,
    levels = c(
      "full-splice_match",
      "incomplete-splice_match",
      "novel_in_catalog",
      "novel_not_in_catalog",
      "genic",
      "antisense",
      "fusion",
      "intergenic",
      "genic_intron"
    )
  )

  structural_category_labels <- c(
    "full-splice_match" = "FSM",
    "incomplete-splice_match" = "ISM",
    "novel_in_catalog" = "NIC",
    "novel_not_in_catalog" = "NNC",
    "genic" = "Genic Genomic",
    "antisense" = "Antisense",
    "fusion" = "Fusion",
    "intergenic" = "Intergenic",
    "genic_intron" = "Genic Intron"
  )
  structural_category_palette <- c(
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
  Classification_file$structural_category_pretty <- structural_category_labels[as.character(Classification_file$structural_category)]
  Classification_file$structural_category_pretty <- factor(
    Classification_file$structural_category_pretty,
    levels = names(structural_category_palette)
  )

  gg_bulk_length_by_category <<- ggplot(Classification_file[Classification_file$length > 0, ], aes(x = length, color = structural_category_pretty)) +
    geom_freqpoly(bins = 80, linewidth = 1.2, na.rm = TRUE) +
    labs(
      title = paste("All", entity_label, "Lengths Distribution by Structural Category"),
      x = "Length, bp",
      y = paste(entity_label_plural, ", count", sep = ""),
      color = NULL
    ) +
    theme_classic(base_size = 16) +
    scale_color_manual(values = structural_category_palette, drop = FALSE) +
    scale_x_log10(breaks = bulk_len_breaks, labels = scales::comma) +
    annotation_logticks(sides = "b") +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(0.8, "cm"),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 40, r = 20, b = 5, l = 5, unit = "pt"),
      axis.title = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16, angle = 35, hjust = 1)
    ) +
    guides(color = guide_legend(nrow = 2))

  # Mono vs multi-exon classification for length
  Classification_file$exons <- as.numeric(Classification_file$exons)

  gg_bulk_length_by_exon_type <<- ggplot(
    Classification_file[Classification_file$length > 0, ],
    aes(x = length, color = ifelse(exons == 1, "Mono-Exon", "Multi-Exon"))
  ) +
    geom_freqpoly(bins = 80, linewidth = 1.2, na.rm = TRUE) +
    labs(
      title = paste("Mono- vs Multi- Exon", entity_label, "Lengths Distribution"),
      x = "Length, bp",
      y = paste(entity_label_plural, ", count", sep = ""),
      color = NULL
    ) +
    theme_classic(base_size = 16) +
    scale_color_manual(
      values = c("Multi-Exon" = "#3B0057", "Mono-Exon" = "#FFE44C")
    ) +
    scale_x_log10(breaks = bulk_len_breaks, labels = scales::comma) +
    annotation_logticks(sides = "b") +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 40, r = 20, b = 5, l = 5, unit = "pt"),
      axis.title = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16, angle = 35, hjust = 1)
    )

  # Cell-level length distributions (all + mono)
  gg_read_distr <<- build_len_violin_for_prefix(
    SQANTI_cell_summary,
    prefix = "Total",
    title = paste(entity_label_plural, "Length Distribution Across Cells"),
    fill_color = "#CC6633",
    box_fill = "#CC6633",
    mono = FALSE,
    violin_alpha = 0.7,
    box_alpha = 0.6,
    box_outline_color = "grey20",
    violin_outline_fill = FALSE
  )

  # Mono-exon length distribution per break
  gg_read_distr_mono <<- build_len_violin_for_prefix(
    SQANTI_cell_summary,
    prefix = "Total",
    title = paste("Mono-exonic", entity_label_plural, "Length Distribution Across Cells"),
    fill_color = "#CC6633",
    box_fill = "#CC6633",
    mono = TRUE,
    violin_alpha = 0.7,
    box_alpha = 0.6,
    box_outline_color = "grey20",
    violin_outline_fill = FALSE
  )

  # Per-category length distributions via loop
  len_specs <- list(
    list(tag = "FSM", pretty = "FSM", color = "#6BAED6"),
    list(tag = "ISM", pretty = "ISM", color = "#FC8D59"),
    list(tag = "NIC", pretty = "NIC", color = "#78C679"),
    list(tag = "NNC", pretty = "NNC", color = "#EE6A50"),
    list(tag = "Genic", pretty = "Genic", color = "#969696"),
    list(tag = "Antisense", pretty = "Antisense", color = "#66C2A4"),
    list(tag = "Fusion", pretty = "Fusion", color = "goldenrod1"),
    list(tag = "Intergenic", pretty = "Intergenic", color = "darksalmon"),
    list(tag = "Genic_intron", pretty = "Genic Intron", color = "#41B6C4")
  )
  for (sp in len_specs) {
    make_len_plot(sp$tag, sp$pretty, sp$color, mono = FALSE)
  }
  for (sp in len_specs) {
    # Mono versions where meaningful (skip NNC and Fusion for PDF)
    if (sp$tag %in% c("NNC", "Fusion")) next
    make_len_plot(sp$tag, sp$pretty, sp$color, mono = TRUE)
  }

  ### Reference coverage across categories ###
  ############################################

  {
    # Only FSM and ISM have a meaningful ref_length association; other categories are excluded.
    cols <- c("FSM_ref_coverage_prop", "ISM_ref_coverage_prop")
    gg_SQANTI_pivot <- pivot_long(SQANTI_cell_summary, cols)
    fill_map <- c("FSM_ref_coverage_prop" = "#6BAED6", "ISM_ref_coverage_prop" = "#FC8D59")
    x_labels <- c("FSM", "ISM")
    # Build dynamic title using cutoff from cell summary
    ref_cov_min_pct <- if ("ref_cov_min_pct" %in% colnames(SQANTI_cell_summary)) {
      vals <- unique(stats::na.omit(SQANTI_cell_summary$ref_cov_min_pct))
      if (length(vals) > 0) as.numeric(vals[1]) else NA_real_
    } else {
      NA_real_
    }
    pct_lbl <- if (is.finite(ref_cov_min_pct)) {
      if (abs(ref_cov_min_pct - round(ref_cov_min_pct)) < 1e-6) sprintf("%.0f", ref_cov_min_pct) else sprintf("%.1f", ref_cov_min_pct)
    } else {
      NULL
    }
    title_txt <- if (!is.null(pct_lbl)) {
      paste0(entity_label_plural, " with Coverage >=", pct_lbl, "% of the Reference Transcript Length\nby Structural Category Across Cells")
    } else {
      "Reference Transcript Length Coverage\nby Structural Category Across Cells"
    }
    gg_ref_coverage_across_category <<- build_violin_plot(
      gg_SQANTI_pivot,
      title = title_txt,
      x_labels = x_labels,
      fill_map = fill_map,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      legend = FALSE,
      violin_outline_fill = TRUE
    )
  }

  # Meta-transcript body coverage profile (bulk, not per-cell)
  gg_meta_transcript_coverage <<- build_meta_coverage_plot(Classification_file)

  # Reference Transcriptome vs Sample length distribution (both isoforms and reads modes)
  gg_isoforms_ref_vs_sample_lengths <<- build_ref_vs_sample_lengths(Classification_file, ref_gtf_path, mode)

  ### Structural categories ###

  category_fill_map <- c(
    "FSM_prop" = "#6BAED6", "ISM_prop" = "#FC8D59", "NIC_prop" = "#78C679", "NNC_prop" = "#EE6A50",
    "Genic_Genomic_prop" = "#969696", "Antisense_prop" = "#66C2A4", "Fusion_prop" = "goldenrod1",
    "Intergenic_prop" = "darksalmon", "Genic_intron_prop" = "#41B6C4"
  )
  pivot_violin(SQANTI_cell_summary, list(
    name = "gg_SQANTI_across_category",
    columns = names(category_fill_map),
    title = "Structural Categories Distribution Across Cells",
    x_labels = cat_labels_pretty,
    y_label = paste(entity_label_plural, ", %", sep = ""),
    fill_map = category_fill_map,
    plot_args = list(override_outline_vars = c("Genic_Genomic_prop"), violin_outline_fill = TRUE)
  ))

  #  Coding/non-coding across structural categories (change it in the future to a combine plot)
  if (include_ORF) {
    # Update to new column naming convention: {tag}_coding_prop
    # Explicitly define columns to match cell_metrics.py output (lowercase for non-canonical categories)
    coding_cols <- c(
      "FSM_coding_prop", "ISM_coding_prop", "NIC_coding_prop", "NNC_coding_prop",
      "genic_coding_prop", "antisense_coding_prop", "fusion_coding_prop",
      "intergenic_coding_prop", "genic_intron_coding_prop"
    )
    coding_fill_map <- c(
      "FSM_coding_prop" = "#6BAED6",
      "ISM_coding_prop" = "#FC8D59",
      "NIC_coding_prop" = "#78C679",
      "NNC_coding_prop" = "#EE6A50",
      "genic_coding_prop" = "#969696",
      "antisense_coding_prop" = "#66C2A4",
      "fusion_coding_prop" = "goldenrod1",
      "intergenic_coding_prop" = "darksalmon",
      "genic_intron_coding_prop" = "#41B6C4"
    )

    pivot_violin(SQANTI_cell_summary, list(
      name = "gg_coding_across_category",
      columns = names(coding_fill_map),
      title = "Coding Proportion of Structural Categories Distribution Across Cells",
      x_labels = cat_labels_pretty,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = coding_fill_map,
      plot_args = list(override_outline_vars = c("genic_coding_prop"), violin_outline_fill = TRUE)
    ))

    # Define colors for non-coding (same as coding but will use alpha)
    noncoding_fill_map <- c(
      "FSM_non_coding_prop" = "#6BAED6",
      "ISM_non_coding_prop" = "#FC8D59",
      "NIC_non_coding_prop" = "#78C679",
      "NNC_non_coding_prop" = "#EE6A50",
      "genic_non_coding_prop" = "#969696",
      "antisense_non_coding_prop" = "#66C2A4",
      "fusion_non_coding_prop" = "goldenrod1",
      "intergenic_non_coding_prop" = "darksalmon",
      "genic_intron_non_coding_prop" = "#41B6C4"
    )

    pivot_violin(SQANTI_cell_summary, list(
      name = "gg_non_coding_across_category",
      columns = names(noncoding_fill_map),
      title = "Non-coding Proportion of Structural Categories Distribution Across Cells",
      x_labels = cat_labels_pretty,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = noncoding_fill_map,
      plot_args = list(
        override_outline_vars = c(),
        violin_outline_fill = TRUE,
        violin_alpha = 0.4,
        box_alpha = 0.1
      )
    ))
  } # End of if (include_ORF)

  subcategory_configs <- list(
    list(
      name = "gg_SQANTI_across_FSM",
      columns = c(
        "FSM_alternative_3end_prop", "FSM_alternative_3end5end_prop", "FSM_alternative_5end_prop",
        "FSM_reference_match_prop", "FSM_mono_exon_prop"
      ),
      title = "FSM Structural Subcategories Distribution Across Cells",
      x_labels = c("Alternative 3'end", "Alternative 3'5'end", "Alternative 5'end", "Reference match", "Mono-exon"),
      fill_map = c(
        "FSM_alternative_3end_prop" = "#02314d", "FSM_alternative_3end5end_prop" = "#0e5a87",
        "FSM_alternative_5end_prop" = "#7ccdfc", "FSM_reference_match_prop" = "#c4e1f2",
        "FSM_mono_exon_prop" = "#cec2d2"
      ),
      plot_args = list(override_outline_vars = c("FSM_alternative_3end_prop", "FSM_alternative_3end5end_prop"), violin_outline_fill = TRUE)
    ),
    list(
      name = "gg_SQANTI_across_ISM",
      columns = c(
        "ISM_3prime_fragment_prop", "ISM_internal_fragment_prop", "ISM_5prime_fragment_prop",
        "ISM_intron_retention_prop", "ISM_mono_exon_prop"
      ),
      title = "ISM Structural Subcategories Distribution Across Cells",
      x_labels = c("3' fragment", "Internal fragment", "5' fragment", "Intron retention", "Mono-exon"),
      fill_map = c(
        "ISM_3prime_fragment_prop" = "#c4531d", "ISM_internal_fragment_prop" = "#e37744",
        "ISM_5prime_fragment_prop" = "#e0936e", "ISM_intron_retention_prop" = "#81eb82",
        "ISM_mono_exon_prop" = "#cec2d2"
      ),
      plot_args = list(override_outline_vars = c("ISM_3prime_fragment_prop", "ISM_internal_fragment_prop"), violin_outline_fill = TRUE)
    ),
    list(
      name = "gg_SQANTI_across_NIC",
      columns = c(
        "NIC_combination_of_known_junctions_prop", "NIC_combination_of_known_splicesites_prop",
        "NIC_intron_retention_prop", "NIC_mono_exon_by_intron_retention_prop", "NIC_mono_exon_prop"
      ),
      title = "NIC Structural Subcategories Distribution Across Cells",
      x_labels = c("Comb. of annot. junctions", "Comb. of annot. splice sites", "Intron retention", "Mono-exon by intron ret.", "Mono-exon"),
      fill_map = c(
        "NIC_combination_of_known_junctions_prop" = "#014d02", "NIC_combination_of_known_splicesites_prop" = "#379637",
        "NIC_intron_retention_prop" = "#81eb82", "NIC_mono_exon_by_intron_retention_prop" = "#4aaa72",
        "NIC_mono_exon_prop" = "#cec2d2"
      ),
      plot_args = list(
        override_outline_vars = c(
          "NIC_combination_of_known_junctions_prop", "NIC_combination_of_known_splicesites_prop",
          "NIC_mono_exon_by_intron_retention_prop", "NIC_mono_exon_prop"
        ),
        violin_outline_fill = TRUE
      )
    ),
    list(
      name = "gg_SQANTI_across_NNC",
      columns = c("NNC_at_least_one_novel_splicesite_prop", "NNC_intron_retention_prop"),
      title = "NNC Structural Subcategories Distribution Across Cells",
      x_labels = c("At least\n1 annot. don./accept.", "Intron retention"),
      fill_map = c("NNC_at_least_one_novel_splicesite_prop" = "#32734d", "NNC_intron_retention_prop" = "#81eb82"),
      plot_args = list(override_outline_vars = c("NNC_at_least_one_novel_splicesite_prop"), violin_outline_fill = TRUE)
    ),
    list(
      name = "gg_SQANTI_across_Fusion",
      columns = c("Fusion_intron_retention_prop", "Fusion_multi_exon_prop"),
      title = "Fusion Structural Subcategories Distribution Across Cells",
      x_labels = c("Intron retention", "Multi-exon"),
      fill_map = c("Fusion_intron_retention_prop" = "#81eb82", "Fusion_multi_exon_prop" = "#876a91"),
      plot_args = list(override_outline_vars = c("Fusion_multi_exon_prop"), violin_outline_fill = TRUE)
    ),
    list(
      name = "gg_SQANTI_across_Genic",
      columns = c("Genic_mono_exon_prop", "Genic_multi_exon_prop"),
      title = "Genic Structural Subcategories Distribution Across Cells",
      x_labels = c("Mono-exon", "Multi-exon"),
      fill_map = c("Genic_mono_exon_prop" = "#81eb82", "Genic_multi_exon_prop" = "#876a91"),
      plot_args = list(override_outline_vars = c("Genic_multi_exon_prop"), violin_outline_fill = TRUE)
    ),
    list(
      name = "gg_SQANTI_across_Genic_Intron",
      columns = c("Genic_intron_mono_exon_prop", "Genic_intron_multi_exon_prop"),
      title = "Genic Intron Structural Subcategories Distribution Across Cells",
      x_labels = c("Mono-exon", "Multi-exon"),
      fill_map = c("Genic_intron_mono_exon_prop" = "#81eb82", "Genic_intron_multi_exon_prop" = "#876a91"),
      plot_args = list(override_outline_vars = c("Genic_intron_multi_exon_prop"), violin_outline_fill = TRUE)
    ),
    list(
      name = "gg_SQANTI_across_Antisense",
      columns = c("Antisense_mono_exon_prop", "Antisense_multi_exon_prop"),
      title = "Antisense Structural Subcategories Distribution Across Cells",
      x_labels = c("Mono-exon", "Multi-exon"),
      fill_map = c("Antisense_mono_exon_prop" = "#81eb82", "Antisense_multi_exon_prop" = "#876a91"),
      plot_args = list(override_outline_vars = c("Antisense_multi_exon_prop"), violin_outline_fill = TRUE)
    ),
    list(
      name = "gg_SQANTI_across_Intergenic",
      columns = c("Intergenic_mono_exon_prop", "Intergenic_multi_exon_prop"),
      title = "Intergenic Structural Subcategories Distribution Across Cells",
      x_labels = c("Mono-exon", "Multi-exon"),
      fill_map = c("Intergenic_mono_exon_prop" = "#81eb82", "Intergenic_multi_exon_prop" = "#876a91"),
      plot_args = list(override_outline_vars = c("Intergenic_multi_exon_prop"), violin_outline_fill = TRUE)
    )
  )
  invisible(lapply(subcategory_configs, function(cfg) pivot_violin(SQANTI_cell_summary, cfg)))

  ### Splice junctions characterization ###
  #########################################

  pivot_violin(SQANTI_cell_summary, list(
    name = "gg_known_novel_canon",
    columns = c("Known_canonical_junctions_prop", "Known_non_canonical_junctions_prop", "Novel_canonical_junctions_prop", "Novel_non_canonical_junctions_prop"),
    title = "Splice Junctions Distribution Across Cells",
    x_labels = c("Known\nCanonical", "Known\nNon-canonical", "Novel\nCanonical", "Novel\nNon-canonical"),
    y_label = "Junctions, %",
    fill_map = c(
      "Known_canonical_junctions_prop" = "#6BAED6",
      "Known_non_canonical_junctions_prop" = "goldenrod1",
      "Novel_canonical_junctions_prop" = "#78C679",
      "Novel_non_canonical_junctions_prop" = "#FC8D59"
    ),
    plot_args = list(violin_outline_fill = TRUE)
  ))
  ### Good features plots (SR & TSS Validation) ###
  #################################################

  # 1. Combined Good Features Plot (All Transcripts)
  all_good_features_map <- list(
    "srjunctions_support_prop" = list(label = "SJs Validated by SRs", color = "#cd4f39"),
    "TSS_ratio_validated_prop" = list(label = "TSS Validated by SRs", color = "#FFC125")
    # Add other good features here if needed (e.g. polyA_motif_found_prop if available)
  )

  # Determine which good feature columns are present
  good_feature_cols_present <- intersect(names(all_good_features_map), colnames(SQANTI_cell_summary))
  good_feature_cols_present <- good_feature_cols_present[sapply(good_feature_cols_present, function(col) any(!is.na(SQANTI_cell_summary[[col]])) && sum(SQANTI_cell_summary[[col]], na.rm = TRUE) > 0)]

  if (length(good_feature_cols_present) > 0) {
    current_good_colors <- sapply(all_good_features_map[good_feature_cols_present], function(x) x$color)
    current_good_labels <- sapply(all_good_features_map[good_feature_cols_present], function(x) x$label)
    names(current_good_colors) <- good_feature_cols_present
    names(current_good_labels) <- good_feature_cols_present

    pivot_violin(SQANTI_cell_summary, list(
      name = "gg_good_feature",
      columns = good_feature_cols_present,
      title = "Validation Features Distribution Across Cells",
      x_labels = current_good_labels,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = current_good_colors,
      plot_args = list(violin_outline_fill = TRUE)
    ))
  }

  # 2. Per-Category Plots
  # Short Read (SJs) Support
  sr_cat_cols <- cat_cols("_srjunctions_support_prop")
  if (all(sr_cat_cols %in% colnames(SQANTI_cell_summary)) && any(colSums(SQANTI_cell_summary[, sr_cat_cols, drop = FALSE], na.rm = TRUE) > 0)) {
    pivot_violin(SQANTI_cell_summary, list(
      name = "gg_sr_support_by_category",
      columns = sr_cat_cols,
      title = "SJs Validated by Short Reads by Structural Category",
      x_labels = cat_labels_pretty,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = setNames(rep("#cd4f39", length(sr_cat_cols)), sr_cat_cols),
      plot_args = list(violin_outline_fill = TRUE)
    ))
  }

  # TSS Validation Support
  tss_cat_cols <- cat_cols("_TSS_ratio_validated_prop")
  if (all(tss_cat_cols %in% colnames(SQANTI_cell_summary)) && any(colSums(SQANTI_cell_summary[, tss_cat_cols, drop = FALSE], na.rm = TRUE) > 0)) {
    pivot_violin(SQANTI_cell_summary, list(
      name = "gg_tss_validation_by_category",
      columns = tss_cat_cols,
      title = "TSS Validated by Short Reads by Structural Category",
      x_labels = cat_labels_pretty,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = setNames(rep("#FFC125", length(tss_cat_cols)), tss_cat_cols),
      plot_args = list(violin_outline_fill = TRUE)
    ))
  }
  ### Bad features plots ###
  ##########################

  bad_specs <- list(
    list(suffix = "_intrapriming_prop", title = "Intrapriming by Structural Category", color = "#78C679", name = "gg_intrapriming_by_category"),
    list(suffix = "_RTS_prop", title = "RT-switching by Structural Category", color = "#FF9933", name = "gg_RTS_by_category"),
    list(suffix = "_noncanon_prop", title = "Non-Canonical Junctions by Structural Category", color = "#41B6C4", name = "gg_noncanon_by_category")
  )
  invisible(lapply(bad_specs, function(sp) {
    cols <- cat_cols(sp$suffix)
    pivot_violin(SQANTI_cell_summary, list(
      name = sp$name,
      columns = cols,
      title = sp$title,
      x_labels = cat_labels_pretty,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = setNames(rep(sp$color, length(cols)), cols),
      plot_args = list(violin_outline_fill = TRUE)
    ))
  }))

  # NMD  (split between categories)
  nmd_cols <- c("FSM_NMD_prop", "ISM_NMD_prop", "NIC_NMD_prop", "NNC_NMD_prop", "Genic_NMD_prop", "Antisense_NMD_prop", "Fusion_NMD_prop", "Intergenic_NMD_prop", "Genic_intron_NMD_prop")
  if (all(nmd_cols %in% colnames(SQANTI_cell_summary))) {
    pivot_violin(SQANTI_cell_summary, list(
      name = "gg_NMD_by_category",
      columns = nmd_cols,
      title = "Nonsense-Mediated Decay by Structural Category",
      x_labels = cat_labels_pretty,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = setNames(rep("#969696", length(nmd_cols)), nmd_cols),
      plot_args = list(violin_outline_fill = TRUE)
    ))
  }

  ## Bad quality features combined figure
  # Define all possible features, their colors, and labels
  all_bad_features_map <- list(
    "Intrapriming_prop_in_cell" = list(label = "Intrapriming", color = "#78C679"),
    "RTS_prop_in_cell" = list(label = "RT-switching", color = "#FF9933"),
    "Non_canonical_prop_in_cell" = list(label = "Non-Canonical Junctions", color = "#41B6C4"),
    "NMD_prop_in_cell" = list(label = "Predicted NMD", color = "#969696")
  )

  # Determine which bad feature columns are actually present in SQANTI_cell_summary.
  # RTS/intrapriming/non-canonical are shown even when all-zero (0 detected is informative).
  # NMD requires at least one non-zero value: all-zero means --include_ORF was not passed
  # and the column carries no real information.
  bad_feature_cols_present <- intersect(names(all_bad_features_map), colnames(SQANTI_cell_summary))
  bad_feature_cols_present <- bad_feature_cols_present[sapply(bad_feature_cols_present, function(col) {
    if (!any(!is.na(SQANTI_cell_summary[[col]]))) return(FALSE)
    if (col == "NMD_prop_in_cell") return(sum(SQANTI_cell_summary[[col]], na.rm = TRUE) > 0)
    return(TRUE)
  })]

  # Order them as originally intended, if present
  ordered_bad_feature_cols <- c("Intrapriming_prop_in_cell", "RTS_prop_in_cell", "Non_canonical_prop_in_cell", "NMD_prop_in_cell")
  bad_feature_cols_present <- intersect(ordered_bad_feature_cols, bad_feature_cols_present)

  if (length(bad_feature_cols_present) > 0) {
    current_colors <- sapply(all_bad_features_map[bad_feature_cols_present], function(x) x$color)
    current_labels <- sapply(all_bad_features_map[bad_feature_cols_present], function(x) x$label)
    # Ensure names are correctly assigned for scales, matching the order in bad_feature_cols_present
    names(current_colors) <- bad_feature_cols_present
    names(current_labels) <- bad_feature_cols_present

    pivot_violin(SQANTI_cell_summary, list(
      name = "gg_bad_feature",
      columns = bad_feature_cols_present,
      title = "Bad Quality Control Attributes Across Cells",
      x_labels = current_labels,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = current_colors,
      plot_args = list(violin_outline_fill = TRUE)
    ))
  } else {
    gg_bad_feature <<- ggplot() +
      labs(title = "Plot not available") +
      theme_minimal()
    layout(
      title = "No bad quality features to display",
      annotations = list(
        text = "No bad quality features to display",
        showarrow = FALSE,
        font = list(size = 18, color = "gray")
      )
    )
  }

  # Good features plots

  good_specs <- list(
    list(cols = cat_cols("_TSSAnnotationSupport"), title = "TSS Annotation Support by Structural Category", color = "#66C2A4", name = "gg_tss_annotation_support", require_all = TRUE),
    list(cols = cat_cols("_CAGE_peak_support_prop"), title = "CAGE Peak Support by Structural Category", color = "#EE6A50", name = "gg_cage_peak_support", require_all = TRUE),
    list(cols = cat_cols("_PolyA_motif_support_prop"), title = "PolyA Support by Structural Category", color = "#78C679", name = "gg_polyA_motif_support", require_all = TRUE),
    list(cols = cat_cols("_canon_prop"), title = "Canonical Junctions by Structural Category", color = "#CC6633", name = "gg_canon_by_category", require_all = TRUE),
    list(cols = cat_cols("_srjunctions_support_prop"), title = "Splice Junctions Support by Structural Category", color = "#cd4f39", name = "gg_sr_support_by_category", require_all = TRUE),
    list(cols = cat_cols("_TSS_ratio_validated_prop"), title = "TSS Support by Structural Category", color = "#ffc125", name = "gg_tss_validation_by_category", require_all = TRUE)
  )
  invisible(lapply(good_specs, function(sp) {
    if (!is.null(sp$require_all) && sp$require_all && !all(sp$cols %in% colnames(SQANTI_cell_summary))) {
      return(NULL)
    }
    # Check if data is not empty (all zeros)
    if (all(sp$cols %in% colnames(SQANTI_cell_summary)) && all(colSums(SQANTI_cell_summary[, sp$cols, drop = FALSE], na.rm = TRUE) == 0)) {
      return(NULL)
    }
    pivot_violin(SQANTI_cell_summary, list(
      name = sp$name,
      columns = sp$cols,
      title = sp$title,
      x_labels = cat_labels_pretty,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = setNames(rep(sp$color, length(sp$cols)), sp$cols),
      plot_args = list(violin_outline_fill = TRUE)
    ))
  }))

  ## Good quality features combined figure
  good_feature_cols <- c("TSSAnnotationSupport_prop")

  if ("CAGE_peak_support_prop" %in% colnames(SQANTI_cell_summary)) {
    good_feature_cols <- c(good_feature_cols, "CAGE_peak_support_prop")
  }
  if ("PolyA_motif_support_prop" %in% colnames(SQANTI_cell_summary)) {
    good_feature_cols <- c(good_feature_cols, "PolyA_motif_support_prop")
  }
  good_feature_cols <- c(good_feature_cols, "Canonical_prop_in_cell")

  if ("srjunctions_support_prop" %in% colnames(SQANTI_cell_summary) && sum(SQANTI_cell_summary$srjunctions_support_prop, na.rm = TRUE) > 0) {
    good_feature_cols <- c(good_feature_cols, "srjunctions_support_prop")
  }
  if ("TSS_ratio_validated_prop" %in% colnames(SQANTI_cell_summary) && sum(SQANTI_cell_summary$TSS_ratio_validated_prop, na.rm = TRUE) > 0) {
    good_feature_cols <- c(good_feature_cols, "TSS_ratio_validated_prop")
  }

  color_map <- c(
    "TSSAnnotationSupport_prop" = "#66C2A4",
    "CAGE_peak_support_prop" = "#EE6A50",
    "PolyA_motif_support_prop" = "#78C679",
    "Canonical_prop_in_cell" = "#CC6633",
    "srjunctions_support_prop" = "#cd4f39",
    "TSS_ratio_validated_prop" = "#ffc125"
  )
  label_map <- c(
    "TSSAnnotationSupport_prop" = "TSS Annotated",
    "CAGE_peak_support_prop" = "Has Coverage CAGE",
    "PolyA_motif_support_prop" = "Has PolyA Motif",
    "Canonical_prop_in_cell" = "Canonical Junctions",
    "srjunctions_support_prop" = "SJs Support by SRs",
    "TSS_ratio_validated_prop" = "TSS Support by SRs"
  )
  color_map <- color_map[good_feature_cols]
  label_map <- label_map[good_feature_cols]

  pivot_violin(SQANTI_cell_summary, list(
    name = "gg_good_feature",
    columns = good_feature_cols,
    title = "Good Quality Control Attributes Across Cells",
    x_labels = label_map,
    y_label = paste(entity_label_plural, ", %", sep = ""),
    fill_map = color_map,
    plot_args = list(violin_outline_fill = TRUE)
  ))

  ### Exon structure across cells by structural category ###
  {
    cat_key_map <- structural_category_map
    all_cats <- structural_category_levels
    x_labels_pretty <- cat_labels_pretty

    fill_map_cat <- cat_fill_map

    # Build per-cell classification with FL-weighted count for ALL exon structure plots.
    # In isoforms mode: explode comma-separated CB/FL so each row = (isoform, cell, count=FL).
    # This means all exon metrics reflect actual expression levels, not just unique isoform counts.
    if (mode == "isoforms" && "FL" %in% colnames(Classification_file) &&
        "CB" %in% colnames(Classification_file)) {
      cls_valid <- Classification_file %>%
        filter(!is.na(CB), CB != "", CB != "unassigned") %>%
        mutate(CB_raw = as.character(CB), FL_raw = as.character(FL)) %>%
        tidyr::separate_rows(CB_raw, FL_raw, sep = ",") %>%
        mutate(CB    = trimws(CB_raw),
               count = suppressWarnings(as.numeric(trimws(FL_raw)))) %>%
        filter(CB != "", CB != "unassigned", !is.na(count), count > 0) %>%
        mutate(cat_key = unname(cat_key_map[structural_category])) %>%
        filter(!is.na(cat_key))
    } else {
      cls_valid <- Classification_file %>%
        dplyr::filter(!is.na(CB) & CB != "" & CB != "unassigned") %>%
        mutate(cat_key = unname(cat_key_map[structural_category]),
               count   = 1) %>%
        filter(!is.na(cat_key))
    }

    # 1) FL-weighted mean exons per cell per category
    exons_mean_by_cell <- cls_valid %>%
      group_by(CB, cat_key) %>%
      summarise(
        mean_exons = sum(as.numeric(exons) * count, na.rm = TRUE) / sum(count, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::complete(CB, cat_key = all_cats, fill = list(mean_exons = NA_real_))

    df_exon_mean_long <- data.frame(
      Variable = factor(exons_mean_by_cell$cat_key, levels = all_cats),
      Value = exons_mean_by_cell$mean_exons
    )

    gg_exon_mean_by_category <<- build_violin_plot(
      df_long = df_exon_mean_long,
      title = paste("Mean Exons per", entity_label, "by Structural Category Across Cells"),
      x_labels = x_labels_pretty,
      fill_map = fill_map_cat,
      y_label = paste("Exons per", entity_label),
      legend = FALSE,
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      override_outline_vars = c("Genic"),
      adjust = 1
    )

    # 2) Percent mono-exonic reads per cell and category (FL-weighted in isoforms mode)
    exons_bin_by_cell <- cls_valid %>%
      mutate(is_mono = as.numeric(exons) == 1) %>%
      group_by(CB, cat_key) %>%
      summarise(total = sum(count, na.rm = TRUE),
                mono  = sum(count[is_mono], na.rm = TRUE), .groups = "drop") %>%
      mutate(perc_mono = ifelse(total > 0, 100 * mono / total, NA_real_)) %>%
      tidyr::complete(CB, cat_key = all_cats)

    df_exon_mono_long <- data.frame(
      Variable = factor(exons_bin_by_cell$cat_key, levels = all_cats),
      Value = exons_bin_by_cell$perc_mono
    )

    gg_exon_mono_by_category <<- build_violin_plot(
      df_long = df_exon_mono_long,
      title = paste("Mono-exonic", entity_label_plural, "by Structural Category Across Cells"),
      x_labels = x_labels_pretty,
      fill_map = fill_map_cat,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      legend = FALSE,
      ylim = c(0, 100),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      override_outline_vars = c("Genic")
    )

    # 3) Exon count bins per structural category across cells (HTML)
    exon_bin_levels <- c("1", "2-3", "4-5", ">=6")
    bin_fill_map <- setNames(c("#6BAED6", "#78C679", "#FC8D59", "#969696"), exon_bin_levels)
    gg_exon_bins_by_category <<- list()
    for (ck in all_cats) {
      cat_df <- cls_valid %>% filter(cat_key == ck)
      if (nrow(cat_df) == 0) {
        next
      }
      bins_by_cell <- cat_df %>%
        mutate(
          exons_n = as.numeric(exons),
          bin = dplyr::case_when(
            exons_n <= 1 ~ "1",
            exons_n <= 3 ~ "2-3",
            exons_n <= 5 ~ "4-5",
            TRUE ~ ">=6"
          )
        ) %>%
        group_by(CB, bin) %>%
        summarise(n = sum(count, na.rm = TRUE), .groups = "drop") %>%
        group_by(CB) %>%
        mutate(perc = ifelse(sum(n) > 0, 100 * n / sum(n), NA_real_)) %>%
        ungroup() %>%
        tidyr::complete(CB, bin = exon_bin_levels)

      df_long_bins <- data.frame(
        Variable = factor(bins_by_cell$bin, levels = exon_bin_levels),
        Value = bins_by_cell$perc
      )

      pretty_name <- switch(ck,
        FSM = "FSM",
        ISM = "ISM",
        NIC = "NIC",
        NNC = "NNC",
        Genic = "Genic Genomic",
        Antisense = "Antisense",
        Fusion = "Fusion",
        Intergenic = "Intergenic",
        Genic_intron = "Genic Intron",
        ck
      )

      gg_exon_bins_by_category[[pretty_name]] <<- build_violin_plot(
        df_long = df_long_bins,
        title = paste0("Exon Count Bins in ", pretty_name, " Across Cells"),
        x_labels = exon_bin_levels,
        fill_map = bin_fill_map,
        y_label = paste(entity_label_plural, ", %", sep = ""),
        legend = FALSE,
        ylim = c(0, 100),
        violin_alpha = 0.7,
        box_alpha = 0.3,
        box_width = 0.05,
        x_tickangle = 45,
        violin_outline_fill = TRUE,
        box_outline_default = "grey20"
      )
    }

    # 4) Exon count profile per category across cells (median + IQR)
    K <- 20
    min_reads <- 5
    gg_exon_profile_by_category <<- list()
    for (ck in all_cats) {
      cat_df <- cls_valid %>%
        mutate(cat_key = unname(cat_key_map[structural_category])) %>%
        dplyr::filter(cat_key == ck)
      if (nrow(cat_df) == 0) next
      # Cells with at least min_reads (FL-weighted) in this category
      cells_ok <- cat_df %>%
        group_by(CB) %>%
        summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>%
        filter(total >= min_reads) %>%
        pull(CB)
      if (length(cells_ok) == 0) {
        cells_ok <- unique(cat_df$CB)
      }
      cat_df2 <- cat_df %>%
        filter(CB %in% cells_ok) %>%
        mutate(exons_n = pmin(as.numeric(exons), K))
      # Per-cell PMF (FL-weighted: sum of FL counts in each exon bin)
      pmf <- cat_df2 %>%
        group_by(CB, exons_n) %>%
        summarise(n = sum(count, na.rm = TRUE), .groups = "drop") %>%
        group_by(CB) %>%
        mutate(perc = ifelse(sum(n) > 0, 100 * n / sum(n), 0)) %>%
        ungroup() %>%
        tidyr::complete(CB, exons_n = seq_len(K), fill = list(n = 0, perc = 0))
      # Aggregate across cells
      prof <- pmf %>%
        group_by(exons_n) %>%
        summarise(
          mean = base::mean(perc, na.rm = TRUE),
          median = stats::median(perc, na.rm = TRUE),
          q1 = stats::quantile(perc, 0.25, na.rm = TRUE, type = 7),
          q3 = stats::quantile(perc, 0.75, na.rm = TRUE, type = 7),
          .groups = "drop"
        ) %>%
        rename(k = exons_n)
      pretty_name <- switch(ck,
        FSM = "FSM",
        ISM = "ISM",
        NIC = "NIC",
        NNC = "NNC",
        Genic = "Genic Genomic",
        Antisense = "Antisense",
        Fusion = "Fusion",
        Intergenic = "Intergenic",
        Genic_intron = "Genic Intron",
        ck
      )
      gg_exon_profile_by_category[[pretty_name]] <<- build_exon_profile_plot(
        df_prof = prof, title = paste0("Exon Count Profile in ", pretty_name, " Across Cells"),
        line_color = fill_map_cat[ck], k_max = K, y_label = paste(entity_label_plural, ", %", sep = ""), n_cells = length(unique(cells_ok))
      )
    }
  }

  ### Presets ###
  ###############

  # t1 <- ttheme_default(core=list(core = list(fg_params = list(cex = 0.6)),
  #                                colhead = list(fg_params = list(cex = 0.7))))

  # Resolve common isoform ID key between Junctions and Classification_file once,
  # shared by all junction blocks below (HTML junc_aug_html, junc_rt, and PDF junc_aug).
  junc_iso_key <- NULL
  for (.k in c("isoform", "readID", "read_id", "ID", "read_name", "read")) {
    if (.k %in% colnames(Junctions) && .k %in% colnames(Classification_file)) {
      junc_iso_key <- .k
      break
    }
  }

  # Build SJ per-type/per-category plots and the all-canonical grouped plot for HTML (and reuse for PDF)
  # This block creates plot objects regardless of generate_pdf so the Rmd can render them.
  {
    # Build junc_aug_html: in isoforms mode, Junctions$CB is comma-separated — must explode.
    if (mode == "isoforms" && !is.null(junc_iso_key) &&
        "FL" %in% colnames(Classification_file) && "CB" %in% colnames(Classification_file) &&
        junc_iso_key %in% colnames(Junctions)) {

      cls_exp_html <- Classification_file %>%
        filter(!is.na(CB), !is.na(structural_category)) %>%
        select(all_of(c(junc_iso_key, "CB", "FL", "structural_category"))) %>%
        mutate(CB_raw = as.character(CB), FL_raw = as.character(FL)) %>%
        tidyr::separate_rows(CB_raw, FL_raw, sep = ",") %>%
        mutate(CB_clean = trimws(CB_raw),
               FL_num   = suppressWarnings(as.numeric(trimws(FL_raw)))) %>%
        filter(CB_clean != "", CB_clean != "unassigned", !is.na(FL_num), FL_num > 0) %>%
        select(all_of(c(junc_iso_key, "CB_clean", "FL_num", "structural_category")))

      junc_aug_html <- Junctions %>%
        select(all_of(c(junc_iso_key, "junction_category", "canonical"))) %>%
        mutate(junction_type = paste(junction_category, canonical, sep = "_")) %>%
        inner_join(cls_exp_html, by = junc_iso_key) %>%
        rename(CB = CB_clean, count = FL_num)

    } else {
      junc_aug_html <- Junctions %>%
        dplyr::filter(!is.na(CB) & CB != "" & CB != "unassigned") %>%
        mutate(junction_type = paste(junction_category, canonical, sep = "_"),
               count = 1)
      if (!("structural_category" %in% colnames(junc_aug_html)) && !is.null(junc_iso_key)) {
        junc_aug_html <- junc_aug_html %>%
          left_join(Classification_file %>% select(all_of(c(junc_iso_key, "structural_category"))),
                    by = junc_iso_key)
      }
    }

    cat_key_map <- structural_category_map
    all_cats <- structural_category_levels
    x_labels_full <- cat_labels_pretty

    junc_summ_html <- junc_aug_html %>%
      filter(!is.na(structural_category)) %>%
      mutate(cat_key = unname(cat_key_map[structural_category])) %>%
      filter(!is.na(cat_key)) %>%
      group_by(CB, cat_key) %>%
      summarise(
        total               = sum(count, na.rm = TRUE),
        known_canonical     = sum(count[junction_type == "known_canonical"], na.rm = TRUE),
        known_non_canonical = sum(count[junction_type == "known_non_canonical"], na.rm = TRUE),
        novel_canonical     = sum(count[junction_type == "novel_canonical"], na.rm = TRUE),
        novel_non_canonical = sum(count[junction_type == "novel_non_canonical"], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        KnownCanonicalPerc    = ifelse(total > 0, 100 * known_canonical    / total, NA_real_),
        KnownNonCanonicalPerc = ifelse(total > 0, 100 * known_non_canonical / total, NA_real_),
        NovelCanonicalPerc    = ifelse(total > 0, 100 * novel_canonical     / total, NA_real_),
        NovelNonCanonicalPerc = ifelse(total > 0, 100 * novel_non_canonical / total, NA_real_)
      ) %>%
      tidyr::complete(CB, cat_key = all_cats) %>%
      ungroup()

    make_df_long_html <- function(col_name) {
      data.frame(Variable = factor(junc_summ_html$cat_key, levels = all_cats), Value = junc_summ_html[[col_name]])
    }

    fill_map_cat <- cat_fill_map

    gg_known_canon_by_category <<- build_violin_plot(
      df_long = make_df_long_html("KnownCanonicalPerc"),
      title = "",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Known Canonical Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )

    gg_known_noncanon_by_category <<- build_violin_plot(
      df_long = make_df_long_html("KnownNonCanonicalPerc"),
      title = "",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Known Non-canonical Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )

    gg_novel_canon_by_category <<- build_violin_plot(
      df_long = make_df_long_html("NovelCanonicalPerc"),
      title = "",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Novel Canonical Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )

    gg_novel_noncanon_by_category <<- build_violin_plot(
      df_long = make_df_long_html("NovelNonCanonicalPerc"),
      title = "",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Novel Non-canonical Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )

    # Stack the four SJ type-by-category plots into one figure
    tick_angle_val <- 45

    p_known_canon_by_category <- build_violin_plot(
      df_long = make_df_long_html("KnownCanonicalPerc"),
      title = "",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Known Canonical Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )

    p_known_noncanon_by_category <- build_violin_plot(
      df_long = make_df_long_html("KnownNonCanonicalPerc"),
      title = "",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Known Non-canonical Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )

    p_novel_canon_by_category <- build_violin_plot(
      df_long = make_df_long_html("NovelCanonicalPerc"),
      title = "",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Novel Canonical Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )

    p_novel_noncanon_by_category <- build_violin_plot(
      df_long = make_df_long_html("NovelNonCanonicalPerc"),
      title = "",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Novel Non-canonical Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )

    # Stack the four SJ type-by-category plots into one static figure using gridExtra
    # Remove x-axis labels/titles for top 3 plots to mimic shared axis
    p1 <- p_known_canon_by_category + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
    p2 <- p_known_noncanon_by_category + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
    p3 <- p_novel_canon_by_category + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
    p4 <- p_novel_noncanon_by_category # Keep x-axis for bottom plot

    gg_sj_type_by_category_stack <<- gridExtra::arrangeGrob(
      p1, p2, p3, p4,
      ncol = 1,
      heights = unit(c(1, 1, 1, 1.25), "null"),
      top = textGrob("Splice Junctions Distribution by Structural Category Across Cells", gp = gpar(fontsize = 18, fontface = "bold"))
    )
  }

  # NEW: RT-switching by splice junction type across cells (all and unique junctions)
  if ("RTS_junction" %in% colnames(Junctions)) {
    # Normalize boolean
    rts_bool_vec <- tolower(as.character(Junctions$RTS_junction)) %in% c("true", "t", "1", "yes")

    # Optional genomic coordinates for collapsing duplicate junctions (same chr + positions)
    rts_coord_cols <- intersect(
      c("chrom", "chr", "genomic_start_coord", "genomic_end_coord", "strand"),
      colnames(Junctions)
    )
    has_rts_junc_coords <- all(c("genomic_start_coord", "genomic_end_coord") %in% rts_coord_cols) &&
      ("chrom" %in% rts_coord_cols || "chr" %in% rts_coord_cols)

    # In isoforms mode, Junctions$CB is a comma-separated list — must explode via classification
    if (mode == "isoforms" && !is.null(junc_iso_key) &&
        "FL" %in% colnames(Classification_file) && "CB" %in% colnames(Classification_file) &&
        junc_iso_key %in% colnames(Junctions)) {

      # Explode Classification_file CB/FL to per-(isoform, cell, FL_count)
      cls_exp_rts <- Classification_file %>%
        filter(!is.na(CB)) %>%
        select(all_of(c(junc_iso_key, "CB", "FL"))) %>%
        mutate(CB_raw = as.character(CB), FL_raw = as.character(FL)) %>%
        tidyr::separate_rows(CB_raw, FL_raw, sep = ",") %>%
        mutate(CB_clean = trimws(CB_raw),
               FL_num   = suppressWarnings(as.numeric(trimws(FL_raw)))) %>%
        filter(CB_clean != "", CB_clean != "unassigned", !is.na(FL_num), FL_num > 0) %>%
        select(all_of(c(junc_iso_key, "CB_clean", "FL_num")))

      junc_rt_cols <- unique(c(
        junc_iso_key, "junction_category", "canonical", "RTS_junction",
        rts_coord_cols
      ))
      junc_rt_cols <- junc_rt_cols[junc_rt_cols %in% colnames(Junctions)]

      # Join junctions (isoform, junction type, RTS flag) to exploded classification
      junc_rt <- Junctions %>%
        select(all_of(junc_rt_cols)) %>%
        mutate(SJ_type  = paste(junction_category, canonical, sep = "_"),
               RTS_bool = tolower(as.character(RTS_junction)) %in% c("true", "t", "1", "yes")) %>%
        inner_join(cls_exp_rts, by = junc_iso_key) %>%
        rename(CB = CB_clean, count = FL_num)

    } else {
      # Reads mode: each junction already has one CB, count = 1
      junc_rt_cols <- unique(c(
        "CB", "junction_category", "canonical", "RTS_junction",
        rts_coord_cols
      ))
      junc_rt_cols <- junc_rt_cols[junc_rt_cols %in% colnames(Junctions)]

      junc_rt <- Junctions %>%
        select(all_of(junc_rt_cols)) %>%
        dplyr::filter(!is.na(CB) & CB != "" & CB != "unassigned") %>%
        mutate(SJ_type  = paste(junction_category, canonical, sep = "_"),
               RTS_bool = rts_bool_vec,
               count    = 1)
    }

    # Ensure consistent SJ type levels and labels
    sj_levels <- c("known_canonical", "known_non_canonical", "novel_canonical", "novel_non_canonical")
    sj_labels <- c("Known\nCanonical", "Known\nNon-canonical", "Novel\nCanonical", "Novel\nNon-canonical")
    junc_rt$SJ_type <- factor(junc_rt$SJ_type, levels = sj_levels)

    # Color map consistent with other SJ type plots
    sj_fill_map <- c(
      known_canonical = "#6BAED6",
      known_non_canonical = "goldenrod1",
      novel_canonical = "#78C679",
      novel_non_canonical = "#FC8D59"
    )

    # Per-cell percentages for ALL junctions (FL-weighted in isoforms mode)
    all_junc_by_cell <- junc_rt %>%
      group_by(CB, SJ_type) %>%
      summarise(total = sum(count, na.rm = TRUE),
                rts   = sum(count[RTS_bool], na.rm = TRUE), .groups = "drop") %>%
      # Only include (CB, SJ_type) that actually have data; NA for absent pairs
      mutate(perc = ifelse(total > 0, 100 * rts / total, NA_real_)) %>%
      tidyr::complete(CB, SJ_type = sj_levels)

    df_long_all <- data.frame(
      Variable = factor(all_junc_by_cell$SJ_type, levels = sj_levels),
      Value = all_junc_by_cell$perc
    )

    gg_rts_all_by_sjtype <<- build_violin_plot(
      df_long = df_long_all,
      title = "RT-switching All Junctions by Splice Junction Type Across Cells",
      x_labels = sj_labels,
      fill_map = sj_fill_map,
      y_label = "Junctions, %",
      legend = FALSE,
      ylim = c(0, 100),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20"
    )

    # Unique junctions per cell: collapse identical genomic introns (reads mode: 1 per junction;
    # isoforms mode: 1 per junction per cell, FL not used). RTS if any supporting row is RTS.
    if (has_rts_junc_coords) {
      uniq_junc_by_cell <- junc_rt %>%
        mutate(
          junc_chr = dplyr::coalesce(
            if ("chrom" %in% names(.)) as.character(.data$chrom) else NA_character_,
            if ("chr" %in% names(.)) as.character(.data$chr) else NA_character_
          ),
          strand_part = if ("strand" %in% names(.)) as.character(.data$strand) else "",
          junc_key = paste(
            junc_chr,
            as.character(.data$genomic_start_coord),
            as.character(.data$genomic_end_coord),
            strand_part,
            sep = ":"
          )
        ) %>%
        dplyr::filter(
          !is.na(junc_chr), nzchar(junc_chr),
          !is.na(.data$genomic_start_coord), !is.na(.data$genomic_end_coord)
        ) %>%
        group_by(CB, SJ_type, junc_key) %>%
        summarise(RTS_any = any(RTS_bool, na.rm = TRUE), .groups = "drop") %>%
        group_by(CB, SJ_type) %>%
        summarise(
          total = dplyr::n(),
          rts   = sum(RTS_any, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(perc = ifelse(total > 0, 100 * rts / total, NA_real_)) %>%
        tidyr::complete(CB, SJ_type = sj_levels)

      df_long_unique <- data.frame(
        Variable = factor(uniq_junc_by_cell$SJ_type, levels = sj_levels),
        Value = uniq_junc_by_cell$perc
      )

      gg_rts_unique_by_sjtype <<- build_violin_plot(
        df_long = df_long_unique,
        title = "RT-switching Unique Junctions by Splice Junction Type Across Cells",
        x_labels = sj_labels,
        fill_map = sj_fill_map,
        y_label = "Junctions, %",
        legend = FALSE,
        ylim = c(0, 100),
        violin_alpha = 0.7,
        box_alpha = 0.3,
        box_width = 0.05,
        x_tickangle = 45,
        violin_outline_fill = TRUE,
        box_outline_default = "grey20"
      )
    } else {
      message(
        "RT-switching Unique Junctions plot skipped: need chrom (or chr), genomic_start_coord, ",
        "and genomic_end_coord in junctions file."
      )
    }

  } else {
    message("RTS_junction column not found in Junctions. Skipping RT-switching by SJ type plots.")
  }

  # Create grouped violins for % reads with all canonical junctions by structural category (HTML)
  if (!exists("gg_allcanon_by_category")) {
    # Build cls2: FL-exploded in isoforms mode so sum(count) = transcript counts
    status_map <- function(x) {
      xch <- tolower(as.character(x))
      ifelse(xch %in% c("true", "canonical", "yes"), "True",
        ifelse(xch %in% c("false", "non_canonical", "no"), "False",
          ifelse(is.logical(x) && x, "True",
            ifelse(is.logical(x) && !x, "False", NA_character_)
          )
        )
      )
    }
    if (mode == "isoforms" && "FL" %in% colnames(Classification_file) &&
        "CB" %in% colnames(Classification_file)) {
      cls2 <- Classification_file %>%
        dplyr::filter(!is.na(CB) & CB != "" & CB != "unassigned") %>%
        mutate(
          CB = strsplit(as.character(CB), ","),
          FL = strsplit(as.character(FL), ",")
        ) %>%
        tidyr::unnest(c(CB, FL)) %>%
        mutate(
          CB    = trimws(CB),
          count = suppressWarnings(as.numeric(trimws(FL)))
        ) %>%
        mutate(count = ifelse(is.na(count) | !is.finite(count), 1, count)) %>%
        dplyr::filter(!is.na(CB) & CB != "")
    } else {
      cls2 <- Classification_file %>%
        dplyr::filter(!is.na(CB) & CB != "" & CB != "unassigned") %>%
        mutate(count = 1)
    }
    cls2 <- cls2 %>% mutate(allcanon_status = status_map(all_canonical))

    cat_key_map    <- structural_category_map
    all_cats       <- structural_category_levels
    bin_pretty_map <- c(FSM = "FSM", ISM = "ISM", NIC = "NIC", NNC = "NNC",
                        Genic = "Genic Genomic", Antisense = "Antisense",
                        Fusion = "Fusion", Intergenic = "Intergenic",
                        Genic_intron = "Genic Intron")
    pretty_levels  <- c("FSM", "ISM", "NIC", "NNC", "Genic Genomic",
                        "Antisense", "Fusion", "Intergenic", "Genic Intron")

    df_allcanon <- cls2 %>%
      mutate(cat_key = unname(cat_key_map[structural_category])) %>%
      filter(!is.na(cat_key), !is.na(allcanon_status)) %>%
      group_by(CB, cat_key, allcanon_status) %>%
      summarise(n = sum(count, na.rm = TRUE), .groups = "drop") %>%
      group_by(CB, cat_key) %>%
      mutate(perc = 100 * n / sum(n)) %>%
      ungroup() %>%
      tidyr::complete(CB, cat_key = all_cats, allcanon_status = c("True", "False"), fill = list(n = 0, perc = 0))

    if (nrow(df_allcanon) > 0) {
      df_allcanon$allcanon_status <- factor(df_allcanon$allcanon_status, levels = c("True", "False"))
      cols_tf <- c("True" = "#6baed6", "False" = "#ffc125")
      gg_allcanon_by_category <<- build_grouped_violin_plot(
        df = df_allcanon %>% transmute(bin = unname(bin_pretty_map[as.character(cat_key)]), group = as.character(allcanon_status), value = perc),
        bin_levels = pretty_levels,
        group_levels = c("True", "False"),
        title = paste(entity_label_plural, "with All Canonical Junctions Distribution by Structural Category Across Cells"),
        fill_map = cols_tf,
        legend_labels = c("True" = "True", "False" = "False"),
        y_label = paste(entity_label_plural, ", %", sep = ""),
        ylim = c(0, 100),
        violin_alpha = 0.7,
        box_alpha = 0.3,
        box_width = 0.08,
        x_tickangle = 45,
        violin_width = 0.45,
        dodge_width = 0.8,
        violangap = 0.05,
        violingroupgap = 0.15,
        legend_title = "all_canonical"
      )
    }
  }

  # Coding / Non-Coding / NMD Plots
  if (exists("SQANTI_cell_summary")) {
    # Check if Coding columns exist
    # Coding: ends with "_coding_prop" but NOT "_non_coding_prop"
    all_coding_like <- grep("_coding_prop$", colnames(SQANTI_cell_summary), value = TRUE)
    non_coding_cols <- grep("_non_coding_prop$", colnames(SQANTI_cell_summary), value = TRUE)
    coding_cols <- setdiff(all_coding_like, non_coding_cols)

    # NMD
    if ("NMD_prop_in_cell" %in% colnames(SQANTI_cell_summary)) {
      nmd_cat_cols <- grep(".*_NMD_prop$", colnames(SQANTI_cell_summary), value = TRUE)
      if (length(nmd_cat_cols) > 0) {
        df_nmd <- SQANTI_cell_summary %>%
          select(CB, all_of(nmd_cat_cols)) %>%
          pivot_longer(cols = all_of(nmd_cat_cols), names_to = "Variable", values_to = "Value") %>%
          mutate(
            Variable = gsub("_NMD_prop$", "", Variable)
          )

        # Helper to map tag to pretty label
        tag_to_pretty <- function(tag) {
          case_map <- c(
            "FSM" = "FSM", "ISM" = "ISM", "NIC" = "NIC", "NNC" = "NNC",
            "genic" = "Genic Genomic", "antisense" = "Antisense", "fusion" = "Fusion",
            "intergenic" = "Intergenic", "genic_intron" = "Genic Intron"
          )
          if (tag %in% names(case_map)) {
            return(case_map[[tag]])
          }
          return(tag)
        }

        df_nmd$PrettyVar <- sapply(df_nmd$Variable, tag_to_pretty)

        # Filter to all known categories (including non-canonical)
        all_nmd_vars <- c("FSM", "ISM", "NIC", "NNC", "genic", "antisense", "fusion", "intergenic", "genic_intron")
        df_nmd <- df_nmd %>% filter(Variable %in% all_nmd_vars)

        # Factor levels
        nmd_levels <- c("FSM", "ISM", "NIC", "NNC", "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic Intron")
        df_nmd$PrettyVar <- factor(df_nmd$PrettyVar, levels = nmd_levels)

        # Use grey color for all NMD plots
        nmd_fill_map <- setNames(rep("#969696", length(nmd_levels)), nmd_levels)

        gg_nmd_by_category <<- build_violin_plot(
          df_long = data.frame(Variable = df_nmd$PrettyVar, Value = df_nmd$Value),
          title = paste("Predicted NMD", entity_label_plural, "Distribution by Structural Category Across Cells"),
          x_labels = levels(df_nmd$PrettyVar),
          fill_map = nmd_fill_map,
          y_label = paste(entity_label_plural, ", %", sep = ""),
          legend = FALSE,
          override_outline_vars = character(0),
          violin_alpha = 0.7,
          box_alpha = 0.3,
          box_width = 0.05,
          x_tickangle = 45,
          violin_outline_fill = TRUE,
          box_outline_default = "grey20",
          ylim = c(0, 100)
        )
      }
    }
  }

  ### Generate PDF report ###
  ###########################

  if (generate_pdf) {
    pdf(file.path(paste0(report_output, ".pdf")), paper = "a4r", width = 14, height = 11, useDingbats = FALSE)
    # Cover page
    grid.newpage()
    title_text <- if (mode == "isoforms") "SQANTI-single cell\nisoforms report" else "SQANTI-single cell\nreads report"
    cover <- textGrob(title_text,
      gp = gpar(fontface = "italic", fontsize = 40, col = "orangered")
    )
    grid.draw(cover)
    # Overview tables
    s <- textGrob("Overview", gp = gpar(fontface = "italic", fontsize = 30), vjust = 0)
    grid.arrange(s)

    # Calculate bulk-level stats
    if (mode == "isoforms") {
      total_reads_count <- sum(Classification_file$count, na.rm = TRUE)
      unique_isoforms <- nrow(Classification_file)
    } else {
      total_reads_count <- nrow(Classification_file)
    }
    unique_genes <- length(unique(Classification_file$associated_gene))
    if (mode == "isoforms") {
      unique_junctions <- 0
    } else {
      unique_junctions <- length(unique(Classification_file$jxn_string))
    }

    annotated_genes <- length(unique(Classification_file$associated_gene[!grepl("^novel", Classification_file$associated_gene)]))
    novel_genes <- length(unique(Classification_file$associated_gene[grepl("^novel", Classification_file$associated_gene)]))
    gene_class_table <- data.frame(
      Category = c("Annotated Genes", "Novel Genes"),
      `Genes, count` = c(annotated_genes, novel_genes),
      `Genes, percent` = c(
        round(100 * annotated_genes / unique_genes, 1),
        round(100 * novel_genes / unique_genes, 1)
      ),
      check.names = FALSE
    )

    # Read Classification table (counts per structural category)
    read_cat_levels <- c(
      "full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog",
      "genic", "antisense", "fusion", "intergenic", "genic_intron"
    )
    read_cat_names <- c(
      "FSM", "ISM", "NIC", "NNC",
      "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic Intron"
    )

    if (mode == "isoforms") {
      read_class_table <- aggregate(Classification_file$count, by = list(Category = factor(Classification_file$structural_category, levels = read_cat_levels)), FUN = sum, na.rm = TRUE)
      colnames(read_class_table) <- c("Category", "Transcripts, count")
      # Ensure all levels are present
      read_class_temp <- data.frame(Category = read_cat_levels)
      read_class_table <- merge(read_class_temp, read_class_table, by = "Category", all.x = TRUE)
      read_class_table <- read_class_table[match(read_cat_levels, read_class_table$Category), ]
      read_class_table[is.na(read_class_table)] <- 0
    } else {
      read_class_table <- as.data.frame(table(factor(Classification_file$structural_category, levels = read_cat_levels)))
      colnames(read_class_table) <- c("Category", "Transcripts, count")
    }
    read_class_table$Category <- read_cat_names
    read_class_table[["Transcripts, percent"]] <- round(100 * read_class_table[["Transcripts, count"]] / sum(read_class_table[["Transcripts, count"]]), 1)

    # Splice Junction Classification table
    Junctions$junction_type <- paste(Junctions$junction_category, Junctions$canonical, sep = "_")

    sj_types <- c("known_canonical", "known_non_canonical", "novel_canonical", "novel_non_canonical")
    if (mode == "isoforms") {
      sj_counts <- sapply(sj_types, function(type) sum(Junctions$count[Junctions$junction_type == type], na.rm = TRUE))
    } else {
      sj_counts <- sapply(sj_types, function(type) sum(Junctions$junction_type == type, na.rm = TRUE))
    }

    # Handle case where there are no junctions
    total_junctions <- sum(sj_counts, na.rm = TRUE)
    sj_perc <- if (total_junctions > 0) {
      round(100 * sj_counts / total_junctions, 2)
    } else {
      rep(0, length(sj_counts))
    }

    SJ_class_table <- data.frame(
      Category = c("Known canonical", "Known Non-canonical", "Novel canonical", "Novel Non-canonical"),
      `SJs, count` = sj_counts,
      `SJs, percent` = sj_perc,
      check.names = FALSE
    )
    rownames(SJ_class_table) <- NULL

    big_table_theme <- ttheme_default(
      core = list(fg_params = list(cex = 1.1)),
      colhead = list(fg_params = list(cex = 1.1, fontface = "bold"))
    )

    title_genes <- textGrob("Gene Classification", gp = gpar(fontface = "italic", fontsize = 20), vjust = -3)
    title_reads <- textGrob(paste(entity_label, "Classification"), gp = gpar(fontface = "italic", fontsize = 20), vjust = -7.7)
    title_sj <- textGrob("Splice Junction Classification", gp = gpar(fontface = "italic", fontsize = 20), vjust = -4.3)

    table_genes <- tableGrob(gene_class_table, rows = NULL, theme = big_table_theme)
    table_reads <- tableGrob(read_class_table, rows = NULL, theme = big_table_theme)
    table_sj <- tableGrob(SJ_class_table, rows = NULL, theme = big_table_theme)

    if (mode == "isoforms") {
      unique_counts_text <- sprintf(
        "Number of %s: %s\nUnique Isoforms: %s\nUnique Genes: %s",
        entity_label_plural,
        format(total_reads_count, big.mark = ","),
        format(unique_isoforms, big.mark = ","),
        format(unique_genes, big.mark = ",")
      )
    } else if (unique_junctions > 0) {
      unique_counts_text <- sprintf(
        "Number of %s: %s\nUnique Genes: %s\nUnique Junction Chains: %s",
        entity_label_plural,
        format(total_reads_count, big.mark = ","),
        format(unique_genes, big.mark = ","),
        format(unique_junctions, big.mark = ",")
      )
    } else {
      unique_counts_text <- sprintf(
        "Number of %s: %s\nUnique Genes: %s",
        entity_label_plural,
        format(total_reads_count, big.mark = ","),
        format(unique_genes, big.mark = ",")
      )
    }
    unique_counts_grob <- textGrob(
      unique_counts_text,
      gp = gpar(fontface = "italic", fontsize = 28), vjust = 0, hjust = 0.5
    )

    # Create gTree objects to overlay titles and tables
    gt_genes <- gTree(children = gList(table_genes, title_genes))
    gt_reads <- gTree(children = gList(table_reads, title_reads))
    gt_sj <- gTree(children = gList(table_sj, title_sj))

    # Arrange left column: Gene Classification + Splice Junction Classification
    left_col <- arrangeGrob(
      gt_genes,
      gt_sj,
      ncol = 1,
      heights = c(0.2, 0.4)
    )

    # Arrange right column: Read Classification
    right_col <- arrangeGrob(
      gt_reads,
      ncol = 1
    )

    # Final page layout
    grid.arrange(
      unique_counts_grob,
      arrangeGrob(left_col, right_col, ncol = 2, widths = c(1.3, 1.3)),
      nrow = 2,
      heights = c(0.8, 1)
    )

    # Single cell tables
    s <- textGrob("Cell summary", gp = gpar(fontface = "italic", fontsize = 30), vjust = 0)
    grid.arrange(s)

    # Number of cells
    num_cells <- nrow(SQANTI_cell_summary)
    num_cells_grob <- textGrob(
      sprintf("Unique Cell Barcodes: %d", num_cells),
      gp = gpar(fontface = "italic", fontsize = 28), vjust = 0.5, hjust = 0.5
    )

    # 1. Unique Genes and Unique Junction Chains summary table
    unique_genes_stats <- c(
      Mean = mean(SQANTI_cell_summary$Genes_in_cell, na.rm = TRUE),
      Median = median(SQANTI_cell_summary$Genes_in_cell, na.rm = TRUE),
      Min = min(SQANTI_cell_summary$Genes_in_cell, na.rm = TRUE),
      Max = max(SQANTI_cell_summary$Genes_in_cell, na.rm = TRUE),
      IQR = IQR(SQANTI_cell_summary$Genes_in_cell, na.rm = TRUE),
      SD = sd(SQANTI_cell_summary$Genes_in_cell, na.rm = TRUE)
    )
    unique_junctions_stats <- c(
      Mean = mean(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE),
      Median = median(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE),
      Min = min(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE),
      Max = max(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE),
      IQR = IQR(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE),
      SD = sd(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE)
    )
    reads_stats <- c(
      Mean = mean(SQANTI_cell_summary[[count_col]], na.rm = TRUE),
      Median = median(SQANTI_cell_summary[[count_col]], na.rm = TRUE),
      Min = min(SQANTI_cell_summary[[count_col]], na.rm = TRUE),
      Max = max(SQANTI_cell_summary[[count_col]], na.rm = TRUE),
      IQR = IQR(SQANTI_cell_summary[[count_col]], na.rm = TRUE),
      SD = sd(SQANTI_cell_summary[[count_col]], na.rm = TRUE)
    )
    umis_stats <- c(
      Mean = mean(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE),
      Median = median(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE),
      Min = min(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE),
      Max = max(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE),
      IQR = IQR(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE),
      SD = sd(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE)
    )
    summary_table1 <- data.frame(
      Feature = c(paste(entity_label_plural, "in cell"), "UMIs in cell", "Unique Genes", "Unique Junction Chains"),
      Mean = c(reads_stats["Mean"], umis_stats["Mean"], unique_genes_stats["Mean"], unique_junctions_stats["Mean"]),
      Median = c(reads_stats["Median"], umis_stats["Median"], unique_genes_stats["Median"], unique_junctions_stats["Median"]),
      Min = c(reads_stats["Min"], umis_stats["Min"], unique_genes_stats["Min"], unique_junctions_stats["Min"]),
      Max = c(reads_stats["Max"], umis_stats["Max"], unique_genes_stats["Max"], unique_junctions_stats["Max"]),
      IQR = c(reads_stats["IQR"], umis_stats["IQR"], unique_genes_stats["IQR"], unique_junctions_stats["IQR"]),
      SD = c(reads_stats["SD"], umis_stats["SD"], unique_genes_stats["SD"], unique_junctions_stats["SD"])
    )
    # If isoforms mode, drop Unique Junction Chains from summary table
    if (mode == "isoforms" || !("UJCs_in_cell" %in% names(SQANTI_cell_summary)) || all(is.na(SQANTI_cell_summary$UJCs_in_cell)) || max(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE) == 0) {
      summary_table1 <- summary_table1[!(summary_table1$Feature %in% c("Unique Junction Chains", "UMIs in cell")), , drop = FALSE]
    }
    summary_table1[, 2:7] <- round(summary_table1[, 2:7], 3)
    table_summary1 <- tableGrob(summary_table1, rows = NULL, theme = big_table_theme)
    gt_summary1 <- gTree(children = gList(table_summary1))

    # 2. Gene Classification summary table (across all cells)
    gene_class_stats <- data.frame(
      Category = c("Annotated Genes", "Novel Genes"),
      Mean = c(mean(SQANTI_cell_summary$Annotated_genes, na.rm = TRUE), mean(SQANTI_cell_summary$Novel_genes, na.rm = TRUE)),
      Median = c(median(SQANTI_cell_summary$Annotated_genes, na.rm = TRUE), median(SQANTI_cell_summary$Novel_genes, na.rm = TRUE)),
      Min = c(min(SQANTI_cell_summary$Annotated_genes, na.rm = TRUE), min(SQANTI_cell_summary$Novel_genes, na.rm = TRUE)),
      Max = c(max(SQANTI_cell_summary$Annotated_genes, na.rm = TRUE), max(SQANTI_cell_summary$Novel_genes, na.rm = TRUE)),
      SD = c(sd(SQANTI_cell_summary$Annotated_genes, na.rm = TRUE), sd(SQANTI_cell_summary$Novel_genes, na.rm = TRUE))
    )
    gene_class_stats[, 2:6] <- round(gene_class_stats[, 2:6], 3)
    table_gene_class_stats <- tableGrob(gene_class_stats, rows = NULL, theme = big_table_theme)
    title_gene_class_stats <- textGrob("Gene Classification (per cell)", gp = gpar(fontface = "italic", fontsize = 22), vjust = -2.9)
    gt_gene_class_stats <- gTree(children = gList(table_gene_class_stats, title_gene_class_stats))

    # 3. Splice Junction Classification summary table (across all cells)

    # Create a junction type column for easier summarization
    Junctions$junction_type <- paste(Junctions$junction_category, Junctions$canonical, sep = "_")

    # Calculate proportions of each junction type per cell
    junction_proportions_per_cell <- Junctions %>%
      filter(CB != "unassigned") %>%
      group_by(CB) %>%
      summarise(
        Known_canonical = sum(junction_type == "known_canonical", na.rm = TRUE) / n() * 100,
        Known_Non_canonical = sum(junction_type == "known_non_canonical", na.rm = TRUE) / n() * 100,
        Novel_canonical = sum(junction_type == "novel_canonical", na.rm = TRUE) / n() * 100,
        Novel_Non_canonical = sum(junction_type == "novel_non_canonical", na.rm = TRUE) / n() * 100,
        .groups = "drop"
      )

    # Calculate summary statistics across all cells
    sj_stats <- junction_proportions_per_cell %>%
      select(-CB) %>%
      summarise(
        across(
          everything(),
          list(
            Mean = ~ mean(.x, na.rm = TRUE),
            Median = ~ median(.x, na.rm = TRUE),
            Min = ~ min(.x, na.rm = TRUE),
            Max = ~ max(.x, na.rm = TRUE),
            SD = ~ sd(.x, na.rm = TRUE)
          )
        )
      )

    # Reshape the data for display
    sj_stats_df <- sj_stats %>%
      pivot_longer(
        cols = everything(),
        names_to = c("Category", ".value"),
        names_pattern = "(.+)_(Mean|Median|Min|Max|SD)$"
      ) %>%
      mutate(Category = gsub("_", " ", Category))

    # Ensure we have the expected number of columns before subsetting
    if (ncol(sj_stats_df) >= 6) {
      sj_stats_df[, 2:6] <- round(sj_stats_df[, 2:6], 3)
    } else {
      # If we have fewer columns, round all numeric columns except the first (Category)
      numeric_cols <- sapply(sj_stats_df[, -1], is.numeric)
      sj_stats_df[, -1][numeric_cols] <- round(sj_stats_df[, -1][numeric_cols], 3)
    }
    table_sj_stats <- tableGrob(sj_stats_df, rows = NULL, theme = big_table_theme)
    title_sj_stats <- textGrob("Splice Junction Classification (per cell, %)", gp = gpar(fontface = "italic", fontsize = 22), vjust = -4.4)
    gt_sj_stats <- gTree(children = gList(table_sj_stats, title_sj_stats))

    grid.arrange(
      num_cells_grob,
      gt_summary1,
      gt_gene_class_stats,
      gt_sj_stats,
      ncol = 1,
      heights = c(0.3, 1, 0.7, 0.9)
    )

    # Cell Summary Statistics Page 2: Read Classification
    title_read_class <- textGrob(paste(entity_label, "Classification"), gp = gpar(fontface = "italic", fontsize = 28), vjust = 0, hjust = 0.5)
    desc_counts <- textGrob(paste("Summary of per cell", entity_label_lower, "counts by structural category"), gp = gpar(fontface = "italic", fontsize = 18), vjust = 0.5)
    desc_props <- textGrob(paste("Summary of per cell", entity_label_lower, "percentages by structural category"), gp = gpar(fontface = "italic", fontsize = 18), vjust = 0.5)
    struct_cat_cols <- c(
      "FSM", "ISM", "NIC", "NNC", "Genic_Genomic", "Antisense", "Fusion", "Intergenic", "Genic_intron"
    )
    struct_cat_names <- c(
      "FSM", "ISM", "NIC", "NNC",
      "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic Intron"
    )

    # Smaller table theme for these two tables
    small_table_theme <- ttheme_default(
      core = list(fg_params = list(cex = 1.2)),
      colhead = list(fg_params = list(cex = 1.2, fontface = "bold"))
    )

    # 1. Counts summary table
    count_stats <- sapply(struct_cat_cols, function(col) {
      vals <- SQANTI_cell_summary[[col]]
      c(
        Mean = mean(vals, na.rm = TRUE),
        Median = median(vals, na.rm = TRUE),
        Min = min(vals, na.rm = TRUE),
        Max = max(vals, na.rm = TRUE),
        IQR = IQR(vals, na.rm = TRUE),
        SD = sd(vals, na.rm = TRUE)
      )
    })
    count_stats_df <- data.frame(
      Category = struct_cat_names,
      t(count_stats)
    )
    colnames(count_stats_df)[2:7] <- c("Mean", "Median", "Min", "Max", "IQR", "SD")
    count_stats_df[, 2:7] <- round(count_stats_df[, 2:7], 3)
    table_count_stats <- tableGrob(count_stats_df, rows = NULL, theme = small_table_theme)

    # 2. Proportions summary table
    prop_cat_cols <- paste0(struct_cat_cols, "_prop")
    prop_stats <- sapply(prop_cat_cols, function(col) {
      vals <- SQANTI_cell_summary[[col]]
      c(
        Mean = mean(vals, na.rm = TRUE),
        Median = median(vals, na.rm = TRUE),
        Min = min(vals, na.rm = TRUE),
        Max = max(vals, na.rm = TRUE),
        IQR = IQR(vals, na.rm = TRUE),
        SD = sd(vals, na.rm = TRUE)
      )
    })
    prop_stats_df <- data.frame(
      Category = struct_cat_names,
      t(prop_stats)
    )
    colnames(prop_stats_df)[2:7] <- c("Mean", "Median", "Min", "Max", "IQR", "SD")
    prop_stats_df[, 2:7] <- round(prop_stats_df[, 2:7], 3)
    table_prop_stats <- tableGrob(prop_stats_df, rows = NULL, theme = small_table_theme)

    grid.arrange(
      title_read_class,
      desc_counts,
      table_count_stats,
      desc_props,
      table_prop_stats,
      ncol = 1,
      heights = c(0.3, 0.12, 1, 0.12, 1)
    )

    # Helper for section title pages
    section_page <- function(title) {
      grid.newpage()
      grid.draw(textGrob(title, gp = gpar(fontface = "italic", fontsize = 30, col = "black")))
    }

    # Per-cell Library Size section
    section_page("Per-cell Library Size")
    render_pdf_plot_centered("gg_reads_in_cells", width_frac = 0.5)
    render_pdf_plot_centered("gg_umis_in_cells", width_frac = 0.5)
    if (exists("gg_JCs_in_cell")) render_pdf_plot_centered("gg_JCs_in_cell", width_frac = 0.5)

    # Gene Characterization section
    section_page("Gene Characterization")
    # Genes Across Cells
    render_pdf_plot_centered("gg_genes_in_cells", width_frac = 0.5)
    render_pdf_plot("gg_annotation_of_genes_in_cell")
    render_pdf_plot("gg_annotation_of_genes_percent_in_cell")
    # Reads per Gene
    render_pdf_plot("gg_annotation_of_reads_in_cell")
    render_pdf_plot("gg_read_bins_all")
    render_pdf_plot("gg_read_bins")
    if (mode == "isoforms" && exists("gg_isoform_bins")) {
      render_pdf_plot("gg_isoform_bins")
    }
    # UJCs per Gene
    if (mode != "isoforms" && exists("gg_ujc_bins_all")) {
      render_pdf_plot("gg_ujc_bins_all")
      render_pdf_plot("gg_ujc_bins")
    }
    # Mitochondrial genes
    render_pdf_plot("gg_MT_perc")

    # Read Length Characterization section
    section_page(paste(entity_label, "Length Characterization"))
    # Bulk Length Distribution
    render_pdf_plot("gg_bulk_all_reads")
    render_pdf_plot("gg_bulk_length_by_category")
    render_pdf_plot("gg_bulk_length_by_exon_type")
    # Overall cell-level distributions: All then Mono on next page
    render_pdf_plot("gg_read_distr")
    render_pdf_plot("gg_read_distr_mono")
    # Category-specific: print All then Mono right after
    for (tag in c("FSM", "ISM", "NIC", "NNC", "genic", "antisense", "fusion", "intergenic", "genic_intron")) {
      all_nm <- paste0("gg_", tag, "_read_distr")
      mono_nm <- paste0("gg_", tag, "_mono_read_distr")
      if (exists(all_nm)) render_pdf_plot(all_nm)
      if (exists(mono_nm)) render_pdf_plot(mono_nm)
    }
    # Reference Transcript Coverage
    render_pdf_plot("gg_ref_coverage_across_category")
    render_pdf_plot("gg_meta_transcript_coverage")
    if (exists("gg_isoforms_ref_vs_sample_lengths") && !is.null(gg_isoforms_ref_vs_sample_lengths)) {
      render_pdf_plot("gg_isoforms_ref_vs_sample_lengths")
    }

    # Structural Read Characterization section
    section_page(paste("Structural", entity_label, "Characterization"))
    # Distribution by Structural Categories
    render_pdf_plot("gg_SQANTI_across_category")
    render_pdf_plot("gg_exon_mono_by_category")
    for (nm in c(
      "gg_SQANTI_across_FSM", "gg_SQANTI_across_ISM", "gg_SQANTI_across_NIC", "gg_SQANTI_across_NNC",
      "gg_SQANTI_across_Genic", "gg_SQANTI_across_Antisense", "gg_SQANTI_across_Fusion", "gg_SQANTI_across_Intergenic",
      "gg_SQANTI_across_Genic_Intron"
    )) {
      render_pdf_plot(nm)
    }
    # Exon Counts
    render_pdf_plot("gg_exon_mean_by_category")
    if (exists("gg_exon_profile_by_category")) {
      prof_order <- c("FSM", "ISM", "NIC", "NNC", "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic Intron")
      for (nm in prof_order) {
        if (!is.null(gg_exon_profile_by_category[[nm]])) {
          print(gg_exon_profile_by_category[[nm]])
        }
      }
    }

    # Coding and Non-Coding Distributions
    if (include_ORF) {
      if (exists("gg_coding_across_category")) render_pdf_plot("gg_coding_across_category")
      if (exists("gg_non_coding_across_category")) render_pdf_plot("gg_non_coding_across_category")
    }

    # Splice Junction Characterization section
    # Compute per-structural-category SJ distributions for PDF pages

    # Find the common isoform ID column between Junctions and Classification_file
    junc_iso_key <- NULL
    for (k in c("isoform", "readID", "read_id", "ID", "read_name", "read")) {
      if (k %in% colnames(Junctions) && k %in% colnames(Classification_file)) {
        junc_iso_key <- k
        break
      }
    }

    if (mode == "isoforms" && !is.null(junc_iso_key) &&
        "FL" %in% colnames(Classification_file) && "CB" %in% colnames(Classification_file)) {
      # Explode Classification_file CB/FL to get one row per (isoform, cell, FL_count)
      cls_exploded <- Classification_file %>%
        filter(!is.na(CB), !is.na(structural_category)) %>%
        select(all_of(c(junc_iso_key, "CB", "FL", "structural_category"))) %>%
        mutate(CB_raw = as.character(CB), FL_raw = as.character(FL)) %>%
        tidyr::separate_rows(CB_raw, FL_raw, sep = ",") %>%
        mutate(
          CB_clean = trimws(CB_raw),
          FL_num   = suppressWarnings(as.numeric(trimws(FL_raw)))
        ) %>%
        filter(CB_clean != "", CB_clean != "unassigned", !is.na(FL_num), FL_num > 0) %>%
        select(all_of(c(junc_iso_key, "CB_clean", "FL_num", "structural_category")))

      # Join junctions to exploded classification by isoform key only.
      # Select only junction-type columns from Junctions first to avoid
      # column name collisions with the CB/count columns in cls_exploded.
      junc_aug <- Junctions %>%
        select(all_of(c(junc_iso_key, "junction_category", "canonical"))) %>%
        mutate(junction_type = paste(junction_category, canonical, sep = "_")) %>%
        inner_join(cls_exploded, by = junc_iso_key) %>%
        rename(CB = CB_clean, count = FL_num)

    } else {
      # Reads mode: each junction row already has a single CB; count = 1
      junc_aug <- Junctions %>%
        dplyr::filter(!is.na(CB) & CB != "" & CB != "unassigned") %>%
        mutate(junction_type = paste(junction_category, canonical, sep = "_"),
               count = 1)

      # Bring in structural_category if not already present
      if (!("structural_category" %in% colnames(junc_aug)) && !is.null(junc_iso_key)) {
        junc_aug <- junc_aug %>%
          left_join(
            Classification_file %>% select(all_of(c(junc_iso_key, "structural_category"))),
            by = junc_iso_key
          )
      }
    }

    cat_key_map <- structural_category_map
    all_cats <- structural_category_levels

    junc_summ <- junc_aug %>%
      filter(!is.na(structural_category)) %>%
      mutate(cat_key = unname(cat_key_map[structural_category])) %>%
      filter(!is.na(cat_key)) %>%
      group_by(CB, cat_key) %>%
      summarise(
        total = sum(count, na.rm = TRUE),
        known_canonical = sum(count[junction_type == "known_canonical"], na.rm = TRUE),
        known_non_canonical = sum(count[junction_type == "known_non_canonical"], na.rm = TRUE),
        novel_canonical = sum(count[junction_type == "novel_canonical"], na.rm = TRUE),
        novel_non_canonical = sum(count[junction_type == "novel_non_canonical"], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      # Do NOT fill absent (CB, cat_key) pairs with 0 — cells with no transcripts
      # of a given category should be NA (excluded from violin), not 0%.
      mutate(
        KnownCanonicalPerc    = ifelse(total > 0, 100 * known_canonical    / total, NA_real_),
        KnownNonCanonicalPerc = ifelse(total > 0, 100 * known_non_canonical / total, NA_real_),
        NovelCanonicalPerc    = ifelse(total > 0, 100 * novel_canonical     / total, NA_real_),
        NovelNonCanonicalPerc = ifelse(total > 0, 100 * novel_non_canonical / total, NA_real_)
      ) %>%
      # Expand to all categories so every category column exists for plotting,
      # but keep NA for cells absent from that category
      tidyr::complete(CB, cat_key = all_cats) %>%
      ungroup()

    # Prepare plotting helpers
    fill_map_cat <- cat_fill_map
    x_labels_full <- cat_labels_pretty
    make_df_long <- function(col_name) {
      data.frame(Variable = factor(junc_summ$cat_key, levels = all_cats), Value = junc_summ[[col_name]])
    }

    p_known_can <- build_violin_plot(
      df_long = make_df_long("KnownCanonicalPerc"),
      title = "Known Canonical Splice Junctions Distribution by Structural Category Across Cells",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )
    p_known_noncan <- build_violin_plot(
      df_long = make_df_long("KnownNonCanonicalPerc"),
      title = "Known Non-canonical Splice Junctions Distribution by Structural Category Across Cells",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )
    p_novel_can <- build_violin_plot(
      df_long = make_df_long("NovelCanonicalPerc"),
      title = "Novel Canonical Splice Junctions Distribution by Structural Category Across Cells",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )
    p_novel_noncan <- build_violin_plot(
      df_long = make_df_long("NovelNonCanonicalPerc"),
      title = "Novel Non-canonical Splice Junctions Distribution by Structural Category Across Cells",
      x_labels = x_labels_full,
      fill_map = fill_map_cat,
      y_label = "Junctions, %",
      legend = FALSE,
      override_outline_vars = c("Genic"),
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20",
      ylim = c(0, 100)
    )

    section_page("Splice Junction Characterization")
    render_pdf_plot("gg_known_novel_canon")
    print(p_known_can)
    print(p_known_noncan)
    print(p_novel_can)
    print(p_novel_noncan)
    render_pdf_plot("gg_allcanon_by_category")
    if (exists("gg_rts_all_by_sjtype")) render_pdf_plot("gg_rts_all_by_sjtype")
    if (exists("gg_rts_unique_by_sjtype")) render_pdf_plot("gg_rts_unique_by_sjtype")

    # Features of Bad Quality section
    section_page("Features of Bad Quality")
    render_pdf_plot("gg_bad_feature")
    render_pdf_plot("gg_intrapriming_by_category")
    render_pdf_plot("gg_RTS_by_category")

    render_pdf_plot("gg_noncanon_by_category")
    if (exists("gg_nmd_by_category")) render_pdf_plot("gg_nmd_by_category")

    # Features of Good Quality section
    section_page("Features of Good Quality")
    render_pdf_plot("gg_good_feature")
    render_pdf_plot("gg_tss_annotation_support")
    if (CAGE_peak) render_pdf_plot("gg_cage_peak_support")
    if (polyA_motif_list) render_pdf_plot("gg_polyA_motif_support")
    render_pdf_plot("gg_canon_by_category")
    if (exists("gg_sr_support_by_category")) render_pdf_plot("gg_sr_support_by_category")
    if (exists("gg_tss_validation_by_category")) render_pdf_plot("gg_tss_validation_by_category")

    # Clustering Analysis
    if (exists("gg_umap") && !is.null(gg_umap)) {
      section_page("Clustering analysis")
      print(gg_umap)

      # Print UMAP by structural category if available (one per page)
      if (exists("gg_umap_by_category") && !is.null(gg_umap_by_category)) {
        for (cat_label in names(gg_umap_by_category)) {
          # Print the UMAP
          print(gg_umap_by_category[[cat_label]])
        }
        
        # Then print out ALL the corresponding Structural Categories Violin Plots
        if (exists("gg_cat_cluster_plots") && !is.null(gg_cat_cluster_plots)) {
           for (cat_label in names(gg_cat_cluster_plots)) {
              print(gg_cat_cluster_plots[[cat_label]])
           }
        }

        # Length Distribution by Cluster (PDF)
        if (exists("gg_len_cluster_plots") && !is.null(gg_len_cluster_plots)) {
          for (len_lbl in names(gg_len_cluster_plots)) {
            if (!is.null(gg_len_cluster_plots[[len_lbl]])) {
              print(gg_len_cluster_plots[[len_lbl]])
            }
          }
        }
      }

      # Print Short Read Support by Cluster (Violin + UMAPs)
      if (exists("gg_sr_cluster_plots") && !is.null(gg_sr_cluster_plots)) {
        # UMAPs first (Global then Category-specific if available)
        if (exists("gg_sr_umap_plots") && !is.null(gg_sr_umap_plots)) {
          if (!is.null(gg_sr_umap_plots[["All Transcripts"]])) print(gg_sr_umap_plots[["All Transcripts"]])
          for (label in setdiff(names(gg_sr_umap_plots), "All Transcripts")) {
            print(gg_sr_umap_plots[[label]])
          }
        }

        # Violin plots
        for (label in names(gg_sr_cluster_plots)) {
          # Print ggplot for PDF
          p_ggplot <- gg_sr_cluster_plots[[label]]
          print(p_ggplot)
        }
      }

      # Print TSS Validation Support by Cluster (Violin + UMAPs)
      if (exists("gg_tss_cluster_plots") && !is.null(gg_tss_cluster_plots)) {
        # UMAPs first
        if (exists("gg_tss_umap_plots") && !is.null(gg_tss_umap_plots)) {
          if (!is.null(gg_tss_umap_plots[["All Transcripts"]])) print(gg_tss_umap_plots[["All Transcripts"]])
          for (label in setdiff(names(gg_tss_umap_plots), "All Transcripts")) {
            print(gg_tss_umap_plots[[label]])
          }
        }

        # Violin plots
        for (label in names(gg_tss_cluster_plots)) {
          # Print ggplot for PDF
          p_ggplot <- gg_tss_cluster_plots[[label]]
          print(p_ggplot)
        }
      }
    }

    dev.off()
  }
}

Classification <- data.table::fread(class.file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)
if (mode == "isoforms" && "FL" %in% colnames(Classification)) {
  Classification$count <- sapply(strsplit(as.character(Classification$FL), ","), function(x) sum(as.numeric(x), na.rm = TRUE))
  Classification$count[is.na(Classification$count) | Classification$count == 0] <- 1
} else {
  Classification$count <- 1
}
Junctions <- data.table::fread(junc.file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)

# Add count column to Junctions for weighted quantification
if (mode == "isoforms") {
  # Try to join by isoform ID
  join_key <- NULL
  for (k in c("isoform", "readID", "read_id", "ID", "read_name", "read")) {
    if (k %in% colnames(Junctions) && k %in% colnames(Classification)) {
      join_key <- k
      break
    }
  }

  if (!is.null(join_key)) {
    # Assign isoform count to each junction row

    # Use match to be faster than merge/join for simple lookup
    Junctions$count <- Classification$count[match(Junctions[[join_key]], Classification[[join_key]])]
    # Handle NAs (shouldn't happen if consistent)
    Junctions$count[is.na(Junctions$count)] <- 1
  } else {
    Junctions$count <- 1
  }
} else {
  Junctions$count <- 1
}

# Require precomputed cell summary produced by sqanti_sc.py
if (!is.null(cell_summary_path) && file.exists(cell_summary_path)) {
  SQANTI_cell_summary <- data.table::fread(cell_summary_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)
} else {
  stop("A precomputed cell summary is required. Pass --cell_summary <path> from sqanti_sc.py.")
}

# Generate reports based on format
if (report.format == "pdf" || report.format == "both") {
  generate_sqantisc_plots(
    SQANTI_cell_summary,
    Classification,
    Junctions,
    report_output
  )
}

if (report.format == "html" || report.format == "both") {
  # Generate plots first so they're available for Rmd
  if (report.format == "html") {
    generate_sqantisc_plots(
      SQANTI_cell_summary,
      Classification,
      Junctions,
      report_output,
      generate_pdf = FALSE
    )
  }
  # Set up HTML report generation
  # Get the directory where this R script is located (utilities folder)
  cmd_args <- commandArgs(trailingOnly = FALSE)
  script_arg <- cmd_args[grep("--file=", cmd_args)]
  if (length(script_arg) > 0) {
    script_path <- gsub("--file=", "", script_arg)
    script_dir <- dirname(normalizePath(script_path))
  } else {
    stop("Cannot determine script location")
  }

  rmd_file <- file.path(script_dir, "SQANTI-sc_report.Rmd")
  css_file <- file.path(script_dir, "style.css")

  # Check if Rmd file exists
  if (!file.exists(rmd_file)) {
    stop(
      "R Markdown file not found: ", rmd_file,
      "\nPlease ensure SQANTI-sc_report.Rmd is in the same directory as this script."
    )
  }

  # Copy CSS file to output directory if it exists
  if (file.exists(css_file)) {
    css_output <- file.path(dirname(report_output), "style.css")
    file.copy(css_file, css_output, overwrite = TRUE)
  }

  # Generate HTML report
  html_output_file <- paste0(report_output, ".html")

  message("Generating HTML report...")

  rmarkdown::render(
    input = rmd_file,
    output_file = html_output_file,
    output_dir = dirname(report_output),
    intermediates_dir = dirname(report_output),
    quiet = FALSE,
    envir = globalenv()
  )

  # Cleanup: remove the copied CSS file
  if (exists("css_output") && file.exists(css_output)) {
    file.remove(css_output)
  }

  message("HTML report generated: ", html_output_file)
}
