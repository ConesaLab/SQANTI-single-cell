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
  library(plotly)
  library(scales)
}))

#********************** Taking arguments from python script

args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1]
junc.file <- args[2]
report.format <- args[3]
outputPathPrefix <- args[4]
mode <- args[5]

# Initialize ignore_cell_summary flag
ignore_cell_summary <- FALSE
skipORF <- FALSE
CAGE_peak <- FALSE
polyA_motif_list <- FALSE
cell_summary_path <- NULL

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
    if (arg == "--skipORF") {
      skipORF <- TRUE
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


generate_sqantisc_plots <- function(SQANTI_cell_summary, Classification_file, Junctions, report_output, generate_pdf = TRUE) {
  # Helper: convert any R color (hex or named) to an rgba() string with alpha without affecting line color
  to_rgba <- function(col, alpha = 1.0) {
    rgb <- grDevices::col2rgb(col)
    sprintf("rgba(%d,%d,%d,%.3f)", rgb[1], rgb[2], rgb[3], alpha)
  }

  # Helper: pivot selected columns to long and return factor-ordered long df
  pivot_long <- function(df, cols) {
    out <- pivot_longer(df, cols = all_of(cols), names_to = "Variable", values_to = "Value") %>%
      select(Variable, Value)
    out$Variable <- factor(out$Variable, levels = cols)
    out
  }

  # Helper: generic violin + box + mean-cross plot with shared theme (ggplot version for PDF)
  build_violin_plot_ggplot <- function(df_long,
                                       title,
                                       x_labels,
                                       fill_map,
                                       color_map = fill_map,
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
                                       bandwidth = NULL) {
    # Determine a robust bandwidth for KDE; floor to avoid bw=0 on constant data
    vals <- df_long$Value
    vals <- vals[is.finite(vals)]
    bw_eff <- bandwidth
    if (is.null(bw_eff) || !is.numeric(bw_eff) || is.na(bw_eff) || bw_eff <= 0) {
      if (length(vals) >= 2) {
        bw_eff <- stats::bw.nrd0(vals)
      } else {
        bw_eff <- NA_real_
      }
    }
    if (is.na(bw_eff) || bw_eff <= 0) bw_eff <- 0.1

    # Create ggplot version for PDF output
    p <- ggplot(df_long, aes(x = Variable, y = Value)) +
      # Violin layer with outline rule
      {
        if (isTRUE(violin_outline_fill)) {
          geom_violin(aes(fill = Variable, color = Variable), alpha = violin_alpha, scale = "width", show.legend = legend, bw = bw_eff, trim = TRUE)
        } else {
          geom_violin(aes(fill = Variable), color = "black", alpha = violin_alpha, scale = "width", show.legend = legend, bw = bw_eff, trim = TRUE)
        }
      } +
      scale_fill_manual(values = fill_map, labels = x_labels) +
      {
        if (isTRUE(violin_outline_fill)) scale_color_manual(values = fill_map, guide = "none") else NULL
      } +
      # Add mean markers on top
      stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1, show.legend = FALSE) +
      scale_x_discrete(labels = x_labels) +
      labs(title = title, x = "", y = y_label) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = x_tickangle, hjust = ifelse(x_tickangle == 0, 0.5, 1)),
        legend.position = if (legend) "bottom" else "none"
      )

    # Add boxplots per variable with correct outline color (grey90 overrides)
    for (var in levels(df_long$Variable)) {
      var_df <- df_long[df_long$Variable == var, , drop = FALSE]
      box_col <- if (var %in% override_outline_vars) "grey90" else box_outline_default
      p <- p + geom_boxplot(
        data = var_df,
        aes(x = Variable, y = Value, fill = Variable),
        width = box_width, outlier.shape = NA, alpha = box_alpha, show.legend = FALSE, color = box_col
      )
    }

    if (!is.null(ylim)) {
      p <- p + coord_cartesian(ylim = ylim)
    }

    return(p)
  }

  # Helper: generic violin + box + mean-cross plot with shared theme
  build_violin_plot <- function(df_long,
                                title,
                                x_labels,
                                fill_map,
                                color_map = fill_map,
                                y_label = paste(entity_label_plural, ", %", sep = ""),
                                legend = FALSE,
                                ylim = NULL,
                                override_outline_vars = character(0),
                                violin_alpha = 0.7,
                                box_alpha = 0.6,
                                box_width = 0.05,
                                x_tickangle = 45,
                                violin_outline_fill = FALSE,
                                box_outline_default = "grey20") {
    # Store data globally for PDF generation
    plot_data_key <- paste0("plot_data_", gsub("[^A-Za-z0-9]", "_", title))
    assign(plot_data_key, list(
      df_long = df_long,
      title = title,
      x_labels = x_labels,
      fill_map = fill_map,
      color_map = color_map,
      y_label = y_label,
      legend = legend,
      ylim = ylim,
      override_outline_vars = override_outline_vars,
      violin_alpha = violin_alpha,
      box_alpha = box_alpha,
      box_width = box_width,
      x_tickangle = x_tickangle,
      violin_outline_fill = violin_outline_fill,
      box_outline_default = box_outline_default
    ), envir = .GlobalEnv)

    # If this is a percentage plot, clamp values to [0,100] so violins don't extend under/over bounds
    df_plot <- df_long
    if (grepl("%", y_label)) {
      df_plot$Value <- pmin(pmax(df_plot$Value, 0), 100)
    } else if (grepl("count", y_label, ignore.case = TRUE)) {
      df_plot$Value <- pmax(df_plot$Value, 0)
    }

    # Compute shared bandwidth for KDE across both HTML and PDF
    valid_vals <- df_plot$Value[is.finite(df_plot$Value)]
    bw_shared <- if (length(valid_vals) >= 2) stats::bw.nrd0(valid_vals) else NULL

    # Store the clamped data for PDF generation as well (keeps parity)
    assign(plot_data_key, list(
      df_long = df_plot,
      title = title,
      x_labels = x_labels,
      fill_map = fill_map,
      color_map = color_map,
      y_label = y_label,
      legend = legend,
      ylim = ylim,
      override_outline_vars = override_outline_vars,
      violin_alpha = violin_alpha,
      box_alpha = box_alpha,
      box_width = box_width,
      x_tickangle = x_tickangle,
      violin_outline_fill = violin_outline_fill,
      box_outline_default = box_outline_default,
      bandwidth = bw_shared
    ), envir = .GlobalEnv)

    # Create plotly plot directly with explicit x-axis positioning
    p <- plot_ly()

    # Get all unique levels and create numeric positions
    all_levels <- levels(df_plot$Variable)
    if (is.null(all_levels)) {
      all_levels <- unique(as.character(df_plot$Variable))
    }

    # Create numeric x positions for each level
    x_positions <- seq_along(all_levels)
    names(x_positions) <- all_levels

    # Determine Plotly tick angle so labels finish at the tick mark (matching PDF hjust = 1)
    tick_angle_plotly <- if (!is.null(x_tickangle) && is.finite(x_tickangle) && x_tickangle != 0) x_tickangle else 0
    tick_label_position <- if (tick_angle_plotly == 0) "outside" else "outside right"

    # Add violin traces first (they will be in the background)
    for (i in seq_along(all_levels)) {
      var <- all_levels[i]
      var_data <- df_plot[df_plot$Variable == var, ]

      # Skip if no data for this level
      if (nrow(var_data) == 0) next

      line_col <- if (isTRUE(violin_outline_fill)) fill_map[var] else "black"
      fill_rgba <- to_rgba(fill_map[var], violin_alpha)

      p <- p %>% add_trace(
        x = rep(x_positions[var], nrow(var_data)),
        y = var_data$Value,
        type = "violin",
        side = "both",
        name = x_labels[i],
        fillcolor = fill_rgba,
        line = list(color = to_rgba(line_col, 1.0), width = 0.6),
        spanmode = "hard",
        bandwidth = bw_shared,
        points = FALSE,
        showlegend = FALSE,
        box = list(visible = FALSE),
        meanline = list(visible = FALSE),
        scalemode = "width",
        width = 0.8
      )
    }

    # Add boxplot traces second (they will be on top)
    for (i in seq_along(all_levels)) {
      var <- all_levels[i]
      var_data <- df_plot[df_plot$Variable == var, ]

      # Determine box outline color
      box_color <- if (var %in% override_outline_vars) "grey90" else box_outline_default
      fill_rgba_box <- to_rgba(fill_map[var], box_alpha)

      p <- p %>% add_trace(
        x = rep(x_positions[var], nrow(var_data)),
        y = var_data$Value,
        type = "box",
        name = paste(x_labels[i], "Box"),
        fillcolor = fill_rgba_box,
        line = list(color = to_rgba(box_color, 1.0), width = 0.8),
        showlegend = FALSE,
        boxpoints = FALSE,
        width = box_width
      )
    }

    # Add mean points last (on top of everything)
    mean_data <- df_plot %>%
      group_by(Variable) %>%
      summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")

    for (i in seq_along(all_levels)) {
      var <- all_levels[i]
      mean_row <- mean_data[mean_data$Variable == var, ]

      # Skip if no data for this level
      if (nrow(mean_row) == 0) next

      p <- p %>% add_trace(
        x = x_positions[var],
        y = mean_row$mean_value,
        type = "scatter",
        mode = "markers",
        name = paste(x_labels[i], "Mean"),
        marker = list(
          symbol = "x-thin-open",
          size = 6,
          color = "red",
          line = list(width = 1.5, color = "red")
        ),
        showlegend = FALSE
      )
    }

    # Store the data key as an attribute for PDF conversion
    attr(p, "plot_data_key") <- plot_data_key

    # Configure layout with explicit tick positions and labels
    html_title <- paste0("<b>", gsub("\n", "<br>", title), "</b>")
    p <- p %>% layout(
      title = list(text = html_title, font = list(size = 18), x = 0.5, xanchor = "center"),
      xaxis = list(
        title = "",
        tickmode = "array",
        tickvals = x_positions,
        ticktext = x_labels,
        tickangle = tick_angle_plotly,
        ticklabelposition = tick_label_position,
        tickfont = list(size = 16),
        showline = TRUE,
        linecolor = "black",
        linewidth = 1,
        zeroline = FALSE,
        range = c(min(x_positions) - 0.5, max(x_positions) + 0.5)
      ),
      yaxis = list(
        title = y_label,
        titlefont = list(size = 16),
        tickfont = list(size = 14),
        showline = TRUE,
        linecolor = "black",
        linewidth = 1,
        zeroline = FALSE
      ),
      showlegend = legend,
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)",
      font = list(family = "Arial", size = 14),
      margin = list(t = 110, l = 80, r = 80, b = ifelse(x_tickangle == 0, 60, 90))
    )

    # Apply y-axis limits if specified
    if (!is.null(ylim)) {
      p <- p %>% layout(yaxis = list(range = ylim))
    } else if (grepl("%", y_label)) {
      p <- p %>% layout(yaxis = list(range = c(0, 100)))
    } else if (grepl("count", y_label, ignore.case = TRUE)) {
      max_y <- suppressWarnings(max(df_plot$Value, na.rm = TRUE))
      if (!is.finite(max_y)) max_y <- 1
      p <- p %>% layout(yaxis = list(range = c(0, max_y * 1.05)))
    }

    return(p)
  }

  # Helper: convert plotly object back to ggplot for PDF output
  plotly_to_ggplot <- function(plotly_obj) {
    if (is.null(plotly_obj) || !inherits(plotly_obj, "plotly")) {
      return(ggplot() +
        labs(title = "Plot not available") +
        theme_minimal())
    }

    # Check if this plotly object has stored data
    plot_data_key <- attr(plotly_obj, "plot_data_key")
    if (!is.null(plot_data_key) && exists(plot_data_key, envir = .GlobalEnv)) {
      plot_data <- get(plot_data_key, envir = .GlobalEnv)

      # Use the stored data to create a proper ggplot
      return(build_violin_plot_ggplot(
        df_long = plot_data$df_long,
        title = plot_data$title,
        x_labels = plot_data$x_labels,
        fill_map = plot_data$fill_map,
        color_map = plot_data$color_map,
        y_label = plot_data$y_label,
        legend = plot_data$legend,
        ylim = plot_data$ylim,
        override_outline_vars = plot_data$override_outline_vars,
        violin_alpha = if (!is.null(plot_data$violin_alpha)) plot_data$violin_alpha else 0.7,
        box_alpha = if (!is.null(plot_data$box_alpha)) plot_data$box_alpha else 0.6,
        box_width = if (!is.null(plot_data$box_width)) plot_data$box_width else 0.05,
        x_tickangle = if (!is.null(plot_data$x_tickangle)) plot_data$x_tickangle else 45,
        violin_outline_fill = isTRUE(plot_data$violin_outline_fill),
        box_outline_default = if (!is.null(plot_data$box_outline_default)) plot_data$box_outline_default else "grey20",
        bandwidth = plot_data$bandwidth
      ))
    }

    # Check if this is a grouped plot built via build_grouped_violin_plot
    grouped_info <- attr(plotly_obj, "grouped_data")
    if (!is.null(grouped_info)) {
      df <- grouped_info$df
      # Ensure factor levels
      df$bin <- factor(df$bin, levels = grouped_info$bin_levels)
      df$group <- factor(df$group, levels = names(grouped_info$fill_map))

      # Compute effective bandwidth for grouped PDF
      vals <- df$value
      vals <- vals[is.finite(vals)]
      bw_eff <- grouped_info$bandwidth
      if (is.null(bw_eff) || !is.numeric(bw_eff) || is.na(bw_eff) || bw_eff <= 0) {
        if (length(vals) >= 2) {
          bw_eff <- stats::bw.nrd0(vals)
        } else {
          bw_eff <- NA_real_
        }
      }
      if (is.na(bw_eff) || bw_eff <= 0) bw_eff <- 0.1

      p <- ggplot(df, aes(x = bin, y = value, fill = group)) +
        # Violin outlines should match fill color
        geom_violin(aes(color = group),
          alpha = if (!is.null(grouped_info$violin_alpha)) grouped_info$violin_alpha else 0.7,
          position = position_dodge(width = if (!is.null(grouped_info$dodge_width)) grouped_info$dodge_width else 0.8), scale = "width", show.legend = TRUE, bw = bw_eff, trim = TRUE
        ) +
        scale_color_manual(values = grouped_info$fill_map, guide = "none") +
        geom_boxplot(
          width = if (!is.null(grouped_info$box_width)) grouped_info$box_width else 0.05,
          outlier.shape = NA,
          alpha = if (!is.null(grouped_info$box_alpha)) grouped_info$box_alpha else 0.6,
          position = position_dodge(width = if (!is.null(grouped_info$dodge_width)) grouped_info$dodge_width else 0.8),
          color = "grey20", show.legend = FALSE
        ) +
        stat_summary(
          fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1,
          position = position_dodge(width = if (!is.null(grouped_info$dodge_width)) grouped_info$dodge_width else 0.8), show.legend = FALSE
        ) +
        scale_fill_manual(values = grouped_info$fill_map, labels = grouped_info$legend_labels) +
        labs(title = grouped_info$title, x = "", y = grouped_info$y_label) +
        theme_classic(base_size = 14) +
        theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(
            size = 12, angle = if (!is.null(grouped_info$x_tickangle)) grouped_info$x_tickangle else 0,
            hjust = ifelse(!is.null(grouped_info$x_tickangle) && grouped_info$x_tickangle == 0, 0.5, 1)
          ),
          legend.position = "bottom",
          legend.title = element_blank()
        )

      if (!is.null(grouped_info$ylim)) {
        p <- p + coord_cartesian(ylim = grouped_info$ylim)
      }
      return(p)
    }

    # Check if this is a faceted plot built via build_violin_plot_facets
    facet_info <- attr(plotly_obj, "facet_data")
    if (!is.null(facet_info)) {
      df <- facet_info$df
      fill_map <- facet_info$fill_map
      x_labels <- facet_info$x_labels
      # Ensure factors and orders
      df$Variable <- factor(df$Variable, levels = names(fill_map))
      df$facet <- factor(df$facet, levels = facet_info$facet_levels)

      p <- ggplot(df, aes(x = Variable, y = Value, fill = Variable)) +
        geom_violin(alpha = if (!is.null(facet_info$violin_alpha)) facet_info$violin_alpha else 0.7, scale = "width", show.legend = FALSE) +
        geom_boxplot(width = if (!is.null(facet_info$box_width)) facet_info$box_width else 0.05, outlier.shape = NA, alpha = if (!is.null(facet_info$box_alpha)) facet_info$box_alpha else 0.6, show.legend = FALSE) +
        stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1, show.legend = FALSE) +
        scale_fill_manual(values = fill_map, labels = x_labels, guide = if (isTRUE(facet_info$show_legend)) guide_legend(override.aes = list(shape = NA)) else "none") +
        scale_x_discrete(labels = x_labels) +
        labs(title = facet_info$title, x = "", y = facet_info$y_label) +
        theme_classic(base_size = 14) +
        theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(
            size = 12, angle = if (!is.null(facet_info$x_tickangle)) facet_info$x_tickangle else 45,
            hjust = ifelse(!is.null(facet_info$x_tickangle) && facet_info$x_tickangle == 0, 0.5, 1)
          ),
          legend.position = if (isTRUE(facet_info$show_legend)) "bottom" else "none",
          strip.placement = "outside",
          strip.text.x = element_text(size = 16),
          strip.background = element_blank()
        ) +
        facet_grid(. ~ facet, scales = "free_x", space = "free", switch = "x")

      if (!is.null(facet_info$ylim)) {
        p <- p + coord_cartesian(ylim = facet_info$ylim)
      }
      return(p)
    }

    # Profile fallback: build ggplot from stored profile_data
    prof_info <- attr(plotly_obj, "profile_data")
    if (!is.null(prof_info)) {
      df <- prof_info$df
      # Canonical Fusion color override
      FUSION_COLOR <- "#F1C40F"
      detect_fusion <- function(df) {
        tryCatch(
          {
            (("category" %in% names(df)) && any(grepl("fusion", df$category, ignore.case = TRUE))) ||
              (("label" %in% names(df)) && any(grepl("fusion", df$label, ignore.case = TRUE)))
          },
          error = function(e) FALSE
        )
      }
      lc <- if (detect_fusion(df)) FUSION_COLOR else prof_info$line_color

      # Compute x-axis break count (1..K-1, ≥K)
      k_max <- if (!is.null(prof_info$k_max) && is.finite(prof_info$k_max)) prof_info$k_max else suppressWarnings(max(df$k[is.finite(df$k)], na.rm = TRUE))
      if (!is.finite(k_max) || is.na(k_max) || k_max < 2) k_max <- 20

      # Helper to lighten HEX colors
      lighten_hex <- function(hex, amount = 0.4) {
        rgb <- grDevices::col2rgb(hex)
        r <- as.integer(round(rgb[1] + (255 - rgb[1]) * amount))
        g <- as.integer(round(rgb[2] + (255 - rgb[2]) * amount))
        b <- as.integer(round(rgb[3] + (255 - rgb[3]) * amount))
        grDevices::rgb(r, g, b, maxColorValue = 255)
      }

      # Prepare summary-line data and aesthetics so PDF mirrors HTML styling
      stat_cols <- intersect(colnames(df), c("mean", "median"))
      line_stats <- if (length(stat_cols)) {
        df %>%
          dplyr::select(k, dplyr::all_of(stat_cols)) %>%
          tidyr::pivot_longer(cols = dplyr::all_of(stat_cols), names_to = "stat", values_to = "value") %>%
          dplyr::filter(!is.na(value)) %>%
          dplyr::mutate(stat = dplyr::recode(stat, mean = "Mean", median = "Median"))
      } else {
        data.frame(k = numeric(0), stat = character(0), value = numeric(0))
      }

      line_levels <- unique(line_stats$stat)
      line_palette <- if (length(line_levels)) setNames(rep(lc, length(line_levels)), line_levels) else character(0)
      linetype_values <- if (length(line_levels)) setNames(rep("solid", length(line_levels)), line_levels) else character(0)
      if ("Median" %in% names(linetype_values)) linetype_values["Median"] <- "dotdash"
      central_stat <- if ("Mean" %in% line_levels) "Mean" else if ("Median" %in% line_levels) "Median" else NULL
      legend_linewidths <- if (length(line_levels)) setNames(ifelse(line_levels == "Median", 1.0, 1.2), line_levels) else numeric(0)

      tick_breaks <- seq_len(k_max)
      label_last <- paste0("\u2265", k_max)
      ticktexts <- c(as.character(seq_len(k_max - 1)), label_last)

      p <- ggplot(df, aes(x = k)) +
        geom_ribbon(aes(ymin = q1, ymax = q3, fill = "IQR"), alpha = 0.25, show.legend = TRUE, key_glyph = "rect") +
        theme(legend.position = "bottom") +
        scale_y_continuous(limits = c(0, 100)) +
        scale_x_continuous(
          breaks = tick_breaks,
          labels = ticktexts,
          expand = expansion(mult = c(0.01, 0.01))
        ) +
        labs(title = prof_info$title, x = paste("Exons per", entity_label), y = prof_info$y_label) +
        theme_classic(base_size = 14) +
        theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.spacing.x = unit(6, "pt")
        )

      if ("Mean" %in% line_levels) {
        p <- p + geom_line(
          data = dplyr::filter(line_stats, stat == "Mean"),
          aes(y = value, color = stat, linetype = stat),
          linewidth = 1.2,
          show.legend = TRUE,
          key_glyph = "path"
        )
      }
      if ("Median" %in% line_levels) {
        p <- p + geom_line(
          data = dplyr::filter(line_stats, stat == "Median"),
          aes(y = value, color = stat, linetype = stat),
          linewidth = 1.0,
          show.legend = TRUE,
          key_glyph = "path"
        )
      }
      if (!is.null(central_stat)) {
        p <- p + geom_point(
          data = dplyr::filter(line_stats, stat == central_stat),
          aes(y = value),
          color = lc,
          size = 1.6,
          show.legend = FALSE
        )
      }

      p <- p +
        scale_fill_manual(
          name = "IQR",
          values = c("IQR" = lighten_hex(lc)),
          guide = guide_legend(order = 1, override.aes = list(alpha = 0.25, color = NA, fill = lighten_hex(lc)))
        )

      if (length(line_levels)) {
        p <- p +
          scale_color_manual(
            name = NULL,
            values = line_palette,
            breaks = line_levels,
            limits = line_levels,
            guide = guide_legend(
              order = 2,
              override.aes = list(
                linetype = unname(linetype_values[line_levels]),
                linewidth = unname(legend_linewidths[line_levels]),
                color = unname(line_palette[line_levels]),
                fill = NA
              )
            )
          ) +
          scale_linetype_manual(values = linetype_values, breaks = line_levels, limits = line_levels, guide = "none")
      }

      return(p)
    }

    # Fallback: create a simple placeholder
    title <- if (!is.null(plotly_obj$x$layout$title$text)) {
      plotly_obj$x$layout$title$text
    } else {
      "Distribution Plot"
    }

    ggplot() +
      labs(title = title) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 12)
      ) +
      annotate("text",
        x = 0.5, y = 0.5,
        label = "Plot data not available for PDF",
        size = 5, hjust = 0.5, vjust = 0.5
      )
  }

  # Helper: grouped violins by bin with legend (Annotated/Novel) using plotly, rebuildable for PDF
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

    # Store metadata for PDF reconstruction
    grouped_info <- list(
      df = data.frame(bin = df$bin, group = df$group, value = df$value),
      bin_levels = bin_levels,
      title = title,
      fill_map = fill_map,
      legend_labels = legend_labels,
      y_label = y_label,
      ylim = ylim,
      violin_alpha = violin_alpha,
      box_alpha = box_alpha,
      box_width = box_width,
      x_tickangle = x_tickangle,
      violin_width = violin_width,
      dodge_width = dodge_width,
      violangap = violangap,
      violingroupgap = violingroupgap
    )

    # Build plotly grouped violins
    # Clamp to the provided ylim to avoid tails outside bounds
    df_clamped <- df
    df_clamped$value <- pmin(pmax(df_clamped$value, ylim[1]), ylim[2])

    # Shared bandwidth across groups for consistent KDE
    valid_vals <- df_clamped$value[is.finite(df_clamped$value)]
    bw_shared <- if (length(valid_vals) >= 2) stats::bw.nrd0(valid_vals) else NULL

    p <- plot_ly()

    # Create combined x-axis labels for proper spacing
    # For each bin, we'll have two positions (one for each group)
    x_positions <- numeric()
    x_tick_positions <- numeric()
    x_tick_labels <- character()
    y_vals <- numeric()
    group_vals <- character()

    base_gap <- 0.3 # gap between groups within a bin
    bin_gap <- 1.0 # gap between bins

    current_x <- 0

    for (i in seq_along(bin_levels)) {
      bin_label <- bin_levels[i]

      # Center position for this bin's tick
      bin_center <- current_x + base_gap * (length(group_levels) - 1) / 2
      x_tick_positions <- c(x_tick_positions, bin_center)
      x_tick_labels <- c(x_tick_labels, bin_label)

      for (j in seq_along(group_levels)) {
        grp <- group_levels[j]
        grp_df <- df_clamped[df_clamped$bin == bin_label & df_clamped$group == grp, , drop = FALSE]

        if (nrow(grp_df) > 0) {
          x_pos <- current_x + (j - 1) * base_gap
          x_positions <- c(x_positions, rep(x_pos, nrow(grp_df)))
          y_vals <- c(y_vals, grp_df$value)
          group_vals <- c(group_vals, rep(grp, nrow(grp_df)))
        }

        current_x_for_trace <- current_x + (j - 1) * base_gap
        fill_rgba <- to_rgba(fill_map[grp], violin_alpha)

        # Only add trace if we have data
        if (nrow(grp_df) > 0) {
          p <- p %>% add_trace(
            x = rep(current_x_for_trace, nrow(grp_df)),
            y = grp_df$value,
            type = "violin",
            name = legend_labels[grp],
            legendgroup = grp,
            scalegroup = grp,
            fillcolor = fill_rgba,
            line = list(color = to_rgba(fill_map[grp], 1.0), width = 0.8),
            spanmode = "hard",
            bandwidth = bw_shared,
            side = "both",
            width = violin_width,
            points = FALSE,
            showlegend = (i == 1), # Only show legend for first bin
            box = list(visible = FALSE),
            meanline = list(visible = FALSE)
          )
        }
      }

      current_x <- current_x + length(group_levels) * base_gap + bin_gap
    }

    # Reset current_x for boxplots and mean markers
    current_x <- 0
    for (i in seq_along(bin_levels)) {
      bin_label <- bin_levels[i]

      for (j in seq_along(group_levels)) {
        grp <- group_levels[j]
        grp_df <- df_clamped[df_clamped$bin == bin_label & df_clamped$group == grp, , drop = FALSE]

        current_x_for_trace <- current_x + (j - 1) * base_gap
        fill_rgba_box <- to_rgba(fill_map[grp], box_alpha)

        # Boxplots
        if (nrow(grp_df) > 0) {
          p <- p %>% add_trace(
            x = rep(current_x_for_trace, nrow(grp_df)),
            y = grp_df$value,
            type = "box",
            name = paste0(legend_labels[grp], " Box"),
            legendgroup = grp,
            fillcolor = fill_rgba_box,
            line = list(color = to_rgba("grey20", 1.0), width = 1),
            showlegend = FALSE,
            boxpoints = FALSE,
            width = box_width
          )

          # Mean markers
          mean_val <- mean(grp_df$value, na.rm = TRUE)
          p <- p %>% add_trace(
            x = current_x_for_trace,
            y = mean_val,
            type = "scatter",
            mode = "markers",
            name = paste0(legend_labels[grp], " Mean"),
            legendgroup = grp,
            marker = list(symbol = "x-thin-open", size = 6, color = "red", line = list(width = 1.5, color = "red")),
            showlegend = FALSE
          )
        }
      }

      current_x <- current_x + length(group_levels) * base_gap + bin_gap
    }

    # Attach metadata for PDF
    attr(p, "grouped_data") <- grouped_info <- c(grouped_info, list(bandwidth = bw_shared))

    tick_angle_plotly <- if (!is.null(x_tickangle) && is.finite(x_tickangle) && x_tickangle != 0) x_tickangle else 0
    tick_label_position <- if (tick_angle_plotly == 0) "outside" else "outside right"

    html_title <- paste0("<b>", gsub("\n", "<br>", title), "</b>")
    # Build legend object with optional title
    legend_obj <- list(orientation = "h", x = 0.5, xanchor = "center", y = -0.15, yanchor = "top")
    if (!is.null(legend_title)) legend_obj$title <- list(text = legend_title)

    p <- p %>% layout(
      title = list(text = html_title, font = list(size = 18), x = 0.5, xanchor = "center"),
      xaxis = list(
        title = "",
        tickmode = "array",
        tickvals = x_tick_positions,
        ticktext = x_tick_labels,
        tickangle = tick_angle_plotly,
        ticklabelposition = tick_label_position,
        tickfont = list(size = 16),
        showline = TRUE,
        linecolor = "black",
        linewidth = 1,
        zeroline = FALSE,
        range = c(-0.5, max(x_tick_positions) + 1)
      ),
      yaxis = list(
        title = y_label, titlefont = list(size = 16), tickfont = list(size = 14), range = ylim,
        showline = TRUE, linecolor = "black", linewidth = 1, zeroline = FALSE
      ),
      violinmode = "overlay",
      legend = legend_obj,
      showlegend = TRUE,
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)",
      font = list(family = "Arial", size = 14),
      margin = list(t = 80, l = 80, r = 80, b = 120)
    )
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
  render_pdf_plot <- function(name, converter = plotly_to_ggplot) {
    if (exists(name)) {
      obj <- get(name)
      print(if (is.null(converter)) obj else converter(obj))
    }
  }

  # Center a ggplot on a page with reduced width
  render_pdf_plot_centered <- function(name, width_frac = 0.45, converter = plotly_to_ggplot) {
    if (!exists(name)) {
      return(invisible(NULL))
    }
    obj <- get(name)
    p <- if (is.null(converter)) obj else converter(obj)
    g <- if (inherits(p, "grob")) p else ggplotGrob(p)
    left_right <- (1 - width_frac) / 2
    grid.arrange(nullGrob(), g, nullGrob(), widths = c(left_right, width_frac, left_right), newpage = TRUE)
  }

  # Helper: build length-distribution violins for given column prefix using native plotly
  # If mono=TRUE, uses *_length_mono_prop columns; otherwise *_length_prop
  build_len_violin_for_prefix <- function(df, prefix, title, fill_color, box_fill = NULL, mono = FALSE, box_outline_color = "grey20", violin_alpha = 0.5, box_alpha = 0.3, violin_outline_fill = FALSE) {
    if (is.null(box_fill)) box_fill <- fill_color
    suffix <- if (mono) "_length_mono_prop" else "_length_prop"
    cols <- c(
      paste0(prefix, "_250b", suffix),
      paste0(prefix, "_500b", suffix),
      paste0(prefix, "_short", suffix),
      paste0(prefix, "_mid", suffix),
      paste0(prefix, "_long", suffix)
    )
    df_long <- pivot_longer(df, cols = all_of(cols), names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
    df_long$Variable <- factor(df_long$Variable, levels = cols)
    # Clamp to [0,100] because these are proportions
    df_long$Value <- pmin(pmax(df_long$Value, 0), 100)

    # Store data globally for PDF generation
    plot_data_key <- paste0("plot_data_", gsub("[^A-Za-z0-9]", "_", title))

    # Create fill and color maps for this plot
    fill_map <- setNames(rep(fill_color, length(cols)), cols)
    color_map <- setNames(rep(fill_color, length(cols)), cols)
    x_labels <- c("0-250bp", "250-500bp", "500-1000bp", "1000-2000bp", ">2000bp")
    names(x_labels) <- cols

    assign(plot_data_key, list(
      df_long = df_long,
      title = title,
      x_labels = x_labels,
      fill_map = fill_map,
      color_map = color_map,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      legend = FALSE,
      ylim = NULL,
      override_outline_vars = character(0),
      violin_alpha = violin_alpha,
      box_alpha = box_alpha,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = violin_outline_fill,
      box_outline_default = box_outline_color
    ), envir = .GlobalEnv)

    # Create plotly plot directly
    p <- plot_ly()

    # Add violin traces first (background)
    for (i in seq_along(levels(df_long$Variable))) {
      var <- levels(df_long$Variable)[i]
      var_data <- df_long[df_long$Variable == var, ]

      line_col <- if (isTRUE(violin_outline_fill)) fill_color else "black"
      fill_rgba <- to_rgba(fill_color, violin_alpha)
      p <- p %>% add_trace(
        data = var_data,
        x = ~Variable,
        y = ~Value,
        type = "violin",
        side = "both",
        name = c("0-250bp", "250-500bp", "500-1000bp", "1000-2000bp", ">2000bp")[i],
        fillcolor = fill_rgba,
        line = list(color = to_rgba(line_col, 1.0), width = 0.6),
        spanmode = "hard",
        points = FALSE,
        showlegend = FALSE,
        box = list(visible = FALSE),
        meanline = list(visible = FALSE)
      )
    }

    # Add boxplot traces second (on top)
    for (i in seq_along(levels(df_long$Variable))) {
      var <- levels(df_long$Variable)[i]
      var_data <- df_long[df_long$Variable == var, ]

      # Use the box_outline_color parameter directly for the line color
      fill_rgba_box <- to_rgba(if (is.null(box_fill)) fill_color else box_fill, box_alpha)

      p <- p %>% add_trace(
        data = var_data,
        x = ~Variable,
        y = ~Value,
        type = "box",
        name = paste(c("0-250bp", "250-500bp", "500-1000bp", "1000-2000bp", ">2000bp")[i], "Box"),
        fillcolor = fill_rgba_box,
        line = list(color = to_rgba(box_outline_color, 1.0), width = 0.8),
        showlegend = FALSE,
        boxpoints = FALSE,
        width = 0.05
      )
    }

    # Add mean points last (on top of everything)
    mean_data <- df_long %>%
      group_by(Variable) %>%
      summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")

    p <- p %>% add_trace(
      data = mean_data,
      x = ~Variable,
      y = ~mean_value,
      type = "scatter",
      mode = "markers",
      name = "Mean",
      marker = list(
        symbol = "x-thin-open",
        size = 6,
        color = "red",
        line = list(width = 1.5, color = "red")
      ),
      showlegend = FALSE
    )

    # Store the data key as an attribute for PDF conversion
    attr(p, "plot_data_key") <- plot_data_key

    # Configure layout
    html_title <- paste0("<b>", gsub("\n", "<br>", title), "</b>")
    tick_angle_plotly <- 45
    tick_label_position <- "outside right"
    p <- p %>% layout(
      title = list(text = html_title, font = list(size = 18), x = 0.5, xanchor = "center"),
      xaxis = list(
        title = "",
        categoryorder = "array",
        categoryarray = cols,
        tickmode = "array",
        tickvals = cols,
        ticktext = x_labels,
        tickangle = tick_angle_plotly,
        ticklabelposition = tick_label_position,
        tickfont = list(size = 16),
        showline = TRUE,
        linecolor = "black",
        linewidth = 1,
        zeroline = FALSE
      ),
      yaxis = list(
        title = paste(entity_label_plural, ", %", sep = ""),
        titlefont = list(size = 16),
        tickfont = list(size = 14),
        range = c(0, 100),
        showline = TRUE,
        linecolor = "black",
        linewidth = 1,
        zeroline = FALSE
      ),
      showlegend = FALSE,
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)",
      font = list(family = "Arial", size = 14),
      margin = list(t = 110, l = 80, r = 80, b = 100)
    )

    return(p)
  }

  # Helper: build a per-category exon count profile (median + IQR across cells)
  build_exon_profile_plot <- function(df_prof, title, line_color, k_max = 20, y_label = paste(entity_label_plural, ", %", sep = ""), n_cells = NULL) {
    # Sanitize title and forbid subtitles
    title <- tryCatch(
      {
        t <- gsub("(?i)cells included:[^|<>]*", "", title, perl = TRUE)
        t <- gsub("(?i)<br[^>]*>.*", "", t, perl = TRUE) # drop any HTML subtitle lines
        trimws(gsub("\n.*", "", t))
      },
      error = function(e) title
    )

    # Canonical Fusion color
    FUSION_COLOR <- "#F1C40F"
    # Helper: detect Fusion profiles using title or data, case-insensitive
    is_fusion_profile <- function(df, ttl) {
      ttl_has <- tryCatch(!is.null(ttl) && grepl("fusion", ttl, ignore.case = TRUE), error = function(e) FALSE)
      in_df <- tryCatch(
        {
          any(grepl("fusion", names(df), ignore.case = TRUE)) ||
            (("category" %in% names(df)) && any(grepl("fusion", df$category, ignore.case = TRUE))) ||
            (("label" %in% names(df)) && any(grepl("fusion", df$label, ignore.case = TRUE)))
        },
        error = function(e) FALSE
      )
      ttl_has || in_df
    }
    # Helpers: color utilities
    hex_to_rgba <- function(hex, alpha = 0.25) {
      rgb <- grDevices::col2rgb(hex)
      sprintf("rgba(%d,%d,%d,%.2f)", rgb[1], rgb[2], rgb[3], alpha)
    }
    lighten_hex <- function(hex, amount = 0.35) {
      # amount in [0,1] toward white
      rgb <- grDevices::col2rgb(hex)
      r <- as.integer(round(rgb[1] + (255 - rgb[1]) * amount))
      g <- as.integer(round(rgb[2] + (255 - rgb[2]) * amount))
      b <- as.integer(round(rgb[3] + (255 - rgb[3]) * amount))
      grDevices::rgb(r, g, b, maxColorValue = 255)
    }

    # Override Fusion color if needed
    if (is_fusion_profile(df_prof, title)) {
      line_color <- FUSION_COLOR
    }

    # df_prof columns: k, median, q1, q3, (optional) mean
    if (is.null(df_prof) || nrow(df_prof) == 0 || all(!is.finite(df_prof$median))) {
      p_empty <- plot_ly() %>% layout(
        title = list(text = paste0("<b>", title, "</b>"), font = list(size = 18), x = 0.5, xanchor = "center"),
        yaxis = list(title = y_label, range = c(0, 100)),
        annotations = list(
          list(
            text = "No data available for this category",
            showarrow = FALSE,
            x = 0.5,
            y = 0.5,
            xref = "paper",
            yref = "paper",
            font = list(size = 14, color = "gray")
          )
        )
      )
      return(p_empty)
    }
    label_last <- paste0("\u2265", k_max) # ≥K
    ticktexts <- c(as.character(seq_len(k_max - 1)), label_last)

    # Choose center line: use mean if present, else median
    y_center <- if (!is.null(df_prof$mean)) df_prof$mean else df_prof$median

    p <- plot_ly()
    # Lower bound (q1)
    p <- p %>% add_trace(
      x = df_prof$k, y = df_prof$q1,
      type = "scatter", mode = "lines",
      line = list(color = line_color, width = 0.0001),
      name = "Q1", showlegend = FALSE
    )
    # Upper bound (q3) with fill to previous trace
    if (!is.null(df_prof$q1) && !is.null(df_prof$q3)) {
      iqr_fill <- hex_to_rgba(lighten_hex(line_color, 0.4), 0.25)
      p <- p %>% plotly::add_ribbons(
        x = df_prof$k, ymin = df_prof$q1, ymax = df_prof$q3,
        fillcolor = iqr_fill,
        line = list(color = "rgba(0,0,0,0)"),
        name = "IQR", showlegend = TRUE
      )
    }
    # Central (mean or median) straight line + markers
    p <- p %>% add_trace(
      x = df_prof$k, y = y_center,
      type = "scatter", mode = "lines+markers",
      line = list(color = line_color, width = 2.5),
      marker = list(color = line_color, size = 6),
      name = if (!is.null(df_prof$mean)) "Mean" else "Median", showlegend = TRUE
    )
    # If mean exists, overlay dashed median
    if (!is.null(df_prof$mean) && !is.null(df_prof$median)) {
      p <- p %>% add_trace(
        x = df_prof$k, y = df_prof$median,
        type = "scatter", mode = "lines",
        line = list(color = line_color, width = 1.5, dash = "dash"),
        name = "Median", showlegend = TRUE
      )
    }
    # Legend at bottom, horizontal (no extra annotations like 'Cells included')
    p <- p %>% layout(
      legend = list(orientation = "h", y = -0.2, x = 0.5, xanchor = "center")
    )
    p <- p %>% layout(
      title = list(text = title, font = list(size = 18), x = 0.5, xanchor = "center"),
      xaxis = list(
        title = paste("Exons per", entity_label),
        tickmode = "array", tickvals = seq_len(k_max), ticktext = ticktexts,
        showline = TRUE, linecolor = "black", linewidth = 1, zeroline = FALSE
      ),
      yaxis = list(
        title = y_label, range = c(0, 100),
        showline = TRUE, linecolor = "black", linewidth = 1, zeroline = FALSE
      ),
      paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
      font = list(family = "Arial", size = 14),
      margin = list(t = 90, l = 80, r = 60, b = 80)
    )

    # Attach metadata for PDF fallback
    attr(p, "profile_data") <- list(
      df = df_prof, title = title, line_color = line_color, k_max = k_max, y_label = y_label
    )
    return(p)
  }

  ### Basic cell informtion ###
  #############################

  single_defaults <- list(
    violin_alpha = 0.5,
    box_alpha = 0.3,
    box_width = 0.05,
    x_tickangle = 45,
    violin_outline_fill = FALSE,
    box_outline_default = "black"
  )
  common_plot_args <- list(
    violin_alpha = 0.5,
    box_alpha = 0.3,
    box_width = 0.05,
    violin_outline_fill = FALSE,
    box_outline_default = "black"
  )

  # 1. Number of Reads Across Cells
  cfg_reads <- list(
    column = count_col,
    name = "gg_reads_in_cells",
    title = paste("Number of", entity_label_plural, "Across Cells"),
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
      title = "Number of UMIs Across Cells",
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
    title = "Number of Genes Across Cells",
    fill = "#CC6633",
    y_label = "Genes, count",
    x_label = "Cells",
    plot_args = common_plot_args
  )
  single_violin(SQANTI_cell_summary, cfg_genes)

  # 4. Number of Unique Junction Chains Across Cells
  if (mode != "isoforms") {
    cfg_ujcs <- list(
      column = "UJCs_in_cell",
      name = "gg_JCs_in_cell",
      title = "Number of Unique Junction Chains Across Cells",
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
    title = "Number of Known/Novel Genes Across Cells",
    x_labels = c("Annotated Genes", "Novel Genes"),
    y_label = paste(entity_label_plural, ", counts", sep = ""),
    fill_map = c("Annotated_genes" = "#CC6633", "Novel_genes" = "#CC6633"),
    plot_args = pivot_defaults
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
      fill_map = c("Annotated_genes_perc" = "#CC6633", "Novel_genes_perc" = "#CC6633"),
      plot_args = pivot_defaults
    ))
  }

  # 5. Percentage of Reads/Transcripts from Known/Novel Genes Across Cells
  # (Enabled for both reads and isoforms modes)
  {
    classification_valid <- Classification_file[Classification_file$CB != "unassigned" & !is.na(Classification_file$CB), ]

    if (nrow(classification_valid) > 0) {
      annotated_reads_per_cell <- classification_valid %>%
        filter(!grepl("^novel", associated_gene)) %>%
        group_by(CB) %>%
        summarise(Annotated_genes_reads = if (mode == "isoforms") sum(count, na.rm = TRUE) else n(), .groups = "drop")

      novel_reads_per_cell <- classification_valid %>%
        filter(grepl("^novel", associated_gene)) %>%
        group_by(CB) %>%
        summarise(Novel_genes_reads = if (mode == "isoforms") sum(count, na.rm = TRUE) else n(), .groups = "drop")

      SQANTI_cell_summary <- SQANTI_cell_summary %>%
        left_join(annotated_reads_per_cell, by = "CB") %>%
        left_join(novel_reads_per_cell, by = "CB")

      SQANTI_cell_summary$Annotated_genes_reads[is.na(SQANTI_cell_summary$Annotated_genes_reads)] <- 0
      SQANTI_cell_summary$Novel_genes_reads[is.na(SQANTI_cell_summary$Novel_genes_reads)] <- 0

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
        fill_map = c("Annotated_reads_perc" = "#CC6633", "Novel_reads_perc" = "#CC6633"),
        plot_args = pivot_defaults
      ))
    } else {
      message("Warning: No valid classification data found. Skipping read expression by gene annotation plot.")
      gg_annotation_of_reads_in_cell <<- plot_ly() %>%
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
    plot_args = list(
      violin_alpha = 0.5,
      box_alpha = 3,
      box_width = 0.05,
      x_tickangle = 45
    )
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

  # Build per-cell per-gene read counts from classification
  genes_by_cb <- Classification_file %>%
    filter(!is.na(CB), CB != "unassigned", !is.na(associated_gene)) %>%
    group_by(CB, associated_gene) %>%
    summarise(reads_per_gene = n(), .groups = "drop") %>%
    mutate(
      gene_type = ifelse(grepl("^novel", associated_gene), "Novel", "Annotated"),
      bin = vapply(reads_per_gene, gene_bin_label, character(1))
    ) %>%
    filter(!is.na(bin))

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

  # Combined (all genes together): one violin per bin
  read_bins_all <- genes_by_cb %>%
    group_by(CB, bin) %>%
    summarise(num_genes = n(), .groups = "drop") %>%
    group_by(CB) %>%
    mutate(percentage = 100 * num_genes / sum(num_genes)) %>%
    ungroup() %>%
    tidyr::complete(CB, bin = gene_bin_levels, fill = list(num_genes = 0, percentage = 0))

  read_bins_all$bin <- factor(read_bins_all$bin, levels = gene_bin_levels)

  {
    df_long <- data.frame(
      Variable = factor(read_bins_all$bin, levels = gene_bin_levels),
      Value = read_bins_all$percentage
    )
    fill_map <- setNames(rep("#CC6633", length(gene_bin_levels)), gene_bin_levels)
    gg_read_bins_all <<- build_violin_plot(
      df_long,
      title = paste("Distribution of Genes by", entity_label, "Count Bins Across Cells"),
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

    ujc_by_cb <- Classification_file %>%
      filter(!is.na(CB), CB != "unassigned", !is.na(associated_gene), exons > 1) %>%
      group_by(CB, associated_gene) %>%
      summarise(ujc_per_gene = dplyr::n_distinct(jxn_string), .groups = "drop") %>%
      mutate(bin = vapply(ujc_per_gene, ujc_bin_label, character(1))) %>%
      filter(!is.na(bin))

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
      title = "Distribution of Genes by UJC Count Bins Across Cells",
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

  # Mitochondrial percentage in cell
  {
    df_long <- data.frame(Variable = "MT_perc", Value = SQANTI_cell_summary$MT_perc)
    df_long$Variable <- factor(df_long$Variable, levels = "MT_perc")
    fill_map <- c("MT_perc" = "#CC6633")
    x_labels <- c("Cell")
    gg_MT_perc <<- build_violin_plot(
      df_long,
      title = paste("Mitochondrial", entity_label_plural, "Across Cells"),
      x_labels = x_labels,
      fill_map = fill_map,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      legend = FALSE,
      violin_alpha = 0.5,
      box_alpha = 3,
      box_width = 0.05,
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
  gg_bulk_all_reads <<- ggplot(Classification_file, aes(x = length)) +
    geom_histogram(binwidth = 50, fill = "#CC6633", color = "black", alpha = 0.5) +
    labs(
      title = paste("All", entity_label, "Lengths Distribution"),
      x = paste(entity_label, "length"),
      y = paste(entity_label_plural, ", counts", sep = "")
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 40, r = 5, b = 5, l = 5, unit = "pt"),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16)
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

  gg_bulk_length_by_category <<- ggplot(Classification_file, aes(x = length, color = structural_category_pretty)) +
    geom_freqpoly(binwidth = 100, linewidth = 1.2, na.rm = TRUE) +
    labs(
      title = paste("All", entity_label, "Lengths Distribution by Structural Category"),
      x = paste(entity_label, "length"),
      y = paste(entity_label_plural, ", counts", sep = ""),
      color = NULL
    ) +
    theme_classic(base_size = 16) +
    scale_color_manual(values = structural_category_palette, drop = FALSE) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 8),
      labels = scales::comma
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(0.8, "cm"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 40, r = 5, b = 5, l = 5, unit = "pt"),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16)
    ) +
    guides(color = guide_legend(nrow = 2))

  # Mono vs multi-exon classification for length
  Classification_file$exons <- as.numeric(Classification_file$exons)

  gg_bulk_length_by_exon_type <<- ggplot(
    Classification_file,
    aes(x = length, color = ifelse(exons == 1, "Mono-Exon", "Multi-Exon"))
  ) +
    geom_freqpoly(binwidth = 100, linewidth = 1.2, na.rm = TRUE) +
    labs(
      title = paste("Mono- vs Multi- Exon", entity_label, "Lengths Distribution"),
      x = paste(entity_label, "length"),
      y = paste(entity_label_plural, ", counts", sep = ""),
      color = NULL
    ) +
    theme_classic(base_size = 16) +
    scale_color_manual(
      values = c("Multi-Exon" = "#3B0057", "Mono-Exon" = "#FFE44C")
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 8),
      labels = scales::comma
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 40, r = 5, b = 5, l = 5, unit = "pt"),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16)
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
    cols <- c(
      "FSM_ref_coverage_prop", "ISM_ref_coverage_prop", "NIC_ref_coverage_prop", "NNC_ref_coverage_prop",
      "Genic_ref_coverage_prop", "Antisense_ref_coverage_prop", "Fusion_ref_coverage_prop", "Intergenic_ref_coverage_prop",
      "Genic_intron_ref_coverage_prop"
    )
    gg_SQANTI_pivot <- pivot_long(SQANTI_cell_summary, cols)
    fill_map <- setNames(unname(cat_fill_map), cols)
    x_labels <- cat_labels_pretty
    # Build dynamic title using cutoff from cell summary when available
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
      override_outline_vars = c("Genic_ref_coverage_prop", "Genic_intron_ref_coverage_prop"),
      violin_outline_fill = TRUE
    )
  }


  ### Structural categories ###
  #############################

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
  if (!skipORF) {
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
      plot_args = list(override_outline_vars = c("genic_coding_prop", "genic_intron_coding_prop"), violin_outline_fill = TRUE)
    ))



    # Define colors for non-coding (lighter versions)
    noncoding_fill_map <- c(
      "FSM_non_coding_prop" = "#D2E6F2",
      "ISM_non_coding_prop" = "#FEDCCD",
      "NIC_non_coding_prop" = "#D6EDD6",
      "NNC_non_coding_prop" = "#F9D2CA",
      "genic_non_coding_prop" = "#DFDFDF",
      "antisense_non_coding_prop" = "#D1ECE3",
      "fusion_non_coding_prop" = "#FFECBD",
      "intergenic_non_coding_prop" = "#F8DFD7",
      "genic_intron_non_coding_prop" = "#C6E9ED"
    )

    pivot_violin(SQANTI_cell_summary, list(
      name = "gg_non_coding_across_category",
      columns = names(noncoding_fill_map),
      title = "Non-coding Proportion of Structural Categories Distribution Across Cells",
      x_labels = cat_labels_pretty,
      y_label = paste(entity_label_plural, ", %", sep = ""),
      fill_map = noncoding_fill_map,
      plot_args = list(override_outline_vars = c("genic_intron_non_coding_prop"), violin_outline_fill = TRUE)
    ))
  } # End of if (!skipORF)

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

  # Determine which bad feature columns are actually present in SQANTI_cell_summary
  # This implicitly handles skipORF, as NMD_prop_in_cell won't be in SQANTI_cell_summary if skipORF is TRUE
  bad_feature_cols_present <- intersect(names(all_bad_features_map), colnames(SQANTI_cell_summary))
  bad_feature_cols_present <- bad_feature_cols_present[sapply(bad_feature_cols_present, function(col) any(!is.na(SQANTI_cell_summary[[col]])) && sum(SQANTI_cell_summary[[col]], na.rm = TRUE) > 0)] # keep only if data exists

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
    gg_bad_feature <<- plot_ly() %>%
      layout(
        title = "No bad quality features to display",
        annotations = list(
          text = "No bad quality features to display",
          showarrow = FALSE,
          font = list(size = 18, color = "gray")
        )
      )
  }

  ### Good features plots ###
  ##########################

  good_specs <- list(
    list(cols = cat_cols("_TSSAnnotationSupport"), title = "TSS Annotation Support by Structural Category", color = "#66C2A4", name = "gg_tss_annotation_support", require_all = TRUE),
    list(cols = cat_cols("_CAGE_peak_support_prop"), title = "CAGE Peak Support by Structural Category", color = "#EE6A50", name = "gg_cage_peak_support", require_all = TRUE),
    list(cols = cat_cols("_PolyA_motif_support_prop"), title = "PolyA Support by Structural Category", color = "#78C679", name = "gg_polyA_motif_support", require_all = TRUE),
    list(cols = cat_cols("_canon_prop"), title = "Canonical Junctions by Structural Category", color = "#CC6633", name = "gg_canon_by_category", require_all = TRUE)
  )
  invisible(lapply(good_specs, function(sp) {
    if (!is.null(sp$require_all) && sp$require_all && !all(sp$cols %in% colnames(SQANTI_cell_summary))) {
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

  color_map <- c(
    "TSSAnnotationSupport_prop" = "#66C2A4",
    "CAGE_peak_support_prop" = "#EE6A50",
    "PolyA_motif_support_prop" = "#78C679",
    "Canonical_prop_in_cell" = "#CC6633"
  )
  label_map <- c(
    "TSSAnnotationSupport_prop" = "TSS Annotated",
    "CAGE_peak_support_prop" = "Has Coverage CAGE",
    "PolyA_motif_support_prop" = "Has PolyA Motif",
    "Canonical_prop_in_cell" = "Canonical Junctions"
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

    # 1) Median exons per read per cell and category
    cls_valid <- Classification_file %>%
      dplyr::filter(CB != "unassigned") %>%
      mutate(cat_key = unname(cat_key_map[structural_category])) %>%
      filter(!is.na(cat_key))
    exons_mean_by_cell <- cls_valid %>%
      group_by(CB, cat_key) %>%
      summarise(median_exons = median(as.numeric(exons), na.rm = TRUE), .groups = "drop") %>%
      tidyr::complete(CB, cat_key = all_cats, fill = list(median_exons = NA_real_))

    df_exon_mean_long <- data.frame(
      Variable = factor(exons_mean_by_cell$cat_key, levels = all_cats),
      Value = exons_mean_by_cell$median_exons
    )

    gg_exon_mean_by_category <<- build_violin_plot(
      df_long = df_exon_mean_long,
      title = paste("Median Exons per", entity_label, "by Structural Category Across Cells"),
      x_labels = x_labels_pretty,
      fill_map = fill_map_cat,
      y_label = paste("Exons per", entity_label),
      legend = FALSE,
      violin_alpha = 0.7,
      box_alpha = 0.3,
      box_width = 0.05,
      x_tickangle = 45,
      violin_outline_fill = TRUE,
      box_outline_default = "grey20"
    )

    # 2) Percent mono-exonic reads per cell and category
    exons_bin_by_cell <- cls_valid %>%
      mutate(is_mono = as.numeric(exons) == 1) %>%
      group_by(CB, cat_key) %>%
      summarise(total = dplyr::n(), mono = sum(is_mono, na.rm = TRUE), .groups = "drop") %>%
      tidyr::complete(CB, cat_key = all_cats, fill = list(total = 0, mono = 0)) %>%
      mutate(perc_mono = ifelse(total > 0, 100 * mono / total, 0))

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
      box_outline_default = "grey20"
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
        summarise(n = dplyr::n(), .groups = "drop") %>%
        group_by(CB) %>%
        mutate(perc = 100 * n / sum(n)) %>%
        ungroup() %>%
        tidyr::complete(CB, bin = exon_bin_levels, fill = list(n = 0, perc = 0))

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
      # Cells with at least min_reads in this category; if none, fallback to all cells with any reads
      cells_ok <- cat_df %>%
        group_by(CB) %>%
        summarise(total = dplyr::n(), .groups = "drop") %>%
        filter(total >= min_reads) %>%
        pull(CB)
      if (length(cells_ok) == 0) {
        cells_ok <- unique(cat_df$CB)
      }
      cat_df2 <- cat_df %>%
        filter(CB %in% cells_ok) %>%
        mutate(exons_n = pmin(as.numeric(exons), K))
      # Per-cell PMF
      pmf <- cat_df2 %>%
        group_by(CB, exons_n) %>%
        summarise(n = dplyr::n(), .groups = "drop") %>%
        group_by(CB) %>%
        mutate(perc = 100 * n / sum(n)) %>%
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

  # Build SJ per-type/per-category plots and the all-canonical grouped plot for HTML (and reuse for PDF)
  # This block creates plot objects regardless of generate_pdf so the Rmd can render them.
  {
    # Junction type percentages by structural category across cells
    junc_aug_html <- Junctions %>%
      dplyr::filter(CB != "unassigned") %>%
      mutate(junction_type = paste(junction_category, canonical, sep = "_"))
    if (!("structural_category" %in% colnames(junc_aug_html))) {
      join_key <- NULL
      for (k in c("isoform", "readID", "read_id", "ID", "read_name", "read")) {
        if (k %in% colnames(Junctions) && k %in% colnames(Classification_file)) {
          join_key <- k
          break
        }
      }
      if (!is.null(join_key)) {
        by_vec <- c(CB = "CB")
        by_vec[[join_key]] <- join_key
        junc_aug_html <- junc_aug_html %>%
          left_join(Classification_file %>% select(all_of(c(join_key, "CB", "structural_category"))), by = by_vec)
      } else {
        junc_aug_html$structural_category <- NA_character_
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
        total = n(),
        known_canonical = sum(junction_type == "known_canonical", na.rm = TRUE),
        known_non_canonical = sum(junction_type == "known_non_canonical", na.rm = TRUE),
        novel_canonical = sum(junction_type == "novel_canonical", na.rm = TRUE),
        novel_non_canonical = sum(junction_type == "novel_non_canonical", na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::complete(CB, cat_key = all_cats, fill = list(total = 0, known_canonical = 0, known_non_canonical = 0, novel_canonical = 0, novel_non_canonical = 0)) %>%
      mutate(
        KnownCanonicalPerc = ifelse(total > 0, 100 * known_canonical / total, 0),
        KnownNonCanonicalPerc = ifelse(total > 0, 100 * known_non_canonical / total, 0),
        NovelCanonicalPerc = ifelse(total > 0, 100 * novel_canonical / total, 0),
        NovelNonCanonicalPerc = ifelse(total > 0, 100 * novel_non_canonical / total, 0)
      ) %>%
      ungroup()

    make_df_long_html <- function(col_name) {
      data.frame(Variable = factor(junc_summ_html$cat_key, levels = all_cats), Value = junc_summ_html[[col_name]])
    }

    fill_map_cat <- cat_fill_map

    # Create plotly versions for HTML
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
    tick_angle_plotly <- 45

    gg_sj_type_by_category_stack <<- subplot(
      gg_known_canon_by_category, gg_known_noncanon_by_category,
      gg_novel_canon_by_category, gg_novel_noncanon_by_category,
      nrows = 4, shareX = TRUE, titleX = FALSE, margin = 0.02, heights = c(0.25, 0.25, 0.25, 0.25)
    ) %>% layout(
      title = list(text = "<b>Splice Junctions Distribution by Structural Category Across Cells</b><br>", x = 0.5, xanchor = "center", font = list(size = 18)),
      showlegend = FALSE,
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)",
      font = list(family = "Arial", size = 14),
      margin = list(t = 100, l = 100, r = 80, b = 130),
      # Configure y-axes
      yaxis = list(
        title = "Known Canonical Junctions, %", titlefont = list(size = 16), tickfont = list(size = 14),
        range = c(0, 100), showline = TRUE, linecolor = "black", linewidth = 1
      ),
      yaxis2 = list(
        title = "Known Non-canonical Junctions, %", titlefont = list(size = 16), tickfont = list(size = 14),
        range = c(0, 100), showline = TRUE, linecolor = "black", linewidth = 1
      ),
      yaxis3 = list(
        title = "Novel Canonical Junctions, %", titlefont = list(size = 16), tickfont = list(size = 14),
        range = c(0, 100), showline = TRUE, linecolor = "black", linewidth = 1
      ),
      yaxis4 = list(
        title = "Novel Non-canonical Junctions, %", titlefont = list(size = 16), tickfont = list(size = 14),
        range = c(0, 100), showline = TRUE, linecolor = "black", linewidth = 1
      ),
      # Configure x-axes: use shared x-axis to display category labels on bottom subplot
      xaxis = list(
        title = "",
        showticklabels = TRUE,
        showline = TRUE,
        linecolor = "black",
        linewidth = 1,
        tickmode = "array",
        tickvals = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
        ticktext = c("FSM", "ISM", "NIC", "NNC", "Genic<br>Genomic", "Antisense", "Fusion", "Intergenic", "Genic<br>Intron"),
        tickfont = list(size = 16),
        tickangle = tick_angle_plotly,
        ticklabelposition = "outside right",
        ticks = "outside",
        ticklen = 8,
        automargin = TRUE,
        range = c(0.5, 9.5)
      ),
      xaxis2 = list(
        showticklabels = FALSE,
        showline = TRUE,
        linecolor = "black",
        linewidth = 1,
        range = c(0.5, 9.5),
        tickangle = tick_angle_plotly,
        ticklabelposition = "outside right"
      ),
      xaxis3 = list(
        showticklabels = FALSE,
        showline = TRUE,
        linecolor = "black",
        linewidth = 1,
        range = c(0.5, 9.5),
        tickangle = tick_angle_plotly,
        ticklabelposition = "outside right"
      ),
      xaxis4 = list(
        showticklabels = FALSE,
        showline = TRUE,
        linecolor = "black",
        linewidth = 1,
        range = c(0.5, 9.5),
        tickangle = tick_angle_plotly,
        ticklabelposition = "outside right"
      )
    )
  }

  # NEW: RT-switching by splice junction type across cells (all and unique junctions)
  if ("RTS_junction" %in% colnames(Junctions)) {
    # Normalize boolean
    rts_bool <- tolower(as.character(Junctions$RTS_junction)) %in% c("true", "t", "1", "yes")
    junc_rt <- Junctions %>%
      dplyr::filter(CB != "unassigned") %>%
      mutate(
        SJ_type = paste(junction_category, canonical, sep = "_"),
        RTS_bool = rts_bool
      )

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

    # Per-cell percentages for ALL junctions
    all_junc_by_cell <- junc_rt %>%
      group_by(CB, SJ_type) %>%
      summarise(total = dplyr::n(), rts = sum(RTS_bool, na.rm = TRUE), .groups = "drop") %>%
      tidyr::complete(CB, SJ_type = sj_levels, fill = list(total = 0, rts = 0)) %>%
      mutate(perc = ifelse(total > 0, 100 * rts / total, 0))

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

    # Build a robust unique junction label if possible
    if (!("junctionLabel" %in% colnames(junc_rt))) {
      if (all(c("chrom", "strand", "genomic_start_coord", "genomic_end_coord") %in% colnames(junc_rt))) {
        junc_rt$junctionLabel <- with(junc_rt, paste(chrom, strand, genomic_start_coord, genomic_end_coord, sep = "_"))
      } else if (all(c("chrom", "strand", "genomic_start", "genomic_end") %in% colnames(junc_rt))) {
        junc_rt$junctionLabel <- with(junc_rt, paste(chrom, strand, genomic_start, genomic_end, sep = "_"))
      } else if ("junction_id" %in% colnames(junc_rt)) {
        junc_rt$junctionLabel <- junc_rt$junction_id
      } else {
        # Fallback to row index within CB as unique proxy
        junc_rt$junctionLabel <- paste0("jl_", seq_len(nrow(junc_rt)))
      }
    }

    # Per-cell percentages for UNIQUE junctions (deduplicate by genomic coordinates per cell & SJ type)
    uniq_junc_by_cell <- junc_rt %>%
      group_by(CB, SJ_type, junctionLabel) %>%
      summarise(rts_any = any(RTS_bool, na.rm = TRUE), .groups = "drop") %>%
      group_by(CB, SJ_type) %>%
      summarise(total = dplyr::n(), rts = sum(rts_any), .groups = "drop") %>%
      tidyr::complete(CB, SJ_type = sj_levels, fill = list(total = 0, rts = 0)) %>%
      mutate(perc = ifelse(total > 0, 100 * rts / total, 0))

    df_long_uniq <- data.frame(
      Variable = factor(uniq_junc_by_cell$SJ_type, levels = sj_levels),
      Value = uniq_junc_by_cell$perc
    )

    gg_rts_unique_by_sjtype <<- build_violin_plot(
      df_long = df_long_uniq,
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
    message("RTS_junction column not found in Junctions. Skipping RT-switching by SJ type plots.")
  }

  # Create grouped violins for % reads with all canonical junctions by structural category (HTML)
  if (!exists("gg_allcanon_by_category")) {
    cls2 <- Classification_file %>% dplyr::filter(CB != "unassigned")
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
    cls2 <- cls2 %>% mutate(allcanon_status = status_map(all_canonical))

    cat_key_map <- structural_category_map
    all_cats <- structural_category_levels
    bin_pretty_map <- c(FSM = "FSM", ISM = "ISM", NIC = "NIC", NNC = "NNC", Genic = "Genic Genomic", Antisense = "Antisense", Fusion = "Fusion", Intergenic = "Intergenic", Genic_intron = "Genic Intron")
    pretty_levels <- c("FSM", "ISM", "NIC", "NNC", "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic Intron")

    df_allcanon <- cls2 %>%
      mutate(cat_key = unname(cat_key_map[structural_category])) %>%
      filter(!is.na(cat_key), !is.na(allcanon_status)) %>%
      group_by(CB, cat_key, allcanon_status) %>%
      summarise(n = dplyr::n(), .groups = "drop") %>%
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

    if (length(coding_cols) > 0) {
      # Coding
      df_coding <- SQANTI_cell_summary %>%
        select(CB, all_of(coding_cols)) %>%
        pivot_longer(cols = all_of(coding_cols), names_to = "Variable", values_to = "Value") %>%
        mutate(
          Variable = gsub("_coding_prop$", "", Variable)
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

      df_coding$PrettyVar <- sapply(df_coding$Variable, tag_to_pretty)

      # Filter to only known categories
      known_vars <- c("FSM", "ISM", "NIC", "NNC", "genic", "antisense", "fusion", "intergenic", "genic_intron")
      df_coding <- df_coding %>% filter(Variable %in% known_vars)

      # Set factor levels
      pretty_levels <- c("FSM", "ISM", "NIC", "NNC", "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic Intron")
      df_coding$PrettyVar <- factor(df_coding$PrettyVar, levels = pretty_levels)

      # Create plot
      gg_coding_by_category <<- build_violin_plot(
        df_long = data.frame(Variable = df_coding$PrettyVar, Value = df_coding$Value),
        title = paste("Coding", entity_label_plural, "Distribution by Structural Category Across Cells"),
        x_labels = levels(df_coding$PrettyVar),
        fill_map = cat_fill_map,
        y_label = paste(entity_label_plural, ", %", sep = ""),
        legend = FALSE,
        override_outline_vars = c("Genic Genomic"),
        violin_alpha = 0.7,
        box_alpha = 0.3,
        box_width = 0.05,
        x_tickangle = 45,
        violin_outline_fill = TRUE,
        box_outline_default = "grey20",
        ylim = c(0, 100)
      )

      # Non-Coding
      # non_coding_cols already defined above
      df_noncoding <- SQANTI_cell_summary %>%
        select(CB, all_of(non_coding_cols)) %>%
        pivot_longer(cols = all_of(non_coding_cols), names_to = "Variable", values_to = "Value") %>%
        mutate(
          Variable = gsub("_non_coding_prop$", "", Variable)
        )

      df_noncoding$PrettyVar <- sapply(df_noncoding$Variable, tag_to_pretty)
      df_noncoding <- df_noncoding %>% filter(Variable %in% known_vars)
      df_noncoding$PrettyVar <- factor(df_noncoding$PrettyVar, levels = pretty_levels)

      gg_noncoding_by_category <<- build_violin_plot(
        df_long = data.frame(Variable = df_noncoding$PrettyVar, Value = df_noncoding$Value),
        title = paste("Non-Coding", entity_label_plural, "Distribution by Structural Category Across Cells"),
        x_labels = levels(df_noncoding$PrettyVar),
        fill_map = cat_fill_map,
        y_label = paste(entity_label_plural, ", %", sep = ""),
        legend = FALSE,
        override_outline_vars = c("Genic Genomic"),
        violin_alpha = 0.7,
        box_alpha = 0.3,
        box_width = 0.05,
        x_tickangle = 45,
        violin_outline_fill = TRUE,
        box_outline_default = "grey20",
        ylim = c(0, 100)
      )
    }

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
    # Add cover page
    grid.newpage()
    title_text <- if (mode == "isoforms") "SQANTI-single cell\nisoforms report" else "SQANTI-single cell\nreads report"
    cover <- textGrob(title_text,
      gp = gpar(fontface = "italic", fontsize = 40, col = "orangered")
    )
    grid.draw(cover)
    # Bulk tables
    s <- textGrob("Bulk summary", gp = gpar(fontface = "italic", fontsize = 30), vjust = 0)
    grid.arrange(s)

    # Calculate bulk-level stats
    if (mode == "isoforms") {
      total_reads_count <- sum(Classification_file$count, na.rm = TRUE)
    } else {
      total_reads_count <- nrow(Classification_file)
    }
    unique_genes <- length(unique(Classification_file$associated_gene))
    if (mode == "isoforms") {
      unique_junctions <- 0
    } else {
      unique_junctions <- length(unique(Classification_file$jxn_string))
    }

    # Gene Classification table
    gene_class_table <- data.frame(
      Category = c("Annotated Genes", "Novel Genes"),
      "Genes, count" = c(
        length(unique(Classification_file$associated_gene[!grepl("^novel", Classification_file$associated_gene)])),
        length(unique(Classification_file$associated_gene[grepl("^novel", Classification_file$associated_gene)]))
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
      colnames(read_class_table) <- c("Category", paste0(entity_label_plural, ", count"))
      # Ensure all levels are present (aggregate might drop empty ones if not careful, but factor levels help)
      # Actually aggregate returns only present levels. Let's use complete.
      read_class_table <- data.frame(Category = read_cat_levels) %>%
        left_join(read_class_table, by = "Category")
      read_class_table[is.na(read_class_table)] <- 0
    } else {
      read_class_table <- as.data.frame(table(factor(Classification_file$structural_category, levels = read_cat_levels)))
      colnames(read_class_table) <- c("Category", paste0(entity_label_plural, ", count"))
    }
    read_class_table$Category <- read_cat_names

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
      Percent = sj_perc,
      check.names = FALSE
    )
    rownames(SJ_class_table) <- NULL

    big_table_theme <- ttheme_default(
      core = list(fg_params = list(cex = 1.5)),
      colhead = list(fg_params = list(cex = 1.5, fontface = "bold"))
    )

    title_genes <- textGrob("Gene Classification", gp = gpar(fontface = "italic", fontsize = 24), vjust = -3)
    title_reads <- textGrob(paste(entity_label, "Classification"), gp = gpar(fontface = "italic", fontsize = 24), vjust = -7.7)
    title_sj <- textGrob("Splice Junction Classification", gp = gpar(fontface = "italic", fontsize = 24), vjust = -4.3)

    table_genes <- tableGrob(gene_class_table, rows = NULL, theme = big_table_theme)
    table_reads <- tableGrob(read_class_table, rows = NULL, theme = big_table_theme)
    table_sj <- tableGrob(SJ_class_table, rows = NULL, theme = big_table_theme)

    unique_counts_text <- if (mode == "isoforms") {
      sprintf(
        "Number of %s: %d\nUnique Genes: %d",
        entity_label_plural, total_reads_count, unique_genes
      )
    } else {
      sprintf(
        "Number of %s: %d\nUnique Genes: %d\nUnique Junction Chains: %d",
        entity_label_plural, total_reads_count, unique_genes, unique_junctions
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
      SD = sd(SQANTI_cell_summary$Genes_in_cell, na.rm = TRUE)
    )
    unique_junctions_stats <- c(
      Mean = mean(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE),
      Median = median(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE),
      Min = min(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE),
      Max = max(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE),
      SD = sd(SQANTI_cell_summary$UJCs_in_cell, na.rm = TRUE)
    )
    reads_stats <- c(
      Mean = mean(SQANTI_cell_summary[[count_col]], na.rm = TRUE),
      Median = median(SQANTI_cell_summary[[count_col]], na.rm = TRUE),
      Min = min(SQANTI_cell_summary[[count_col]], na.rm = TRUE),
      Max = max(SQANTI_cell_summary[[count_col]], na.rm = TRUE),
      SD = sd(SQANTI_cell_summary[[count_col]], na.rm = TRUE)
    )
    umis_stats <- c(
      Mean = mean(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE),
      Median = median(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE),
      Min = min(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE),
      Max = max(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE),
      SD = sd(SQANTI_cell_summary$UMIs_in_cell, na.rm = TRUE)
    )
    summary_table1 <- data.frame(
      Feature = c(paste(entity_label_plural, "in cell"), "UMIs in cell", "Unique Genes", "Unique Junction Chains"),
      Mean = c(reads_stats["Mean"], umis_stats["Mean"], unique_genes_stats["Mean"], unique_junctions_stats["Mean"]),
      Median = c(reads_stats["Median"], umis_stats["Median"], unique_genes_stats["Median"], unique_junctions_stats["Median"]),
      Min = c(reads_stats["Min"], umis_stats["Min"], unique_genes_stats["Min"], unique_junctions_stats["Min"]),
      Max = c(reads_stats["Max"], umis_stats["Max"], unique_genes_stats["Max"], unique_junctions_stats["Max"]),
      SD = c(reads_stats["SD"], umis_stats["SD"], unique_genes_stats["SD"], unique_junctions_stats["SD"])
    )
    # If isoforms mode, drop Unique Junction Chains from summary table
    if (mode == "isoforms") {
      summary_table1 <- summary_table1[!(summary_table1$Feature %in% c("Unique Junction Chains", "UMIs in cell")), , drop = FALSE]
    }
    summary_table1[, 2:6] <- round(summary_table1[, 2:6], 3)
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
        SD = sd(vals, na.rm = TRUE)
      )
    })
    count_stats_df <- data.frame(
      Category = struct_cat_names,
      t(count_stats)
    )
    colnames(count_stats_df)[2:6] <- c("Mean", "Median", "Min", "Max", "SD")
    count_stats_df[, 2:6] <- round(count_stats_df[, 2:6], 3)
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
        SD = sd(vals, na.rm = TRUE)
      )
    })
    prop_stats_df <- data.frame(
      Category = struct_cat_names,
      t(prop_stats)
    )
    colnames(prop_stats_df)[2:6] <- c("Mean", "Median", "Min", "Max", "SD")
    prop_stats_df[, 2:6] <- round(prop_stats_df[, 2:6], 3)
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
    render_pdf_plot_centered("gg_JCs_in_cell", width_frac = 0.5)

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
    # UJCs per Gene
    if (mode != "isoforms") {
      render_pdf_plot("gg_ujc_bins_all")
      render_pdf_plot("gg_ujc_bins")
    }
    # Mitochondrial genes
    render_pdf_plot("gg_MT_perc")

    # Read Length Characterization section
    section_page(paste(entity_label, "Length Characterization"))
    # Bulk Length Distribution
    render_pdf_plot("gg_bulk_all_reads", converter = NULL)
    render_pdf_plot("gg_bulk_length_by_category", converter = NULL)
    render_pdf_plot("gg_bulk_length_by_exon_type", converter = NULL)
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
          print(plotly_to_ggplot(gg_exon_profile_by_category[[nm]]))
        }
      }
    }

    # Coding and Non-Coding Distributions
    if (exists("gg_coding_by_category")) render_pdf_plot("gg_coding_by_category")
    if (exists("gg_noncoding_by_category")) render_pdf_plot("gg_noncoding_by_category")

    # Splice Junction Characterization section
    # Compute per-structural-category SJ distributions for PDF pages
    junc_aug <- Junctions %>%
      dplyr::filter(CB != "unassigned") %>%
      mutate(junction_type = paste(junction_category, canonical, sep = "_"))

    # Add count column for weighted quantification
    if (mode == "isoforms" && "FL" %in% colnames(Classification_file)) {
      Classification_file$count <- Classification_file$FL
    } else {
      Classification_file$count <- 1
    }

    if (!("structural_category" %in% colnames(junc_aug))) {
      # Try to bring structural_category from Classification_file by common key + CB
      join_key <- NULL
      for (k in c("isoform", "readID", "read_id", "ID", "read_name", "read")) {
        if (k %in% colnames(Junctions) && k %in% colnames(Classification_file)) {
          join_key <- k
          break
        }
      }
      if (!is.null(join_key)) {
        by_vec <- c(CB = "CB")
        by_vec[[join_key]] <- join_key
        junc_aug <- junc_aug %>%
          left_join(Classification_file %>% select(all_of(c(join_key, "CB", "structural_category"))), by = by_vec)
      } else {
        junc_aug$structural_category <- NA_character_
        junc_aug$count <- 1
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
      tidyr::complete(CB, cat_key = all_cats, fill = list(total = 0, known_canonical = 0, known_non_canonical = 0, novel_canonical = 0, novel_non_canonical = 0)) %>%
      mutate(
        KnownCanonicalPerc = ifelse(total > 0, 100 * known_canonical / total, 0),
        KnownNonCanonicalPerc = ifelse(total > 0, 100 * known_non_canonical / total, 0),
        NovelCanonicalPerc = ifelse(total > 0, 100 * novel_canonical / total, 0),
        NovelNonCanonicalPerc = ifelse(total > 0, 100 * novel_non_canonical / total, 0)
      ) %>%
      ungroup()

    # Prepare plotting helpers
    fill_map_cat <- cat_fill_map
    x_labels_full <- cat_labels_pretty
    make_df_long <- function(col_name) {
      data.frame(Variable = factor(junc_summ$cat_key, levels = all_cats), Value = junc_summ[[col_name]])
    }

    p_known_can <- build_violin_plot_ggplot(
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
    p_known_noncan <- build_violin_plot_ggplot(
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
    p_novel_can <- build_violin_plot_ggplot(
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
    p_novel_noncan <- build_violin_plot_ggplot(
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

    dev.off()
  }
}

Classification <- read.table(class.file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (mode == "isoforms" && "FL" %in% colnames(Classification)) {
  Classification$count <- sapply(strsplit(as.character(Classification$FL), ","), function(x) sum(as.numeric(x), na.rm = TRUE))
  Classification$count[is.na(Classification$count) | Classification$count == 0] <- 1
} else {
  Classification$count <- 1
}
Junctions <- read.table(junc.file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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
    # We need to be careful about duplicates if Junctions has multiple rows per isoform (it does, one per junction)
    # We want to assign the isoform's count to each junction row?
    # Yes, because when we sum junctions, we want (count * junctions_per_isoform).
    # Wait, if we sum *types* of junctions (e.g. canonical), we want sum(count) for all junctions of that type.
    # So yes, each junction row should have the isoform's count.

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
  message("Using precomputed cell summary: ", cell_summary_path)
  SQANTI_cell_summary <- read.table(cell_summary_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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

  rmd_file <- file.path(script_dir, "SQANTI-sc_reads_report.Rmd")
  css_file <- file.path(script_dir, "style.css")

  # Check if Rmd file exists
  if (!file.exists(rmd_file)) {
    stop(
      "R Markdown file not found: ", rmd_file,
      "\nPlease ensure SQANTI-sc_reads_report.Rmd is in the same directory as this script."
    )
  }

  # Copy CSS file to output directory if it exists
  if (file.exists(css_file)) {
    css_output <- file.path(dirname(report_output), "style.css")
    file.copy(css_file, css_output, overwrite = TRUE)
    message("CSS file copied to output directory: ", css_output)
  }

  # Generate HTML report
  html_output_file <- paste0(report_output, ".html")

  message("Generating HTML report...")
  message("Rmd file: ", rmd_file)
  message("Output file: ", html_output_file)
  message("Output directory: ", dirname(report_output))

  # Check if plot objects exist
  plot_objects <- c("gg_reads_in_cells", "gg_umis_in_cells", "gg_genes_in_cells", "gg_JCs_in_cell")
  for (obj in plot_objects) {
    if (exists(obj)) {
      message("Plot object '", obj, "' exists")
    } else {
      message("Plot object '", obj, "' does not exist")
    }
  }

  rmarkdown::render(
    input = rmd_file,
    output_file = html_output_file,
    output_dir = dirname(report_output),
    intermediates_dir = dirname(report_output),
    quiet = FALSE,
    envir = globalenv()
  )

  message("HTML report generated: ", html_output_file)
}
