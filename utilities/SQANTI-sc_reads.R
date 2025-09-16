#!/usr/env/bin Rscript

######################################################
##### SQANTI single-cell reads report generation #####
######################################################



### Author: Juan Francisco Cervilla & Carlos Blanco

#********************** Packages 

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(grid)
library(gridExtra)

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

# Check for optional arguments
if (length(args) > 5) {
  for (arg in args[6:length(args)]) {
    if (arg == "--ignore_cell_summary") {
      ignore_cell_summary <- TRUE
    }
    if (arg == "--skipORF") {
      skipORF <- TRUE
    }
    if (arg == "--CAGE_peak") {
      CAGE_peak <- TRUE
    }
    if (arg == "--polyA_motif_list") {
      polyA_motif_list <- TRUE
    }
  }
}

# Validate arguments
if (length(args) < 5) {
  stop("Incorrect number of arguments! Required: [classification file] [junc file] [report format] [outputPathPrefix] [mode]. Abort!")
}

if (!(report.format %in% c('pdf', 'html', 'both'))) {
  stop("Report format needs to be: pdf, html, or both. Abort!")
}

# Validate mode argument
if (!(mode %in% c('reads', 'isoforms'))) {
  stop("Mode needs to be: reads or isoforms. Abort!")
}

# Set labels based on mode
if (mode == "isoforms") {
  entity_label <- "Isoform"
  entity_label_plural <- "Isoforms"
} else {
  entity_label <- "Read"
  entity_label_plural <- "Reads"
}

# Print cell summary saving status
if (ignore_cell_summary) {
  print("Cell summary table will not be saved (--ignore_cell_summary flag is active).")
} else {
  print("Cell summary table will be saved.")
}

# Call the function with the appropriate Save parameter
save_option <- ifelse(ignore_cell_summary, "N", "Y")

# Generate output file names with full paths
cell_summary_output <- file.path(paste0(outputPathPrefix, "_SQANTI_cell_summary"))
report_output <- file.path(paste0(outputPathPrefix, "_SQANTI_sc_report_reads"))

calculate_metrics_per_cell <- function(Classification, Junctions, cell_summary_output, Save){

  # Helper function to safely calculate proportions
  safe_prop <- function(numerator, denominator, default_val = 0) {
    prop <- ifelse(denominator > 0, (numerator / denominator) * 100, default_val)
    prop[is.na(prop) | is.infinite(prop)] <- default_val
    return(prop)
  }

  if (mode == "isoforms") {
  print("Isoform mode detected: expanding comma-separated cell barcodes.")
  Classification <- Classification %>%
    filter(!is.na(CB) & CB != "") %>%
    separate_rows(CB, sep = ",")
  }

  

  # Filter out reads with no CB and print a warning if any exist
  if (any(Classification$CB == "")) {
    print("There are molecules with no cell barcode assigned and they will not be considered. Check your classification file.")
    Classification <- Classification %>% filter(CB != "")
  }


  # Initial summary of reads/isoforms and UMIs per cell
  if (mode == "isoforms") {
      SQANTI_cell_summary_base <- Classification %>%
        group_by(CB) %>%
        summarise(
          total_reads = n(), # This is total isoforms per cell
          total_UMI = 0, # UMI is not applicable in isoform mode
          .groups = 'drop'
        )
  } else {
      SQANTI_cell_summary_base <- Classification %>%
        group_by(CB) %>%
        summarise(
          total_reads = n(),
          total_UMI = n_distinct(UMI),
          .groups = 'drop'
        )
  }

  # Reads no monoexon
  reads_no_monoexon_per_cell <- Classification %>%
    filter(!exons == 1) %>%
    group_by(CB) %>%
    summarise(total_reads_no_monoexon = n(), .groups = 'drop')

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary_base, reads_no_monoexon_per_cell, by = "CB") %>%
    mutate(total_reads_no_monoexon = ifelse(is.na(total_reads_no_monoexon), 0, total_reads_no_monoexon))


  # Structural category counts
  structural_categories <- c('full-splice_match', 'incomplete-splice_match', 'novel_in_catalog',
                             'novel_not_in_catalog', 'genic', 'antisense', 'fusion',
                             'intergenic', 'genic_intron')

  category_counts <- Classification %>%
    group_by(CB, structural_category) %>%
    summarise(count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = structural_category, values_from = count, values_fill = 0)

  # Ensure all category columns exist
  for (cat_col in structural_categories) {
    if (!cat_col %in% names(category_counts)) {
      category_counts[[cat_col]] <- 0
    }
  }
  # Rename columns for clarity, e.g. FSM_count
  category_counts <- category_counts %>%
    rename_with(~ paste0(toupper(gsub("[-_]", "", .)), "_count"), all_of(structural_categories))


  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, category_counts, by = "CB")
  # Fill NA counts with 0 for cells that might not have any reads in certain categories
  for (col_name in names(category_counts %>% select(-CB))) { # Iterate over new count columns
    SQANTI_cell_summary[[col_name]] <- ifelse(is.na(SQANTI_cell_summary[[col_name]]), 0, SQANTI_cell_summary[[col_name]])
  }


  # Genes per cell
  genes_per_cell <- Classification %>%
    group_by(CB) %>%
    summarise(genes_in_cell = n_distinct(associated_gene), .groups = 'drop')
  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, genes_per_cell, by = "CB")

  # Models (unique junction chains) in multiexonic transcripts per cell
  models_per_cell <- Classification %>%
    filter(exons > 1) %>%
    group_by(CB) %>%
    summarise(models_in_cell = n_distinct(jxn_string), .groups = 'drop') # Sum of distinct jxn_string per gene is not needed if we count across all genes in a cell

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, models_per_cell, by = "CB") %>%
    mutate(models_in_cell = ifelse(is.na(models_in_cell), 0, models_in_cell))


  # Mitochondrial reads
  mt_reads_per_cell <- Classification %>%
    filter(chrom == "MT") %>%
    group_by(CB) %>%
    summarise(MT_reads_count = n(), .groups = 'drop')

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, mt_reads_per_cell, by = "CB") %>%
    mutate(
      MT_reads_count = ifelse(is.na(MT_reads_count), 0, MT_reads_count),
      MT_perc = safe_prop(MT_reads_count, total_reads)
    )

  # Novel vs annotated genes
  annotated_novel_genes <- Classification %>%
    group_by(CB) %>%
    summarise(
      annotated_genes = n_distinct(associated_gene[!grepl("^novel", associated_gene)]),
      novel_genes = n_distinct(associated_gene[grepl("^novel", associated_gene)]),
      .groups = 'drop'
    )
  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, annotated_novel_genes, by = "CB")


  # Binned read count metrics for annotated and novel genes
  gene_counts_per_cell <- Classification %>%
    group_by(CB, associated_gene) %>%
    summarise(read_count = n(), .groups = 'drop')

  binned_read_counts <- gene_counts_per_cell %>%
    mutate(gene_type = ifelse(grepl("^novel", associated_gene), "novel", "annotated")) %>%
    group_by(CB, gene_type) %>%
    summarise(
      bin1_count = sum(read_count == 1),
      bin2_3_count = sum(read_count >= 2 & read_count <= 3),
      bin4_5_count = sum(read_count >= 4 & read_count <= 5),
      bin6plus_count = sum(read_count >= 6),
      total_genes_in_type = n_distinct(associated_gene),
      .groups = 'drop'
    ) %>%
    mutate(
      anno_bin1_perc = ifelse(gene_type == "annotated", safe_prop(bin1_count, total_genes_in_type), NA),
      anno_bin2_3_perc = ifelse(gene_type == "annotated", safe_prop(bin2_3_count, total_genes_in_type), NA),
      anno_bin4_5_perc = ifelse(gene_type == "annotated", safe_prop(bin4_5_count, total_genes_in_type), NA),
      anno_bin6plus_perc = ifelse(gene_type == "annotated", safe_prop(bin6plus_count, total_genes_in_type), NA),
      novel_bin1_perc = ifelse(gene_type == "novel", safe_prop(bin1_count, total_genes_in_type), NA),
      novel_bin2_3_perc = ifelse(gene_type == "novel", safe_prop(bin2_3_count, total_genes_in_type), NA),
      novel_bin4_5_perc = ifelse(gene_type == "novel", safe_prop(bin4_5_count, total_genes_in_type), NA),
      novel_bin6plus_perc = ifelse(gene_type == "novel", safe_prop(bin6plus_count, total_genes_in_type), NA)
    ) %>%
    group_by(CB) %>% # Group again to consolidate novel and annotated into one row per CB
    summarise(
      anno_bin1_perc = sum(anno_bin1_perc, na.rm = TRUE), # sum will pick the non-NA value
      anno_bin2_3_perc = sum(anno_bin2_3_perc, na.rm = TRUE),
      anno_bin4_5_perc = sum(anno_bin4_5_perc, na.rm = TRUE),
      anno_bin6plus_perc = sum(anno_bin6plus_perc, na.rm = TRUE),
      novel_bin1_perc = sum(novel_bin1_perc, na.rm = TRUE),
      novel_bin2_3_perc = sum(novel_bin2_3_perc, na.rm = TRUE),
      novel_bin4_5_perc = sum(novel_bin4_5_perc, na.rm = TRUE),
      novel_bin6plus_perc = sum(novel_bin6plus_perc, na.rm = TRUE),
      .groups = 'drop'
    )

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, binned_read_counts, by = "CB")
  # Fill NA percentages with 0
  binned_cols <- grep("bin\\d.*_perc$", names(SQANTI_cell_summary), value = TRUE)
  for (col_name in binned_cols) {
    SQANTI_cell_summary[[col_name]] <- ifelse(is.na(SQANTI_cell_summary[[col_name]]), 0, SQANTI_cell_summary[[col_name]])
  }

  # UJC bins for multiexonic reads
  gene_ujc_counts_per_cell <- Classification %>%
    filter(exons > 1) %>%
    group_by(CB, associated_gene) %>%
    summarise(ujc_count = n_distinct(jxn_string), .groups = 'drop')

  binned_ujc_counts <- gene_ujc_counts_per_cell %>%
    mutate(gene_type = ifelse(grepl("^novel", associated_gene), "novel", "annotated")) %>%
    group_by(CB, gene_type) %>%
    summarise(
      ujc_bin1_count = sum(ujc_count == 1),
      ujc_bin2_3_count = sum(ujc_count >= 2 & ujc_count <= 3),
      ujc_bin4_5_count = sum(ujc_count >= 4 & ujc_count <= 5),
      ujc_bin6plus_count = sum(ujc_count >= 6),
      total_genes_in_type_ujc = n_distinct(associated_gene),
      .groups = 'drop'
    ) %>%
    mutate(
      anno_ujc_bin1_perc = ifelse(gene_type == "annotated", safe_prop(ujc_bin1_count, total_genes_in_type_ujc), NA),
      anno_ujc_bin2_3_perc = ifelse(gene_type == "annotated", safe_prop(ujc_bin2_3_count, total_genes_in_type_ujc), NA),
      anno_ujc_bin4_5_perc = ifelse(gene_type == "annotated", safe_prop(ujc_bin4_5_count, total_genes_in_type_ujc), NA),
      anno_ujc_bin6plus_perc = ifelse(gene_type == "annotated", safe_prop(ujc_bin6plus_count, total_genes_in_type_ujc), NA),
      novel_ujc_bin1_perc = ifelse(gene_type == "novel", safe_prop(ujc_bin1_count, total_genes_in_type_ujc), NA),
      novel_ujc_bin2_3_perc = ifelse(gene_type == "novel", safe_prop(ujc_bin2_3_count, total_genes_in_type_ujc), NA),
      novel_ujc_bin4_5_perc = ifelse(gene_type == "novel", safe_prop(ujc_bin4_5_count, total_genes_in_type_ujc), NA),
      novel_ujc_bin6plus_perc = ifelse(gene_type == "novel", safe_prop(ujc_bin6plus_count, total_genes_in_type_ujc), NA)
    ) %>%
    group_by(CB) %>%
    summarise(
      anno_ujc_bin1_perc = sum(anno_ujc_bin1_perc, na.rm = TRUE),
      anno_ujc_bin2_3_perc = sum(anno_ujc_bin2_3_perc, na.rm = TRUE),
      anno_ujc_bin4_5_perc = sum(anno_ujc_bin4_5_perc, na.rm = TRUE),
      anno_ujc_bin6plus_perc = sum(anno_ujc_bin6plus_perc, na.rm = TRUE),
      novel_ujc_bin1_perc = sum(novel_ujc_bin1_perc, na.rm = TRUE),
      novel_ujc_bin2_3_perc = sum(novel_ujc_bin2_3_perc, na.rm = TRUE),
      novel_ujc_bin4_5_perc = sum(novel_ujc_bin4_5_perc, na.rm = TRUE),
      novel_ujc_bin6plus_perc = sum(novel_ujc_bin6plus_perc, na.rm = TRUE),
      .groups = 'drop'
    )
  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, binned_ujc_counts, by = "CB")
  binned_ujc_cols <- grep("ujc_bin\\d.*_perc$", names(SQANTI_cell_summary), value = TRUE)
  for (col_name in binned_ujc_cols) {
    SQANTI_cell_summary[[col_name]] <- ifelse(is.na(SQANTI_cell_summary[[col_name]]), 0, SQANTI_cell_summary[[col_name]])
  }


  # Known/novel canonical/non-canonical - calculated from Junctions file
  # Only keep junctions with a valid CB
  Junctions_valid <- Junctions[!is.na(Junctions$CB) & Junctions$CB != '', ]
  
  if (nrow(Junctions_valid) > 0) {
    # Create junction type by combining junction_category and canonical
    Junctions_valid$junction_type <- paste(Junctions_valid$junction_category, 
                                         Junctions_valid$canonical, 
                                         sep = "_")
    
    # Count junctions by CB and junction type
    junction_counts <- Junctions_valid %>%
      group_by(CB, junction_type) %>%
      summarise(count = n(), .groups = 'drop') %>%
      pivot_wider(names_from = junction_type, values_from = count, values_fill = 0)
    
    # Ensure all required columns exist (in case some junction types are missing)
    required_cols <- c("known_canonical", "known_non_canonical", "novel_canonical", "novel_non_canonical")
    for (col in required_cols) {
      if (!col %in% colnames(junction_counts)) {
        junction_counts[[col]] <- 0
      }
    }
    
    # Rename columns to match expected naming convention
    junction_counts <- junction_counts %>%
      rename(
        known_canonical_count = known_canonical,
        known_non_canonical_count = known_non_canonical,
        novel_canonical_count = novel_canonical,
        novel_non_canonical_count = novel_non_canonical
      )
    
    # Calculate total junctions per cell
    junction_counts$total_junctions <- rowSums(junction_counts[, c("known_canonical_count", "known_non_canonical_count", 
                                                                  "novel_canonical_count", "novel_non_canonical_count")])
    
    # Calculate proportions
    junction_counts <- junction_counts %>%
      mutate(
        known_canonical_prop = safe_prop(known_canonical_count, total_junctions),
        known_non_canonical_prop = safe_prop(known_non_canonical_count, total_junctions),
        novel_canonical_prop = safe_prop(novel_canonical_count, total_junctions),
        novel_non_canonical_prop = safe_prop(novel_non_canonical_count, total_junctions)
      )
    
    # Join with cell summary
    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, junction_counts, by = "CB")
  } else {
    # If no valid junctions, create empty columns
    SQANTI_cell_summary$known_canonical_count <- 0
    SQANTI_cell_summary$known_non_canonical_count <- 0
    SQANTI_cell_summary$novel_canonical_count <- 0
    SQANTI_cell_summary$novel_non_canonical_count <- 0
    SQANTI_cell_summary$total_junctions <- 0
    SQANTI_cell_summary$known_canonical_prop <- 0
    SQANTI_cell_summary$known_non_canonical_prop <- 0
    SQANTI_cell_summary$novel_canonical_prop <- 0
    SQANTI_cell_summary$novel_non_canonical_prop <- 0
  }
  
  # Fill NA values with 0 for any cells that didn't have junctions
  SQANTI_cell_summary$known_canonical_count[is.na(SQANTI_cell_summary$known_canonical_count)] <- 0
  SQANTI_cell_summary$known_non_canonical_count[is.na(SQANTI_cell_summary$known_non_canonical_count)] <- 0
  SQANTI_cell_summary$novel_canonical_count[is.na(SQANTI_cell_summary$novel_canonical_count)] <- 0
  SQANTI_cell_summary$novel_non_canonical_count[is.na(SQANTI_cell_summary$novel_non_canonical_count)] <- 0
  SQANTI_cell_summary$total_junctions[is.na(SQANTI_cell_summary$total_junctions)] <- 0
  SQANTI_cell_summary$known_canonical_prop[is.na(SQANTI_cell_summary$known_canonical_prop)] <- 0
  SQANTI_cell_summary$known_non_canonical_prop[is.na(SQANTI_cell_summary$known_non_canonical_prop)] <- 0
  SQANTI_cell_summary$novel_canonical_prop[is.na(SQANTI_cell_summary$novel_canonical_prop)] <- 0
  SQANTI_cell_summary$novel_non_canonical_prop[is.na(SQANTI_cell_summary$novel_non_canonical_prop)] <- 0


  # Sqanti category proportions
  # Denominator is total_reads
  for (cat_name in structural_categories) {
    count_col <- paste0(toupper(gsub("[-_]", "", cat_name)), "_count")
    prop_col <- paste0(toupper(gsub("[-_]", "", cat_name)), "_prop")
    SQANTI_cell_summary[[prop_col]] <- safe_prop(SQANTI_cell_summary[[count_col]], SQANTI_cell_summary$total_reads)
  }


  # Subcategory proportions
  # FSM subcategories
  fsm_subcat_levels <- c("alternative_3end", "alternative_3end5end", "alternative_5end", "reference_match", "mono-exon")
  fsm_subcat_props <- Classification %>%
    filter(structural_category == "full-splice_match") %>%
    group_by(CB, subcategory) %>%
    summarise(reads_in_subcat = n(), .groups = 'drop') %>%
    complete(CB = unique(SQANTI_cell_summary$CB), subcategory = factor(fsm_subcat_levels, levels = fsm_subcat_levels), fill = list(reads_in_subcat = 0)) %>% # Ensure all CBs and subcats
    left_join(SQANTI_cell_summary %>% select(CB, FSM_count = FULLSPLICEMATCH_count), by = "CB") %>% # FSM_count comes from category_counts
    mutate(prop_subcat = safe_prop(reads_in_subcat, FSM_count)) %>%
    select(CB, subcategory, prop_subcat) %>%
    pivot_wider(names_from = subcategory, values_from = prop_subcat, names_prefix = "sub_FSM_", values_fill = 0)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, fsm_subcat_props, by = "CB")
  for(col in names(fsm_subcat_props %>% select(-CB))) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])


  # ISM subcategories
  ism_subcat_levels <- c("3prime_fragment", "internal_fragment", "5prime_fragment", "intron_retention", "mono-exon")
  ism_subcat_props <- Classification %>%
    filter(structural_category == "incomplete-splice_match") %>%
    group_by(CB, subcategory) %>%
    summarise(reads_in_subcat = n(), .groups = 'drop') %>%
    complete(CB = unique(SQANTI_cell_summary$CB), subcategory = factor(ism_subcat_levels, levels = ism_subcat_levels), fill = list(reads_in_subcat = 0)) %>%
    left_join(SQANTI_cell_summary %>% select(CB, ISM_count = INCOMPLETESPLICEMATCH_count), by = "CB") %>%
    mutate(prop_subcat = safe_prop(reads_in_subcat, ISM_count)) %>%
    select(CB, subcategory, prop_subcat) %>%
    pivot_wider(names_from = subcategory, values_from = prop_subcat, names_prefix = "sub_ISM_", values_fill = 0)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, ism_subcat_props, by = "CB")
  for(col in names(ism_subcat_props %>% select(-CB))) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])

  # NIC subcategories
  nic_subcat_levels <- c("combination_of_known_junctions", "combination_of_known_splicesites", "intron_retention", "mono-exon_by_intron_retention", "mono-exon")
  nic_subcat_props <- Classification %>%
    filter(structural_category == "novel_in_catalog") %>%
    group_by(CB, subcategory) %>%
    summarise(reads_in_subcat = n(), .groups = 'drop') %>%
    complete(CB = unique(SQANTI_cell_summary$CB), subcategory = factor(nic_subcat_levels, levels = nic_subcat_levels), fill = list(reads_in_subcat = 0)) %>%
    left_join(SQANTI_cell_summary %>% select(CB, NIC_count = NOVELINCATALOG_count), by = "CB") %>%
    mutate(prop_subcat = safe_prop(reads_in_subcat, NIC_count)) %>%
    select(CB, subcategory, prop_subcat) %>%
    pivot_wider(names_from = subcategory, values_from = prop_subcat, names_prefix = "sub_NIC_", values_fill = 0)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, nic_subcat_props, by = "CB")
  for(col in names(nic_subcat_props %>% select(-CB))) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])

  # NNC subcategories
  nnc_subcat_levels <- c("at_least_one_novel_splicesite", "intron_retention")
  nnc_subcat_props <- Classification %>%
    filter(structural_category == "novel_not_in_catalog") %>%
    group_by(CB, subcategory) %>%
    summarise(reads_in_subcat = n(), .groups = 'drop') %>%
    complete(CB = unique(SQANTI_cell_summary$CB), subcategory = factor(nnc_subcat_levels, levels = nnc_subcat_levels), fill = list(reads_in_subcat = 0)) %>%
    left_join(SQANTI_cell_summary %>% select(CB, NNC_count = NOVELNOTINCATALOG_count), by = "CB") %>%
    mutate(prop_subcat = safe_prop(reads_in_subcat, NNC_count)) %>%
    select(CB, subcategory, prop_subcat) %>%
    pivot_wider(names_from = subcategory, values_from = prop_subcat, names_prefix = "sub_NNC_", values_fill = 0)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, nnc_subcat_props, by = "CB")
  for(col in names(nnc_subcat_props %>% select(-CB))) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])


  # Genic subcategories
  genic_subcat_levels <- c("mono-exon", "multi-exon")
  genic_subcat_props <- Classification %>%
    filter(structural_category == "genic") %>%
    group_by(CB, subcategory) %>%
    summarise(reads_in_subcat = n(), .groups = 'drop') %>%
    complete(CB = unique(SQANTI_cell_summary$CB), subcategory = factor(genic_subcat_levels, levels = genic_subcat_levels), fill = list(reads_in_subcat = 0)) %>%
    left_join(SQANTI_cell_summary %>% select(CB, Genic_count = GENIC_count), by = "CB") %>% # Name from category_counts
    mutate(prop_subcat = safe_prop(reads_in_subcat, Genic_count)) %>%
    select(CB, subcategory, prop_subcat) %>%
    pivot_wider(names_from = subcategory, values_from = prop_subcat, names_prefix = "sub_genic_", values_fill = 0)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, genic_subcat_props, by = "CB")
  for(col in names(genic_subcat_props %>% select(-CB))) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])


  # Antisense subcategories
  antisense_subcat_levels <- c("mono-exon", "multi-exon")
  antisense_subcat_props <- Classification %>%
    filter(structural_category == "antisense") %>%
    group_by(CB, subcategory) %>%
    summarise(reads_in_subcat = n(), .groups = 'drop') %>%
    complete(CB = unique(SQANTI_cell_summary$CB), subcategory = factor(antisense_subcat_levels, levels = antisense_subcat_levels), fill = list(reads_in_subcat = 0)) %>%
    left_join(SQANTI_cell_summary %>% select(CB, Antisense_count = ANTISENSE_count), by = "CB") %>%
    mutate(prop_subcat = safe_prop(reads_in_subcat, Antisense_count)) %>%
    select(CB, subcategory, prop_subcat) %>%
    pivot_wider(names_from = subcategory, values_from = prop_subcat, names_prefix = "sub_antisense_", values_fill = 0)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, antisense_subcat_props, by = "CB")
  for(col in names(antisense_subcat_props %>% select(-CB))) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])

  # Fusion subcategories
  fusion_subcat_levels <- c("intron_retention", "multi-exon")
  fusion_subcat_props <- Classification %>%
    filter(structural_category == "fusion") %>%
    group_by(CB, subcategory) %>%
    summarise(reads_in_subcat = n(), .groups = 'drop') %>%
    complete(CB = unique(SQANTI_cell_summary$CB), subcategory = factor(fusion_subcat_levels, levels = fusion_subcat_levels), fill = list(reads_in_subcat = 0)) %>%
    left_join(SQANTI_cell_summary %>% select(CB, Fusion_count = FUSION_count), by = "CB") %>%
    mutate(prop_subcat = safe_prop(reads_in_subcat, Fusion_count)) %>%
    select(CB, subcategory, prop_subcat) %>%
    pivot_wider(names_from = subcategory, values_from = prop_subcat, names_prefix = "sub_fusion_", values_fill = 0)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, fusion_subcat_props, by = "CB")
  for(col in names(fusion_subcat_props %>% select(-CB))) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])

  # Intergenic subcategories
  intergenic_subcat_levels <- c("mono-exon", "multi-exon")
  intergenic_subcat_props <- Classification %>%
    filter(structural_category == "intergenic") %>%
    group_by(CB, subcategory) %>%
    summarise(reads_in_subcat = n(), .groups = 'drop') %>%
    complete(CB = unique(SQANTI_cell_summary$CB), subcategory = factor(intergenic_subcat_levels, levels = intergenic_subcat_levels), fill = list(reads_in_subcat = 0)) %>%
    left_join(SQANTI_cell_summary %>% select(CB, Intergenic_count = INTERGENIC_count), by = "CB") %>%
    mutate(prop_subcat = safe_prop(reads_in_subcat, Intergenic_count)) %>%
    select(CB, subcategory, prop_subcat) %>%
    pivot_wider(names_from = subcategory, values_from = prop_subcat, names_prefix = "sub_intergenic_", values_fill = 0)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, intergenic_subcat_props, by = "CB")
  for(col in names(intergenic_subcat_props %>% select(-CB))) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])


  # Genic_intron subcategories
  genic_intron_subcat_levels <- c("mono-exon", "multi-exon")
  genic_intron_subcat_props <- Classification %>%
    filter(structural_category == "genic_intron") %>%
    group_by(CB, subcategory) %>%
    summarise(reads_in_subcat = n(), .groups = 'drop') %>%
    complete(CB = unique(SQANTI_cell_summary$CB), subcategory = factor(genic_intron_subcat_levels, levels = genic_intron_subcat_levels), fill = list(reads_in_subcat = 0)) %>%
    left_join(SQANTI_cell_summary %>% select(CB, Genic_intron_count = GENICINTRON_count), by = "CB") %>%
    mutate(prop_subcat = safe_prop(reads_in_subcat, Genic_intron_count)) %>%
    select(CB, subcategory, prop_subcat) %>%
    pivot_wider(names_from = subcategory, values_from = prop_subcat, names_prefix = "sub_genic_intron_", values_fill = 0)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, genic_intron_subcat_props, by = "CB")
  for(col in names(genic_intron_subcat_props %>% select(-CB))) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])


  # Read lengths general
  read_lengths_general <- Classification %>%
    group_by(CB) %>%
    summarise(
      two_fifty_length_reads_count = sum(length <= 250, na.rm = TRUE),
      five_hund_length_reads_count = sum(length > 250 & length <= 500, na.rm = TRUE),
      short_length_reads_count = sum(length > 500 & length <= 1000, na.rm = TRUE),
      mid_length_reads_count = sum(length > 1000 & length <= 2000, na.rm = TRUE),
      long_length_reads_count = sum(length > 2000, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    left_join(SQANTI_cell_summary %>% select(CB, total_reads), by = "CB") %>%
    mutate(
      two_fifty_length_reads = safe_prop(two_fifty_length_reads_count, total_reads),
      five_hund_length_reads = safe_prop(five_hund_length_reads_count, total_reads),
      short_length_reads = safe_prop(short_length_reads_count, total_reads),
      mid_length_reads = safe_prop(mid_length_reads_count, total_reads),
      long_length_reads = safe_prop(long_length_reads_count, total_reads)
    ) %>% select(CB, starts_with(c("two_", "five_", "short_", "mid_", "long_")) & ends_with("length_reads"))

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, read_lengths_general, by = "CB")

  # Monoexons read length general
  mono_read_lengths_general <- Classification %>%
    filter(exons == 1) %>%
    group_by(CB) %>%
    summarise(
      mono_two_fifty_length_reads_count = sum(length <= 250, na.rm = TRUE),
      mono_five_hund_length_reads_count = sum(length > 250 & length <= 500, na.rm = TRUE),
      mono_short_length_reads_count = sum(length > 500 & length <= 1000, na.rm = TRUE),
      mono_mid_length_reads_count = sum(length > 1000 & length <= 2000, na.rm = TRUE),
      mono_long_length_reads_count = sum(length > 2000, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    left_join(SQANTI_cell_summary %>% select(CB, total_reads), by = "CB") %>%
    mutate(
      mono_two_fifty_length_reads = safe_prop(mono_two_fifty_length_reads_count, total_reads),
      mono_five_hund_length_reads = safe_prop(mono_five_hund_length_reads_count, total_reads),
      mono_short_length_reads = safe_prop(mono_short_length_reads_count, total_reads),
      mono_mid_length_reads = safe_prop(mono_mid_length_reads_count, total_reads),
      mono_long_length_reads = safe_prop(mono_long_length_reads_count, total_reads)
    ) %>% select(CB, starts_with("mono_") & ends_with("length_reads"))

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, mono_read_lengths_general, by = "CB")
  mono_len_cols <- grep("^mono_.*length_reads$", names(SQANTI_cell_summary), value = TRUE) # select only the general mono
  # Fill NAs for these general mono length columns specifically
    for(col_name in mono_len_cols) {
        # Check if the column exists before trying to modify it
        if(col_name %in% names(SQANTI_cell_summary)) {
            SQANTI_cell_summary[[col_name]] <- ifelse(is.na(SQANTI_cell_summary[[col_name]]), 0, SQANTI_cell_summary[[col_name]])
        }
    }


  # Read length per category and monoexon per length break per category
  for (cat_name_full in structural_categories) {
    cat_name_short <- toupper(gsub("[-_]", "", cat_name_full))
    cat_count_col_name <- paste0(cat_name_short, "_count") # This is the actual column name in SQANTI_cell_summary

    # Overall length breaks for the category
    cat_lengths <- Classification %>%
      filter(structural_category == cat_name_full) %>%
      group_by(CB) %>%
      summarise(
        cat_two_fifty_count = sum(length <= 250, na.rm = TRUE),
        cat_five_hund_count = sum(length > 250 & length <= 500, na.rm = TRUE),
        cat_short_count = sum(length > 500 & length <= 1000, na.rm = TRUE),
        cat_mid_count = sum(length > 1000 & length <= 2000, na.rm = TRUE),
        cat_long_count = sum(length > 2000, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      left_join(SQANTI_cell_summary %>% select(CB, !!sym(cat_count_col_name)), by = "CB") %>%
      mutate(
        !!paste0("two_fifty_length_reads_", cat_name_short) := safe_prop(cat_two_fifty_count, !!sym(cat_count_col_name)),
        !!paste0("five_hund_length_reads_", cat_name_short) := safe_prop(cat_five_hund_count, !!sym(cat_count_col_name)),
        !!paste0("short_length_reads_", cat_name_short) := safe_prop(cat_short_count, !!sym(cat_count_col_name)),
        !!paste0("mid_length_reads_", cat_name_short) := safe_prop(cat_mid_count, !!sym(cat_count_col_name)),
        !!paste0("long_length_reads_", cat_name_short) := safe_prop(cat_long_count, !!sym(cat_count_col_name))
      ) %>% select(CB, starts_with(c("two_", "five_", "short_", "mid_", "long_")) & ends_with(cat_name_short))

    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, cat_lengths, by = "CB")
    new_cat_len_cols <- names(cat_lengths %>% select(-CB))
    for(col in new_cat_len_cols) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])


    # Monoexon length breaks for the category
    mono_cat_lengths <- Classification %>%
      filter(structural_category == cat_name_full & exons == 1) %>%
      group_by(CB) %>%
      summarise(
        mono_cat_two_fifty_count = sum(length <= 250, na.rm = TRUE),
        mono_cat_five_hund_count = sum(length > 250 & length <= 500, na.rm = TRUE),
        mono_cat_short_count = sum(length > 500 & length <= 1000, na.rm = TRUE),
        mono_cat_mid_count = sum(length > 1000 & length <= 2000, na.rm = TRUE),
        mono_cat_long_count = sum(length > 2000, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      left_join(SQANTI_cell_summary %>% select(CB, !!sym(cat_count_col_name)), by = "CB") %>%
      mutate(
        !!paste0("mono_two_fifty_length_reads_", cat_name_short) := safe_prop(mono_cat_two_fifty_count, !!sym(cat_count_col_name)),
        !!paste0("mono_five_hund_length_reads_", cat_name_short) := safe_prop(mono_cat_five_hund_count, !!sym(cat_count_col_name)),
        !!paste0("mono_short_length_reads_", cat_name_short) := safe_prop(mono_cat_short_count, !!sym(cat_count_col_name)),
        !!paste0("mono_mid_length_reads_", cat_name_short) := safe_prop(mono_cat_mid_count, !!sym(cat_count_col_name)),
        !!paste0("mono_long_length_reads_", cat_name_short) := safe_prop(mono_cat_long_count, !!sym(cat_count_col_name))
      ) %>% select(CB, starts_with("mono_") & ends_with(cat_name_short))

    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, mono_cat_lengths, by = "CB")
    new_mono_cat_len_cols <- names(mono_cat_lengths %>% select(-CB))
    for(col in new_mono_cat_len_cols) SQANTI_cell_summary[[col]] <- ifelse(is.na(SQANTI_cell_summary[[col]]), 0, SQANTI_cell_summary[[col]])

  }


  # Reference body coverage per category
  # ref_body_cover_in_cell was not used in the final row_data, so skipping direct equivalent for now
  for (cat_name_full in structural_categories) {
    cat_name_short <- toupper(gsub("[-_]", "", cat_name_full))
    cat_count_col_name <- paste0(cat_name_short, "_count")
    prop_col_name <- paste0("ref_body_cover_", cat_name_short)

    ref_body_cover_cat <- Classification %>%
      filter(structural_category == cat_name_full & (length / ref_length * 100) >= 45) %>%
      group_by(CB) %>%
      summarise(covered_reads_cat = n(), .groups = 'drop') %>%
      left_join(SQANTI_cell_summary %>% select(CB, !!sym(cat_count_col_name)), by = "CB") %>%
      mutate(!!prop_col_name := safe_prop(covered_reads_cat, !!sym(cat_count_col_name))) %>%
      select(CB, !!prop_col_name)

    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, ref_body_cover_cat, by = "CB")
    SQANTI_cell_summary[[prop_col_name]] <- ifelse(is.na(SQANTI_cell_summary[[prop_col_name]]), 0, SQANTI_cell_summary[[prop_col_name]])
  }


  # RTS metrics
  RTS_overall <- Classification %>%
    filter(RTS_stage == TRUE) %>%
    group_by(CB) %>%
    summarise(RTS_count = n(), .groups = 'drop') %>%
    left_join(SQANTI_cell_summary %>% select(CB, total_reads), by = "CB") %>%
    mutate(RTS_in_cell_prop = safe_prop(RTS_count, total_reads)) %>%
    select(CB, RTS_in_cell_prop)
  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, RTS_overall, by = "CB")
  SQANTI_cell_summary$RTS_in_cell_prop[is.na(SQANTI_cell_summary$RTS_in_cell_prop)] <- 0

  for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) { # Renamed to avoid conflict with cat_name_short
    cat_name_full_iter <- case_when( # Renamed to avoid conflict
      cat_name_short_iter == "FSM" ~ "full-splice_match",
      cat_name_short_iter == "ISM" ~ "incomplete-splice_match",
      cat_name_short_iter == "NIC" ~ "novel_in_catalog",
      cat_name_short_iter == "NNC" ~ "novel_not_in_catalog",
      TRUE ~ NA_character_
    )
    # cat_count_col <- paste0(cat_name_short_iter, "_count") # Incorrect: would be FSM_count
    actual_cat_count_col_name <- paste0(toupper(gsub("[-_]", "", cat_name_full_iter)), "_count") # Correct: e.g. FULLSPLICEMATCH_count
    prop_col <- paste0(cat_name_short_iter, "_RTS_prop")

    RTS_cat <- Classification %>%
      filter(structural_category == cat_name_full_iter & RTS_stage == TRUE) %>%
      group_by(CB) %>%
      summarise(RTS_cat_count = n(), .groups = 'drop') %>%
      left_join(SQANTI_cell_summary %>% select(CB, !!sym(actual_cat_count_col_name)), by = "CB") %>%
      mutate(!!prop_col := safe_prop(RTS_cat_count, !!sym(actual_cat_count_col_name))) %>%
      select(CB, !!prop_col)

    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, RTS_cat, by = "CB")
    SQANTI_cell_summary[[prop_col]][is.na(SQANTI_cell_summary[[prop_col]])] <- 0
  }


  # Non-canonical metrics
  non_canon_overall <- Classification %>%
    filter(all_canonical == "non_canonical") %>%
    group_by(CB) %>%
    summarise(non_canon_count = n(), .groups = 'drop') %>%
    left_join(SQANTI_cell_summary %>% select(CB, total_reads_no_monoexon), by = "CB") %>%
    mutate(non_canonical_in_cell_prop = safe_prop(non_canon_count, total_reads_no_monoexon)) %>%
    select(CB, non_canonical_in_cell_prop)
  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, non_canon_overall, by = "CB")
  SQANTI_cell_summary$non_canonical_in_cell_prop[is.na(SQANTI_cell_summary$non_canonical_in_cell_prop)] <- 0


  for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
    cat_name_full_iter <- case_when(
      cat_name_short_iter == "FSM" ~ "full-splice_match",
      cat_name_short_iter == "ISM" ~ "incomplete-splice_match",
      cat_name_short_iter == "NIC" ~ "novel_in_catalog",
      cat_name_short_iter == "NNC" ~ "novel_not_in_catalog",
      TRUE ~ NA_character_
    )
    prop_col <- paste0(cat_name_short_iter, "_noncanon_prop")

    # Denominator needs to be count of reads in category with exons > 1
    denom_noncanon_cat <- Classification %>%
        filter(structural_category == cat_name_full_iter & exons > 1) %>%
        group_by(CB) %>%
        summarise(denom_count = n(), .groups = 'drop')


    noncanon_cat <- Classification %>%
      filter(structural_category == cat_name_full_iter & exons > 1 & all_canonical == "non_canonical") %>%
      group_by(CB) %>%
      summarise(noncanon_cat_count = n(), .groups = 'drop') %>%
      left_join(denom_noncanon_cat, by = "CB") %>% # Make sure denom_count is never NA before division
      mutate(denom_count = ifelse(is.na(denom_count), 0, denom_count)) %>%
      mutate(!!prop_col := safe_prop(noncanon_cat_count, denom_count)) %>% 
      select(CB, !!prop_col)

    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, noncanon_cat, by = "CB")
    SQANTI_cell_summary[[prop_col]][is.na(SQANTI_cell_summary[[prop_col]])] <- 0
  }


  # Intrapriming metrics
  intraprim_overall <- Classification %>%
    filter(perc_A_downstream_TTS >= 60) %>%
    group_by(CB) %>%
    summarise(intraprim_count = n(), .groups = 'drop') %>%
    left_join(SQANTI_cell_summary %>% select(CB, total_reads), by = "CB") %>%
    mutate(intrapriming_in_cell_prop = safe_prop(intraprim_count, total_reads)) %>%
    select(CB, intrapriming_in_cell_prop)
  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, intraprim_overall, by = "CB")
  SQANTI_cell_summary$intrapriming_in_cell_prop[is.na(SQANTI_cell_summary$intrapriming_in_cell_prop)] <- 0

  for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
     cat_name_full_iter <- case_when(
      cat_name_short_iter == "FSM" ~ "full-splice_match",
      cat_name_short_iter == "ISM" ~ "incomplete-splice_match",
      cat_name_short_iter == "NIC" ~ "novel_in_catalog",
      cat_name_short_iter == "NNC" ~ "novel_not_in_catalog",
      TRUE ~ NA_character_
    )
    # cat_count_col <- paste0(cat_name_short_iter, "_count") # Incorrect
    actual_cat_count_col_name <- paste0(toupper(gsub("[-_]", "", cat_name_full_iter)), "_count") # Correct
    prop_col <- paste0(cat_name_short_iter, "_intrapriming_prop")

    intraprim_cat <- Classification %>%
      filter(structural_category == cat_name_full_iter & perc_A_downstream_TTS >= 60) %>%
      group_by(CB) %>%
      summarise(intraprim_cat_count = n(), .groups = 'drop') %>%
      left_join(SQANTI_cell_summary %>% select(CB, !!sym(actual_cat_count_col_name)), by = "CB") %>%
      mutate(!!prop_col := safe_prop(intraprim_cat_count, !!sym(actual_cat_count_col_name))) %>%
      select(CB, !!prop_col)

    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, intraprim_cat, by = "CB")
    SQANTI_cell_summary[[prop_col]][is.na(SQANTI_cell_summary[[prop_col]])] <- 0
  }

  # TSS Annotation Support
  if ("diff_to_gene_TSS" %in% colnames(Classification)) {
    tss_support_overall <- Classification %>%
        mutate(is_supported_tss = abs(diff_to_gene_TSS) <= 50) %>%
        group_by(CB) %>%
        summarise(supported_tss_count = sum(is_supported_tss, na.rm=TRUE), .groups = 'drop') %>%
        left_join(SQANTI_cell_summary %>% select(CB, total_reads), by="CB") %>%
        mutate(tss_annotation_support_in_cell_prop = safe_prop(supported_tss_count, total_reads)) %>%
        select(CB, tss_annotation_support_in_cell_prop)
    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, tss_support_overall, by="CB")
    SQANTI_cell_summary$tss_annotation_support_in_cell_prop[is.na(SQANTI_cell_summary$tss_annotation_support_in_cell_prop)] <- 0

    for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
        cat_name_full_iter <- case_when(
            cat_name_short_iter == "FSM" ~ "full-splice_match",
            cat_name_short_iter == "ISM" ~ "incomplete-splice_match",
            cat_name_short_iter == "NIC" ~ "novel_in_catalog",
            cat_name_short_iter == "NNC" ~ "novel_not_in_catalog",
            TRUE ~ NA_character_
        )
        prop_col <- paste0(cat_name_short_iter, "_tss_annotation_support") # Matches original implicit naming for row_data

        # Denom: total reads in category
        denom_cat_reads <- Classification %>%
            filter(structural_category == cat_name_full_iter) %>% 
            group_by(CB) %>%
            summarise(denom_cat_read_count = n(), .groups = 'drop')


        tss_cat <- Classification %>%
            filter(structural_category == cat_name_full_iter & abs(diff_to_gene_TSS) <= 50) %>%
            group_by(CB) %>%
            summarise(supported_tss_cat_count = n(), .groups = 'drop') %>%
            left_join(denom_cat_reads, by="CB") %>%
             mutate(denom_cat_read_count = ifelse(is.na(denom_cat_read_count), 0, denom_cat_read_count)) %>%
            mutate(!!prop_col := safe_prop(supported_tss_cat_count, denom_cat_read_count)) %>%
            select(CB, !!prop_col)

        SQANTI_cell_summary <- left_join(SQANTI_cell_summary, tss_cat, by = "CB")
        SQANTI_cell_summary[[prop_col]][is.na(SQANTI_cell_summary[[prop_col]])] <- 0
    }
  } else { # diff_to_gene_TSS not in input
      SQANTI_cell_summary$tss_annotation_support_in_cell_prop <- 0
      for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
          SQANTI_cell_summary[[paste0(cat_name_short_iter, "_tss_annotation_support")]] <- 0
      }
  }


  # Annotated genes in cell prop
  SQANTI_cell_summary <- SQANTI_cell_summary %>%
    mutate(annotated_genes_in_cell_prop = safe_prop(annotated_genes, genes_in_cell))

  # Annotated genes by category
  for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
    cat_name_full_iter <- case_when(
      cat_name_short_iter == "FSM" ~ "full-splice_match",
      cat_name_short_iter == "ISM" ~ "incomplete-splice_match",
      cat_name_short_iter == "NIC" ~ "novel_in_catalog",
      cat_name_short_iter == "NNC" ~ "novel_not_in_catalog",
      TRUE ~ NA_character_
    )
    prop_col <- paste0(cat_name_short_iter, "_anno_genes_prop")

    # Denominator is total distinct genes in that category for that cell
    total_genes_cat <- Classification %>%
        filter(structural_category == cat_name_full_iter) %>%
        group_by(CB) %>%
        summarise(total_genes_in_cat = n_distinct(associated_gene), .groups = 'drop')

    anno_genes_cat <- Classification %>%
      filter(structural_category == cat_name_full_iter & !grepl("^novel", associated_gene)) %>%
      group_by(CB) %>%
      summarise(anno_genes_in_cat = n_distinct(associated_gene), .groups = 'drop') %>%
      left_join(total_genes_cat, by = "CB") %>%
      mutate(total_genes_in_cat = ifelse(is.na(total_genes_in_cat), 0, total_genes_in_cat)) %>%
      mutate(!!prop_col := safe_prop(anno_genes_in_cat, total_genes_in_cat)) %>%
      select(CB, !!prop_col)

    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, anno_genes_cat, by = "CB")
    SQANTI_cell_summary[[prop_col]][is.na(SQANTI_cell_summary[[prop_col]])] <- 0
  }

  # Annotated junction strings prop
  # Denominator is models_in_cell (unique jxn_string for multiexonic across all genes in cell)
  # Numerator is sum of distinct jxn_string for multiexonic, non-novel transcripts
  anno_models_in_cell <- Classification %>%
    filter(!grepl("^novel", associated_transcript) & exons > 1) %>% # Consider only multiexonic non-novel transcripts
    group_by(CB) %>%
    summarise(anno_jxn_count = n_distinct(jxn_string), .groups = 'drop') %>% # Total distinct annotated jxns
    left_join(SQANTI_cell_summary %>% select(CB, models_in_cell), by = "CB") %>%
    mutate(anno_models_in_cell_prop = safe_prop(anno_jxn_count, models_in_cell)) %>%
    select(CB, anno_models_in_cell_prop)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, anno_models_in_cell, by = "CB")
  SQANTI_cell_summary$anno_models_in_cell_prop[is.na(SQANTI_cell_summary$anno_models_in_cell_prop)] <- 0

  # Percentage of canonical (overall and by category)
  # Denominator for overall: total_reads_no_monoexon
  # Numerator for overall: all canonical reads (mono or multi)
  canon_overall_prop_calc <- Classification %>%
    filter(all_canonical == "canonical") %>% # Numerator: all canonical reads
    group_by(CB) %>%
    summarise(num_canon_overall = n(), .groups = 'drop') %>%
    left_join(SQANTI_cell_summary %>% select(CB, total_reads_no_monoexon), by = "CB") %>%
    mutate(canonical_in_cell_prop = safe_prop(num_canon_overall, total_reads_no_monoexon)) %>%
    select(CB, canonical_in_cell_prop)

  SQANTI_cell_summary <- left_join(SQANTI_cell_summary, canon_overall_prop_calc, by = "CB")
  SQANTI_cell_summary$canonical_in_cell_prop[is.na(SQANTI_cell_summary$canonical_in_cell_prop)] <- 0


  for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
    cat_name_full_iter <- case_when(
      cat_name_short_iter == "FSM" ~ "full-splice_match",
      cat_name_short_iter == "ISM" ~ "incomplete-splice_match",
      cat_name_short_iter == "NIC" ~ "novel_in_catalog",
      cat_name_short_iter == "NNC" ~ "novel_not_in_catalog",
      TRUE ~ NA_character_
    )
    prop_col <- paste0(cat_name_short_iter, "_canon_prop")

    # Denom: count of reads in category with exons > 1
    denom_canon_cat <- Classification %>%
        filter(structural_category == cat_name_full_iter & exons > 1) %>%
        group_by(CB) %>%
        summarise(denom_count_canon = n(), .groups = 'drop')

    canon_cat <- Classification %>%
      filter(structural_category == cat_name_full_iter & exons > 1 & all_canonical == "canonical") %>%
      group_by(CB) %>%
      summarise(canon_cat_count = n(), .groups = 'drop') %>%
      left_join(denom_canon_cat, by = "CB") %>%
      mutate(denom_count_canon = ifelse(is.na(denom_count_canon),0,denom_count_canon)) %>%
      mutate(!!prop_col := safe_prop(canon_cat_count, denom_count_canon)) %>%
      select(CB, !!prop_col)

    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, canon_cat, by = "CB")
    SQANTI_cell_summary[[prop_col]][is.na(SQANTI_cell_summary[[prop_col]])] <- 0
  }


  # Conditional metrics: skipORF, CAGE_peak, polyA_motif_list
  if (!skipORF) {
    # NMD metrics
    nmd_overall <- Classification %>%
      filter(predicted_NMD == TRUE) %>%
      group_by(CB) %>%
      summarise(nmd_count = n(), .groups = 'drop') %>%
      left_join(SQANTI_cell_summary %>% select(CB, total_reads), by = "CB") %>%
      mutate(NMD_in_cell_prop = safe_prop(nmd_count, total_reads)) %>%
      select(CB, NMD_in_cell_prop)
    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, nmd_overall, by = "CB")
    SQANTI_cell_summary$NMD_in_cell_prop[is.na(SQANTI_cell_summary$NMD_in_cell_prop)] <- 0

    for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
      cat_name_full_iter <- case_when(
        cat_name_short_iter == "FSM" ~ "full-splice_match",
        cat_name_short_iter == "ISM" ~ "incomplete-splice_match",
        cat_name_short_iter == "NIC" ~ "novel_in_catalog",
        cat_name_short_iter == "NNC" ~ "novel_not_in_catalog",
        TRUE ~ NA_character_
      )
      # cat_count_col <- paste0(cat_name_short_iter, "_count") # Incorrect
      actual_cat_count_col_name <- paste0(toupper(gsub("[-_]", "", cat_name_full_iter)), "_count") # Correct
      prop_col <- paste0(cat_name_short_iter, "_NMD_prop")

      nmd_cat <- Classification %>%
        filter(structural_category == cat_name_full_iter & predicted_NMD == TRUE) %>%
        group_by(CB) %>%
        summarise(nmd_cat_count = n(), .groups = 'drop') %>%
        left_join(SQANTI_cell_summary %>% select(CB, !!sym(actual_cat_count_col_name)), by = "CB") %>%
        mutate(!!prop_col := safe_prop(nmd_cat_count, !!sym(actual_cat_count_col_name))) %>%
        select(CB, !!prop_col)

      SQANTI_cell_summary <- left_join(SQANTI_cell_summary, nmd_cat, by = "CB")
      SQANTI_cell_summary[[prop_col]][is.na(SQANTI_cell_summary[[prop_col]])] <- 0
    }

    # Coding/non-coding by category
    for (cat_name_full_iter in structural_categories) {
      cat_name_short_upper <- toupper(gsub("[-_]", "", cat_name_full_iter))
      cat_count_col <- paste0(cat_name_short_upper, "_count")
      cod_prop_col <- paste0("cod_", cat_name_short_upper)
      ncod_prop_col <- paste0("ncod_", cat_name_short_upper)

      coding_status_cat <- Classification %>%
        filter(structural_category == cat_name_full_iter) %>%
        group_by(CB) %>%
        summarise(
          coding_cat_count = sum(coding == "coding", na.rm = TRUE),
          non_coding_cat_count = sum(coding == "non_coding", na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        left_join(SQANTI_cell_summary %>% select(CB, !!sym(cat_count_col)), by = "CB") %>%
        mutate(
          !!cod_prop_col := safe_prop(coding_cat_count, !!sym(cat_count_col)),
          !!ncod_prop_col := safe_prop(non_coding_cat_count, !!sym(cat_count_col), default_val = 0) # Changed default_val to 0
        ) %>%
        select(CB, !!cod_prop_col, !!ncod_prop_col)

      SQANTI_cell_summary <- left_join(SQANTI_cell_summary, coding_status_cat, by = "CB")
      SQANTI_cell_summary[[cod_prop_col]][is.na(SQANTI_cell_summary[[cod_prop_col]])] <- 0
      # Special handling for ncod default if category count is 0 (original script set to 100)
      # Check if the count column exists and is not NA before using it
      # Simplified NA handling for ncod_prop_col
      SQANTI_cell_summary[[ncod_prop_col]][is.na(SQANTI_cell_summary[[ncod_prop_col]])] <- 0
    }
  } else { # skipORF is TRUE
      # Set ORF related columns to 0 or default
      SQANTI_cell_summary$NMD_in_cell_prop <- 0
      for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
          SQANTI_cell_summary[[paste0(cat_name_short_iter, "_NMD_prop")]] <- 0
      }
      for (cat_name_full_iter in structural_categories) {
          cat_name_short_upper <- toupper(gsub("[-_]", "", cat_name_full_iter))
          SQANTI_cell_summary[[paste0("cod_", cat_name_short_upper)]] <- 0
          SQANTI_cell_summary[[paste0("ncod_", cat_name_short_upper)]] <- 100 # Default from original
      }
  }


  if (CAGE_peak) {
    cage_support_overall <- Classification %>%
        filter(within_CAGE_peak == TRUE) %>%
        group_by(CB) %>%
        summarise(cage_supported_count = n(), .groups = 'drop') %>%
        left_join(SQANTI_cell_summary %>% select(CB, total_reads), by = "CB") %>%
        mutate(CAGE_peak_support_in_cell_prop = safe_prop(cage_supported_count, total_reads)) %>%
        select(CB, CAGE_peak_support_in_cell_prop)
    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, cage_support_overall, by = "CB")
    SQANTI_cell_summary$CAGE_peak_support_in_cell_prop[is.na(SQANTI_cell_summary$CAGE_peak_support_in_cell_prop)] <- 0

    for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
        cat_name_full_iter <- case_when(
            cat_name_short_iter == "FSM" ~ "full-splice_match",
            cat_name_short_iter == "ISM" ~ "incomplete-splice_match",
            cat_name_short_iter == "NIC" ~ "novel_in_catalog",
            cat_name_short_iter == "NNC" ~ "novel_not_in_catalog",
            TRUE ~ NA_character_
        )
        prop_col <- paste0(cat_name_short_iter, "_cage_peak_support")

        denom_cage_cat <- Classification %>%
            filter(structural_category == cat_name_full_iter) %>%
            group_by(CB) %>%
            summarise(denom_count = n(), .groups = 'drop')

        cage_cat <- Classification %>%
            filter(structural_category == cat_name_full_iter & within_CAGE_peak == TRUE) %>%
            group_by(CB) %>%
            summarise(cage_cat_count = n(), .groups = 'drop') %>%
            left_join(denom_cage_cat, by = "CB") %>%
            mutate(denom_count = ifelse(is.na(denom_count), 0, denom_count)) %>%
            mutate(!!prop_col := safe_prop(cage_cat_count, denom_count)) %>%
            select(CB, !!prop_col)

        SQANTI_cell_summary <- left_join(SQANTI_cell_summary, cage_cat, by = "CB")
        SQANTI_cell_summary[[prop_col]][is.na(SQANTI_cell_summary[[prop_col]])] <- 0
    }
  } else { # CAGE_peak is FALSE
      SQANTI_cell_summary$CAGE_peak_support_in_cell_prop <- 0
      for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
          SQANTI_cell_summary[[paste0(cat_name_short_iter, "_cage_peak_support")]] <- 0
      }
  }

  if (polyA_motif_list) {
    polya_support_overall <- Classification %>%
        filter(polyA_motif_found == TRUE) %>% 
        group_by(CB) %>%
        summarise(polya_supported_count = n(), .groups = 'drop') %>%
        left_join(SQANTI_cell_summary %>% select(CB, total_reads), by = "CB") %>%
        mutate(polyA_motif_support_in_cell_prop = safe_prop(polya_supported_count, total_reads)) %>%
        select(CB, polyA_motif_support_in_cell_prop)
    SQANTI_cell_summary <- left_join(SQANTI_cell_summary, polya_support_overall, by = "CB")
    SQANTI_cell_summary$polyA_motif_support_in_cell_prop[is.na(SQANTI_cell_summary$polyA_motif_support_in_cell_prop)] <- 0

    for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
        cat_name_full_iter <- case_when(
            cat_name_short_iter == "FSM" ~ "full-splice_match",
            cat_name_short_iter == "ISM" ~ "incomplete-splice_match",
            cat_name_short_iter == "NIC" ~ "novel_in_catalog",
            cat_name_short_iter == "NNC" ~ "novel_not_in_catalog",
            TRUE ~ NA_character_
        )
        prop_col <- paste0(cat_name_short_iter, "_polyA_motif_support")

        denom_polya_cat <- Classification %>%
            filter(structural_category == cat_name_full_iter) %>%
            group_by(CB) %>%
            summarise(denom_count = n(), .groups = 'drop')

        polya_cat <- Classification %>%
            filter(structural_category == cat_name_full_iter & polyA_motif_found == TRUE) %>% 
            group_by(CB) %>%
            summarise(polya_cat_count = n(), .groups = 'drop') %>%
            left_join(denom_polya_cat, by = "CB") %>%
             mutate(denom_count = ifelse(is.na(denom_count), 0, denom_count)) %>%
            mutate(!!prop_col := safe_prop(polya_cat_count, denom_count)) %>%
            select(CB, !!prop_col)

        SQANTI_cell_summary <- left_join(SQANTI_cell_summary, polya_cat, by = "CB")
        SQANTI_cell_summary[[prop_col]][is.na(SQANTI_cell_summary[[prop_col]])] <- 0
    }
  } else { # polyA_motif_list is FALSE
      SQANTI_cell_summary$polyA_motif_support_in_cell_prop <- 0
      for (cat_name_short_iter in c("FSM", "ISM", "NIC", "NNC")) {
          SQANTI_cell_summary[[paste0(cat_name_short_iter, "_polyA_motif_support")]] <- 0
      }
  }

  # Define the final column names in the desired order from the original script
  # This is crucial for matching the exact output format.
  # Base columns
    col_names_final <- c(
    "CB", "Reads_in_cell", "UMIs_in_cell", "Genes_in_cell", "UJCs_in_cell",
    "Annotated_genes", "Novel_genes", # Std cell counts
    "MT_perc",
    "Total_junctions",
    "Known_canonical_junctions", "Known_non_canonical_junctions", "Novel_canonical_junctions", "Novel_non_canonical_junctions",
    "Known_canonical_junctions_prop", "Known_non_canonical_junctions_prop", "Novel_canonical_junctions_prop", "Novel_non_canonical_junctions_prop", # Junction-based metrics
    "FSM", "ISM", "NIC", "NNC",
    "Genic_Genomic", "Antisense", "Fusion", "Intergenic", "Genic_intron", # Structural categories counts (reads)
    "FSM_prop", "ISM_prop", "NIC_prop", "NNC_prop",
    "Genic_Genomic_prop", "Antisense_prop", "Fusion_prop", "Intergenic_prop", "Genic_intron_prop", # Structural categories props (reads)
    "FSM_alternative_3.end_prop", "FSM_alternative_3.5.end_prop", "FSM_alternative_5.end_prop", "FSM_reference_match_prop", "FSM_mono.exon_prop",
    "ISM_3._fragment_prop", "ISM_internal_fragment_prop", "ISM_5._fragment_prop", "ISM_intron_retention_prop", "ISM_mono.exon_prop",
    "NIC_comb_annot_junctions_prop", "NIC_comb_annot_splice_sites_prop", "NIC_intron_retention_prop", "NIC_mono.exon_by_intron_retention_prop", "NIC_mono.exon_prop",
    "NNC_at_least_1_don_accept_prop", "NNC_intron_retention_prop",
    "Genic_mono.exon_prop", "Genic_multi.exon_prop",
    "Antisense_mono.exon_prop", "Antisense_multi.exon_prop",
    "Fusion_intron_retention_prop", "Fusion_multi.exon_prop",
    "Intergenic_mono.exon_prop", "Intergenic_multi.exon_prop",
    "Genic_intron_mono.exon_prop", "Genic_intron_multi.exon_prop", # Structural subcategories props (reads)
    "Total_250b_length_prop", "Total_250b_length_mono_prop",
    "Total_500b_length_prop", "Total_500b_length_mono_prop",
    "Total_short_length_prop", "Total_short_length_mono_prop",
    "Total_mid_length_prop", "Total_mid_length_mono_prop",
    "Total_long_length_prop", "Total_long_length_mono_prop", # Read lengths breaks general + monoexons (reads)
    "FSM_250b_length_prop", "FSM_250b_length_mono_prop",
    "FSM_500b_length_prop", "FSM_500b_length_mono_prop",
    "FSM_short_length_prop", "FSM_short_length_mono_prop",
    "FSM_mid_length_prop", "FSM_mid_length_mono_prop",
    "FSM_long_length_prop", "FSM_long_length_mono_prop",
    "ISM_250b_length_prop", "ISM_250b_length_mono_prop",
    "ISM_500b_length_prop", "ISM_500b_length_mono_prop",
    "ISM_short_length_prop", "ISM_short_length_mono_prop",
    "ISM_mid_length_prop", "ISM_mid_length_mono_prop",
    "ISM_long_length_prop", "ISM_long_length_mono_prop",
    "NIC_250b_length_prop", "NIC_250b_length_mono_prop",
    "NIC_500b_length_prop", "NIC_500b_length_mono_prop",
    "NIC_short_length_prop", "NIC_short_length_mono_prop",
    "NIC_mid_length_prop", "NIC_mid_length_mono_prop",
    "NIC_long_length_prop", "NIC_long_length_mono_prop",
    "NNC_250b_length_prop", "NNC_250b_length_mono_prop",
    "NNC_500b_length_prop", "NNC_500b_length_mono_prop",
    "NNC_short_length_prop", "NNC_short_length_mono_prop",
    "NNC_mid_length_prop", "NNC_mid_length_mono_prop",
    "NNC_long_length_prop", "NNC_long_length_mono_prop",
    "Genic_250b_length_prop", "Genic_250b_length_mono_prop",
    "Genic_500b_length_prop", "Genic_500b_length_mono_prop",
    "Genic_short_length_prop", "Genic_short_length_mono_prop",
    "Genic_mid_length_prop", "Genic_mid_length_mono_prop",
    "Genic_long_length_prop", "Genic_long_length_mono_prop",
    "Antisense_250b_length_prop", "Antisense_250b_length_mono_prop",
    "Antisense_500b_length_prop", "Antisense_500b_length_mono_prop",
    "Antisense_short_length_prop", "Antisense_short_length_mono_prop",
    "Antisense_mid_length_prop", "Antisense_mid_length_mono_prop",
    "Antisense_long_length_prop", "Antisense_long_length_mono_prop",
    "Fusion_250b_length_prop", "Fusion_250b_length_mono_prop",
    "Fusion_500b_length_prop", "Fusion_500b_length_mono_prop",
    "Fusion_short_length_prop", "Fusion_short_length_mono_prop",
    "Fusion_mid_length_prop", "Fusion_mid_length_mono_prop",
    "Fusion_long_length_prop", "Fusion_long_length_mono_prop",
    "Intergenic_250b_length_prop", "Intergenic_250b_length_mono_prop",
    "Intergenic_500b_length_prop", "Intergenic_500b_length_mono_prop",
    "Intergenic_short_length_prop", "Intergenic_short_length_mono_prop",
    "Intergenic_mid_length_prop", "Intergenic_mid_length_mono_prop",
    "Intergenic_long_length_prop", "Intergenic_long_length_mono_prop",
    "Genic_intron_250b_length_prop", "Genic_intron_250b_length_mono_prop",
    "Genic_intron_500b_length_prop", "Genic_intron_500b_length_mono_prop",
    "Genic_intron_short_length_prop", "Genic_intron_short_length_mono_prop",
    "Genic_intron_mid_length_prop", "Genic_intron_mid_length_mono_prop",
    "Genic_intron_long_length_prop", "Genic_intron_long_length_mono_prop",
    "FSM_ref_coverage_prop", "ISM_ref_coverage_prop", "NIC_ref_coverage_prop", "NNC_ref_coverage_prop",
    "Genic_ref_coverage_prop", "Antisense_ref_coverage_prop", "Fusion_ref_coverage_prop", "Intergenic_ref_coverage_prop", "Genic_intron_ref_coverage_prop", # Coverage of reference length (set at 45% default)
    "anno_bin1_perc", "anno_bin2_3_perc", "anno_bin4_5_perc", "anno_bin6plus_perc", # Gene count bins
    "novel_bin1_perc", "novel_bin2_3_perc", "novel_bin4_5_perc", "novel_bin6plus_perc",
    "anno_ujc_bin1_perc", "anno_ujc_bin2_3_perc", "anno_ujc_bin4_5_perc", "anno_ujc_bin6plus_perc", # Junction count bins
    "novel_ujc_bin1_perc", "novel_ujc_bin2_3_perc", "novel_ujc_bin4_5_perc", "novel_ujc_bin6plus_perc",
    "RTS_prop_in_cell", "FSM_RTS_prop", "ISM_RTS_prop", "NIC_RTS_prop", "NNC_RTS_prop", # Features of bad quality
    "Non_canonical_prop_in_cell", "FSM_noncanon_prop", "ISM_noncanon_prop", "NIC_noncanon_prop", "NNC_noncanon_prop",
    "Intrapriming_prop_in_cell", "FSM_intrapriming_prop", "ISM_intrapriming_prop", "NIC_intrapriming_prop", "NNC_intrapriming_prop", # Features of bad quality
    "TSSAnnotationSupport_prop", "FSM_TSSAnnotationSupport", "ISM_TSSAnnotationSupport", "NIC_TSSAnnotationSupport", "NNC_TSSAnnotationSupport",
    "Annotated_genes_prop_in_cell", "FSM_anno_genes_prop", "ISM_anno_genes_prop", "NIC_anno_genes_prop", "NNC_anno_genes_prop",
    "Annotated_juction_strings_prop_in_cell",
    "Canonical_prop_in_cell", "FSM_canon_prop", "ISM_canon_prop", "NIC_canon_prop", "NNC_canon_prop" # Features of good quality
    )

  if (!skipORF) {
  col_names_final <- c(col_names_final,
    "Coding_FSM_prop", "Non_coding_FSM_prop", "Coding_ISM_prop", "Non_coding_ISM_prop",
    "Coding_NIC_prop", "Non_coding_NIC_prop", "Coding_NNC_prop", "Non_coding_NNC_prop",
    "Coding_genic_prop", "Non_coding_genic_prop", "Coding_antisense_prop", "Non_coding_antisense_prop",
    "Coding_fusion_prop", "Non_coding_fusion_prop", "Coding_intergenic_prop", "Non_coding_intergenic_prop",
    "Coding_genic_intron_prop", "Non_coding_genic_intron_prop",
    "NMD_prop_in_cell", "FSM_NMD_prop", "ISM_NMD_prop", "NIC_NMD_prop", "NNC_NMD_prop"
    )
  }

  if (CAGE_peak) {
  col_names_final <- c(col_names_final,
  "CAGE_peak_support_prop", "FSM_CAGE_peak_support_prop", "ISM_CAGE_peak_support_prop", "NIC_CAGE_peak_support_prop", "NNC_CAGE_peak_support_prop"
    )
  }

  if (polyA_motif_list) {
    col_names_final <- c(col_names_final,
    "PolyA_motif_support_prop", "FSM_PolyA_motif_support_prop", "ISM_PolyA_motif_support_prop", "NIC_PolyA_motif_support_prop", "NNC_PolyA_motif_support_prop"
    )
  }

  # Create a data frame with all unique CBs from the input Classification file.
  # This ensures that all cells are present in the final output, even if some had no reads
  # or were filtered out during intermediate steps.
  all_CBs_df <- data.frame(CB = unique(Classification$CB))

  # Left join the calculated SQANTI_cell_summary with all_CBs_df.
  # This will add rows for CBs that might be missing from SQANTI_cell_summary, filling with NAs.
  SQANTI_cell_summary_merged <- left_join(all_CBs_df, SQANTI_cell_summary, by = "CB")


  # Map generated column names to the final desired column names
  # This is a critical step to ensure the output matches the original script's format.
  # current_to_final_map defines the mapping.
  # (Using the same extensive current_to_final_map as defined in the previous attempt)
  
  current_to_final_map <- c(
      "total_reads" = "Reads_in_cell", "total_UMI" = "UMIs_in_cell", "genes_in_cell" = "Genes_in_cell", "models_in_cell" = "UJCs_in_cell",
      "annotated_genes" = "Annotated_genes", "novel_genes" = "Novel_genes",
      "total_junctions" = "Total_junctions",
      "known_canonical_count" = "Known_canonical_junctions",
      "known_non_canonical_count" = "Known_non_canonical_junctions",
      "novel_canonical_count" = "Novel_canonical_junctions",
      "novel_non_canonical_count" = "Novel_non_canonical_junctions",
      "known_canonical_prop" = "Known_canonical_junctions_prop",
      "known_non_canonical_prop" = "Known_non_canonical_junctions_prop",
      "novel_canonical_prop" = "Novel_canonical_junctions_prop",
      "novel_non_canonical_prop" = "Novel_non_canonical_junctions_prop",
      # Structural category counts
      "FULLSPLICEMATCH_count" = "FSM", "INCOMPLETESPLICEMATCH_count" = "ISM", "NOVELINCATALOG_count" = "NIC", "NOVELNOTINCATALOG_count" = "NNC",
      "GENIC_count" = "Genic_Genomic", "ANTISENSE_count" = "Antisense", "FUSION_count" = "Fusion", "INTERGENIC_count" = "Intergenic", "GENICINTRON_count" = "Genic_intron",
      # Structural category proportions
      "FULLSPLICEMATCH_prop" = "FSM_prop", "INCOMPLETESPLICEMATCH_prop" = "ISM_prop", "NOVELINCATALOG_prop" = "NIC_prop", "NOVELNOTINCATALOG_prop" = "NNC_prop",
      "GENIC_prop" = "Genic_Genomic_prop", "ANTISENSE_prop" = "Antisense_prop", "FUSION_prop" = "Fusion_prop", "INTERGENIC_prop" = "Intergenic_prop", "GENICINTRON_prop" = "Genic_intron_prop",
      # Subcategory proportions (FSM)
      "sub_FSM_alternative_3end" = "FSM_alternative_3.end_prop", "sub_FSM_alternative_3end5end" = "FSM_alternative_3.5.end_prop", "sub_FSM_alternative_5end" = "FSM_alternative_5.end_prop",
      "sub_FSM_reference_match" = "FSM_reference_match_prop", "sub_FSM_mono-exon" = "FSM_mono.exon_prop",
      # Subcategory proportions (ISM)
      "sub_ISM_3prime_fragment" = "ISM_3._fragment_prop", "sub_ISM_internal_fragment" = "ISM_internal_fragment_prop", "sub_ISM_5prime_fragment" = "ISM_5._fragment_prop",
      "sub_ISM_intron_retention" = "ISM_intron_retention_prop", "sub_ISM_mono-exon" = "ISM_mono.exon_prop",
      # Subcategory proportions (NIC)
      "sub_NIC_combination_of_known_junctions" = "NIC_comb_annot_junctions_prop", "sub_NIC_combination_of_known_splicesites" = "NIC_comb_annot_splice_sites_prop",
      "sub_NIC_intron_retention" = "NIC_intron_retention_prop", "sub_NIC_mono-exon_by_intron_retention" = "NIC_mono.exon_by_intron_retention_prop", "sub_NIC_mono-exon" = "NIC_mono.exon_prop",
      # Subcategory proportions (NNC)
      "sub_NNC_at_least_one_novel_splicesite" = "NNC_at_least_1_don_accept_prop", "sub_NNC_intron_retention" = "NNC_intron_retention_prop",
      # Subcategory proportions (Genic, Antisense, Fusion, Intergenic, Genic_intron)
      "sub_genic_mono-exon" = "Genic_mono.exon_prop", "sub_genic_multi-exon" = "Genic_multi.exon_prop",
      "sub_antisense_mono-exon" = "Antisense_mono.exon_prop", "sub_antisense_multi-exon" = "Antisense_multi.exon_prop",
      "sub_fusion_intron_retention" = "Fusion_intron_retention_prop", "sub_fusion_multi-exon" = "Fusion_multi.exon_prop",
      "sub_intergenic_mono-exon" = "Intergenic_mono.exon_prop", "sub_intergenic_multi-exon" = "Intergenic_multi.exon_prop",
      "sub_genic_intron_mono-exon" = "Genic_intron_mono.exon_prop", "sub_genic_intron_multi-exon" = "Genic_intron_multi.exon_prop",
      # General Read Lengths
      "two_fifty_length_reads" = "Total_250b_length_prop", "mono_two_fifty_length_reads" = "Total_250b_length_mono_prop",
      "five_hund_length_reads" = "Total_500b_length_prop", "mono_five_hund_length_reads" = "Total_500b_length_mono_prop",
      "short_length_reads" = "Total_short_length_prop", "mono_short_length_reads" = "Total_short_length_mono_prop",
      "mid_length_reads" = "Total_mid_length_prop", "mono_mid_length_reads" = "Total_mid_length_mono_prop",
      "long_length_reads" = "Total_long_length_prop", "mono_long_length_reads" = "Total_long_length_mono_prop",
      # Category Read Lengths (FSM example, repeat for all categories)
      "two_fifty_length_reads_FULLSPLICEMATCH" = "FSM_250b_length_prop", "mono_two_fifty_length_reads_FULLSPLICEMATCH" = "FSM_250b_length_mono_prop",
        "five_hund_length_reads_FULLSPLICEMATCH" = "FSM_500b_length_prop", "mono_five_hund_length_reads_FULLSPLICEMATCH" = "FSM_500b_length_mono_prop",
        "short_length_reads_FULLSPLICEMATCH" = "FSM_short_length_prop", "mono_short_length_reads_FULLSPLICEMATCH" = "FSM_short_length_mono_prop",
        "mid_length_reads_FULLSPLICEMATCH" = "FSM_mid_length_prop", "mono_mid_length_reads_FULLSPLICEMATCH" = "FSM_mid_length_mono_prop",
        "long_length_reads_FULLSPLICEMATCH" = "FSM_long_length_prop", "mono_long_length_reads_FULLSPLICEMATCH" = "FSM_long_length_mono_prop",
      "two_fifty_length_reads_INCOMPLETESPLICEMATCH" = "ISM_250b_length_prop", "mono_two_fifty_length_reads_INCOMPLETESPLICEMATCH" = "ISM_250b_length_mono_prop",
        "five_hund_length_reads_INCOMPLETESPLICEMATCH" = "ISM_500b_length_prop", "mono_five_hund_length_reads_INCOMPLETESPLICEMATCH" = "ISM_500b_length_mono_prop",
        "short_length_reads_INCOMPLETESPLICEMATCH" = "ISM_short_length_prop", "mono_short_length_reads_INCOMPLETESPLICEMATCH" = "ISM_short_length_mono_prop",
        "mid_length_reads_INCOMPLETESPLICEMATCH" = "ISM_mid_length_prop", "mono_mid_length_reads_INCOMPLETESPLICEMATCH" = "ISM_mid_length_mono_prop",
        "long_length_reads_INCOMPLETESPLICEMATCH" = "ISM_long_length_prop", "mono_long_length_reads_INCOMPLETESPLICEMATCH" = "ISM_long_length_mono_prop",
      "two_fifty_length_reads_NOVELINCATALOG" = "NIC_250b_length_prop", "mono_two_fifty_length_reads_NOVELINCATALOG" = "NIC_250b_length_mono_prop",
        "five_hund_length_reads_NOVELINCATALOG" = "NIC_500b_length_prop", "mono_five_hund_length_reads_NOVELINCATALOG" = "NIC_500b_length_mono_prop",
        "short_length_reads_NOVELINCATALOG" = "NIC_short_length_prop", "mono_short_length_reads_NOVELINCATALOG" = "NIC_short_length_mono_prop",
        "mid_length_reads_NOVELINCATALOG" = "NIC_mid_length_prop", "mono_mid_length_reads_NOVELINCATALOG" = "NIC_mid_length_mono_prop",
        "long_length_reads_NOVELINCATALOG" = "NIC_long_length_prop", "mono_long_length_reads_NOVELINCATALOG" = "NIC_long_length_mono_prop",
      "two_fifty_length_reads_NOVELNOTINCATALOG" = "NNC_250b_length_prop", "mono_two_fifty_length_reads_NOVELNOTINCATALOG" = "NNC_250b_length_mono_prop",
        "five_hund_length_reads_NOVELNOTINCATALOG" = "NNC_500b_length_prop", "mono_five_hund_length_reads_NOVELNOTINCATALOG" = "NNC_500b_length_mono_prop",
        "short_length_reads_NOVELNOTINCATALOG" = "NNC_short_length_prop", "mono_short_length_reads_NOVELNOTINCATALOG" = "NNC_short_length_mono_prop",
        "mid_length_reads_NOVELNOTINCATALOG" = "NNC_mid_length_prop", "mono_mid_length_reads_NOVELNOTINCATALOG" = "NNC_mid_length_mono_prop",
        "long_length_reads_NOVELNOTINCATALOG" = "NNC_long_length_prop", "mono_long_length_reads_NOVELNOTINCATALOG" = "NNC_long_length_mono_prop",
      "two_fifty_length_reads_GENIC" = "Genic_250b_length_prop", "mono_two_fifty_length_reads_GENIC" = "Genic_250b_length_mono_prop",
        "five_hund_length_reads_GENIC" = "Genic_500b_length_prop", "mono_five_hund_length_reads_GENIC" = "Genic_500b_length_mono_prop",
        "short_length_reads_GENIC" = "Genic_short_length_prop", "mono_short_length_reads_GENIC" = "Genic_short_length_mono_prop",
        "mid_length_reads_GENIC" = "Genic_mid_length_prop", "mono_mid_length_reads_GENIC" = "Genic_mid_length_mono_prop",
        "long_length_reads_GENIC" = "Genic_long_length_prop", "mono_long_length_reads_GENIC" = "Genic_long_length_mono_prop",
      "two_fifty_length_reads_ANTISENSE" = "Antisense_250b_length_prop", "mono_two_fifty_length_reads_ANTISENSE" = "Antisense_250b_length_mono_prop",
        "five_hund_length_reads_ANTISENSE" = "Antisense_500b_length_prop", "mono_five_hund_length_reads_ANTISENSE" = "Antisense_500b_length_mono_prop",
        "short_length_reads_ANTISENSE" = "Antisense_short_length_prop", "mono_short_length_reads_ANTISENSE" = "Antisense_short_length_mono_prop",
        "mid_length_reads_ANTISENSE" = "Antisense_mid_length_prop", "mono_mid_length_reads_ANTISENSE" = "Antisense_mid_length_mono_prop",
        "long_length_reads_ANTISENSE" = "Antisense_long_length_prop", "mono_long_length_reads_ANTISENSE" = "Antisense_long_length_mono_prop",
      "two_fifty_length_reads_FUSION" = "Fusion_250b_length_prop", "mono_two_fifty_length_reads_FUSION" = "Fusion_250b_length_mono_prop",
        "five_hund_length_reads_FUSION" = "Fusion_500b_length_prop", "mono_five_hund_length_reads_FUSION" = "Fusion_500b_length_mono_prop",
        "short_length_reads_FUSION" = "Fusion_short_length_prop", "mono_short_length_reads_FUSION" = "Fusion_short_length_mono_prop",
        "mid_length_reads_FUSION" = "Fusion_mid_length_prop", "mono_mid_length_reads_FUSION" = "Fusion_mid_length_mono_prop",
        "long_length_reads_FUSION" = "Fusion_long_length_prop", "mono_long_length_reads_FUSION" = "Fusion_long_length_mono_prop",
      "two_fifty_length_reads_INTERGENIC" = "Intergenic_250b_length_prop", "mono_two_fifty_length_reads_INTERGENIC" = "Intergenic_250b_length_mono_prop",
        "five_hund_length_reads_INTERGENIC" = "Intergenic_500b_length_prop", "mono_five_hund_length_reads_INTERGENIC" = "Intergenic_500b_length_mono_prop",
        "short_length_reads_INTERGENIC" = "Intergenic_short_length_prop", "mono_short_length_reads_INTERGENIC" = "Intergenic_short_length_mono_prop",
        "mid_length_reads_INTERGENIC" = "Intergenic_mid_length_prop", "mono_mid_length_reads_INTERGENIC" = "Intergenic_mid_length_mono_prop",
        "long_length_reads_INTERGENIC" = "Intergenic_long_length_prop", "mono_long_length_reads_INTERGENIC" = "Intergenic_long_length_mono_prop",
      "two_fifty_length_reads_GENICINTRON" = "Genic_intron_250b_length_prop", "mono_two_fifty_length_reads_GENICINTRON" = "Genic_intron_250b_length_mono_prop",
        "five_hund_length_reads_GENICINTRON" = "Genic_intron_500b_length_prop", "mono_five_hund_length_reads_GENICINTRON" = "Genic_intron_500b_length_mono_prop",
        "short_length_reads_GENICINTRON" = "Genic_intron_short_length_prop", "mono_short_length_reads_GENICINTRON" = "Genic_intron_short_length_mono_prop",
        "mid_length_reads_GENICINTRON" = "Genic_intron_mid_length_prop", "mono_mid_length_reads_GENICINTRON" = "Genic_intron_mid_length_mono_prop",
        "long_length_reads_GENICINTRON" = "Genic_intron_long_length_prop", "mono_long_length_reads_GENICINTRON" = "Genic_intron_long_length_mono_prop",
      # Ref Coverage
      "ref_body_cover_FULLSPLICEMATCH" = "FSM_ref_coverage_prop", "ref_body_cover_INCOMPLETESPLICEMATCH" = "ISM_ref_coverage_prop", "ref_body_cover_NOVELINCATALOG" = "NIC_ref_coverage_prop", "ref_body_cover_NOVELNOTINCATALOG" = "NNC_ref_coverage_prop",
      "ref_body_cover_GENIC" = "Genic_ref_coverage_prop", "ref_body_cover_ANTISENSE" = "Antisense_ref_coverage_prop", "ref_body_cover_FUSION" = "Fusion_ref_coverage_prop", "ref_body_cover_INTERGENIC" = "Intergenic_ref_coverage_prop", "ref_body_cover_GENICINTRON" = "Genic_intron_ref_coverage_prop",
      # Bad Quality Features
      "RTS_in_cell_prop" = "RTS_prop_in_cell", "FSM_RTS_prop" = "FSM_RTS_prop", "ISM_RTS_prop" = "ISM_RTS_prop", "NIC_RTS_prop" = "NIC_RTS_prop", "NNC_RTS_prop" = "NNC_RTS_prop",
      "non_canonical_in_cell_prop" = "Non_canonical_prop_in_cell", "FSM_noncanon_prop" = "FSM_noncanon_prop", "ISM_noncanon_prop" = "ISM_noncanon_prop", "NIC_noncanon_prop" = "NIC_noncanon_prop", "NNC_noncanon_prop" = "NNC_noncanon_prop",
      "intrapriming_in_cell_prop" = "Intrapriming_prop_in_cell", "FSM_intrapriming_prop" = "FSM_intrapriming_prop", "ISM_intrapriming_prop" = "ISM_intrapriming_prop", "NIC_intrapriming_prop" = "NIC_intrapriming_prop", "NNC_intrapriming_prop" = "NNC_intrapriming_prop",
      # Good Quality Features
      "tss_annotation_support_in_cell_prop" = "TSSAnnotationSupport_prop", "FSM_tss_annotation_support" = "FSM_TSSAnnotationSupport", "ISM_tss_annotation_support" = "ISM_TSSAnnotationSupport", "NIC_tss_annotation_support" = "NIC_TSSAnnotationSupport", "NNC_tss_annotation_support" = "NNC_TSSAnnotationSupport",
      "annotated_genes_in_cell_prop" = "Annotated_genes_prop_in_cell", "FSM_anno_genes_prop" = "FSM_anno_genes_prop", "ISM_anno_genes_prop" = "ISM_anno_genes_prop", "NIC_anno_genes_prop" = "NIC_anno_genes_prop", "NNC_anno_genes_prop" = "NNC_anno_genes_prop",
      "anno_models_in_cell_prop" = "Annotated_juction_strings_prop_in_cell",
      "canonical_in_cell_prop" = "Canonical_prop_in_cell", "FSM_canon_prop" = "FSM_canon_prop", "ISM_canon_prop" = "ISM_canon_prop", "NIC_canon_prop" = "NIC_canon_prop", "NNC_canon_prop" = "NNC_canon_prop",
       # ORF related columns
      "cod_FULLSPLICEMATCH" = "Coding_FSM_prop", "ncod_FULLSPLICEMATCH" = "Non_coding_FSM_prop", "cod_INCOMPLETESPLICEMATCH" = "Coding_ISM_prop", "ncod_INCOMPLETESPLICEMATCH" = "Non_coding_ISM_prop",
      "cod_NOVELINCATALOG" = "Coding_NIC_prop", "ncod_NOVELINCATALOG" = "Non_coding_NIC_prop", "cod_NOVELNOTINCATALOG" = "Coding_NNC_prop", "ncod_NOVELNOTINCATALOG" = "Non_coding_NNC_prop",
      "cod_GENIC" = "Coding_genic_prop", "ncod_GENIC" = "Non_coding_genic_prop", "cod_ANTISENSE" = "Coding_antisense_prop", "ncod_ANTISENSE" = "Non_coding_antisense_prop",
      "cod_FUSION" = "Coding_fusion_prop", "ncod_FUSION" = "Non_coding_fusion_prop", "cod_INTERGENIC" = "Coding_intergenic_prop", "ncod_INTERGENIC" = "Non_coding_intergenic_prop",
      "cod_GENICINTRON" = "Coding_genic_intron_prop", "ncod_GENICINTRON" = "Non_coding_genic_intron_prop",
      "NMD_in_cell_prop" = "NMD_prop_in_cell", "FSM_NMD_prop" = "FSM_NMD_prop", "ISM_NMD_prop" = "ISM_NMD_prop", "NIC_NMD_prop" = "NIC_NMD_prop", "NNC_NMD_prop" = "NNC_NMD_prop",
      # CAGE related
      "CAGE_peak_support_in_cell_prop" = "CAGE_peak_support_prop", "FSM_cage_peak_support" = "FSM_CAGE_peak_support_prop", "ISM_cage_peak_support" = "ISM_CAGE_peak_support_prop", "NIC_cage_peak_support" = "NIC_CAGE_peak_support_prop", "NNC_cage_peak_support" = "NNC_CAGE_peak_support_prop",
      # PolyA related
      "polyA_motif_support_in_cell_prop" = "PolyA_motif_support_prop", "FSM_polyA_motif_support" = "FSM_PolyA_motif_support_prop", "ISM_polyA_motif_support" = "ISM_PolyA_motif_support_prop", "NIC_polyA_motif_support" = "NIC_PolyA_motif_support_prop", "NNC_polyA_motif_support" = "NNC_PolyA_motif_support_prop"
    )

  # Rename columns in SQANTI_cell_summary_merged
  current_names_in_merged_df <- names(SQANTI_cell_summary_merged)
  new_names_for_merged <- current_names_in_merged_df
  for(i in 1:length(current_names_in_merged_df)){
      gen_name <- current_names_in_merged_df[i]
      if(gen_name %in% names(current_to_final_map)){
          new_names_for_merged[i] <- current_to_final_map[[gen_name]]
      }
  }
  names(SQANTI_cell_summary_merged) <- new_names_for_merged
  
  # Ensure all columns from col_names_final are present, adding them with 0 if missing
  # This is important if some metrics were not calculated for any cell (e.g. optional CAGE/PolyA)
  for (final_col_name in col_names_final) {
    if (!final_col_name %in% names(SQANTI_cell_summary_merged)) {
      SQANTI_cell_summary_merged[[final_col_name]] <- 0
    }
  }
  
  # Select columns in the final specified order
  SQANTI_cell_summary_final_ordered <- SQANTI_cell_summary_merged[, col_names_final]

  # Change data type of columns from 2nd to last to numeric and fill any remaining NAs with 0
  SQANTI_cell_summary_final_ordered <- SQANTI_cell_summary_final_ordered %>%
    mutate(across(2:ncol(.), as.numeric)) %>%
    mutate(across(everything(), ~replace_na(., 0)))


  if (Save == "Y"){
    print(paste0("Saving cell summary table. Starting at ", Sys.time()))
    write.table(SQANTI_cell_summary_final_ordered, file = gzfile(paste0(cell_summary_output, ".txt.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
    print(paste0("Cell summary table saved as ", cell_summary_output, ".txt.gz"))
  }
  print(paste0("All steps completed at ",Sys.time()))
  return(SQANTI_cell_summary_final_ordered)
}

generate_sqantisc_plots <- function(SQANTI_cell_summary, Classification_file, Junctions, report_output){
  
  ### Basic cell informtion ###
  #############################
  
  # Reads in cell
  gg_reads_in_cells <- ggplot(SQANTI_cell_summary, aes(x = "", y = Reads_in_cell)) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") +  
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) + 
    labs(title = "Number of Reads\nAcross Cells",
         x = "Cell",
         y = "Reads, count") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14))
  
  # UMIs in cell
  gg_umis_in_cells <- ggplot(SQANTI_cell_summary, aes(x = "", y = UMIs_in_cell)) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") +  
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) + 
    labs(title = "Number of UMIs\nAcross Cells",
         x = "Cell",
         y = "UMI, count") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14))
  
  # Genes in cell
  gg_genes_in_cells <- ggplot(SQANTI_cell_summary, aes(x = "", y = Genes_in_cell)) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha = 0.3) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) + 
    labs(title = "Number of Genes\nAcross Cells",
         x = "Cell",
         y = "Genes, count") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14))
  
  # Junctions strings in cell
  gg_JCs_in_cell <- ggplot(SQANTI_cell_summary, aes(x = "", y = UJCs_in_cell)) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") + 
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha=0.3) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) + 
    labs(title = "Number of Unique Junction\nChains Across Cells",
         x = "Cell",
         y = "Unique Junction Chains, count") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14))
  
  # Composite plot of all cell info plots
  # gg_cell_report1 <- grid.arrange(gg_reads_in_cells, gg_umis_in_cells, ncol=2)
  # gg_cell_report2 <- grid.arrange(gg_genes_in_cells, gg_JCs_in_cell, ncol=2)
  
  # Anno/novel genes in cell
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Annotated_genes", "Novel_genes"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Annotated_genes, Novel_genes)))
  
  gg_annotation_of_genes_in_cell <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") + 
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha=0.3) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("Annotated Genes", "Novel Genes")) +
    labs(title = "Number of Known/Novel Genes Across Cells",
         x = "",
         y = "Genes, counts") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  ### Gene Distribution by Read Count Bins ###
  ###########################################
  
  read_bins_data <- data.frame(
    CB = rep(SQANTI_cell_summary$CB, 8),
    bin = rep(c("1", "2-3", "4-5", ">=6", "1", "2-3", "4-5", ">=6"), each = nrow(SQANTI_cell_summary)),
    gene_type = rep(c("Annotated", "Annotated", "Annotated", "Annotated", "Novel", "Novel", "Novel", "Novel"), each = nrow(SQANTI_cell_summary)),
    percentage = c(
      SQANTI_cell_summary$anno_bin1_perc, 
      SQANTI_cell_summary$anno_bin2_3_perc, 
      SQANTI_cell_summary$anno_bin4_5_perc, 
      SQANTI_cell_summary$anno_bin6plus_perc,
      SQANTI_cell_summary$novel_bin1_perc, 
      SQANTI_cell_summary$novel_bin2_3_perc, 
      SQANTI_cell_summary$novel_bin4_5_perc, 
      SQANTI_cell_summary$novel_bin6plus_perc
    )
  )
  
  read_bins_data$bin <- factor(read_bins_data$bin, levels = c("1", "2-3", "4-5", ">=6"))
  read_bins_data$gene_type <- factor(read_bins_data$gene_type, levels = c("Annotated", "Novel"))
  
  gg_read_bins <- ggplot(read_bins_data, aes(x = gene_type, y = percentage, fill = gene_type)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6, show.legend = FALSE) +  # Suppress legend here
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1, show.legend = FALSE) +  # Suppress legend
    scale_fill_manual(values = c("Annotated" = "#e37744", "Novel" = "#78C679")) +
    scale_color_manual(values = c("Annotated" = "#e37744", "Novel" = "#78C679"), guide = "none") +  # Remove color legend
    facet_grid(. ~ bin, scales = "free_x", space = "free", switch = "x") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_classic(base_size = 14) +
    labs(
      title = "Distribution of Genes by Read Count Bins Across Cells",
      x = "",
      y = "Genes, %"
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(0.8, "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_blank(),
      strip.placement = "outside", 
      strip.text.x = element_text(size = 16),
      strip.background = element_blank(),
      legend.text = element_text(size = 14)
    ) +
    guides(
      fill = guide_legend(override.aes = list(
        alpha = 0.7,
        color = "black"
      )),
      color = "none"
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
  
  gg_ujc_bins <- ggplot(ujc_bins_data, aes(x = gene_type, y = percentage, fill = gene_type)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6, show.legend = FALSE) +  # Suppress legend here
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1, show.legend = FALSE) +  # Suppress legend
    scale_fill_manual(values = c("Annotated" = "#e37744", "Novel" = "#78C679")) +
    scale_color_manual(values = c("Annotated" = "#e37744", "Novel" = "#78C679"), guide = "none") +  # Remove color legend
    facet_grid(. ~ bin, scales = "free_x", space = "free", switch = "x") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_classic(base_size = 14) +
    labs(
      title = "Distribution of Genes by UJC Count Bins Across Cells",
      x = "",
      y = "Genes, %"
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(0.8, "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_blank(),
      strip.placement = "outside", 
      strip.text.x = element_text(size = 16),
      strip.background = element_blank(),
      legend.text = element_text(size = 14)
    ) +
    guides(
      fill = guide_legend(override.aes = list(
        alpha = 0.7,
        color = "black"
      )),
      color = "none"
    )

  # Mitochondrial percentage in cell
  gg_MT_perc <- ggplot(SQANTI_cell_summary, aes(x = "", y = MT_perc)) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") +  
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) + 
    labs(title = "Mitochondrial Reads Across Cells",
         x = "Cell",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14))

  #  Mono/multi-exon prop novel vs annotated genes
  
  ### Length distribution ###
  ###########################
  
  # All reads length distribution
  gg_bulk_all_reads <- ggplot(Classification_file, aes(x=length)) +
    geom_histogram(binwidth=50, fill="#CC6633", color="black", alpha=0.5) +
    labs(title = "All Read Lengths Distribution",
         x = "Read length",
         y = "Reads, counts") +
    theme_classic() +
    theme(
     legend.position = "none",
     plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
     axis.title = element_text(size = 16),
     axis.text.y = element_text(size = 14),
     axis.text.x = element_text(size = 16))

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

  gg_bulk_length_by_category <- ggplot(Classification_file, aes(x = length, color = structural_category)) +
    geom_freqpoly(binwidth = 100, linewidth = 1.2, na.rm = TRUE) +
    labs(
      title = "All Read Lengths Distribution by Structural Category",
      x = "Read length",
      y = "Reads, counts"
    ) +
    theme_classic(base_size = 16) +
    scale_color_manual(
      values = c(
        "full-splice_match" = "#6BAED6",
        "incomplete-splice_match" = "#FC8D59",
        "novel_in_catalog" = "#78C679",
        "novel_not_in_catalog" = "#EE6A50",
        "genic" = "#969696",
        "antisense" = "#66C2A4",
        "fusion" = "goldenrod1",
        "intergenic" = "darksalmon",
        "genic_intron" = "#41B6C4"
      ),
      labels = c(
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
    ) +
    scale_x_continuous(
      breaks = seq(0, max(Classification_file$length, na.rm = TRUE), by = 500)
    ) + 
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16)
    )

  # Bulk read length distribution by exonic structure
  Classification_file$exon_type <- ifelse(
    Classification_file$exons == 1, "Mono-Exon", "Multi-Exon"
  )
  Classification_file$exon_type <- factor(
    Classification_file$exon_type, levels = c("Multi-Exon", "Mono-Exon")
  )

  gg_bulk_length_by_exon_type <- ggplot(Classification_file, aes(x = length, color = exon_type)) +
    geom_freqpoly(binwidth = 100, linewidth = 1.2, na.rm = TRUE) +
    labs(
      title = "Mono- vs Multi- Exon Read Lengths Distribution",
      x = "Read length",
      y = "Reads, counts"
    ) +
    theme_classic(base_size = 16) +
    scale_color_manual(
      values = c("Multi-Exon" = "#3B0057", "Mono-Exon" = "#FFE44C")
    ) +
    scale_x_continuous(
      breaks = seq(0, max(Classification_file$length, na.rm = TRUE), by = 500)
    ) + 
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16)
    )


  # Length distribution per break (cells)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Total_250b_length_prop", "Total_500b_length_prop",
                                                                "Total_short_length_prop", "Total_mid_length_prop",
                                                                "Total_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Total_250b_length_prop, Total_500b_length_prop,
                                                                                                       Total_short_length_prop, Total_mid_length_prop,
                                                                                                       Total_long_length_prop)))
  gg_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha=0.3) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Mono-exon length distribution per break
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Total_250b_length_mono_prop", "Total_500b_length_mono_prop",
                                                                "Total_short_length_mono_prop", "Total_mid_length_mono_prop",
                                                                "Total_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Total_250b_length_mono_prop, Total_500b_length_mono_prop,
                                                                                                       Total_short_length_mono_prop, Total_mid_length_mono_prop,
                                                                                                       Total_long_length_mono_prop)))
  gg_read_distr_mono <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha=0.3) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Distribution of Mono-Exonic Read Length Proportions Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Length distribution per break per structural category (cell)
  # FSM
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("FSM_250b_length_prop", "FSM_500b_length_prop",
                                                                "FSM_short_length_prop", "FSM_mid_length_prop",
                                                                "FSM_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(FSM_250b_length_prop, FSM_500b_length_prop,
                                                                                                       FSM_short_length_prop, FSM_mid_length_prop,
                                                                                                       FSM_long_length_prop)))
  gg_FSM_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#6BAED6", color = "#6BAED6", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#6BAED6", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "FSM Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # ISM
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("ISM_250b_length_prop", "ISM_500b_length_prop",
                                                                "ISM_short_length_prop", "ISM_mid_length_prop",
                                                                "ISM_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(ISM_250b_length_prop, ISM_500b_length_prop,
                                                                                                       ISM_short_length_prop, ISM_mid_length_prop,
                                                                                                       ISM_long_length_prop)))
  gg_ISM_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#FC8D59", color = "#FC8D59", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#FC8D59", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "ISM Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # NIC
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("NIC_250b_length_prop", "NIC_500b_length_prop",
                                                                "NIC_short_length_prop", "NIC_mid_length_prop",
                                                                "NIC_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(NIC_250b_length_prop, NIC_500b_length_prop,
                                                                                                       NIC_short_length_prop, NIC_mid_length_prop,
                                                                                                       NIC_long_length_prop)))
  gg_NIC_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#78C679", color = "#78C679", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#78C679", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "NIC Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # NNC
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("NNC_250b_length_prop", "NNC_500b_length_prop",
                                                                "NNC_short_length_prop", "NNC_mid_length_prop",
                                                                "NNC_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(NNC_250b_length_prop, NNC_500b_length_prop,
                                                                                                       NNC_short_length_prop, NNC_mid_length_prop,
                                                                                                       NNC_long_length_prop)))
  gg_NNC_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#EE6A50", color = "#EE6A50", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#EE6A50", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "NNC Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Genic
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Genic_250b_length_prop", "Genic_500b_length_prop",
                                                                "Genic_short_length_prop", "Genic_mid_length_prop",
                                                                "Genic_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Genic_250b_length_prop, Genic_500b_length_prop,
                                                                                                       Genic_short_length_prop, Genic_mid_length_prop,
                                                                                                       Genic_long_length_prop)))
  gg_genic_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#969696", color = "#969696", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#969696", color = "grey90", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Genic Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Antisense
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Antisense_250b_length_prop", "Antisense_500b_length_prop",
                                                                "Antisense_short_length_prop", "Antisense_mid_length_prop",
                                                                "Antisense_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Antisense_250b_length_prop, Antisense_500b_length_prop,
                                                                                                       Antisense_short_length_prop, Antisense_mid_length_prop,
                                                                                                       Antisense_long_length_prop)))
  gg_antisense_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#66C2A4", color = "#66C2A4", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#66C2A4", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Antisense Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Fusion
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Fusion_250b_length_prop", "Fusion_500b_length_prop",
                                                                "Fusion_short_length_prop", "Fusion_mid_length_prop",
                                                                "Fusion_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Fusion_250b_length_prop, Fusion_500b_length_prop,
                                                                                                       Fusion_short_length_prop, Fusion_mid_length_prop,
                                                                                                       Fusion_long_length_prop)))
  gg_fusion_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "goldenrod1", color = "goldenrod1", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "goldenrod1", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Fusion Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))

  # Intergenic
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Intergenic_250b_length_prop", "Intergenic_500b_length_prop",
                                                                "Intergenic_short_length_prop", "Intergenic_mid_length_prop",
                                                                "Intergenic_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Intergenic_250b_length_prop, Intergenic_500b_length_prop,
                                                                                                       Intergenic_short_length_prop, Intergenic_mid_length_prop,
                                                                                                       Intergenic_long_length_prop)))
  gg_intergenic_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "darksalmon", color = "darksalmon", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "darksalmon", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Intergenic Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Genic intron
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Genic_intron_250b_length_prop", "Genic_intron_500b_length_prop",
                                                                "Genic_intron_short_length_prop", "Genic_intron_mid_length_prop",
                                                                "Genic_intron_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Genic_intron_250b_length_prop, Genic_intron_500b_length_prop,
                                                                                                       Genic_intron_short_length_prop, Genic_intron_mid_length_prop,
                                                                                                       Genic_intron_long_length_prop)))
  gg_genic_intron_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#41B6C4", color = "#41B6C4", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#41B6C4", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Genic Intron Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Mono-exon length distribution across categories
  # FSM
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("FSM_250b_length_mono_prop", "FSM_500b_length_mono_prop",
                                                                "FSM_short_length_mono_prop", "FSM_mid_length_mono_prop",
                                                                "FSM_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(FSM_250b_length_mono_prop, FSM_500b_length_mono_prop,
                                                                                                       FSM_short_length_mono_prop, FSM_mid_length_mono_prop,
                                                                                                       FSM_long_length_mono_prop)))
  gg_FSM_mono_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#6BAED6", color = "#6BAED6", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#6BAED6", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "FSM Mono-exonic Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # ISM
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("ISM_250b_length_mono_prop", "ISM_500b_length_mono_prop",
                                                                "ISM_short_length_mono_prop", "ISM_mid_length_mono_prop",
                                                                "ISM_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(ISM_250b_length_mono_prop, ISM_500b_length_mono_prop,
                                                                                                       ISM_short_length_mono_prop, ISM_mid_length_mono_prop,
                                                                                                       ISM_long_length_mono_prop)))
  gg_ISM_mono_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#FC8D59", color = "#FC8D59", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#FC8D59", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "ISM Mono-exonic Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # NIC
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("NIC_250b_length_mono_prop", "NIC_500b_length_mono_prop",
                                                                "NIC_short_length_mono_prop", "NIC_mid_length_mono_prop",
                                                                "NIC_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(NIC_250b_length_mono_prop, NIC_500b_length_mono_prop,
                                                                                                       NIC_short_length_mono_prop, NIC_mid_length_mono_prop,
                                                                                                       NIC_long_length_mono_prop)))
  gg_NIC_mono_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#78C679", color = "#78C679", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#78C679", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "NIC Mono-exonic Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Genic
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Genic_250b_length_mono_prop", "Genic_500b_length_mono_prop",
                                                                "Genic_short_length_mono_prop", "Genic_mid_length_mono_prop",
                                                                "Genic_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Genic_250b_length_mono_prop, Genic_500b_length_mono_prop,
                                                                                                       Genic_short_length_mono_prop, Genic_mid_length_mono_prop,
                                                                                                       Genic_long_length_mono_prop)))
  gg_genic_mono_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#969696", color = "#969696", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#969696", color = "grey90", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Genic Mono-exonic Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Antisense
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Antisense_250b_length_mono_prop", "Antisense_500b_length_mono_prop",
                                                                "Antisense_short_length_mono_prop", "Antisense_mid_length_mono_prop",
                                                                "Antisense_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Antisense_250b_length_mono_prop, Antisense_500b_length_mono_prop,
                                                                                                       Antisense_short_length_mono_prop, Antisense_mid_length_mono_prop,
                                                                                                       Antisense_long_length_mono_prop)))
  gg_antisense_mono_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#66C2A4", color = "#66C2A4", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#66C2A4", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Antisense Mono-exonic Reads Length Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Intergenic
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Intergenic_250b_length_mono_prop", "Intergenic_500b_length_mono_prop",
                                                                "Intergenic_short_length_mono_prop", "Intergenic_mid_length_mono_prop",
                                                                "Intergenic_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Intergenic_250b_length_mono_prop, Intergenic_500b_length_mono_prop,
                                                                                                       Intergenic_short_length_mono_prop, Intergenic_mid_length_mono_prop,
                                                                                                       Intergenic_long_length_mono_prop)))
  gg_intergenic_mono_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "darksalmon", color = "darksalmon", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "darksalmon", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Intergenic Mono-exonic Read Lengths Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Genic intron
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Genic_intron_250b_length_mono_prop", "Genic_intron_500b_length_mono_prop",
                                                                "Genic_intron_short_length_mono_prop", "Genic_intron_mid_length_mono_prop",
                                                                "Genic_intron_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Genic_intron_250b_length_mono_prop, Genic_intron_500b_length_mono_prop,
                                                                                                       Genic_intron_short_length_mono_prop, Genic_intron_mid_length_mono_prop,
                                                                                                       Genic_intron_long_length_mono_prop)))
  gg_genic_intron_mono_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(fill = "#41B6C4", color = "#41B6C4", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#41B6C4", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Genic Intron Mono-exonic Read Lengths Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  ### Reference coverage across categories ###
  ############################################
  
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("FSM_ref_coverage_prop",
                                                                "ISM_ref_coverage_prop",
                                                                "NIC_ref_coverage_prop",
                                                                "NNC_ref_coverage_prop",
                                                                "Genic_ref_coverage_prop",
                                                                "Antisense_ref_coverage_prop",
                                                                "Fusion_ref_coverage_prop",
                                                                "Intergenic_ref_coverage_prop",
                                                                "Genic_intron_ref_coverage_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(FSM_ref_coverage_prop,
                                                                                                       ISM_ref_coverage_prop,
                                                                                                       NIC_ref_coverage_prop,
                                                                                                       NNC_ref_coverage_prop,
                                                                                                       Genic_ref_coverage_prop,
                                                                                                       Antisense_ref_coverage_prop,
                                                                                                       Fusion_ref_coverage_prop,
                                                                                                       Intergenic_ref_coverage_prop,
                                                                                                       Genic_intron_ref_coverage_prop)))
  
  gg_ref_coverage_across_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +

    geom_violin(aes(color = Variable, 
                    fill = Variable), scale = "width", alpha = 0.7) +  
    geom_boxplot(aes(fill = Variable), color = c(rep("grey20",4),"grey90",rep("grey20",4)),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_fill_manual(values = c("FSM_ref_coverage_prop"="#6BAED6", "ISM_ref_coverage_prop"="#FC8D59",
                                 "NIC_ref_coverage_prop"="#78C679", "NNC_ref_coverage_prop"="#EE6A50", 
                                 "Genic_ref_coverage_prop"="#969696", "Antisense_ref_coverage_prop"="#66C2A4",
                                 "Fusion_ref_coverage_prop"="goldenrod1", "Intergenic_ref_coverage_prop"="darksalmon",
                                 "Genic_intron_ref_coverage_prop"="#41B6C4")) + 
    scale_color_manual(values = c("FSM_ref_coverage_prop"="#6BAED6", "ISM_ref_coverage_prop"="#FC8D59",
                                  "NIC_ref_coverage_prop"="#78C679", "NNC_ref_coverage_prop"="#EE6A50", 
                                  "Genic_ref_coverage_prop"="#969696", "Antisense_ref_coverage_prop"="#66C2A4",
                                  "Fusion_ref_coverage_prop"="goldenrod1", "Intergenic_ref_coverage_prop"="darksalmon",
                                  "Genic_intron_ref_coverage_prop"="#41B6C4")) +
    scale_x_discrete(labels = c("FSM","ISM","NIC","NNC","Genic\nGenomic",
                                "Antisense","Fusion","Intergenic","Genic\nintron")) +
    theme_classic(base_size = 14) +  
    labs(title = "Reference Length Coverage by Structural Category Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  
  ### Structural categories ###
  #############################
  
  # Isoform Distribution Across Structural Categories
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("FSM_prop", "ISM_prop", "NIC_prop", "NNC_prop", "Genic_Genomic_prop",
                                                                "Antisense_prop", "Fusion_prop", "Intergenic_prop", "Genic_intron_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(FSM_prop,ISM_prop,NIC_prop,NNC_prop,
                                                                                                       Genic_Genomic_prop,Antisense_prop,
                                                                                                       Fusion_prop,Intergenic_prop,Genic_intron_prop)))
  
  gg_SQANTI_across_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c(rep("grey20",4),"grey90",rep("grey20",4)),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_fill_manual(values = c("FSM_prop"="#6BAED6", "ISM_prop"="#FC8D59", "NIC_prop"="#78C679", "NNC_prop"="#EE6A50", 
                                 "Genic_Genomic_prop"="#969696", "Antisense_prop"="#66C2A4", "Fusion_prop"="goldenrod1",
                                 "Intergenic_prop"="darksalmon", "Genic_intron_prop"="#41B6C4")) + 
    scale_color_manual(values = c("FSM_prop"="#6BAED6", "ISM_prop"="#FC8D59", "NIC_prop"="#78C679", "NNC_prop"="#EE6A50",
                                  "Genic_Genomic_prop"="#969696", "Antisense_prop"="#66C2A4", "Fusion_prop"="goldenrod1",
                                  "Intergenic_prop"="darksalmon", "Genic_intron_prop"="#41B6C4")) +
    scale_x_discrete(labels = c("FSM","ISM","NIC","NNC","Genic\nGenomic",
                                "Antisense","Fusion","Intergenic","Genic\nintron")) +
    theme_classic(base_size = 14) +  
    labs(title = "Structural Categories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  #  Coding/non-coding across structural categories (change it in the future to a combine plot)
  if (!skipORF) {
    gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Coding_FSM_prop",
                                                                  "Coding_ISM_prop",
                                                                  "Coding_NIC_prop",
                                                                  "Coding_NNC_prop",
                                                                  "Coding_genic_prop",
                                                                  "Coding_antisense_prop",
                                                                  "Coding_fusion_prop",
                                                                  "Coding_intergenic_prop",
                                                                  "Coding_genic_intron_prop"), 
                                    names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
    gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Coding_FSM_prop,
                                                                                                         Coding_ISM_prop,
                                                                                                         Coding_NIC_prop,
                                                                                                         Coding_NNC_prop,
                                                                                                         Coding_genic_prop,
                                                                                                         Coding_antisense_prop,
                                                                                                         Coding_fusion_prop,
                                                                                                         Coding_intergenic_prop,
                                                                                                         Coding_genic_intron_prop)))
  
    gg_coding_across_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
      geom_violin(aes(color = Variable, 
                      fill = Variable),
                      alpha = 0.7, scale = "width") +  
      geom_boxplot(aes(fill = Variable), color = c(rep("grey20",4),"grey90",rep("grey20",4)),
                   width = 0.08, outlier.shape = NA, alpha = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
      scale_fill_manual(values = c("Coding_FSM_prop"="#6BAED6", "Coding_ISM_prop"="#FC8D59", "Coding_NIC_prop"="#78C679", "Coding_NNC_prop"="#EE6A50", 
                                   "Coding_genic_prop"="#969696", "Coding_antisense_prop"="#66C2A4", "Coding_fusion_prop"="goldenrod1",
                                   "Coding_intergenic_prop"="darksalmon", "Coding_genic_intron_prop"="#41B6C4")) + 
      scale_color_manual(values = c("Coding_FSM_prop"="#6BAED6", "Coding_ISM_prop"="#FC8D59", "Coding_NIC_prop"="#78C679", "Coding_NNC_prop"="#EE6A50", 
                                    "Coding_genic_prop"="#969696", "Coding_antisense_prop"="#66C2A4", "Coding_fusion_prop"="goldenrod1",
                                    "Coding_intergenic_prop"="darksalmon", "Coding_genic_intron_prop"="#41B6C4")) +
      scale_x_discrete(labels = c("FSM","ISM","NIC","NNC","Genic\nGenomic",
                                  "Antisense","Fusion","Intergenic","Genic\nintron")) +
      theme_classic(base_size = 14) +  
      labs(title = "Coding Proportion of Structural Categories Distribution Across Cells",
           x = "",
           y = "Reads, %") +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        axis.title = element_text(size = 16), 
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
      )
  
  
    gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Non_coding_FSM_prop",
                                                                  "Non_coding_ISM_prop",
                                                                  "Non_coding_NIC_prop",
                                                                  "Non_coding_NNC_prop",
                                                                  "Non_coding_genic_prop",
                                                                  "Non_coding_antisense_prop",
                                                                  "Non_coding_fusion_prop",
                                                                  "Non_coding_intergenic_prop",
                                                                  "Non_coding_genic_intron_prop"), 
                                    names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
    gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Non_coding_FSM_prop,
                                                                                                         Non_coding_ISM_prop,
                                                                                                         Non_coding_NIC_prop,
                                                                                                         Non_coding_NNC_prop,
                                                                                                         Non_coding_genic_prop,
                                                                                                         Non_coding_antisense_prop,
                                                                                                         Non_coding_fusion_prop,
                                                                                                         Non_coding_intergenic_prop,
                                                                                                         Non_coding_genic_intron_prop)))
  
    gg_non_coding_across_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
      geom_violin(aes(color = Variable, 
                      fill = Variable),
                      alpha = 0.7, scale = "width") +  
      geom_boxplot(aes(fill = Variable), color = "grey20",
                   width = 0.08, outlier.shape = NA, alpha = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
      scale_fill_manual(values = c("Non_coding_FSM_prop"="#D2E6F2", "Non_coding_ISM_prop"="#FEDCCD", "Non_coding_NIC_prop"="#D6EDD6", "Non_coding_NNC_prop"="#F9D2CA", 
                                   "Non_coding_genic_prop"="#DFDFDF", "Non_coding_antisense_prop"="#D1ECE3", "Non_coding_fusion_prop"="#FFECBD",
                                   "Non_coding_intergenic_prop"="#F8DFD7", "Non_coding_genic_intron_prop"="#C6E9ED")) + 
      scale_color_manual(values = c("Non_coding_FSM_prop"="#D2E6F2", "Non_coding_ISM_prop"="#FEDCCD", "Non_coding_NIC_prop"="#D6EDD6", "Non_coding_NNC_prop"="#F9D2CA", 
                                    "Non_coding_genic_prop"="#DFDFDF", "Non_coding_antisense_prop"="#D1ECE3", "Non_coding_fusion_prop"="#FFECBD",
                                    "Non_coding_intergenic_prop"="#F8DFD7", "Non_coding_genic_intron_prop"="#C6E9ED")) +
      scale_x_discrete(labels = c("FSM","ISM","NIC","NNC","Genic\nGenomic",
                                  "Antisense","Fusion","Intergenic","Genic\nintron")) +
      theme_classic(base_size = 14) +  
      labs(title = "Non-coding Proportion of Structural Categories Distribution Across Cells",
           x = "",
           y = "Reads, %") +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        axis.title = element_text(size = 16), 
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
      )
  } # End of if (!skipORF)

  # Isoform Distribution Across FSM 
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("FSM_alternative_3.end_prop","FSM_alternative_3.5.end_prop",
                                                                "FSM_alternative_5.end_prop","FSM_reference_match_prop",
                                                                "FSM_mono.exon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(FSM_alternative_3.end_prop,
                                                                                                       FSM_alternative_3.5.end_prop,
                                                                                                       FSM_alternative_5.end_prop,
                                                                                                       FSM_reference_match_prop,
                                                                                                       FSM_mono.exon_prop)))
  
  gg_SQANTI_across_FSM <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c(rep("grey90",2),rep("grey20",3)),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_color_manual(values = c("FSM_alternative_3.end_prop"='#02314d', "FSM_alternative_3.5.end_prop"='#0e5a87',
                                  "FSM_alternative_5.end_prop"='#7ccdfc', 'FSM_reference_match_prop'='#c4e1f2',
                                  "FSM_mono.exon_prop"='#cec2d2')) + 
    scale_fill_manual(values = c("FSM_alternative_3.end_prop"='#02314d', "FSM_alternative_3.5.end_prop"='#0e5a87',
                                 "FSM_alternative_5.end_prop"='#7ccdfc', 'FSM_reference_match_prop'='#c4e1f2',
                                 "FSM_mono.exon_prop"='#cec2d2')) +  
    scale_x_discrete(labels = c("Alternative 3'end", "Alternative 3'5'end", "Alternative 5'end",
                                'Reference match', "Mono-exon")) +
    theme_classic(base_size = 14) +  
    labs(title = "FSM Structural Subcategories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  # Isoform Distribution Across ISM
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("ISM_3._fragment_prop", "ISM_internal_fragment_prop",
                                                                "ISM_5._fragment_prop", "ISM_intron_retention_prop", "ISM_mono.exon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(ISM_3._fragment_prop,
                                                                                                       ISM_internal_fragment_prop,
                                                                                                       ISM_5._fragment_prop,
                                                                                                       ISM_intron_retention_prop,
                                                                                                       ISM_mono.exon_prop))) 
  
  gg_SQANTI_across_ISM <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c(rep("grey90",2),rep("grey20",3)),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_color_manual(values = c("ISM_3._fragment_prop"='#c4531d', "ISM_internal_fragment_prop"='#e37744',  
                                  "ISM_5._fragment_prop"='#e0936e', "ISM_intron_retention_prop"='#81eb82',
                                  "ISM_mono.exon_prop"='#cec2d2')) + 
    scale_fill_manual(values = c("ISM_3._fragment_prop"='#c4531d', "ISM_internal_fragment_prop"='#e37744',  
                                 "ISM_5._fragment_prop"='#e0936e', "ISM_intron_retention_prop"='#81eb82',
                                 "ISM_mono.exon_prop"='#cec2d2')) +  
    scale_x_discrete(labels = c("3' fragment", "Internal fragment", "5' fragment",
                                'Intron retention', "Mono-exon")) +
    theme_classic(base_size = 14) +  
    labs(title = "ISM Structural Subcategories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  # Isoform Distribution Across NIC
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("NIC_comb_annot_junctions_prop", 
                                                                "NIC_comb_annot_splice_sites_prop",
                                                                "NIC_intron_retention_prop",
                                                                "NIC_mono.exon_by_intron_retention_prop",
                                                                "NIC_mono.exon_prop"), 
                                  names_to = "Variable", values_to = "Value")
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(NIC_comb_annot_junctions_prop, 
                                                                                                       NIC_comb_annot_splice_sites_prop,
                                                                                                       NIC_intron_retention_prop,
                                                                                                       NIC_mono.exon_by_intron_retention_prop,
                                                                                                       NIC_mono.exon_prop)))
  gg_SQANTI_across_NIC <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c(rep("grey90",2),"grey20",rep("grey90",2)),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_color_manual(values = c("NIC_comb_annot_junctions_prop"='#014d02', "NIC_comb_annot_splice_sites_prop"='#379637',  
                                  "NIC_intron_retention_prop"='#81eb82', "NIC_mono.exon_by_intron_retention_prop"="#4aaa72",
                                  "NIC_mono-exon_prop"="#cec2d2")) + 
    scale_fill_manual(values = c("NIC_comb_annot_junctions_prop"='#014d02', "NIC_comb_annot_splice_sites_prop"='#379637',  
                                 "NIC_intron_retention_prop"='#81eb82', "NIC_mono.exon_by_intron_retention_prop"="#4aaa72",
                                 "NIC_mono-exon_prop"="#cec2d2")) +  
    scale_x_discrete(labels = c("Comb. of annot. junctions", "Comb. of annot. splice sites",  
                                "Intron retention", "Mono-exon by intron ret.", "Mono-exon")) +
    theme_classic(base_size = 14) +  
    labs(title = "NIC Structural Subcategories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  # Isoform Distribution Across NNC
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("NNC_at_least_1_don_accept_prop",
                                                                "NNC_intron_retention_prop"), 
                                  names_to = "Variable", values_to = "Value")
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(NNC_at_least_1_don_accept_prop,
                                                                                                       NNC_intron_retention_prop)))
  
  gg_SQANTI_across_NNC <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c("grey90","grey20"),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_color_manual(values = c("NNC_at_least_1_don_accept_prop"="#32734d", "NNC_intron_retention_prop"="#81eb82")) + 
    scale_fill_manual(values = c("NNC_at_least_1_don_accept_prop"="#32734d", "NNC_intron_retention_prop"="#81eb82")) +  
    scale_x_discrete(labels = c("At least\n1 annot. don./accept.", "Intron retention")) +
    theme_classic(base_size = 14) +  
    labs(title = "NNC Structural Subcategories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  # Isoform Distribution Across Fusion
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary,
                                  cols = c("Fusion_intron_retention_prop", 
                                           "Fusion_multi.exon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable,
                                     colnames(SQANTI_cell_summary %>% select(Fusion_intron_retention_prop, 
                                                                             Fusion_multi.exon_prop)))
  
  gg_SQANTI_across_Fusion <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c("grey20","grey90"),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_color_manual(values = c("Fusion_intron_retention_prop"="#81eb82", "Fusion_multi.exon_prop"="#876a91")) +
    scale_fill_manual(values = c("Fusion_intron_retention_prop"="#81eb82", "Fusion_multi.exon_prop"="#876a91")) +  
    scale_x_discrete(labels = c("Intron retention", "Multi-exon")) +
    theme_classic(base_size = 14) +  
    labs(title = "Fusion Structural Subcategories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  # Isoform Distribution Across Genic
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary,
                                  cols = c("Genic_mono.exon_prop", 
                                           "Genic_multi.exon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable,
                                     colnames(SQANTI_cell_summary %>% select(Genic_mono.exon_prop, 
                                                                             Genic_multi.exon_prop)))
  
  gg_SQANTI_across_Genic <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c("grey20","grey90"),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_color_manual(values = c("Genic_mono.exon_prop"="#81eb82", "Genic_multi.exon_prop"="#876a91")) +
    scale_fill_manual(values = c("Genic_mono.exon_prop"="#81eb82", "Genic_multi.exon_prop"="#876a91")) +  
    scale_x_discrete(labels = c("Mono-exon", "Multi-exon")) +
    theme_classic(base_size = 14) +  
    labs(title = "Genic Structural Subcategories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  # Isoform Distribution Across Genic Intron
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary,
                                  cols = c("Genic_intron_mono.exon_prop", 
                                           "Genic_intron_multi.exon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable,
                                     colnames(SQANTI_cell_summary %>% select(Genic_intron_mono.exon_prop, 
                                                                             Genic_intron_multi.exon_prop)))
  
  gg_SQANTI_across_Genic_Intron <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c("grey20","grey90"),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_color_manual(values = c("Genic_intron_mono.exon_prop"="#81eb82", "Genic_intron_multi.exon_prop"="#876a91")) +
    scale_fill_manual(values = c("Genic_intron_mono.exon_prop"="#81eb82", "Genic_intron_multi.exon_prop"="#876a91")) +  
    scale_x_discrete(labels = c("Mono-exon", "Multi-exon")) +
    theme_classic(base_size = 14) +  
    labs(title = "Genic Intron Structural Subcategories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  # Isoform Distribution Across Antisense
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary,
                                  cols = c("Antisense_mono.exon_prop", 
                                           "Antisense_multi.exon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable,
                                     colnames(SQANTI_cell_summary %>% select(Antisense_mono.exon_prop, 
                                                                             Antisense_multi.exon_prop)))
  
  gg_SQANTI_across_Antisense <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c("grey20","grey90"),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_color_manual(values = c("Antisense_mono.exon_prop"="#81eb82", "Antisense_multi.exon_prop"="#876a91")) +
    scale_fill_manual(values = c("Antisense_mono.exon_prop"="#81eb82", "Antisense_multi.exon_prop"="#876a91")) +  
    scale_x_discrete(labels = c("Mono-exon", "Multi-exon")) +
    theme_classic(base_size = 14) +  
    labs(title = "Antisense Structural Subcategories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  # Isoform Distribution Across Intergenic
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary,
                                  cols = c("Intergenic_mono.exon_prop", 
                                           "Intergenic_multi.exon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable,
                                     colnames(SQANTI_cell_summary %>% select(Intergenic_mono.exon_prop, 
                                                                             Intergenic_multi.exon_prop)))
  
  gg_SQANTI_across_Intergenic <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = c("grey20","grey90"),
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    scale_color_manual(values = c("Intergenic_mono.exon_prop"="#81eb82", "Intergenic_multi.exon_prop"="#876a91")) +
    scale_fill_manual(values = c("Intergenic_mono.exon_prop"="#81eb82", "Intergenic_multi.exon_prop"="#876a91")) +  
    scale_x_discrete(labels = c("Mono-exon", "Multi-exon")) +
    theme_classic(base_size = 14) +  
    labs(title = "Intergenic Structural Subcategories Distribution Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
    )
  
  ### Splice junctions characterization ###
  #########################################
  
  # Known/novel canonical/non-canonical SJs
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Known_canonical_junctions_prop", "Known_non_canonical_junctions_prop",
                                                                "Novel_canonical_junctions_prop", "Novel_non_canonical_junctions_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Known_canonical_junctions_prop, Known_non_canonical_junctions_prop,
                                                                                                       Novel_canonical_junctions_prop, Novel_non_canonical_junctions_prop)))
  gg_known_novel_canon <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    theme_classic(base_size = 14) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = c("Known_canonical_junctions_prop" = "#6BAED6",
                                  "Known_non_canonical_junctions_prop" = "goldenrod1",
                                  "Novel_canonical_junctions_prop" = "#78C679",
                                  "Novel_non_canonical_junctions_prop" = "#FC8D59")) +
    scale_fill_manual(values = c("Known_canonical_junctions_prop" = "#6BAED6",
                                  "Known_non_canonical_junctions_prop" = "goldenrod1",
                                  "Novel_canonical_junctions_prop" = "#78C679",
                                  "Novel_non_canonical_junctions_prop" = "#FC8D59")) +
    scale_x_discrete(labels = c("Known\nCanonical", "Known\nNon-canonical",
                                "Novel\nCanonical", "Novel\nNon-canonical")) +
    labs(title = "Distribution of Splice Junctions Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16))
  
  ### Bad features plots ###
  ##########################
  
  # Intrapriming  (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_intrapriming_prop", "ISM_intrapriming_prop", 
                                           "NIC_intrapriming_prop", "NNC_intrapriming_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_intrapriming_prop", "ISM_intrapriming_prop", 
                                               "NIC_intrapriming_prop", "NNC_intrapriming_prop"))
  
  gg_intrapriming_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    # Use the same intrapriming green color from the original bad_feature plot
    scale_color_manual(values = rep("#78C679", 4)) +
    scale_fill_manual(values = rep("#78C679", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "Intrapriming by Structural Category",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))

  # RTS  (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_RTS_prop", "ISM_RTS_prop", 
                                           "NIC_RTS_prop", "NNC_RTS_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_RTS_prop", "ISM_RTS_prop", 
                                               "NIC_RTS_prop", "NNC_RTS_prop"))
  
  gg_RTS_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = rep("#FF9933", 4)) +
    scale_fill_manual(values = rep("#FF9933", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "RT-switching by Structural Category",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))

  # Non-canonical  (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_noncanon_prop", "ISM_noncanon_prop", 
                                           "NIC_noncanon_prop", "NNC_noncanon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_noncanon_prop", "ISM_noncanon_prop", 
                                               "NIC_noncanon_prop", "NNC_noncanon_prop"))
  
  gg_noncanon_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    # Use the same non-canonical blue color from the original bad_feature plot
    scale_color_manual(values = rep("#41B6C4", 4)) +
    scale_fill_manual(values = rep("#41B6C4", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "Non-Canonical Junctions by Structural Category",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))

  # NMD  (split between categories)
  nmd_cols <- c("FSM_NMD_prop", "ISM_NMD_prop", "NIC_NMD_prop", "NNC_NMD_prop")
  if (all  (nmd_cols %in% colnames(SQANTI_cell_summary))) {
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_NMD_prop", "ISM_NMD_prop", 
                                           "NIC_NMD_prop", "NNC_NMD_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_NMD_prop", "ISM_NMD_prop", 
                                               "NIC_NMD_prop", "NNC_NMD_prop"))
  
  gg_NMD_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = rep("#969696", 4)) +
    scale_fill_manual(values = rep("#969696", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "Nonsense-Mediated Decay by Structural Category",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
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
  bad_feature_cols_present <- bad_feature_cols_present[sapply(bad_feature_cols_present, function(col) any(!is.na(SQANTI_cell_summary[[col]])) && sum(SQANTI_cell_summary[[col]], na.rm=TRUE) > 0 )] # keep only if data exists

  # Order them as originally intended, if present
  ordered_bad_feature_cols <- c("Intrapriming_prop_in_cell", "RTS_prop_in_cell", "Non_canonical_prop_in_cell", "NMD_prop_in_cell")
  bad_feature_cols_present <- intersect(ordered_bad_feature_cols, bad_feature_cols_present)


  if (length(bad_feature_cols_present) > 0) {
    current_colors <- sapply(all_bad_features_map[bad_feature_cols_present], function(x) x$color)
    current_labels <- sapply(all_bad_features_map[bad_feature_cols_present], function(x) x$label)
    # Ensure names are correctly assigned for scales, matching the order in bad_feature_cols_present
    names(current_colors) <- bad_feature_cols_present
    names(current_labels) <- bad_feature_cols_present

    gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = all_of(bad_feature_cols_present),
                                    names_to = "Variable", values_to = "Value") %>% select(Variable, Value)

    gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, levels = bad_feature_cols_present)

    gg_bad_feature <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
      geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +
      geom_boxplot(aes(fill = Variable), color = "grey20", width = 0.08, outlier.shape = NA, alpha = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
      theme_classic(base_size = 14) +
      scale_color_manual(values = current_colors, labels = current_labels, name = NULL) +
      scale_fill_manual(values = current_colors, labels = current_labels, name = NULL) +
      scale_x_discrete(labels = current_labels) +
      labs(title = "Bad Quality Control Attributes Across Cells", x = "", y = "Reads, %") +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16)
      )
  } else {
    gg_bad_feature <- ggplot() + theme_void() + ggtitle("No bad quality features to display") +
                       theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  }

  ### Good features plots ###
  ##########################
  # Junction strings mapped to annotated transcripts  (split between categories) # maybe not
  # TSS annotation support (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_TSSAnnotationSupport", "ISM_TSSAnnotationSupport", 
                                  "NIC_TSSAnnotationSupport", "NNC_TSSAnnotationSupport"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_TSSAnnotationSupport", "ISM_TSSAnnotationSupport", 
                                    "NIC_TSSAnnotationSupport", "NNC_TSSAnnotationSupport"))

  gg_tss_annotation_support <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = rep("#66C2A4", 4)) +
    scale_fill_manual(values = rep("#66C2A4", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "TSS Annotation Support by Structural Category",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))

  # CAGE peak support  (split between categories)
  cage_peak_cols <- c("FSM_CAGE_peak_support_prop", "ISM_CAGE_peak_support_prop", 
                      "NIC_CAGE_peak_support_prop", "NNC_CAGE_peak_support_prop")
  if (all(cage_peak_cols %in% colnames(SQANTI_cell_summary))) {
    gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_CAGE_peak_support_prop", "ISM_CAGE_peak_support_prop", 
                                           "NIC_CAGE_peak_support_prop", "NNC_CAGE_peak_support_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
    gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                      levels = c("FSM_CAGE_peak_support_prop", "ISM_CAGE_peak_support_prop", 
                                                "NIC_CAGE_peak_support_prop", "NNC_CAGE_peak_support_prop")) 

    gg_cage_peak_support <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
      geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
      geom_boxplot(aes(fill = Variable), color = "grey20",
                  width = 0.08, outlier.shape = NA, alpha = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
      theme_classic(base_size = 14) +
      scale_color_manual(values = rep("#EE6A50", 4)) +
      scale_fill_manual(values = rep("#EE6A50", 4)) + 
      scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
      labs(title = "CAGE Peak Support by Structural Category",
          x = "",
          y = "Reads, %") +
      theme(
        legend.position = "none", 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        axis.title = element_text(size = 16), 
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 16))
  }

  # PolyA motif support  (split between categories)
  polyA_motif_cols <- c("FSM_PolyA_motif_support_prop", "ISM_PolyA_motif_support_prop", 
                       "NIC_PolyA_motif_support_prop", "NNC_PolyA_motif_support_prop")
  if (all(polyA_motif_cols %in% colnames(SQANTI_cell_summary))) {
    gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_PolyA_motif_support_prop", "ISM_PolyA_motif_support_prop", 
                                           "NIC_PolyA_motif_support_prop", "NNC_PolyA_motif_support_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value) 

    gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                      levels = c("FSM_PolyA_motif_support_prop", "ISM_PolyA_motif_support_prop", 
                                                "NIC_PolyA_motif_support_prop", "NNC_PolyA_motif_support_prop"))

    gg_polyA_motif_support <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
      geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
      geom_boxplot(aes(fill = Variable), color = "grey20",
                  width = 0.08, outlier.shape = NA, alpha = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
      theme_classic(base_size = 14) +
      scale_color_manual(values = rep("#78C679", 4)) +
      scale_fill_manual(values = rep("#78C679", 4)) + 
      scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
      labs(title = "PolyA Support by Structural Category",
          x = "",
          y = "Reads, %") +
      theme(
        legend.position = "none", 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),   
        axis.title = element_text(size = 16), 
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 16))
  } 

  # Canonical  (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_canon_prop", "ISM_canon_prop", 
                                           "NIC_canon_prop", "NNC_canon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_canon_prop", "ISM_canon_prop", 
                                               "NIC_canon_prop", "NNC_canon_prop"))
  
  gg_canon_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = rep("#CC6633", 4)) +
    scale_fill_manual(values = rep("#CC6633", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "Canonical Junctions by Structural Category",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))

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

  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = all_of(good_feature_cols), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, good_feature_cols)
  gg_good_feature <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map) +
    scale_x_discrete(labels = label_map) +
    labs(title = "Good Quality Control Attributes Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16))
  
  ### Presets ###
  ###############
  
  # t1 <- ttheme_default(core=list(core = list(fg_params = list(cex = 0.6)),
  #                                colhead = list(fg_params = list(cex = 0.7))))
  
  ### Generate PDF report ###
  ###########################
  
  pdf(file.path(paste0(report_output,".pdf")), paper = "a4r", width = 14, height = 11)
  # Add cover page
  grid.newpage()
  cover <- textGrob("SQANTI-single cell\nreads report",
                    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
  grid.draw(cover)
  # Bulk tables
  s <- textGrob("Bulk summary", gp=gpar(fontface="italic", fontsize=30), vjust = 0)
  grid.arrange(s)

  # Calculate bulk-level stats
  total_reads_count <- nrow(Classification_file)
  unique_genes <- length(unique(Classification_file$associated_gene))
  unique_junctions <- length(unique(Classification_file$jxn_string))

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

  read_class_table <- as.data.frame(table(factor(Classification_file$structural_category, levels = read_cat_levels)))
  colnames(read_class_table) <- c("Category", "Reads, count")
  read_class_table$Category <- read_cat_names

  # Splice Junction Classification table
  Junctions$junction_type <- paste(Junctions$junction_category, Junctions$canonical, sep = "_")
  
  sj_types <- c("known_canonical", "known_non_canonical", "novel_canonical", "novel_non_canonical")
  sj_counts <- sapply(sj_types, function(type) sum(Junctions$junction_type == type, na.rm = TRUE))
  
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

  # Table theme with larger font
  big_table_theme <- ttheme_default(
    core = list(fg_params = list(cex = 1.5)),
    colhead = list(fg_params = list(cex = 1.5, fontface = "bold"))
  )

  # Titles with larger font and negative vjust to bring them closer to tables
  title_genes <- textGrob("Gene Classification", gp=gpar(fontface="italic", fontsize=24), vjust = -3)
  title_reads <- textGrob("Read Classification", gp=gpar(fontface="italic", fontsize=24), vjust = -7.7)
  title_sj <- textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=24), vjust = -4.3)

  # Table grobs with bigger font
  table_genes <- tableGrob(gene_class_table, rows = NULL, theme = big_table_theme)
  table_reads <- tableGrob(read_class_table, rows = NULL, theme = big_table_theme)
  table_sj <- tableGrob(SJ_class_table, rows = NULL, theme = big_table_theme)
  
  # Unique counts grob, bigger font
  unique_counts_grob <- textGrob(
    sprintf("Number of Reads: %d\nUnique Genes: %d\nUnique Junction Chains: %d", total_reads_count, unique_genes, unique_junctions),
    gp=gpar(fontface="italic", fontsize=28), vjust = 0, hjust = 0.5
  )

  # Create gTree objects to overlay titles and tables
  gt_genes <- gTree(children = gList(table_genes, title_genes))
  gt_reads <- gTree(children = gList(table_reads, title_reads))
  gt_sj <- gTree(children = gList(table_sj, title_sj))

  # Arrange left column: Gene Classification + Splice Junction Classification
  left_col <- arrangeGrob(
    gt_genes,
    gt_sj,
    ncol=1,
    heights=c(0.2, 0.4)
  )

  # Arrange right column: Read Classification
  right_col <- arrangeGrob(
    gt_reads,
    ncol=1
  )

  # Final page layout
  grid.arrange(
    unique_counts_grob,
    arrangeGrob(left_col, right_col, ncol=2, widths=c(1.3, 1.3)),
    nrow=2,
    heights=c(0.8, 1)
  )

  # Single cell tables
  s <- textGrob("Cell summary", gp=gpar(fontface="italic", fontsize=30), vjust = 0)
  grid.arrange(s)

  # Number of cells
  num_cells <- nrow(SQANTI_cell_summary)
  num_cells_grob <- textGrob(
    sprintf("Unique Cell Barcodes: %d", num_cells),
    gp=gpar(fontface="italic", fontsize=28), vjust=0.5, hjust=0.5
  )
  
  # 1. Unique Genes and Unique Junction Chains summary table
  unique_genes_stats <- c(
    Mean = mean(SQANTI_cell_summary$Genes_in_cell, na.rm=TRUE),
    Median = median(SQANTI_cell_summary$Genes_in_cell, na.rm=TRUE),
    Min = min(SQANTI_cell_summary$Genes_in_cell, na.rm=TRUE),
    Max = max(SQANTI_cell_summary$Genes_in_cell, na.rm=TRUE),
    SD = sd(SQANTI_cell_summary$Genes_in_cell, na.rm=TRUE)
  )
  unique_junctions_stats <- c(
    Mean = mean(SQANTI_cell_summary$UJCs_in_cell, na.rm=TRUE),
    Median = median(SQANTI_cell_summary$UJCs_in_cell, na.rm=TRUE),
    Min = min(SQANTI_cell_summary$UJCs_in_cell, na.rm=TRUE),
    Max = max(SQANTI_cell_summary$UJCs_in_cell, na.rm=TRUE),
    SD = sd(SQANTI_cell_summary$UJCs_in_cell, na.rm=TRUE)
  )
  reads_stats <- c(
    Mean = mean(SQANTI_cell_summary$Reads_in_cell, na.rm=TRUE),
    Median = median(SQANTI_cell_summary$Reads_in_cell, na.rm=TRUE),
    Min = min(SQANTI_cell_summary$Reads_in_cell, na.rm=TRUE),
    Max = max(SQANTI_cell_summary$Reads_in_cell, na.rm=TRUE),
    SD = sd(SQANTI_cell_summary$Reads_in_cell, na.rm=TRUE)
  )
  umis_stats <- c(
    Mean = mean(SQANTI_cell_summary$UMIs_in_cell, na.rm=TRUE),
    Median = median(SQANTI_cell_summary$UMIs_in_cell, na.rm=TRUE),
    Min = min(SQANTI_cell_summary$UMIs_in_cell, na.rm=TRUE),
    Max = max(SQANTI_cell_summary$UMIs_in_cell, na.rm=TRUE),
    SD = sd(SQANTI_cell_summary$UMIs_in_cell, na.rm=TRUE)
  )
  summary_table1 <- data.frame(
    Feature = c("Reads in cell", "UMIs in cell", "Unique Genes", "Unique Junction Chains"),
    Mean = c(reads_stats["Mean"], umis_stats["Mean"], unique_genes_stats["Mean"], unique_junctions_stats["Mean"]),
    Median = c(reads_stats["Median"], umis_stats["Median"], unique_genes_stats["Median"], unique_junctions_stats["Median"]),
    Min = c(reads_stats["Min"], umis_stats["Min"], unique_genes_stats["Min"], unique_junctions_stats["Min"]),
    Max = c(reads_stats["Max"], umis_stats["Max"], unique_genes_stats["Max"], unique_junctions_stats["Max"]),
    SD = c(reads_stats["SD"], umis_stats["SD"], unique_genes_stats["SD"], unique_junctions_stats["SD"])
  )
  summary_table1[, 2:6] <- round(summary_table1[, 2:6], 3)
  table_summary1 <- tableGrob(summary_table1, rows = NULL, theme = big_table_theme)
  gt_summary1 <- gTree(children = gList(table_summary1))

  # 2. Gene Classification summary table (across all cells)
  gene_class_stats <- data.frame(
    Category = c("Annotated Genes", "Novel Genes"),
    Mean = c(mean(SQANTI_cell_summary$Annotated_genes, na.rm=TRUE), mean(SQANTI_cell_summary$Novel_genes, na.rm=TRUE)),
    Median = c(median(SQANTI_cell_summary$Annotated_genes, na.rm=TRUE), median(SQANTI_cell_summary$Novel_genes, na.rm=TRUE)),
    Min = c(min(SQANTI_cell_summary$Annotated_genes, na.rm=TRUE), min(SQANTI_cell_summary$Novel_genes, na.rm=TRUE)),
    Max = c(max(SQANTI_cell_summary$Annotated_genes, na.rm=TRUE), max(SQANTI_cell_summary$Novel_genes, na.rm=TRUE)),
    SD = c(sd(SQANTI_cell_summary$Annotated_genes, na.rm=TRUE), sd(SQANTI_cell_summary$Novel_genes, na.rm=TRUE))
  )
  gene_class_stats[, 2:6] <- round(gene_class_stats[, 2:6], 3)
  table_gene_class_stats <- tableGrob(gene_class_stats, rows = NULL, theme = big_table_theme)
  title_gene_class_stats <- textGrob("Gene Classification (per cell)", gp=gpar(fontface="italic", fontsize=22), vjust = -2.9)
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
      .groups = 'drop'
    )
  
  # Calculate summary statistics across all cells
  sj_stats <- junction_proportions_per_cell %>%
    select(-CB) %>%
    summarise(
      across(
        everything(),
        list(
          Mean = ~mean(.x, na.rm = TRUE),
          Median = ~median(.x, na.rm = TRUE),
          Min = ~min(.x, na.rm = TRUE),
          Max = ~max(.x, na.rm = TRUE),
          SD = ~sd(.x, na.rm = TRUE)
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
    numeric_cols <- sapply(sj_stats_df[,-1], is.numeric)
    sj_stats_df[,-1][numeric_cols] <- round(sj_stats_df[,-1][numeric_cols], 3)
  }
  table_sj_stats <- tableGrob(sj_stats_df, rows = NULL, theme = big_table_theme)
  title_sj_stats <- textGrob("Splice Junction Classification (per cell, %)", gp=gpar(fontface="italic", fontsize=22), vjust = -4.4)
  gt_sj_stats <- gTree(children = gList(table_sj_stats, title_sj_stats))

  grid.arrange(
    num_cells_grob,
    gt_summary1,
    gt_gene_class_stats,
    gt_sj_stats,
    ncol=1,
    heights=c(0.3, 1, 0.7, 0.9)
  )

  #Cell Summary Statistics Page 2: Read Classification
  title_read_class <- textGrob("Read Classification", gp=gpar(fontface="italic", fontsize=28), vjust = 0, hjust = 0.5)
  desc_counts <- textGrob("Summary of per cell read counts by structural category", gp=gpar(fontface="italic", fontsize=18), vjust = 0.5)
  desc_props <- textGrob("Summary of per cell read percentages by structural category", gp=gpar(fontface="italic", fontsize=18), vjust = 0.5)
  struct_cat_cols <- c(
    "FSM", "ISM", "NIC", "NNC", "Genic_Genomic", "Antisense", "Fusion", "Intergenic", "Genic_intron"
  )
  struct_cat_names <- c(
    "FSM", "ISM", "NIC", "NNC", "Genic Genomic", "Antisense", "Fusion", "Intergenic", "Genic Intron"
  )

  # Smaller table theme for these two tables
  small_table_theme <- ttheme_default(
    core = list(fg_params = list(cex = 1.2)),
    colhead = list(fg_params = list(cex = 1.2, fontface = "bold"))
  )

  # 1. Counts summary table
  count_stats <- sapply(struct_cat_cols, function(col) {
    vals <- SQANTI_cell_summary[[col]]
    c(Mean = mean(vals, na.rm=TRUE),
      Median = median(vals, na.rm=TRUE),
      Min = min(vals, na.rm=TRUE),
      Max = max(vals, na.rm=TRUE),
      SD = sd(vals, na.rm=TRUE))
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
    c(Mean = mean(vals, na.rm=TRUE),
      Median = median(vals, na.rm=TRUE),
      Min = min(vals, na.rm=TRUE),
      Max = max(vals, na.rm=TRUE),
      SD = sd(vals, na.rm=TRUE))
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
    ncol=1,
    heights=c(0.3, 0.12, 1, 0.12, 1)
  )

  grid.arrange(gg_reads_in_cells, gg_umis_in_cells, ncol=2)
  grid.arrange(gg_genes_in_cells, gg_JCs_in_cell, ncol=2)
  print(gg_annotation_of_genes_in_cell)
  print(gg_read_bins)
  print(gg_ujc_bins)
  print(gg_MT_perc)
  ### Read lengths ###
  print(gg_bulk_all_reads)
  print(gg_bulk_length_by_category)
  print(gg_bulk_length_by_exon_type)
  print(gg_read_distr)
  print(gg_read_distr_mono)
  print(gg_FSM_read_distr)
  print(gg_FSM_mono_read_distr)
  print(gg_ISM_read_distr)
  print(gg_ISM_mono_read_distr)
  print(gg_NIC_read_distr)
  print(gg_NIC_mono_read_distr)
  print(gg_NNC_read_distr)
  print(gg_genic_read_distr)
  print(gg_genic_mono_read_distr)
  print(gg_antisense_read_distr)
  print(gg_antisense_mono_read_distr)
  print(gg_fusion_read_distr)
  print(gg_intergenic_read_distr)
  print(gg_intergenic_mono_read_distr)
  print(gg_genic_intron_read_distr)
  print(gg_genic_intron_mono_read_distr)
  ### SQANTI structural categories ###
  print(gg_SQANTI_across_category)
  print(gg_SQANTI_across_FSM)
  print(gg_SQANTI_across_ISM)
  print(gg_SQANTI_across_NIC)
  print(gg_SQANTI_across_NNC)
  print(gg_SQANTI_across_Genic)
  print(gg_SQANTI_across_Antisense)
  print(gg_SQANTI_across_Fusion)
  print(gg_SQANTI_across_Intergenic)
  print(gg_SQANTI_across_Genic_Intron)
  ### Coding/non-coding ###
  if (!skipORF) {
    # Create a temporary copy of SQANTI_cell_summary for plotting coding/non-coding
    # to avoid altering the main data table that might be saved.
    # Here, set proportions to NA if the category count is 0 for that cell,
    # so violin plots focus on cells with actual reads in that category.
    plot_data_temp_coding <- SQANTI_cell_summary

    prop_count_pairs_for_coding_ncoding <- list(
      c("Coding_FSM_prop", "FSM"), c("Non_coding_FSM_prop", "FSM"),
      c("Coding_ISM_prop", "ISM"), c("Non_coding_ISM_prop", "ISM"),
      c("Coding_NIC_prop", "NIC"), c("Non_coding_NIC_prop", "NIC"),
      c("Coding_NNC_prop", "NNC"), c("Non_coding_NNC_prop", "NNC"),
      c("Coding_genic_prop", "Genic_Genomic"), c("Non_coding_genic_prop", "Genic_Genomic"),
      c("Coding_antisense_prop", "Antisense"), c("Non_coding_antisense_prop", "Antisense"),
      c("Coding_fusion_prop", "Fusion"), c("Non_coding_fusion_prop", "Fusion"),
      c("Coding_intergenic_prop", "Intergenic"), c("Non_coding_intergenic_prop", "Intergenic"),
      c("Coding_genic_intron_prop", "Genic_intron"), c("Non_coding_genic_intron_prop", "Genic_intron")
    )

    for (pair in prop_count_pairs_for_coding_ncoding) {
      prop_col <- pair[1]
      count_col <- pair[2]
      if (prop_col %in% names(plot_data_temp_coding) && count_col %in% names(plot_data_temp_coding)) {
        plot_data_temp_coding[[prop_col]] <- ifelse(plot_data_temp_coding[[count_col]] == 0, NA_real_, plot_data_temp_coding[[prop_col]])
      }
    }

    gg_SQANTI_pivot_coding <- pivot_longer(plot_data_temp_coding, cols = c("Coding_FSM_prop",
                                                                  "Coding_ISM_prop",
                                                                  "Coding_NIC_prop",
                                                                  "Coding_NNC_prop",
                                                                  "Coding_genic_prop",
                                                                  "Coding_antisense_prop",
                                                                  "Coding_fusion_prop",
                                                                  "Coding_intergenic_prop",
                                                                  "Coding_genic_intron_prop"), 
                                    names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
    
    gg_SQANTI_pivot_coding$Variable <- factor(gg_SQANTI_pivot_coding$Variable, colnames(plot_data_temp_coding %>% select(Coding_FSM_prop,
                                                                                                         Coding_ISM_prop,
                                                                                                         Coding_NIC_prop,
                                                                                                         Coding_NNC_prop,
                                                                                                         Coding_genic_prop,
                                                                                                         Coding_antisense_prop,
                                                                                                         Coding_fusion_prop,
                                                                                                         Coding_intergenic_prop,
                                                                                                         Coding_genic_intron_prop)))
    
    gg_coding_across_category <- ggplot(gg_SQANTI_pivot_coding, aes(x = Variable, y = Value)) + # Use pivoted temp data
      geom_violin(aes(color = Variable, 
                      fill = Variable),
                      alpha = 0.7, scale = "width", na.rm = TRUE) +  
      geom_boxplot(aes(fill = Variable), color = c(rep("grey20",4),"grey90",rep("grey20",4)),
                   width = 0.08, outlier.shape = NA, alpha = 0.6, na.rm = TRUE) +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1, na.rm = TRUE) +
      scale_fill_manual(values = c("Coding_FSM_prop"="#6BAED6", "Coding_ISM_prop"="#FC8D59", "Coding_NIC_prop"="#78C679", "Coding_NNC_prop"="#EE6A50", 
                                   "Coding_genic_prop"="#969696", "Coding_antisense_prop"="#66C2A4", "Coding_fusion_prop"="goldenrod1",
                                   "Coding_intergenic_prop"="darksalmon", "Coding_genic_intron_prop"="#41B6C4")) + 
      scale_color_manual(values = c("Coding_FSM_prop"="#6BAED6", "Coding_ISM_prop"="#FC8D59", "Coding_NIC_prop"="#78C679", "Coding_NNC_prop"="#EE6A50", 
                                    "Coding_genic_prop"="#969696", "Coding_antisense_prop"="#66C2A4", "Coding_fusion_prop"="goldenrod1",
                                    "Coding_intergenic_prop"="darksalmon", "Coding_genic_intron_prop"="#41B6C4")) +
      scale_x_discrete(labels = c("FSM","ISM","NIC","NNC","Genic\nGenomic",
                                  "Antisense","Fusion","Intergenic","Genic\nintron")) +
      theme_classic(base_size = 14) +  
      labs(title = "Coding Proportion of Structural Categories Distribution Across Cells",
           x = "",
           y = "Reads, %") +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        axis.title = element_text(size = 16), 
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
      )
    print(gg_coding_across_category)
    
    gg_SQANTI_pivot_non_coding <- pivot_longer(plot_data_temp_coding, cols = c("Non_coding_FSM_prop",
                                                                  "Non_coding_ISM_prop",
                                                                  "Non_coding_NIC_prop",
                                                                  "Non_coding_NNC_prop",
                                                                  "Non_coding_genic_prop",
                                                                  "Non_coding_antisense_prop",
                                                                  "Non_coding_fusion_prop",
                                                                  "Non_coding_intergenic_prop",
                                                                  "Non_coding_genic_intron_prop"), 
                                    names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
    
    gg_SQANTI_pivot_non_coding$Variable <- factor(gg_SQANTI_pivot_non_coding$Variable, colnames(plot_data_temp_coding %>% select(Non_coding_FSM_prop,
                                                                                                         Non_coding_ISM_prop,
                                                                                                         Non_coding_NIC_prop,
                                                                                                         Non_coding_NNC_prop,
                                                                                                         Non_coding_genic_prop,
                                                                                                         Non_coding_antisense_prop,
                                                                                                         Non_coding_fusion_prop,
                                                                                                         Non_coding_intergenic_prop,
                                                                                                         Non_coding_genic_intron_prop)))
    
    gg_non_coding_across_category <- ggplot(gg_SQANTI_pivot_non_coding, aes(x = Variable, y = Value)) + # Use pivoted temp data
      geom_violin(aes(color = Variable, 
                      fill = Variable),
                      alpha = 0.7, scale = "width", na.rm = TRUE) +  
      geom_boxplot(aes(fill = Variable), color = "grey20",
                   width = 0.08, outlier.shape = NA, alpha = 0.6, na.rm = TRUE) +
      stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1, na.rm = TRUE) +
      scale_fill_manual(values = c("Non_coding_FSM_prop"="#D2E6F2", "Non_coding_ISM_prop"="#FEDCCD", "Non_coding_NIC_prop"="#D6EDD6", "Non_coding_NNC_prop"="#F9D2CA", 
                                   "Non_coding_genic_prop"="#DFDFDF", "Non_coding_antisense_prop"="#D1ECE3", "Non_coding_fusion_prop"="#FFECBD",
                                   "Non_coding_intergenic_prop"="#F8DFD7", "Non_coding_genic_intron_prop"="#C6E9ED")) + 
      scale_color_manual(values = c("Non_coding_FSM_prop"="#D2E6F2", "Non_coding_ISM_prop"="#FEDCCD", "Non_coding_NIC_prop"="#D6EDD6", "Non_coding_NNC_prop"="#F9D2CA", 
                                    "Non_coding_genic_prop"="#DFDFDF", "Non_coding_antisense_prop"="#D1ECE3", "Non_coding_fusion_prop"="#FFECBD",
                                    "Non_coding_intergenic_prop"="#F8DFD7", "Non_coding_genic_intron_prop"="#C6E9ED")) +
      scale_x_discrete(labels = c("FSM","ISM","NIC","NNC","Genic\nGenomic",
                                  "Antisense","Fusion","Intergenic","Genic\nintron")) +
      theme_classic(base_size = 14) +  
      labs(title = "Non-coding Proportion of Structural Categories Distribution Across Cells",
           x = "",
           y = "Reads, %") +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        axis.title = element_text(size = 16), 
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16) 
      )
    print(gg_non_coding_across_category)
  }
  ### Coverage (TSS/TTS in the future) ###
  print(gg_ref_coverage_across_category)
  ### Splice Junctions Categories ###
  print(gg_known_novel_canon)
  ### Bad features ###
  print(gg_bad_feature)
  ### Bad features by structural category ###
  print(gg_intrapriming_by_category)
  print(gg_RTS_by_category)
  print(gg_noncanon_by_category)
  if (!skipORF) {
    if (exists("gg_NMD_by_category") && !is.null(gg_NMD_by_category)){
      print(gg_NMD_by_category)
    }
  }
  ### Good features ###
  print(gg_good_feature)
  ### Good features by structural category ###
  print(gg_tss_annotation_support)
  if (CAGE_peak) {
    print(gg_cage_peak_support)
  }
  if (polyA_motif_list) {
    print(gg_polyA_motif_support)
  }
  print(gg_canon_by_category)

  dev.off()
}

Classification <- read.table(class.file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
Junctions <- read.table(junc.file, header=TRUE, sep="\t", stringsAsFactors=FALSE) 
SQANTI_cell_summary <- calculate_metrics_per_cell(
  Classification, 
  Junctions,
  cell_summary_output, 
  Save = save_option
)

generate_sqantisc_plots(
  SQANTI_cell_summary, 
  Classification, 
  Junctions,
  report_output
)