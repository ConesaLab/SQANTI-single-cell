
# Load dependencies
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(grid)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1]
report.format <- args[2]
outputPathPrefix <- args[3]

# Initialize ignore_cell_summary flag
ignore_cell_summary <- FALSE

# Check for optional arguments
if (length(args) > 3) {
  for (arg in args[4:length(args)]) {
    if (arg == "--ignore_cell_summary") {
      ignore_cell_summary <- TRUE
    }
  }
}

# Validate arguments
if (length(args) < 3) {
  stop("Incorrect number of arguments! Required: [classification file] [report format] [outputPathPrefix]. Abort!")
}

if (!(report.format %in% c('pdf', 'html', 'both'))) {
  stop("Report format needs to be: pdf, html, or both. Abort!")
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

calculate_metrics_per_cell <- function(Classification, cell_summary_output, Save){
  CBs <- Classification %>% select(CB) %>% distinct() %>% .$CB
  maxCB <- length(CBs)
  for (CB_id in CBs){
    if (CB_id==""){
      print("There are reads with no cell barcode assigned and they will not be considered. Check your classification file.")
      next
    }
    print(paste0("Calculating metrics for CB ",CB_id,". (",which(CBs==CB_id),"/",maxCB,")"))
    sorted_classification <- Classification %>% filter(CB==CB_id)
  
    ##### READS ########
    # Reads per cell 
    total_reads <- nrow(sorted_classification)
    
    total_UMI <- sorted_classification %>% select(UMI) %>% n_distinct()
    
    total_reads_no_monoexon <- sorted_classification %>%
                               filter(!exons==1) %>%
                               nrow()
    
    FSM_count <- sorted_classification %>%
                 filter(structural_category=="full-splice_match") %>%
                 nrow()
    
    ISM_count <- sorted_classification %>%
                 filter(structural_category=="incomplete-splice_match") %>%
                 nrow()
    
    NIC_count <- sorted_classification %>%
                 filter(structural_category=="novel_in_catalog") %>%
                 nrow()
    
    NNC_count <- sorted_classification %>%
                 filter(structural_category=="novel_not_in_catalog") %>%
                 nrow()
    
    Genic_count <- sorted_classification %>%
                   filter(structural_category=="genic") %>%
                   nrow()
    
    Antisense_count <- sorted_classification %>%
                       filter(structural_category=="antisense") %>%
                       nrow()
    
    Fusion_count <- sorted_classification %>%
                    filter(structural_category=="fusion") %>%
                    nrow()
    
    Intergenic_count <- sorted_classification %>%
                        filter(structural_category=="intergenic") %>%
                        nrow()
    
    Genic_intron_count <- sorted_classification %>%
                          filter(structural_category=="genic_intron") %>%
                          nrow()
    
    # Genes per cell
    genes_in_cell <- sorted_classification %>%
                     select(associated_gene) %>%
                     n_distinct()
    
    models_in_cell <- sorted_classification %>%
                      group_by(associated_gene) %>%
                      summarise(t_chains=n_distinct(jxn_string)) %>%
                      .$t_chains %>% sum()
    # Take into consideration that monoexons of the same gene could not be counted correctly. Need to update jxn_string to include UTRs to distinguish them
    
    # Mitochondrial reads
    MT_reads_count <- sorted_classification %>%
                      filter(chrom == "MT") %>%
                      nrow()

    MT_perc <- (MT_reads_count / total_reads) * 100
            
    # Novel vs annotated metrics
    annotated_genes <- sorted_classification %>%
                       select(associated_gene) %>%
                       filter(!grepl("^novel", associated_gene)) %>%
                       n_distinct()
    
    novel_genes <- sorted_classification %>%
                   select(associated_gene) %>%
                   filter(grepl("^novel", associated_gene)) %>%
                   n_distinct()
    
    # Known/novel canonical/non-canonical
    if (total_reads_no_monoexon==0){
      known_canonical_prop <- 0
      known_non_canonical_prop <- 0
      novel_canonical_prop <- 0
      novel_non_canonical_prop <- 0
    } else {
    known_canonical_prop <- sorted_classification %>%
                            filter(!grepl("^novel", associated_transcript) & all_canonical=="canonical") %>%
                            nrow()/total_reads_no_monoexon*100
    
    known_non_canonical_prop <- sorted_classification %>%
                                filter(!grepl("^novel", associated_transcript) & all_canonical=="non_canonical") %>%
                                nrow()/total_reads_no_monoexon*100
    
    novel_canonical_prop <- sorted_classification %>%
                            filter(grepl("^novel", associated_transcript) & all_canonical=="canonical") %>%
                            nrow()/total_reads_no_monoexon*100
    
    novel_non_canonical_prop <- sorted_classification %>%
                                filter(grepl("^novel", associated_transcript) & all_canonical=="non_canonical") %>%
                                nrow()/total_reads_no_monoexon*100
    }
    
    # Sqanti category
    sqanti_props <- sorted_classification %>%
                    count(structural_category, name = "reads_in_cat") %>%
                    mutate(prop_cat=reads_in_cat/total_reads*100) %>%
                    complete(structural_category = factor(structural_category,
                                                   levels=c('full-splice_match',
                                                            'incomplete-splice_match',
                                                            'novel_in_catalog',
                                                            'novel_not_in_catalog',
                                                            'genic',
                                                            'antisense',
                                                            'fusion',
                                                            'intergenic',
                                                            'genic_intron')),  # Ensure all categories exist
                             fill = list(prop_cat = 0)) %>% # Assign 0 to missing cat if not present
                    arrange(fct_relevel(structural_category,
                                        'full-splice_match',
                                        'incomplete-splice_match',
                                        'novel_in_catalog',
                                        'novel_not_in_catalog',
                                        'genic',
                                        'antisense',
                                        'fusion',
                                        'intergenic',
                                        'genic_intron')) %>% .$prop_cat
    
    ### COVERAGE OF GENE BODY ###
    # Ref length corresponds to transcript length, not gene. Maybe can be added.
    # Fusion genes can be excluded?
    ref_body_cover_in_cell <- sorted_classification %>%
                              filter((length/ref_length*100)>=45) %>%
                              nrow()/total_reads*100
    
    # Read lengths general
    two_fifty_length_reads <- sorted_classification %>%
                              filter(length<=250) %>%
                              nrow()/total_reads*100
    
    five_hund_length_reads <- sorted_classification %>%
                              filter(length>250 & length<=500) %>%
                              nrow()/total_reads*100
    
    short_length_reads <- sorted_classification %>%
                          filter(length>500 & length<=1000) %>%
                          nrow()/total_reads*100
    
    mid_length_reads <- sorted_classification %>%
                        filter(length>1000 & length<=2000) %>%
                        nrow()/total_reads*100
    
    long_length_reads <- sorted_classification %>%
                         filter(length>2000) %>%
                         nrow()/total_reads*100
     
    # Monoexons read length general
    mono_two_fifty_length_reads <- sorted_classification %>%
                                   filter(exons==1 & length<=250) %>%
                                   nrow()/total_reads*100
    
    mono_five_hund_length_reads <- sorted_classification %>%
                                   filter(exons==1 & length>250 & length<=500) %>%
                                   nrow()/total_reads*100
    
    mono_short_length_reads <- sorted_classification %>%
                               filter(exons==1 & length>500 & length<=1000) %>%
                               nrow()/total_reads*100
    
    mono_mid_length_reads <- sorted_classification %>%
                             filter(exons==1 & length>1000 & length<=2000) %>%
                             nrow()/total_reads*100
    
    mono_long_length_reads <- sorted_classification %>%
                              filter(exons==1 & length>2000) %>%
                              nrow()/total_reads*100
    
    ### Sqanti subcategories ###
    if (FSM_count==0){
      sub_FSM_sqanti_props <- c(Alternative_3end=0,
                                Alternative_3end5end=0,
                                Alternative_5end=0,
                                Reference_match=0,
                                Mono_exon_FSM=0) %>% as.numeric()
      
      two_fifty_length_reads_FSM <- 0
      five_hund_length_reads_FSM <- 0
      short_length_reads_FSM <- 0
      mid_length_reads_FSM <- 0
      long_length_reads_FSM <- 0
      mono_two_fifty_length_reads_FSM <- 0
      mono_five_hund_length_reads_FSM <- 0
      mono_short_length_reads_FSM <- 0
      mono_mid_length_reads_FSM <- 0
      mono_long_length_reads_FSM <- 0
      ref_body_cover_FSM <- 0
      cod_FSM <- 0
      ncod_FSM <- 100 #    
    } else {
      sub_FSM_sqanti_props <- sorted_classification %>% 
                              filter(structural_category=="full-splice_match") %>%
                              count(subcategory, name = "reads_in_subcat") %>%
                              mutate(prop_subcat=reads_in_subcat/FSM_count*100) %>%
                              complete(subcategory = factor(subcategory,
                                                            levels = c("alternative_3end",
                                                                       "alternative_3end5end",
                                                                       "alternative_5end",
                                                                       "reference_match",
                                                                       "mono-exon")), # Ensure all subcategories exist
                                       fill = list(prop_subcat = 0)) %>% # Assign 0 to missing subcat if not present
                              arrange(fct_relevel(subcategory,
                                                  "alternative_3end",
                                                  "alternative_3end5end",
                                                  "alternative_5end",
                                                  "reference_match",
                                                  "mono-exon")) %>% .$prop_subcat
      
      # Read length per category
      two_fifty_length_reads_FSM <- sorted_classification %>%
                                    filter(structural_category=="full-splice_match" & length<=250) %>%
                                    nrow()/FSM_count*100
     
      five_hund_length_reads_FSM <- sorted_classification %>%
                                     filter(structural_category=="full-splice_match" & length>250 & length<=500) %>%
                                     nrow()/FSM_count*100
      
      short_length_reads_FSM <- sorted_classification %>%
                                filter(structural_category=="full-splice_match" & length>500 & length<=1000) %>%
                                nrow()/FSM_count*100
      
      mid_length_reads_FSM <- sorted_classification %>%
                              filter(structural_category=="full-splice_match" & length>1000 & length<=2000) %>%
                              nrow()/FSM_count*100
      
      long_length_reads_FSM <- sorted_classification %>%
                               filter(structural_category=="full-splice_match" & length>2000) %>%
                               nrow()/FSM_count*100
      
      # Monoexon per length break
      mono_two_fifty_length_reads_FSM <- sorted_classification %>%
                                         filter(structural_category=="full-splice_match" & exons==1 & length<=250) %>%
                                         nrow()/FSM_count*100
      
      mono_five_hund_length_reads_FSM <- sorted_classification %>%
                                         filter(structural_category=="full-splice_match" & exons==1 & length>250 & length<=500) %>%
                                         nrow()/FSM_count*100
      
      mono_short_length_reads_FSM <- sorted_classification %>%
                                     filter(structural_category=="full-splice_match" & exons==1 & length>500 & length<=1000) %>%
                                     nrow()/FSM_count*100
      
      mono_mid_length_reads_FSM <- sorted_classification %>%
                                   filter(structural_category=="full-splice_match" & exons==1 & length>1000 & length<=2000) %>%
                                   nrow()/FSM_count*100
      
      mono_long_length_reads_FSM <- sorted_classification %>%
                                    filter(structural_category=="full-splice_match" & exons==1 & length>2000) %>%
                                    nrow()/FSM_count*100
      
      # Reference body coverage
      ref_body_cover_FSM <- sorted_classification %>%
                            filter(structural_category=="full-splice_match" & length/ref_length*100>=45) %>%
                            nrow()/FSM_count*100
      
      # Coding/non-coding
      cod_FSM <- sorted_classification %>%
                 filter(structural_category=="full-splice_match" & coding=="coding") %>%
                 nrow()/FSM_count*100
      
      ncod_FSM <- sorted_classification %>%
                  filter(structural_category=="full-splice_match" & coding=="non_coding") %>%
                  nrow()/FSM_count*100
    }
    
    if (ISM_count==0){
      sub_ISM_sqanti_props <- c(Fragment_3prime=0,
                                Internal_fragment=0,
                                Fragment_5prime=0,
                                Intron_retention_ISM=0,
                                Mono_exon_ISM=0) %>% as.numeric()
      
      two_fifty_length_reads_ISM <-0
      five_hund_length_reads_ISM <- 0
      short_length_reads_ISM <- 0
      mid_length_reads_ISM <- 0
      long_length_reads_ISM <- 0
      mono_two_fifty_length_reads_ISM <-0
      mono_five_hund_length_reads_ISM <- 0
      mono_short_length_reads_ISM <- 0
      mono_mid_length_reads_ISM <- 0
      mono_long_length_reads_ISM <- 0
      ref_body_cover_ISM <- 0
      cod_ISM <- 0
      ncod_ISM <- 100 #     
    } else {
      sub_ISM_sqanti_props <- sorted_classification %>%
                              filter(structural_category=="incomplete-splice_match") %>%
                              count(subcategory, name = "reads_in_subcat") %>%
                              mutate(prop_subcat=reads_in_subcat/ISM_count*100) %>%
                              complete(subcategory = factor(subcategory,
                                                            levels = c("3prime_fragment",
                                                                       "internal_fragment",
                                                                       "5prime_fragment",
                                                                       "intron_retention",
                                                                       "mono-exon")),  
                                      fill = list(prop_subcat = 0)) %>% 
                              arrange(fct_relevel(subcategory,
                                                  "3prime_fragment",
                                                  "internal_fragment",
                                                  "5prime_fragment",
                                                  "intron_retention",
                                                  "mono-exon")) %>% .$prop_subcat
      
      # Read length per category
      two_fifty_length_reads_ISM <- sorted_classification %>%
                                    filter(structural_category=="incomplete-splice_match") %>%
                                    filter(length<=250) %>%
                                    nrow()/ISM_count*100
      
      five_hund_length_reads_ISM <- sorted_classification %>%
                                    filter(structural_category=="incomplete-splice_match") %>%
                                    filter(length>250 & length<=500) %>%
                                    nrow()/ISM_count*100
      
      short_length_reads_ISM <- sorted_classification %>%
                                filter(structural_category=="incomplete-splice_match") %>%
                                filter(length>500 & length<=1000) %>%
                                nrow()/ISM_count*100
      
      mid_length_reads_ISM <- sorted_classification %>%
                              filter(structural_category=="incomplete-splice_match") %>%
                              filter(length>1000 & length<=2000) %>%
                              nrow()/ISM_count*100
      
      long_length_reads_ISM <- sorted_classification %>%
                               filter(structural_category=="incomplete-splice_match") %>%
                               filter(length>2000) %>%
                               nrow()/ISM_count*100
      
      # Monoexons
      mono_two_fifty_length_reads_ISM <- sorted_classification %>%
                                         filter(structural_category=="incomplete-splice_match" & exons==1 & length<=250) %>%
                                         nrow()/ISM_count*100
      
      mono_five_hund_length_reads_ISM <- sorted_classification %>%
                                         filter(structural_category=="incomplete-splice_match" & exons==1 & length>250 & length<=500) %>%
                                         nrow()/ISM_count*100
      
      mono_short_length_reads_ISM <- sorted_classification %>%
                                     filter(structural_category=="incomplete-splice_match" & exons==1 & length>500 & length<=1000) %>%
                                     nrow()/ISM_count*100
      
      mono_mid_length_reads_ISM <- sorted_classification %>%
                                   filter(structural_category=="incomplete-splice_match" & exons==1 & length>1000 & length<=2000) %>%
                                   nrow()/ISM_count*100
      
      mono_long_length_reads_ISM <- sorted_classification %>%
                                    filter(structural_category=="incomplete-splice_match" & exons==1 & length>2000) %>%
                                    nrow()/ISM_count*100
      
      # Reference body coverage
      ref_body_cover_ISM <- sorted_classification %>%
                            filter(structural_category=="incomplete-splice_match" & length/ref_length*100>=45) %>%
                            nrow()/ISM_count*100
      
      # Coding/non_coding
      cod_ISM <- sorted_classification %>%
                 filter(structural_category=="incomplete-splice_match" & coding=="coding") %>%
                 nrow()/ISM_count*100
      
      ncod_ISM <- sorted_classification %>%
                  filter(structural_category=="incomplete-splice_match" & coding=="non_coding") %>%
                  nrow()/ISM_count*100
    }
    
    if (NIC_count==0){
      sub_NIC_sqanti_probs <- c(Combination_of_known_junctions=0,
                                Combination_of_known_splicesites=0,
                                Intron_retention_NIC=0,
                                Mono_exon_by_intron_retention=0,
                                Mono_exon_NIC=0) %>% as.numeric()
      
      two_fifty_length_reads_NIC <- 0
      five_hund_length_reads_NIC <- 0
      short_length_reads_NIC <- 0
      mid_length_reads_NIC <- 0
      long_length_reads_NIC <- 0
      mono_two_fifty_length_reads_NIC <- 0
      mono_five_hund_length_reads_NIC <- 0
      mono_short_length_reads_NIC <- 0
      mono_mid_length_reads_NIC <- 0
      mono_long_length_reads_NIC <- 0
      ref_body_cover_NIC <- 0
      cod_NIC <- 0
      ncod_NIC <- 100 #    
    } else {
      sub_NIC_sqanti_probs <- sorted_classification %>%
                              filter(structural_category=="novel_in_catalog") %>%
                              count(subcategory, name = "reads_in_subcat") %>%
                              mutate(prop_subcat=reads_in_subcat/NIC_count*100) %>%
                              complete(subcategory = factor(subcategory,
                                                            levels = c("combination_of_known_junctions",
                                                                       "combination_of_known_splicesites",
                                                                       "intron_retention",
                                                                       "mono-exon_by_intron_retention",
                                                                       "mono-exon")),  
                                       fill = list(prop_subcat = 0)) %>% 
                              arrange(fct_relevel(subcategory,
                                                  "combination_of_known_junctions",
                                                  "combination_of_known_splicesites",
                                                  "intron_retention",
                                                  "mono-exon_by_intron_retention",
                                                  "mono-exon")) %>% .$prop_subcat
      
      # Read length per category
      two_fifty_length_reads_NIC <- sorted_classification %>%
                                    filter(structural_category=="novel_in_catalog" & length<=250) %>%
                                    nrow()/NIC_count*100
      
      five_hund_length_reads_NIC <- sorted_classification %>%
                                    filter(structural_category=="novel_in_catalog" & length>250 & length<=500) %>%
                                    nrow()/NIC_count*100
      
      short_length_reads_NIC <- sorted_classification %>%
                                filter(structural_category=="novel_in_catalog" & length>500 & length<=1000) %>%
                                nrow()/NIC_count*100
      
      mid_length_reads_NIC <- sorted_classification %>%
                              filter(structural_category=="novel_in_catalog" & length>1000 & length<=2000 ) %>%
                              nrow()/NIC_count*100
      
      long_length_reads_NIC <- sorted_classification %>%
                               filter(structural_category=="novel_in_catalog" & length>2000) %>%
                               nrow()/NIC_count*100
      
      # Monoexons
      mono_two_fifty_length_reads_NIC <- sorted_classification %>%
                                         filter(structural_category=="novel_in_catalog" & exons==1 & length<=250) %>%
                                         nrow()/NIC_count*100
      
      mono_five_hund_length_reads_NIC <- sorted_classification %>%
                                         filter(structural_category=="novel_in_catalog" & exons==1 & length>250 & length<=500) %>%
                                         nrow()/NIC_count*100
      
      mono_short_length_reads_NIC <- sorted_classification %>%
                                     filter(structural_category=="novel_in_catalog" & exons==1 & length>500 & length<=1000) %>%
                                     nrow()/NIC_count*100
      
      mono_mid_length_reads_NIC <- sorted_classification %>%
                                   filter(structural_category=="novel_in_catalog" & exons==1 & length>1000 & length<=2000) %>%
                                   nrow()/NIC_count*100
      
      mono_long_length_reads_NIC <- sorted_classification %>%
                                    filter(structural_category=="novel_in_catalog" & exons==1 & length>2000) %>%
                                    nrow()/NIC_count*100
      
      # Reference body coverage
      ref_body_cover_NIC <- sorted_classification %>%
                            filter(structural_category=="novel_in_catalog" & length/ref_length*100>=45) %>%
                            nrow()/NIC_count*100
      
      # Coding/non-coding
      cod_NIC <- sorted_classification %>%
                 filter(structural_category=="novel_in_catalog" & coding=="coding") %>%
                 nrow()/NIC_count*100
      
      ncod_NIC <- sorted_classification %>%
                  filter(structural_category=="novel_in_catalog" & coding=="non_coding") %>%
                  nrow()/NIC_count*100
    }
    
    if (NNC_count==0){
      sub_NNC_sqanti_probs <- c(At_least_one_novel_splicesite=0,
                                Intron_retention_NNC=0) %>% as.numeric()
      
      two_fifty_length_reads_NNC <- 0
      five_hund_length_reads_NNC <- 0
      short_length_reads_NNC <- 0
      mid_length_reads_NNC <- 0
      long_length_reads_NNC <- 0
      mono_two_fifty_length_reads_NNC <- 0
      mono_five_hund_length_reads_NNC <- 0
      mono_short_length_reads_NNC <- 0
      mono_mid_length_reads_NNC <- 0
      mono_long_length_reads_NNC <- 0
      ref_body_cover_NNC <- 0
      cod_NNC <- 0
      ncod_NNC <- 100
    } else {
      sub_NNC_sqanti_probs <- sorted_classification %>% 
                              filter(structural_category=="novel_not_in_catalog") %>%
                              count(subcategory, name = "reads_in_subcat") %>%
                              mutate(prop_subcat=reads_in_subcat/NNC_count*100) %>%
                              complete(subcategory = factor(subcategory,
                                                            levels = c("at_least_one_novel_splicesite",
                                                                       "intron_retention")),  
                                       fill = list(prop_subcat = 0)) %>% 
                              arrange(fct_relevel(subcategory,
                                                  "at_least_one_novel_splicesite",
                                                  "intron_retention")) %>% .$prop_subcat
      
      # Read length per category
      two_fifty_length_reads_NNC <- sorted_classification %>%
                                    filter(structural_category=="novel_not_in_catalog" & length<=250) %>%
                                    nrow()/NNC_count*100
      
      five_hund_length_reads_NNC <- sorted_classification %>%
                                    filter(structural_category=="novel_not_in_catalog" & length>250 & length<=500) %>%
                                    nrow()/NNC_count*100
      
      short_length_reads_NNC <- sorted_classification %>%
                                filter(structural_category=="novel_not_in_catalog" & length>500 & length<=1000) %>%
                                nrow()/NNC_count*100
      
      mid_length_reads_NNC <- sorted_classification %>%
                              filter(structural_category=="novel_not_in_catalog" & length>1000 & length<=2000) %>%
                              nrow()/NNC_count*100
      
      long_length_reads_NNC <- sorted_classification %>%
                               filter(structural_category=="novel_not_in_catalog" & length>2000) %>%
                               nrow()/NNC_count*100
      
      # Monoexons
      mono_two_fifty_length_reads_NNC <- sorted_classification %>%
                                         filter(structural_category=="novel_not_in_catalog" & exons==1 & length<=250) %>%
                                         nrow()/NNC_count*100
      
      mono_five_hund_length_reads_NNC <- sorted_classification %>%
                                         filter(structural_category=="novel_not_in_catalog" & exons==1 & length>250 & length<=500) %>%
                                         nrow()/NNC_count*100
      
      mono_short_length_reads_NNC <- sorted_classification %>%
                                     filter(structural_category=="novel_not_in_catalog" & exons==1 & length>500 & length<=1000) %>%
                                     nrow()/NNC_count*100
      
      mono_mid_length_reads_NNC <- sorted_classification %>%
                                   filter(structural_category=="novel_not_in_catalog" & exons==1 & length>1000 & length<=2000) %>%
                                   nrow()/NNC_count*100
      
      mono_long_length_reads_NNC <- sorted_classification %>%
                                    filter(structural_category=="novel_not_in_catalog" & exons==1 & length>2000) %>%
                                    nrow()/NNC_count*100
      
      # Reference body coverage
      ref_body_cover_NNC <- sorted_classification %>%
                            filter(structural_category=="novel_not_in_catalog" & length/ref_length*100>=45) %>%
                            nrow()/NNC_count*100
      
      # Coding/non-coding
      cod_NNC <- sorted_classification %>%
                 filter(structural_category=="novel_not_in_catalog" & coding=="coding") %>%
                 nrow()/NNC_count*100
      
      ncod_NNC <- sorted_classification %>%
                  filter(structural_category=="novel_not_in_catalog" & coding=="non_coding") %>%
                  nrow()/NNC_count*100
    }
    
    if (Fusion_count==0){
      sub_fusion_sqanti_props <- c(Intron_retention_fusion=0,
                                   Multi_exon_fusion=0) %>% as.numeric()
      
      two_fifty_length_reads_fusion <- 0
      five_hund_length_reads_fusion <- 0
      short_length_reads_fusion <- 0
      mid_length_reads_fusion <- 0 
      long_length_reads_fusion <- 0
      mono_two_fifty_length_reads_fusion <- 0
      mono_five_hund_length_reads_fusion <- 0
      mono_short_length_reads_fusion <- 0
      mono_mid_length_reads_fusion <- 0 
      mono_long_length_reads_fusion <- 0
      ref_body_cover_fusion <- 0
      cod_fusion <- 0
      ncod_fusion <- 100 #    
    } else {
      sub_fusion_sqanti_props <- sorted_classification %>%
                                 filter(structural_category=="fusion") %>%
                                 count(subcategory, name = "reads_in_subcat") %>%
                                 mutate(prop_subcat=reads_in_subcat/Fusion_count*100) %>%
                                 complete(subcategory = factor(subcategory,
                                                               levels = c("intron_retention",
                                                                          "multi-exon")),  
                                          fill = list(prop_subcat = 0)) %>%
                                 arrange(fct_relevel(subcategory,
                                                     "intron_retention",
                                                     "multi-exon")) %>% .$prop_subcat
      
      # Read length per category
      two_fifty_length_reads_fusion <- sorted_classification %>%
                                       filter(structural_category=="fusion" & length<=250) %>%
                                       nrow()/Fusion_count*100
      
      five_hund_length_reads_fusion <- sorted_classification %>%
                                       filter(structural_category=="fusion" & length>250 & length<=500) %>%
                                       nrow()/Fusion_count*100
      
      short_length_reads_fusion <- sorted_classification %>%
                                   filter(structural_category=="fusion" & length>500 & length<=1000) %>%
                                   nrow()/Fusion_count*100
      
      mid_length_reads_fusion <- sorted_classification %>%
                                 filter(structural_category=="fusion" & length>1000 & length<=2000) %>%
                                 nrow()/Fusion_count*100
      
      long_length_reads_fusion <- sorted_classification %>%
                                  filter(structural_category=="fusion" & length>2000) %>%
                                  nrow()/Fusion_count*100
      
      # Monoexon
      mono_two_fifty_length_reads_fusion <- sorted_classification %>%
                                            filter(structural_category=="fusion" & exons==1 & length<=250) %>%
                                            nrow()/Fusion_count*100
      
      mono_five_hund_length_reads_fusion <- sorted_classification %>%
                                            filter(structural_category=="fusion" & exons==1 & length>250 & length<=500) %>%
                                            nrow()/Fusion_count*100
      
      mono_short_length_reads_fusion <- sorted_classification %>%
                                        filter(structural_category=="fusion" & exons==1 & length>500 & length<=1000) %>%
                                        nrow()/Fusion_count*100
      
      mono_mid_length_reads_fusion <- sorted_classification %>%
                                      filter(structural_category=="fusion" & exons==1 & length>1000 & length<=2000) %>%
                                      nrow()/Fusion_count*100
      
      mono_long_length_reads_fusion <- sorted_classification %>%
                                       filter(structural_category=="fusion" & exons==1 & length>2000) %>%
                                       nrow()/Fusion_count*100
      
      # Reference body coverage
      ref_body_cover_fusion <- sorted_classification %>%
                               filter(structural_category=="fusion" & length/ref_length*100>=45) %>% ####  Take a decision about this
                               nrow()/Fusion_count*100
      
      # Coding/non-coding
      cod_fusion <- sorted_classification %>%
                    filter(structural_category=="fusion" & coding=="coding") %>%
                    nrow()/Fusion_count*100
      
      ncod_fusion <- sorted_classification %>%
                     filter(structural_category=="fusion" & coding=="non_coding") %>%
                     nrow()/Fusion_count*100
    }
    
    if (Genic_count==0){
      sub_genic_sqanti_props <- c(Mono_exon_genic=0,
                                  Multi_exon_genic=0) %>% as.numeric()
      
      two_fifty_length_reads_genic <- 0
      five_hund_length_reads_genic <- 0
      short_length_reads_genic <- 0
      mid_length_reads_genic <- 0
      long_length_reads_genic <- 0
      mono_two_fifty_length_reads_genic <- 0
      mono_five_hund_length_reads_genic <- 0
      mono_short_length_reads_genic <- 0
      mono_mid_length_reads_genic <- 0
      mono_long_length_reads_genic <- 0
      ref_body_cover_genic <- 0 
      cod_genic <- 0
      ncod_genic <- 100
    } else {
      sub_genic_sqanti_props <- sorted_classification %>%
                                filter(structural_category=="genic") %>%
                                count(subcategory, name = "reads_in_subcat") %>%
                                mutate(prop_subcat=reads_in_subcat/Genic_count*100) %>%
                                complete(subcategory = factor(subcategory,
                                                              levels = c("mono-exon",
                                                                         "multi-exon")),  
                                         fill = list(prop_subcat = 0)) %>%
                                arrange(fct_relevel(subcategory,
                                                    "mono-exon",
                                                    "multi-exon")) %>% .$prop_subcat
      
      # Read length per category
      two_fifty_length_reads_genic <- sorted_classification %>%
                                      filter(structural_category=="genic" & length<=250) %>%
                                      nrow()/Genic_count*100
      
      five_hund_length_reads_genic <- sorted_classification %>%
                                      filter(structural_category=="genic" & length>250 & length<=500) %>%
                                      nrow()/Genic_count*100
      
      short_length_reads_genic <- sorted_classification %>%
                                  filter(structural_category=="genic" & length>500 & length<=1000) %>%
                                  nrow()/Genic_count*100
      
      mid_length_reads_genic <- sorted_classification %>%
                                filter(structural_category=="genic" & length>1000 & length<=2000) %>%
                                nrow()/Genic_count*100
      
      long_length_reads_genic <- sorted_classification %>%
                                 filter(structural_category=="genic" & length>2000) %>%
                                 nrow()/Genic_count*100
      
      # Monoexons
      mono_two_fifty_length_reads_genic <- sorted_classification %>%
                                           filter(structural_category=="genic" & exons==1 & length<=250) %>%
                                           nrow()/Genic_count*100
       
      mono_five_hund_length_reads_genic <- sorted_classification %>%
                                           filter(structural_category=="genic" & exons==1 & length>250 & length<=500) %>%
                                           nrow()/Genic_count*100
      
      mono_short_length_reads_genic <- sorted_classification %>%
                                       filter(structural_category=="genic" & exons==1 & length>500 & length<=1000) %>%
                                       nrow()/Genic_count*100
      
      mono_mid_length_reads_genic <- sorted_classification %>%
                                     filter(structural_category=="genic" & exons==1 & length>1000 & length<=2000) %>%
                                     nrow()/Genic_count*100
      
      mono_long_length_reads_genic <- sorted_classification %>%
                                      filter(structural_category=="genic" & exons==1 & length>2000) %>%
                                      nrow()/Genic_count*100
      
      # Reference body coverage
      ref_body_cover_genic <- sorted_classification %>%
                              filter(structural_category=="genic" & length/ref_length*100>=45) %>%
                              nrow()/Genic_count*100
      
      # Coding/non-coding
      cod_genic <- sorted_classification %>%
                   filter(structural_category=="genic" & coding=="coding") %>%
                   nrow()/Genic_count*100
      
      ncod_genic <- sorted_classification %>%
                    filter(structural_category=="genic" & coding=="non_coding") %>%
                    nrow()/Genic_count*100
    } 
    
    if (Genic_intron_count==0){
      sub_genic_intron_sqanti_props <- c(Mono_exon_genic_intron=0,
                                         Multi_exon_genic_intron=0) %>% as.numeric()
      
      two_fifty_length_reads_genic_intron <- 0
      five_hund_length_reads_genic_intron <- 0
      short_length_reads_genic_intron <- 0 
      mid_length_reads_genic_intron <- 0
      long_length_reads_genic_intron <- 0
      mono_two_fifty_length_reads_genic_intron <- 0
      mono_five_hund_length_reads_genic_intron <- 0
      mono_short_length_reads_genic_intron <- 0 
      mono_mid_length_reads_genic_intron <- 0
      mono_long_length_reads_genic_intron <- 0
      ref_body_cover_genic_intron <- 0
      cod_genic_intron <- 0
      ncod_genic_intron <- 100
    } else {
      sub_genic_intron_sqanti_props <- sorted_classification %>%
                                       filter(structural_category=="genic_intron") %>%
                                       count(subcategory, name = "reads_in_subcat") %>%
                                       mutate(prop_subcat=reads_in_subcat/Genic_intron_count*100) %>%
                                       complete(subcategory = factor(subcategory,
                                                                     levels = c("mono-exon",
                                                                                "multi-exon")),  
                                                fill = list(prop_subcat = 0)) %>% 
                                       arrange(fct_relevel(subcategory,
                                                           "mono-exon",
                                                           "multi-exon")) %>% .$prop_subcat
      
      # Read length per category
      two_fifty_length_reads_genic_intron <- sorted_classification %>%
                                             filter(structural_category=="genic_intron" & length<=250) %>%
                                             nrow()/Genic_intron_count*100
      
      five_hund_length_reads_genic_intron <- sorted_classification %>%
                                             filter(structural_category=="genic_intron" & length>250 & length<=500) %>%
                                             nrow()/Genic_intron_count*100
      
      short_length_reads_genic_intron <- sorted_classification %>%
                                         filter(structural_category=="genic_intron" & length>500 & length<=1000) %>%
                                         nrow()/Genic_intron_count*100
      
      mid_length_reads_genic_intron <- sorted_classification %>%
                                       filter(structural_category=="genic_intron" & length>1000 & length<=2000) %>%
                                       nrow()/Genic_intron_count*100
      
      long_length_reads_genic_intron <- sorted_classification %>%
                                        filter(structural_category=="genic_intron" & length>2000) %>%
                                        nrow()/Genic_intron_count*100
      
      # Monoexons
      mono_two_fifty_length_reads_genic_intron <- sorted_classification %>%
                                                  filter(structural_category=="genic_intron" & exons==1 & length<=250) %>%
                                                  nrow()/Genic_intron_count*100
      
      mono_five_hund_length_reads_genic_intron <- sorted_classification %>%
                                                  filter(structural_category=="genic_intron" & exons==1 & length>250 & length<=500) %>%
                                                  nrow()/Genic_intron_count*100
      
      mono_short_length_reads_genic_intron <- sorted_classification %>%
                                              filter(structural_category=="genic_intron" & exons==1 & length>500 & length<=1000) %>%
                                              nrow()/Genic_intron_count*100
      
      mono_mid_length_reads_genic_intron <- sorted_classification %>%
                                            filter(structural_category=="genic_intron" & exons==1 & length>1000 & length<=2000) %>%
                                            nrow()/Genic_intron_count*100
      
      mono_long_length_reads_genic_intron <- sorted_classification %>%
                                             filter(structural_category=="genic_intron" & exons==1 & length>2000) %>%
                                             nrow()/Genic_intron_count*100
      
      # Reference body coverage
      ref_body_cover_genic_intron <- sorted_classification %>%
                                     filter(structural_category=="genic_intron" & (length/ref_length*100)>=45) %>%
                                     nrow()/Genic_intron_count*100
      
      # Coding/non-coding
      cod_genic_intron <- sorted_classification %>%
                          filter(structural_category=="genic_intron" & coding=="coding") %>%
                          nrow()/Genic_intron_count*100
      
      ncod_genic_intron <- sorted_classification %>%
                           filter(structural_category=="genic_intron" & coding=="non_coding") %>%
                           nrow()/Genic_intron_count*100
    } 
    
    if (Antisense_count==0){
      sub_antisense_sqanti_props <- c(Mono_exon_antisense=0,
                                      Multi_exon_antisense=0) %>% as.numeric()
      
      two_fifty_length_reads_antisense <- 0
      five_hund_length_reads_antisense <- 0
      short_length_reads_antisense <- 0
      mid_length_reads_antisense <- 0 
      long_length_reads_antisense <- 0 
      mono_two_fifty_length_reads_antisense <- 0
      mono_five_hund_length_reads_antisense <- 0
      mono_short_length_reads_antisense <- 0
      mono_mid_length_reads_antisense <- 0 
      mono_long_length_reads_antisense <- 0 
      ref_body_cover_antisense <- 0
      cod_antisense <- 0
      ncod_antisense <- 100
    } else {
      sub_antisense_sqanti_props <- sorted_classification %>%
                                    filter(structural_category=="antisense") %>%
                                    count(subcategory, name = "reads_in_subcat") %>%
                                    mutate(prop_subcat=reads_in_subcat/Antisense_count*100) %>%
                                    complete(subcategory = factor(subcategory,
                                                                  levels = c("mono-exon",
                                                                             "multi-exon")),  
                                             fill = list(prop_subcat = 0)) %>%
                                    arrange(fct_relevel(subcategory,
                                                        "mono-exon",
                                                        "multi-exon")) %>% .$prop_subcat
      
      # Read length per category
      two_fifty_length_reads_antisense <- sorted_classification %>%
                                          filter(structural_category=="antisense" & length<=250) %>%
                                          nrow()/Antisense_count*100
      
      five_hund_length_reads_antisense <- sorted_classification %>%
                                          filter(structural_category=="antisense" & length>250 & length<=500) %>%
                                          nrow()/Antisense_count*100
      
      short_length_reads_antisense <- sorted_classification %>%
                                      filter(structural_category=="antisense" & length>500 & length<=1000) %>%
                                      nrow()/Antisense_count*100
      
      mid_length_reads_antisense <- sorted_classification %>%
                                    filter(structural_category=="antisense" & length>1000 & length<=2000) %>%
                                    nrow()/Antisense_count*100
      
      long_length_reads_antisense <- sorted_classification %>%
                                     filter(structural_category=="antisense" & length>2000) %>%
                                     nrow()/Antisense_count*100
      
      # Monoexons
      mono_two_fifty_length_reads_antisense <- sorted_classification %>%
                                               filter(structural_category=="antisense" & exons==1 & length<=250) %>%
                                               nrow()/Antisense_count*100
      
      mono_five_hund_length_reads_antisense <- sorted_classification %>%
                                               filter(structural_category=="antisense" & exons==1 & length>250 & length<=500) %>%
                                               nrow()/Antisense_count*100
      
      mono_short_length_reads_antisense <- sorted_classification %>%
                                           filter(structural_category=="antisense" & exons==1 & length>500 & length<=1000) %>%
                                           nrow()/Antisense_count*100
      
      mono_mid_length_reads_antisense <- sorted_classification %>%
                                         filter(structural_category=="antisense" & exons==1 & length>1000 & length<=2000) %>%
                                         nrow()/Antisense_count*100
      
      mono_long_length_reads_antisense <- sorted_classification %>%
                                          filter(structural_category=="antisense" & exons==1 & length>2000) %>%
                                          nrow()/Antisense_count*100
      
      # Reference body coverage
      ref_body_cover_antisense <- sorted_classification %>%
                                  filter(structural_category=="antisense" & length/ref_length*100>=45) %>%
                                  nrow()/Antisense_count*100
      
      # Coding/non-coding
      cod_antisense <- sorted_classification %>%
                       filter(structural_category=="antisense" & coding=="coding") %>%
                       nrow()/Antisense_count*100
      
      ncod_antisense <- sorted_classification %>%
                        filter(structural_category=="antisense" & coding=="non_coding") %>%
                        nrow()/Antisense_count*100
      }
    
    if (Intergenic_count==0){
      sub_intergenic_sqanti_props <- c(Mono_exon_antisense=0,
                                       Multi_exon_antisense=0) %>% as.numeric()
      
      two_fifty_length_reads_intergenic <- 0
      five_hund_length_reads_intergenic <- 0
      short_length_reads_intergenic <- 0
      mid_length_reads_intergenic <- 0
      long_length_reads_intergenic <- 0
      mono_two_fifty_length_reads_intergenic <- 0
      mono_five_hund_length_reads_intergenic <- 0
      mono_short_length_reads_intergenic <- 0
      mono_mid_length_reads_intergenic <- 0
      mono_long_length_reads_intergenic <- 0
      ref_body_cover_intergenic <- 0
      cod_intergenic <- 0
      ncod_intergenic <- 100
    } else {
      sub_intergenic_sqanti_props <- sorted_classification %>%
                                     filter(structural_category=="intergenic") %>%
                                     count(subcategory, name = "reads_in_subcat") %>%
                                     mutate(prop_subcat=reads_in_subcat/Intergenic_count*100) %>%
                                     complete(subcategory = factor(subcategory,
                                                                   levels = c("mono-exon",
                                                                              "multi-exon")),  
                                              fill = list(prop_subcat = 0)) %>% 
                                     arrange(fct_relevel(subcategory,
                                                         "mono-exon",
                                                         "multi-exon")) %>% .$prop_subcat
      
      # Read length per category
      two_fifty_length_reads_intergenic <- sorted_classification %>%
                                           filter(structural_category=="intergenic" & length<=250) %>%
                                           nrow()/Intergenic_count*100
      
      five_hund_length_reads_intergenic <- sorted_classification %>%
                                           filter(structural_category=="intergenic" & length>250 & length<=500) %>%
                                           nrow()/Intergenic_count*100
      
      short_length_reads_intergenic <- sorted_classification %>%
                                       filter(structural_category=="intergenic" & length>500 & length<=1000) %>%
                                       nrow()/Intergenic_count*100
      
      mid_length_reads_intergenic <- sorted_classification %>%
                                     filter(structural_category=="intergenic" & length>1000 & length<=2000) %>%
                                     nrow()/Intergenic_count*100
      
      long_length_reads_intergenic <- sorted_classification %>%
                                      filter(structural_category=="intergenic" & length>2000) %>%
                                      nrow()/Intergenic_count*100
      
      # Monoexons
      mono_two_fifty_length_reads_intergenic <- sorted_classification %>%
                                                filter(structural_category=="intergenic" & exons==1 & length<=250) %>%
                                                nrow()/Intergenic_count*100
      
      mono_five_hund_length_reads_intergenic <- sorted_classification %>%
                                                filter(structural_category=="intergenic" & exons==1 & length>250 & length<=500) %>%
                                                nrow()/Intergenic_count*100
      
      mono_short_length_reads_intergenic <- sorted_classification %>%
                                            filter(structural_category=="intergenic" & exons==1 & length>500 & length<=1000) %>%
                                            nrow()/Intergenic_count*100
      
      mono_mid_length_reads_intergenic <- sorted_classification %>%
                                          filter(structural_category=="intergenic" & exons==1 & length>1000 & length<=2000) %>%
                                          nrow()/Intergenic_count*100
      
      mono_long_length_reads_intergenic <- sorted_classification %>%
                                           filter(structural_category=="intergenic" & exons==1 & length>2000) %>%
                                           nrow()/Intergenic_count*100
      
      # Reference body coverage
      ref_body_cover_intergenic <- sorted_classification %>%
                                   filter(structural_category=="intergenic" & length/ref_length*100>=45) %>%
                                   nrow()/Intergenic_count*100
      
      # Coding/non-coding
      cod_intergenic <- sorted_classification %>%
                        filter(structural_category=="intergenic" & coding=="coding") %>%
                        nrow()/Intergenic_count*100
      
      ncod_intergenic <- sorted_classification %>%
                         filter(structural_category=="intergenic" & coding=="non_coding") %>%
                         nrow()/Intergenic_count*100
    }
    
    ### BAD QUALITY METRICS ###
    # Percentage of RTS
    RTS_in_cell_prop <- sorted_classification %>%
                        filter(RTS_stage==TRUE) %>%
                        nrow()/total_reads*100

    # RTS by category
    FSM_RTS_prop <- if(FSM_count == 0) 0 else sorted_classification %>%
                    filter(structural_category == "full-splice_match" & RTS_stage == TRUE) %>%
                    nrow() / FSM_count * 100
  
    ISM_RTS_prop <- if(ISM_count == 0) 0 else sorted_classification %>%
                    filter(structural_category == "incomplete-splice_match" & RTS_stage == TRUE) %>%
                    nrow() / ISM_count * 100
    
    NIC_RTS_prop <- if(NIC_count == 0) 0 else sorted_classification %>%
                    filter(structural_category == "novel_in_catalog" & RTS_stage == TRUE) %>%
                    nrow() / NIC_count * 100
    
    NNC_RTS_prop <- if(NNC_count == 0) 0 else sorted_classification %>%
                    filter(structural_category == "novel_not_in_catalog" & RTS_stage == TRUE) %>%
                    nrow() / NNC_count * 100

    # Percentage of non_canonical
    if (total_reads_no_monoexon==0){
      non_canonical_in_cell_prop <- 0
    } else {
    non_canonical_in_cell_prop <- sorted_classification %>%
                                  filter(all_canonical=="non_canonical") %>%
                                  nrow()/total_reads_no_monoexon*100
    }
    
    # Non-canonical by category (exclude monoexons)
    FSM_noncanon_prop <- if(FSM_count == 0) 0 else sorted_classification %>%
                         filter(structural_category == "full-splice_match" & exons > 1 & all_canonical == "non_canonical") %>%
                         nrow() / max(1, nrow(sorted_classification %>% filter(structural_category == "full-splice_match" & exons > 1))) * 100
    
    ISM_noncanon_prop <- if(ISM_count == 0) 0 else sorted_classification %>%
                         filter(structural_category == "incomplete-splice_match" & exons > 1 & all_canonical == "non_canonical") %>%
                         nrow() / max(1, nrow(sorted_classification %>% filter(structural_category == "incomplete-splice_match" & exons > 1))) * 100
    
    NIC_noncanon_prop <- if(NIC_count == 0) 0 else sorted_classification %>%
                         filter(structural_category == "novel_in_catalog" & exons > 1 & all_canonical == "non_canonical") %>%
                         nrow() / max(1, nrow(sorted_classification %>% filter(structural_category == "novel_in_catalog" & exons > 1))) * 100
    
    NNC_noncanon_prop <- if(NNC_count == 0) 0 else sorted_classification %>%
                         filter(structural_category == "novel_not_in_catalog" & exons > 1 & all_canonical == "non_canonical") %>%
                         nrow() / max(1, nrow(sorted_classification %>% filter(structural_category == "novel_not_in_catalog" & exons > 1))) * 100

    # Percentage of intrapriming
    intrapriming_in_cell_prop <- sorted_classification %>%
                                 filter(perc_A_downstream_TTS>=60) %>%
                                 nrow()/total_reads*100
    

    # Intrapriming by category
    FSM_intrapriming_prop <- if(FSM_count == 0) 0 else sorted_classification %>%
                             filter(structural_category == "full-splice_match" & perc_A_downstream_TTS >= 60) %>%
                             nrow() / FSM_count * 100
    
    ISM_intrapriming_prop <- if(ISM_count == 0) 0 else sorted_classification %>%
                             filter(structural_category == "incomplete-splice_match" & perc_A_downstream_TTS >= 60) %>%
                             nrow() / ISM_count * 100
    
    NIC_intrapriming_prop <- if(NIC_count == 0) 0 else sorted_classification %>%
                             filter(structural_category == "novel_in_catalog" & perc_A_downstream_TTS >= 60) %>%
                             nrow() / NIC_count * 100
    
    NNC_intrapriming_prop <- if(NNC_count == 0) 0 else sorted_classification %>%
                             filter(structural_category == "novel_not_in_catalog" & perc_A_downstream_TTS >= 60) %>%
                             nrow() / NNC_count * 100

    ##   Adding NMD in the future. 
    ##   Maybe also dist to TTS/TES
    
    ### GOOD QUALITY METRICS ###
    # Annotated genes -- done
    annotated_genes_in_cell_prop<- annotated_genes/genes_in_cell*100

    # Annotated genes by category
    total_genes_FSM <- sorted_classification %>% 
                       filter(structural_category == "full-splice_match") %>%
                       select(associated_gene) %>% n_distinct()
    
    total_genes_ISM <- sorted_classification %>% 
                       filter(structural_category == "incomplete-splice_match") %>%
                       select(associated_gene) %>% n_distinct()
    
    total_genes_NIC <- sorted_classification %>% 
                       filter(structural_category == "novel_in_catalog") %>%
                       select(associated_gene) %>% n_distinct()
    
    total_genes_NNC <- sorted_classification %>% 
                       filter(structural_category == "novel_not_in_catalog") %>%
                       select(associated_gene) %>% n_distinct()
    
    FSM_anno_genes_prop <- if(total_genes_FSM == 0) 0 else sorted_classification %>%
                           filter(structural_category == "full-splice_match") %>%
                           select(associated_gene) %>%
                           filter(!grepl("^novel", associated_gene)) %>%
                           n_distinct() / total_genes_FSM * 100
    
    ISM_anno_genes_prop <- if(total_genes_ISM == 0) 0 else sorted_classification %>%
                           filter(structural_category == "incomplete-splice_match") %>%
                           select(associated_gene) %>%
                           filter(!grepl("^novel", associated_gene)) %>%
                           n_distinct() / total_genes_ISM * 100
    
    NIC_anno_genes_prop <- if(total_genes_NIC == 0) 0 else sorted_classification %>%
                           filter(structural_category == "novel_in_catalog") %>%
                           select(associated_gene) %>%
                           filter(!grepl("^novel", associated_gene)) %>%
                           n_distinct() / total_genes_NIC * 100
    
    NNC_anno_genes_prop <- if(total_genes_NNC == 0) 0 else sorted_classification %>%
                           filter(structural_category == "novel_not_in_catalog") %>%
                           select(associated_gene) %>%
                           filter(!grepl("^novel", associated_gene)) %>%
                           n_distinct() / total_genes_NNC * 100

    # Annotated junction strings
    anno_models_in_cell_prop <- sorted_classification %>%
                                group_by(associated_gene) %>%
                                filter(!grepl("^novel", associated_transcript)) %>%
                                summarise(t_chains=n_distinct(jxn_string)) %>%
                                .$t_chains %>% sum()/models_in_cell*100

    # Percentage of canonical
    if (total_reads_no_monoexon==0){
      canonical_in_cell_prop <- 0
    } else {
    canonical_in_cell_prop <- sorted_classification %>%
                              filter(all_canonical=="canonical") %>%
                              nrow()/total_reads_no_monoexon*100
    }
    
    # Canonical by category (exclude monoexons)
    FSM_canon_prop <- if(FSM_count == 0) 0 else sorted_classification %>%
                      filter(structural_category == "full-splice_match" & exons > 1 & all_canonical == "canonical") %>%
                      nrow() / max(1, nrow(sorted_classification %>% filter(structural_category == "full-splice_match" & exons > 1))) * 100
    
    ISM_canon_prop <- if(ISM_count == 0) 0 else sorted_classification %>%
                      filter(structural_category == "incomplete-splice_match" & exons > 1 & all_canonical == "canonical") %>%
                      nrow() / max(1, nrow(sorted_classification %>% filter(structural_category == "incomplete-splice_match" & exons > 1))) * 100
    
    NIC_canon_prop <- if(NIC_count == 0) 0 else sorted_classification %>%
                      filter(structural_category == "novel_in_catalog" & exons > 1 & all_canonical == "canonical") %>%
                      nrow() / max(1, nrow(sorted_classification %>% filter(structural_category == "novel_in_catalog" & exons > 1))) * 100
    
    NNC_canon_prop <- if(NNC_count == 0) 0 else sorted_classification %>%
                      filter(structural_category == "novel_not_in_catalog" & exons > 1 & all_canonical == "canonical") %>%
                      nrow() / max(1, nrow(sorted_classification %>% filter(structural_category == "novel_not_in_catalog" & exons > 1))) * 100

    # Hacer tabla intermadia con los datos por celula (guardar en temp o al final del report?)
    if (exists("SQANTI_cell_summary")==FALSE){
      SQANTI_cell_summary <- c(CB_id,
                              total_reads,
                              total_UMI,
                              genes_in_cell,
                              models_in_cell,
                              annotated_genes,
                              novel_genes, # Std cell counts
                              MT_perc,
                              known_canonical_prop,
                              known_non_canonical_prop,
                              novel_canonical_prop,
                              novel_non_canonical_prop, # Add canonical/non-canonical per (sub)structural category 
                              FSM_count,
                              ISM_count,
                              NIC_count,
                              NNC_count,
                              Genic_count,
                              Antisense_count,
                              Fusion_count,
                              Intergenic_count,
                              Genic_intron_count, # Structural categories counts (reads)
                              sqanti_props, # Structural categories props (reads)
                              cod_FSM,
                              ncod_FSM,
                              cod_ISM,
                              ncod_ISM,
                              cod_NIC,
                              ncod_NIC,
                              cod_NNC,
                              ncod_NNC,
                              cod_genic,
                              ncod_genic,
                              cod_antisense,
                              ncod_antisense,
                              cod_fusion,
                              ncod_fusion,
                              cod_intergenic,
                              ncod_intergenic,
                              cod_genic_intron,
                              ncod_genic_intron, # Coding and non-coding per structural category (reads)
                              sub_FSM_sqanti_props,
                              sub_ISM_sqanti_props,
                              sub_NIC_sqanti_probs,
                              sub_NNC_sqanti_probs,
                              sub_genic_sqanti_props,
                              sub_antisense_sqanti_props,
                              sub_fusion_sqanti_props,
                              sub_intergenic_sqanti_props,
                              sub_genic_intron_sqanti_props, # Structural subcategories props (reads)
                              two_fifty_length_reads,
                              mono_two_fifty_length_reads,
                              five_hund_length_reads,
                              mono_five_hund_length_reads,
                              short_length_reads,
                              mono_short_length_reads,
                              mid_length_reads,
                              mono_mid_length_reads,
                              long_length_reads,
                              mono_long_length_reads, # Read lengths breaks general (reads)
                              two_fifty_length_reads_FSM,
                              mono_two_fifty_length_reads_FSM,
                              five_hund_length_reads_FSM,
                              mono_five_hund_length_reads_FSM,
                              short_length_reads_FSM,
                              mono_short_length_reads_FSM,
                              mid_length_reads_FSM,
                              mono_mid_length_reads_FSM,
                              long_length_reads_FSM,
                              mono_long_length_reads_FSM,
                              two_fifty_length_reads_ISM,
                              mono_two_fifty_length_reads_ISM,
                              five_hund_length_reads_ISM,
                              mono_five_hund_length_reads_ISM,
                              short_length_reads_ISM,
                              mono_short_length_reads_ISM,
                              mid_length_reads_ISM,
                              mono_mid_length_reads_ISM,
                              long_length_reads_ISM,
                              mono_long_length_reads_ISM,
                              two_fifty_length_reads_NIC,
                              mono_two_fifty_length_reads_NIC,
                              five_hund_length_reads_NIC,
                              mono_five_hund_length_reads_NIC,
                              short_length_reads_NIC,
                              mono_short_length_reads_NIC,
                              mid_length_reads_NIC,
                              mono_mid_length_reads_NIC,
                              long_length_reads_NIC,
                              mono_long_length_reads_NIC,
                              two_fifty_length_reads_NNC,
                              mono_two_fifty_length_reads_NNC,
                              five_hund_length_reads_NNC,
                              mono_five_hund_length_reads_NNC,
                              short_length_reads_NNC,
                              mono_short_length_reads_NNC,
                              mid_length_reads_NNC,
                              mono_mid_length_reads_NNC,
                              long_length_reads_NNC,
                              mono_long_length_reads_NNC,
                              two_fifty_length_reads_genic,
                              mono_two_fifty_length_reads_genic,
                              five_hund_length_reads_genic,
                              mono_five_hund_length_reads_genic,
                              short_length_reads_genic,
                              mono_short_length_reads_genic,
                              mid_length_reads_genic,
                              mono_mid_length_reads_genic,
                              long_length_reads_genic,
                              mono_long_length_reads_genic,
                              two_fifty_length_reads_antisense,
                              mono_two_fifty_length_reads_antisense,
                              five_hund_length_reads_antisense,
                              mono_five_hund_length_reads_antisense,
                              short_length_reads_antisense,
                              mono_short_length_reads_antisense,
                              mid_length_reads_antisense,
                              mono_mid_length_reads_antisense,
                              long_length_reads_antisense,
                              mono_long_length_reads_antisense, 
                              two_fifty_length_reads_fusion,
                              mono_two_fifty_length_reads_fusion,
                              five_hund_length_reads_fusion,
                              mono_five_hund_length_reads_fusion,
                              short_length_reads_fusion,
                              mono_short_length_reads_fusion,
                              mid_length_reads_fusion,
                              mono_mid_length_reads_fusion,
                              long_length_reads_fusion,
                              mono_long_length_reads_fusion,
                              two_fifty_length_reads_intergenic,
                              mono_two_fifty_length_reads_intergenic,
                              five_hund_length_reads_intergenic,
                              mono_five_hund_length_reads_intergenic,
                              short_length_reads_intergenic,
                              mono_short_length_reads_intergenic,
                              mid_length_reads_intergenic,
                              mono_mid_length_reads_intergenic,
                              long_length_reads_intergenic,
                              mono_long_length_reads_intergenic,
                              two_fifty_length_reads_genic_intron,
                              mono_two_fifty_length_reads_genic_intron,
                              five_hund_length_reads_genic_intron,
                              mono_five_hund_length_reads_genic_intron,
                              short_length_reads_genic_intron,
                              mono_short_length_reads_genic_intron,
                              mid_length_reads_genic_intron,
                              mono_mid_length_reads_genic_intron,
                              long_length_reads_genic_intron,
                              mono_long_length_reads_genic_intron, # Reads length breaks per structural category (including monoexon breaks)
                              ref_body_cover_FSM,
                              ref_body_cover_ISM,
                              ref_body_cover_NIC,
                              ref_body_cover_NNC,
                              ref_body_cover_genic,
                              ref_body_cover_antisense,
                              ref_body_cover_fusion,
                              ref_body_cover_intergenic,
                              ref_body_cover_genic_intron, # Coverage of reference length (set at 45% default)
                              RTS_in_cell_prop,
                              non_canonical_in_cell_prop,
                              intrapriming_in_cell_prop, 
                              FSM_RTS_prop, ISM_RTS_prop, NIC_RTS_prop, NNC_RTS_prop,
                              FSM_noncanon_prop, ISM_noncanon_prop, NIC_noncanon_prop, NNC_noncanon_prop,
                              FSM_intrapriming_prop, ISM_intrapriming_prop, NIC_intrapriming_prop, NNC_intrapriming_prop, # Features of bad quality
                              annotated_genes_in_cell_prop,
                              anno_models_in_cell_prop,
                              canonical_in_cell_prop,
                              FSM_anno_genes_prop, ISM_anno_genes_prop, NIC_anno_genes_prop, NNC_anno_genes_prop,
                              FSM_canon_prop, ISM_canon_prop, NIC_canon_prop, NNC_canon_prop) # Features of good quality
    } else {
      SQANTI_cell_summary <- rbind(SQANTI_cell_summary, c(CB_id,
                                                        total_reads,
                                                        total_UMI,
                                                        genes_in_cell,
                                                        models_in_cell,
                                                        annotated_genes,
                                                        novel_genes, # Std cell counts
                                                        MT_perc,
                                                        known_canonical_prop,
                                                        known_non_canonical_prop,
                                                        novel_canonical_prop,
                                                        novel_non_canonical_prop, # Add canonical/non-canonical per (sub)structural category 
                                                        FSM_count,
                                                        ISM_count,
                                                        NIC_count,
                                                        NNC_count,
                                                        Genic_count,
                                                        Antisense_count,
                                                        Fusion_count,
                                                        Intergenic_count,
                                                        Genic_intron_count, # Structural categories counts (reads)
                                                        sqanti_props, # Structural categories props (reads)
                                                        cod_FSM,
                                                        ncod_FSM,
                                                        cod_ISM,
                                                        ncod_ISM,
                                                        cod_NIC,
                                                        ncod_NIC,
                                                        cod_NNC,
                                                        ncod_NNC,
                                                        cod_genic,
                                                        ncod_genic,
                                                        cod_antisense,
                                                        ncod_antisense,
                                                        cod_fusion,
                                                        ncod_fusion,
                                                        cod_intergenic,
                                                        ncod_intergenic,
                                                        cod_genic_intron,
                                                        ncod_genic_intron, # Coding and non-coding per structural category (reads)
                                                        sub_FSM_sqanti_props,
                                                        sub_ISM_sqanti_props,
                                                        sub_NIC_sqanti_probs,
                                                        sub_NNC_sqanti_probs,
                                                        sub_genic_sqanti_props,
                                                        sub_antisense_sqanti_props,
                                                        sub_fusion_sqanti_props,
                                                        sub_intergenic_sqanti_props,
                                                        sub_genic_intron_sqanti_props, # Structural subcategories props (reads)
                                                        two_fifty_length_reads,
                                                        mono_two_fifty_length_reads,
                                                        five_hund_length_reads,
                                                        mono_five_hund_length_reads,
                                                        short_length_reads,
                                                        mono_short_length_reads,
                                                        mid_length_reads,
                                                        mono_mid_length_reads,
                                                        long_length_reads,
                                                        mono_long_length_reads, # Read lengths breaks general (reads)
                                                        two_fifty_length_reads_FSM,
                                                        mono_two_fifty_length_reads_FSM,
                                                        five_hund_length_reads_FSM,
                                                        mono_five_hund_length_reads_FSM,
                                                        short_length_reads_FSM,
                                                        mono_short_length_reads_FSM,
                                                        mid_length_reads_FSM,
                                                        mono_mid_length_reads_FSM,
                                                        long_length_reads_FSM,
                                                        mono_long_length_reads_FSM,
                                                        two_fifty_length_reads_ISM,
                                                        mono_two_fifty_length_reads_ISM,
                                                        five_hund_length_reads_ISM,
                                                        mono_five_hund_length_reads_ISM,
                                                        short_length_reads_ISM,
                                                        mono_short_length_reads_ISM,
                                                        mid_length_reads_ISM,
                                                        mono_mid_length_reads_ISM,
                                                        long_length_reads_ISM,
                                                        mono_long_length_reads_ISM,
                                                        two_fifty_length_reads_NIC,
                                                        mono_two_fifty_length_reads_NIC,
                                                        five_hund_length_reads_NIC,
                                                        mono_five_hund_length_reads_NIC,
                                                        short_length_reads_NIC,
                                                        mono_short_length_reads_NIC,
                                                        mid_length_reads_NIC,
                                                        mono_mid_length_reads_NIC,
                                                        long_length_reads_NIC,
                                                        mono_long_length_reads_NIC,
                                                        two_fifty_length_reads_NNC,
                                                        mono_two_fifty_length_reads_NNC,
                                                        five_hund_length_reads_NNC,
                                                        mono_five_hund_length_reads_NNC,
                                                        short_length_reads_NNC,
                                                        mono_short_length_reads_NNC,
                                                        mid_length_reads_NNC,
                                                        mono_mid_length_reads_NNC,
                                                        long_length_reads_NNC,
                                                        mono_long_length_reads_NNC,
                                                        two_fifty_length_reads_genic,
                                                        mono_two_fifty_length_reads_genic,
                                                        five_hund_length_reads_genic,
                                                        mono_five_hund_length_reads_genic,
                                                        short_length_reads_genic,
                                                        mono_short_length_reads_genic,
                                                        mid_length_reads_genic,
                                                        mono_mid_length_reads_genic,
                                                        long_length_reads_genic,
                                                        mono_long_length_reads_genic,
                                                        two_fifty_length_reads_antisense,
                                                        mono_two_fifty_length_reads_antisense,
                                                        five_hund_length_reads_antisense,
                                                        mono_five_hund_length_reads_antisense,
                                                        short_length_reads_antisense,
                                                        mono_short_length_reads_antisense,
                                                        mid_length_reads_antisense,
                                                        mono_mid_length_reads_antisense,
                                                        long_length_reads_antisense,
                                                        mono_long_length_reads_antisense, 
                                                        two_fifty_length_reads_fusion,
                                                        mono_two_fifty_length_reads_fusion,
                                                        five_hund_length_reads_fusion,
                                                        mono_five_hund_length_reads_fusion,
                                                        short_length_reads_fusion,
                                                        mono_short_length_reads_fusion,
                                                        mid_length_reads_fusion,
                                                        mono_mid_length_reads_fusion,
                                                        long_length_reads_fusion,
                                                        mono_long_length_reads_fusion,
                                                        two_fifty_length_reads_intergenic,
                                                        mono_two_fifty_length_reads_intergenic,
                                                        five_hund_length_reads_intergenic,
                                                        mono_five_hund_length_reads_intergenic,
                                                        short_length_reads_intergenic,
                                                        mono_short_length_reads_intergenic,
                                                        mid_length_reads_intergenic,
                                                        mono_mid_length_reads_intergenic,
                                                        long_length_reads_intergenic,
                                                        mono_long_length_reads_intergenic,
                                                        two_fifty_length_reads_genic_intron,
                                                        mono_two_fifty_length_reads_genic_intron,
                                                        five_hund_length_reads_genic_intron,
                                                        mono_five_hund_length_reads_genic_intron,
                                                        short_length_reads_genic_intron,
                                                        mono_short_length_reads_genic_intron,
                                                        mid_length_reads_genic_intron,
                                                        mono_mid_length_reads_genic_intron,
                                                        long_length_reads_genic_intron,
                                                        mono_long_length_reads_genic_intron, # Reads length breaks per structural category (including monoexon breaks)
                                                        ref_body_cover_FSM,
                                                        ref_body_cover_ISM,
                                                        ref_body_cover_NIC,
                                                        ref_body_cover_NNC,
                                                        ref_body_cover_genic,
                                                        ref_body_cover_antisense,
                                                        ref_body_cover_fusion,
                                                        ref_body_cover_intergenic,
                                                        ref_body_cover_genic_intron, # Coverage of reference length (set at 45% default)
                                                        RTS_in_cell_prop,
                                                        non_canonical_in_cell_prop,
                                                        intrapriming_in_cell_prop, 
                                                        FSM_RTS_prop, ISM_RTS_prop, NIC_RTS_prop, NNC_RTS_prop,
                                                        FSM_noncanon_prop, ISM_noncanon_prop, NIC_noncanon_prop, NNC_noncanon_prop,
                                                        FSM_intrapriming_prop, ISM_intrapriming_prop, NIC_intrapriming_prop, NNC_intrapriming_prop, # Features of bad quality
                                                        annotated_genes_in_cell_prop,
                                                        anno_models_in_cell_prop,
                                                        canonical_in_cell_prop,
                                                        FSM_anno_genes_prop, ISM_anno_genes_prop, NIC_anno_genes_prop, NNC_anno_genes_prop,
                                                        FSM_canon_prop, ISM_canon_prop, NIC_canon_prop, NNC_canon_prop)) # Features of good quality
    }
  }
  SQANTI_cell_summary <- as.data.frame(SQANTI_cell_summary)
  rownames(SQANTI_cell_summary) <- NULL
  colnames(SQANTI_cell_summary) <- c("CB",
                                    "Reads_in_cell",
                                    "UMIs_in_cell",
                                    "Genes_in_cell",
                                    "UJCs_in_cell",
                                    "Annotated_genes",
                                    "Novel_genes", # Std cell counts
                                    "MT_perc",
                                    "Known_canonical_prop",
                                    "Known_non_canonical_prop",
                                    "Novel_canonical_prop",
                                    "Novel_non_canonical_prop", # Canonical/non-canonical (reads without monoexons)
                                    "FSM",
                                    "ISM",
                                    "NIC",
                                    "NNC",
                                    "Genic_Genomic",
                                    "Antisense",
                                    "Fusion",
                                    "Intergenic",
                                    "Genic_intron", # Structural categories counts (reads) 
                                    "FSM_prop", 
                                    "ISM_prop", 
                                    "NIC_prop", 
                                    "NNC_prop", 
                                    "Genic_Genomic_prop",
                                    "Antisense_prop",
                                    "Fusion_prop",
                                    "Intergenic_prop", 
                                    "Genic_intron_prop", # Structural categories props (reads)
                                    "Coding_FSM_prop",
                                    "Non_coding_FSM_prop",
                                    "Coding_ISM_prop",
                                    "Non_coding_ISM_prop",
                                    "Coding_NIC_prop",
                                    "Non_coding_NIC_prop",
                                    "Coding_NNC_prop",
                                    "Non_coding_NNC_prop",
                                    "Coding_genic_prop",
                                    "Non_coding_genic_prop",
                                    "Coding_antisense_prop",
                                    "Non_coding_antisense_prop",
                                    "Coding_fusion_prop",
                                    "Non_coding_fusion_prop",
                                    "Coding_intergenic_prop",
                                    "Non_coding_intergenic_prop",
                                    "Coding_genic_intron_prop",
                                    "Non_coding_genic_intron_prop", # Coding and non-coding per structural category (reads)
                                    "FSM_alternative_3.end_prop",
                                    "FSM_alternative_3.5.end_prop",
                                    "FSM_alternative_5.end_prop",
                                    "FSM_reference_match_prop",
                                    "FSM_mono.exon_prop",
                                    "ISM_3._fragment_prop",
                                    "ISM_internal_fragment_prop",
                                    "ISM_5._fragment_prop",
                                    "ISM_intron_retention_prop",
                                    "ISM_mono.exon_prop",
                                    "NIC_comb_annot_junctions_prop", 
                                    "NIC_comb_annot_splice_sites_prop",
                                    "NIC_intron_retention_prop",
                                    "NIC_mono.exon_prop",
                                    "NIC_mono.exon_by_intron_retention_prop",
                                    "NNC_at_least_1_don_accept_prop",
                                    "NNC_intron_retention_prop",
                                    "Genic_mono.exon_prop", 
                                    "Genic_multi.exon_prop", 
                                    "Antisense_mono.exon_prop", 
                                    "Antisense_multi.exon_prop",
                                    "Fusion_intron_retention_prop",
                                    "Fusion_multi.exon_prop", 
                                    "Intergenic_mono.exon_prop",
                                    "Intergenic_multi.exon_prop",
                                    "Genic_intron_mono.exon_prop",
                                    "Genic_intron_multi.exon_prop", # Structural subcategories props (reads)
                                    "Total_250b_length_prop",
                                    "Total_250b_length_mono_prop",
                                    "Total_500b_length_prop",
                                    "Total_500b_length_mono_prop",
                                    "Total_short_length_prop",
                                    "Total_short_length_mono_prop",
                                    "Total_mid_length_prop",
                                    "Total_mid_length_mono_prop",
                                    "Total_long_length_prop",
                                    "Total_long_length_mono_prop", # Read lengths breaks general + monoexons (reads)
                                    "FSM_250b_length_prop",
                                    "FSM_250b_length_mono_prop",
                                    "FSM_500b_length_prop",
                                    "FSM_500b_length_mono_prop",
                                    "FSM_short_length_prop",
                                    "FSM_short_length_mono_prop",
                                    "FSM_mid_length_prop",
                                    "FSM_mid_length_mono_prop",
                                    "FSM_long_length_prop",
                                    "FSM_long_length_mono_prop",
                                    "ISM_250b_length_prop",
                                    "ISM_250b_length_mono_prop",
                                    "ISM_500b_length_prop",
                                    "ISM_500b_length_mono_prop",
                                    "ISM_short_length_prop",
                                    "ISM_short_length_mono_prop",
                                    "ISM_mid_length_prop",
                                    "ISM_mid_length_mono_prop",
                                    "ISM_long_length_prop",
                                    "ISM_long_length_mono_prop",
                                    "NIC_250b_length_prop",
                                    "NIC_250b_length_mono_prop",
                                    "NIC_500b_length_prop",
                                    "NIC_500b_length_mono_prop",
                                    "NIC_short_length_prop",
                                    "NIC_short_length_mono_prop",
                                    "NIC_mid_length_prop",
                                    "NIC_mid_length_mono_prop",
                                    "NIC_long_length_prop",
                                    "NIC_long_length_mono_prop",
                                    "NNC_250b_length_prop",
                                    "NNC_250b_length_mono_prop",
                                    "NNC_500b_length_prop",
                                    "NNC_500b_length_mono_prop",
                                    "NNC_short_length_prop",
                                    "NNC_short_length_mono_prop",
                                    "NNC_mid_length_prop",
                                    "NNC_mid_length_mono_prop",
                                    "NNC_long_length_prop",
                                    "NNC_long_length_mono_prop",
                                    "Genic_250b_length_prop",
                                    "Genic_250b_length_mono_prop",
                                    "Genic_500b_length_prop",
                                    "Genic_500b_length_mono_prop",
                                    "Genic_short_length_prop",
                                    "Genic_short_length_mono_prop",
                                    "Genic_mid_length_prop",
                                    "Genic_mid_length_mono_prop",
                                    "Genic_long_length_prop",
                                    "Genic_long_length_mono_prop",
                                    "Antisense_250b_length_prop",
                                    "Antisense_250b_length_mono_prop",
                                    "Antisense_500b_length_prop",
                                    "Antisense_500b_length_mono_prop",
                                    "Antisense_short_length_prop",
                                    "Antisense_short_length_mono_prop",
                                    "Antisense_mid_length_prop",
                                    "Antisense_mid_length_mono_prop",
                                    "Antisense_long_length_prop",
                                    "Antisense_long_length_mono_prop",
                                    "Fusion_250b_length_prop",
                                    "Fusion_250b_length_mono_prop",
                                    "Fusion_500b_length_prop",
                                    "Fusion_500b_length_mono_prop",
                                    "Fusion_short_length_prop",
                                    "Fusion_short_length_mono_prop",
                                    "Fusion_mid_length_prop",
                                    "Fusion_mid_length_mono_prop",
                                    "Fusion_long_length_prop",
                                    "Fusion_long_length_mono_prop",
                                    "Intergenic_250b_length_prop",
                                    "Intergenic_250b_length_mono_prop",
                                    "Intergenic_500b_length_prop",
                                    "Intergenic_500b_length_mono_prop",
                                    "Intergenic_short_length_prop",
                                    "Intergenic_short_length_mono_prop",
                                    "Intergenic_mid_length_prop",
                                    "Intergenic_mid_length_mono_prop",
                                    "Intergenic_long_length_prop",
                                    "Intergenic_long_length_mono_prop",
                                    "Genic_intron_250b_length_prop",
                                    "Genic_intron_250b_length_mono_prop",
                                    "Genic_intron_500b_length_prop",
                                    "Genic_intron_500b_length_mono_prop",
                                    "Genic_intron_short_length_prop",
                                    "Genic_intron_short_length_mono_prop",
                                    "Genic_intron_mid_length_prop",
                                    "Genic_intron_mid_length_mono_prop",
                                    "Genic_intron_long_length_prop",
                                    "Genic_intron_long_length_mono_prop", # Reads length breaks per structural category
                                    "FSM_ref_coverage_prop",
                                    "ISM_ref_coverage_prop",
                                    "NIC_ref_coverage_prop",
                                    "NNC_ref_coverage_prop",
                                    "Genic_ref_coverage_prop",
                                    "Antisense_ref_coverage_prop",
                                    "Fusion_ref_coverage_prop",
                                    "Intergenic_ref_coverage_prop",
                                    "Genic_intron_ref_coverage_prop", # Coverage of reference length (set at 45% default)
                                    "RTS_prop_in_cell",
                                    "Non_canonical_prop_in_cell",
                                    "Intrapriming_prop_in_cell", 
                                    "FSM_RTS_prop", "ISM_RTS_prop", "NIC_RTS_prop", "NNC_RTS_prop",
                                    "FSM_noncanon_prop", "ISM_noncanon_prop", "NIC_noncanon_prop", "NNC_noncanon_prop",
                                    "FSM_intrapriming_prop", "ISM_intrapriming_prop", "NIC_intrapriming_prop", "NNC_intrapriming_prop", # Features of bad quality
                                    "Annotated_genes_prop_in_cell",
                                    "Annotated_juction_strings_prop_in_cell",
                                    "Canonical_prop_in_cell",
                                    "FSM_anno_genes_prop", "ISM_anno_genes_prop", "NIC_anno_genes_prop", "NNC_anno_genes_prop",
                                    "FSM_canon_prop", "ISM_canon_prop", "NIC_canon_prop", "NNC_canon_prop") # Features of good quality. Add annotated genes
  
  # Change data type of columns
  SQANTI_cell_summary <- SQANTI_cell_summary %>%
    mutate(across(2:ncol(.), as.numeric))  

  if (Save == "Y"){
    print(paste0("Saving cell summary table. Starting at ", Sys.time()))
    write.table(SQANTI_cell_summary, file = gzfile(paste0(cell_summary_output, ".txt.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
    print(paste0("Cell summary table saved as ", cell_summary_output, ".txt.gz"))
  }
  print(paste0("All steps completed at ",Sys.time()))
  return(SQANTI_cell_summary)
}

generate_sqantisc_plots <- function(SQANTI_cell_summary, Classification_file, report_output){
  
  ### Basic cell informtion ###
  #############################
  
  # Reads in cell
  gg_reads_in_cells <- ggplot(SQANTI_cell_summary, aes(x = "", y = Reads_in_cell)) +
    geom_point(color = "#CC6633", position = position_dodge2(width = 0.8), size = 0.5, alpha = 0.2) +
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
    geom_point(color = "#CC6633", position = position_dodge2(width = 0.8), size = 0.5, alpha = 0.2) +
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
    geom_point(color = "#CC6633", position = position_dodge2(width = 0.8), size = 0.5, alpha = 0.2) + 
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
    geom_point(color = "#CC6633", position = position_dodge2(width = 0.8), size = 0.5, alpha = 0.2) +
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
    geom_point(color = "#CC6633", position = position_dodge2(width = 0.8), size = 0.5, alpha = 0.2) +
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
  
  # Mitochondrial percentage in cell
  gg_MT_perc <- ggplot(SQANTI_cell_summary, aes(x = "", y = MT_perc)) +
    geom_point(color = "#CC6633", position = position_dodge2(width = 0.8), size = 0.5, alpha = 0.2) +
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

  #  Mono/multi-exon prop novel vs annotated genes
  
  ### Length distribution ###
  ###########################
  
  # All reads length distribution
  gg_bulk_all_reads <- ggplot(Classification_file, aes(x=length)) +
    geom_histogram(binwidth=50, fill="#CC6633", color="black", alpha=0.5) +
    labs(title = "All Read Lengths Distribution",
         x = "",
         y = "Reads, counts") +
    theme_classic() +
    theme(
     legend.position = "none",
     plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
     axis.title = element_text(size = 16),
     axis.text.y = element_text(size = 14),
     axis.text.x = element_text(size = 16))
    
  # Length distribution per break (cells)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Total_250b_length_prop", "Total_500b_length_prop",
                                                                "Total_short_length_prop", "Total_mid_length_prop",
                                                                "Total_long_length_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Total_250b_length_prop, Total_500b_length_prop,
                                                                                                       Total_short_length_prop, Total_mid_length_prop,
                                                                                                       Total_long_length_prop)))
  gg_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_point(color = "#CC6633", position = position_dodge2(width = 0.8), size = 0.5, alpha = 0.2) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha=0.3) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Distribution of Read Length Proportions per Cell",
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
    geom_point(color = "#CC6633", position = position_dodge2(width = 0.8), size = 0.5, alpha = 0.2) +
    geom_violin(fill = "#CC6633", color = "black", alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.05, fill = "#CC6633", color = "grey20", outlier.shape = NA, alpha=0.3) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Distribution of Mono-Exonic Read Length Proportions per Cell",
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
    geom_point(alpha = 0.2, color = "#6BAED6", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#6BAED6", color = "#6BAED6", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#6BAED6", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Distribution of FSM Read Length Proportions per Cell",
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
    geom_point(alpha = 0.5, color = "#FC8D59", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#FC8D59", color = "#FC8D59", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#FC8D59", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "ISM Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#78C679", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#78C679", color = "#78C679", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#78C679", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "NIC Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#EE6A50", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#EE6A50", color = "#EE6A50", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#EE6A50", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "NNC Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#969696", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#969696", color = "#969696", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#969696", color = "grey90", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Genic Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#66C2A4", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#66C2A4", color = "#66C2A4", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#66C2A4", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Antisense Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "goldenrod1", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "goldenrod1", color = "goldenrod1", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "goldenrod1", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Fusion Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "darksalmon", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "darksalmon", color = "darksalmon", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "darksalmon", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Intergenic Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#41B6C4", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#41B6C4", color = "#41B6C4", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#41B6C4", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Genic Intron Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#6BAED6", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#6BAED6", color = "#6BAED6", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#6BAED6", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "FSM Mono-exonic Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#FC8D59", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#FC8D59", color = "#FC8D59", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#FC8D59", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "ISM Mono-exonic Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#78C679", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#78C679", color = "#78C679", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#78C679", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "NIC Mono-exonic Reads Length Distribution per Cell",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # NNC
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("NNC_250b_length_mono_prop", "NNC_500b_length_mono_prop",
                                                                "NNC_short_length_mono_prop", "NNC_mid_length_mono_prop",
                                                                "NNC_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(NNC_250b_length_mono_prop, NNC_500b_length_mono_prop,
                                                                                                       NNC_short_length_mono_prop, NNC_mid_length_mono_prop,
                                                                                                       NNC_long_length_mono_prop)))
  gg_NNC_mono_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_point(alpha = 0.5, color = "#EE6A50", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#EE6A50", color = "#EE6A50", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#EE6A50", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "NNC Mono-exonic Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#969696", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#969696", color = "#969696", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#969696", color = "grey90", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Genic Mono-exonic Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#66C2A4", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "#66C2A4", color = "#66C2A4", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "#66C2A4", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Antisense Mono-exonic Reads Length Distribution per Cell",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 16))
  
  # Fusion
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Fusion_250b_length_mono_prop", "Fusion_500b_length_mono_prop",
                                                                "Fusion_short_length_mono_prop", "Fusion_mid_length_mono_prop",
                                                                "Fusion_long_length_mono_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Fusion_250b_length_mono_prop, Fusion_500b_length_mono_prop,
                                                                                                       Fusion_short_length_mono_prop, Fusion_mid_length_mono_prop,
                                                                                                       Fusion_long_length_mono_prop)))
  gg_fusion_mono_read_distr <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_point(alpha = 0.5, color = "goldenrod1", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "goldenrod1", color = "goldenrod1", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "goldenrod1", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Fusion Mono-exonic Reads Length Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "darksalmon", size = 0.5, position = position_dodge2(width = 0.8)) +
    geom_violin(fill = "darksalmon", color = "darksalmon", alpha = 0.7, scale = "width") +
    geom_boxplot(fill = "darksalmon", color = "grey20", alpha=0.3, outlier.shape = NA, width = 0.05) + 
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_x_discrete(labels = c("0-250bp", "250-500bp",
                                "500-1000bp", "1000-2000bp",
                                ">2000bp")) +
    labs(title = "Intergenic Mono-exonic Read Lengths Distribution per Cell",
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
    geom_point(alpha = 0.5, color = "#41B6C4", size = 0.5, position = position_dodge2(width = 0.8)) +
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) +
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
  
  #  Coding/non-coding across structural categories (change it in the future to a combine plot)
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
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
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Known_canonical_prop", "Known_non_canonical_prop",
                                                                "Novel_canonical_prop", "Novel_non_canonical_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Known_canonical_prop, Known_non_canonical_prop,
                                                                                                       Novel_canonical_prop, Novel_non_canonical_prop)))
  gg_known_novel_canon <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    theme_classic(base_size = 14) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = c("Known_canonical_prop" = "#6BAED6",
                                  "Known_non_canonical_prop" = "goldenrod1",
                                  "Novel_canonical_prop" = "#78C679",
                                  "Novel_non_canonical_prop" = "#FC8D59")) +
    scale_fill_manual(values = c("Known_canonical_prop" = "#6BAED6",
                                  "Known_non_canonical_prop" = "goldenrod1",
                                  "Novel_canonical_prop" = "#78C679",
                                  "Novel_non_canonical_prop" = "#FC8D59")) +
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
  
  # Intrapriming  (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_intrapriming_prop", "ISM_intrapriming_prop", 
                                           "NIC_intrapriming_prop", "NNC_intrapriming_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_intrapriming_prop", "ISM_intrapriming_prop", 
                                               "NIC_intrapriming_prop", "NNC_intrapriming_prop"))
  
  gg_intrapriming_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_point(aes(color = Variable), position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    # Use the same intrapriming green color from the original bad_feature plot
    scale_color_manual(values = rep("#78C679", 4)) +
    scale_fill_manual(values = rep("#78C679", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "Intrapriming",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16))

  # RTS  (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_RTS_prop", "ISM_RTS_prop", 
                                           "NIC_RTS_prop", "NNC_RTS_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_RTS_prop", "ISM_RTS_prop", 
                                               "NIC_RTS_prop", "NNC_RTS_prop"))
  
  gg_RTS_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_point(aes(color = Variable), position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = rep("#FF9933", 4)) +
    scale_fill_manual(values = rep("#FF9933", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "RT-switching",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16))

  # Non-canonical  (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_noncanon_prop", "ISM_noncanon_prop", 
                                           "NIC_noncanon_prop", "NNC_noncanon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_noncanon_prop", "ISM_noncanon_prop", 
                                               "NIC_noncanon_prop", "NNC_noncanon_prop"))
  
  gg_noncanon_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_point(aes(color = Variable), position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    # Use the same non-canonical blue color from the original bad_feature plot
    scale_color_manual(values = rep("#41B6C4", 4)) +
    scale_fill_manual(values = rep("#41B6C4", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "Non-Canonical Junctions",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16))



  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Intrapriming_prop_in_cell", "RTS_prop_in_cell",
                                                                "Non_canonical_prop_in_cell"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Intrapriming_prop_in_cell, RTS_prop_in_cell,
                                                                                                       Non_canonical_prop_in_cell)))
  gg_bad_feature <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = c("Intrapriming_prop_in_cell" = "#78C679",
                                  "RTS_prop_in_cell" = "#FF9933",
                                  "Non_canonical_prop_in_cell" = "#41B6C4")) +
    scale_fill_manual(values = c("Intrapriming_prop_in_cell" = "#78C679",
                                 "RTS_prop_in_cell" = "#FF9933",
                                 "Non_canonical_prop_in_cell" = "#41B6C4")) +
    scale_x_discrete(labels = c("Intrapriming", "RT-switching",
                                "Non-Canonical Junctions")) +
    labs(title = "Bad Quality Control Attributes Across Cells",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16))
  
  #  NMD  (split between categories)
  
  ### Good features plot ###
  ##########################
  
  # Annotated genes  (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_anno_genes_prop", "ISM_anno_genes_prop", 
                                           "NIC_anno_genes_prop", "NNC_anno_genes_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_anno_genes_prop", "ISM_anno_genes_prop", 
                                               "NIC_anno_genes_prop", "NNC_anno_genes_prop"))
  
  gg_anno_genes_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_point(aes(color = Variable), position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = rep("#0e5a87", 4)) +
    scale_fill_manual(values = rep("#0e5a87", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "Annotated Genes",
         x = "",
         y = "Genes, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16))

  # Junction strings mapped to annotated transcripts  (split between categories) # maybe not
  # Canonical  (split between categories)
  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, 
                                  cols = c("FSM_canon_prop", "ISM_canon_prop", 
                                           "NIC_canon_prop", "NNC_canon_prop"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, 
                                    levels = c("FSM_canon_prop", "ISM_canon_prop", 
                                               "NIC_canon_prop", "NNC_canon_prop"))
  
  gg_canon_by_category <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, fill = Variable), alpha = 0.7, scale = "width") +  
    geom_point(aes(color = Variable), position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
    geom_boxplot(aes(fill = Variable), color = "grey20",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    # Use the same canonical color from the original good_feature plot
    scale_color_manual(values = rep("#CC6633", 4)) +
    scale_fill_manual(values = rep("#CC6633", 4)) +
    scale_x_discrete(labels = c("FSM", "ISM", "NIC", "NNC")) +
    labs(title = "Canonical Junctions",
         x = "",
         y = "Reads, %") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      axis.title = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = 16))


  gg_SQANTI_pivot <- pivot_longer(SQANTI_cell_summary, cols = c("Annotated_genes_prop_in_cell",
                                                                "Annotated_juction_strings_prop_in_cell",
                                                                "Canonical_prop_in_cell"), 
                                  names_to = "Variable", values_to = "Value") %>% select(Variable, Value)
  
  gg_SQANTI_pivot$Variable <- factor(gg_SQANTI_pivot$Variable, colnames(SQANTI_cell_summary %>% select(Annotated_genes_prop_in_cell,
                                                                                                       Annotated_juction_strings_prop_in_cell,
                                                                                                       Canonical_prop_in_cell)))
  gg_good_feature <- ggplot(gg_SQANTI_pivot, aes(x = Variable, y = Value)) +
    geom_violin(aes(color = Variable, 
                    fill = Variable),
                    alpha = 0.7, scale = "width") +  
    geom_point(aes(color = Variable), 
               position = position_dodge2(width = 0.8), 
               size = 0.5, alpha = 0.8) + 
    geom_boxplot(aes(fill = Variable), color = "grey90",
                 width = 0.08, outlier.shape = NA, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 1, color = "red", stroke = 1) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = c("Annotated_genes_prop_in_cell" = "#0e5a87",
                                  "Annotated_juction_strings_prop_in_cell" = "#6699CC",
                                  "Canonical_prop_in_cell" = "#CC6633")) +
    scale_fill_manual(values = c("Annotated_genes_prop_in_cell" = "#0e5a87",
                                 "Annotated_juction_strings_prop_in_cell" = "#6699CC",
                                 "Canonical_prop_in_cell" = "#CC6633")) +
    scale_x_discrete(labels = c("Annotated genes", "Annotated UJCs",
                                "Reads with\ncanonical splice\njunctions")) +
    labs(title = "Good Quality Control Attributes Across Cells",
         x = "",
         y = "Percentage, %") + #       Maybe we need to do two/three plots. Y axis is different 
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
  ### Basic cell informtion ###
  grid.arrange(tableGrob(summary(SQANTI_cell_summary[,12:20]), rows = NULL),
               tableGrob(summary(SQANTI_cell_summary[,21:29]), rows = NULL),
               nrow=2)
  grid.arrange(gg_reads_in_cells, gg_umis_in_cells, ncol=2)
  grid.arrange(gg_genes_in_cells, gg_JCs_in_cell, ncol=2)
  print(gg_annotation_of_genes_in_cell)
  print(gg_MT_perc)
  ### Read lengths ###
  print(gg_bulk_all_reads)
  print(gg_read_distr)
  print(gg_read_distr_mono)
  grid.arrange(gg_read_distr, gg_read_distr_mono, nrow=2)
  print(gg_FSM_read_distr)
  print(gg_FSM_mono_read_distr)
  grid.arrange(gg_FSM_read_distr, gg_FSM_mono_read_distr, nrow=2)
  print(gg_ISM_read_distr)
  print(gg_ISM_mono_read_distr)
  grid.arrange(gg_ISM_read_distr, gg_ISM_mono_read_distr, nrow=2)
  print(gg_NIC_read_distr)
  print(gg_NIC_mono_read_distr)
  grid.arrange(gg_NIC_read_distr, gg_NIC_mono_read_distr, nrow=2)
  print(gg_NNC_read_distr)
  print(gg_NNC_mono_read_distr)
  grid.arrange(gg_NNC_read_distr, gg_NNC_mono_read_distr, nrow=2)
  print(gg_genic_read_distr)
  print(gg_genic_mono_read_distr)
  grid.arrange(gg_genic_read_distr, gg_genic_mono_read_distr, nrow=2)
  print(gg_antisense_read_distr)
  print(gg_antisense_mono_read_distr)
  grid.arrange(gg_antisense_read_distr, gg_antisense_mono_read_distr, nrow=2)
  print(gg_fusion_read_distr)
  print(gg_fusion_mono_read_distr)
  grid.arrange(gg_fusion_read_distr, gg_fusion_mono_read_distr, nrow=2)
  print(gg_intergenic_read_distr)
  print(gg_intergenic_mono_read_distr)
  grid.arrange(gg_intergenic_read_distr, gg_intergenic_mono_read_distr, nrow=2)
  print(gg_genic_intron_read_distr)
  print(gg_genic_intron_mono_read_distr)
  grid.arrange(gg_genic_intron_read_distr, gg_genic_intron_mono_read_distr, nrow=2)
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
  print(gg_coding_across_category)
  print(gg_non_coding_across_category)
  grid.arrange(gg_coding_across_category, gg_non_coding_across_category, nrow=2)
  ### Coverage (TSS/TES in the future) ###
  print(gg_ref_coverage_across_category)
  ### Unique Splice Junctions ###
  print(gg_known_novel_canon)
  ### Bad features ###
  print(gg_bad_feature)
  ### Bad features by structural category ###
  print(gg_intrapriming_by_category)
  print(gg_RTS_by_category)
  print(gg_noncanon_by_category)
  ### Good features ###
  print(gg_good_feature)
  ### Good features by structural category ###
  print(gg_anno_genes_by_category)
  print(gg_canon_by_category)
  dev.off()
}

Classification <- read.table(class.file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
SQANTI_cell_summary <- calculate_metrics_per_cell(
  Classification, 
  cell_summary_output, 
  Save = save_option
)

generate_sqantisc_plots(
  SQANTI_cell_summary, 
  Classification, 
  report_output
)