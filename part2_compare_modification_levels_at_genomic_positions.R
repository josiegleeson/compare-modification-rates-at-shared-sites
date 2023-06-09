
# execute as:
# Rscript part2_compare_modification_levels_at_genomic_positions.R annotated_modification_sites.csv brain_region

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(reshape2)
  library(rstatix)
  library(purrr)
  library(lme4)
  library(multcomp)
})
options(dplyr.summarise.inform = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# command line variables
annotated_modification_sites <- fread(args[1])
brain_region <- args[2]

# filter for at least 3 obs per site within a transcript and multiple transcripts per gene
mod_sites_for_testing <- annotated_modification_sites %>% 
  dplyr::filter(group_id == paste0(brain_region)) %>% 
  group_by(gene_id, genome_pos, transcript_id) %>% 
  dplyr::filter(n() >= 2) %>% 
  dplyr::filter(length(unique(mod_ratio)) > 1) %>% 
  group_by(gene_id, genome_pos) %>% 
  dplyr::filter(length(unique(transcript_id)) > 1) %>% 
  ungroup() %>% 
  mutate(genome_locus = paste0(gene_id, "_", genome_pos),
         modified_reads = mod_ratio * n_reads) %>% 
  mutate(across(where(is.character), as.factor))

# function to perform the analysis comparing mod_ratios of different txs per genome position
perform_glmm_analysis <- function(df) {
  # fit GLMM
  model <- lme4::glmer(cbind(modified_reads, n_reads - modified_reads) ~ transcript_id + (1|sample_id), 
                       data = df, family = binomial)
  
  # perform pairwise comps
  comp <- glht(model, mcp(transcript_id = "Tukey"))
  sum_comp <- summary(comp)
  
  # information to dataframe
  df_comp <- data.frame(
    estimate = sum_comp$test$coefficients,
    SE = sum_comp$test$sigma,
    z_value = sum_comp$test$tstat,
    p_value = sum_comp$test$pvalues
  )
  
  # add row names as a new column
  df_comp$row_names <- rownames(df_comp)
  #message(rownames(df_comp))
  rownames(df_comp) <- NULL
  # separate row_names into group1 and group2
  df_comp <- df_comp %>% separate(row_names, into = c("group1", "group2"), sep = " - ")
  
  # make a metadata df of the desired extra info
  other_data <- df %>% 
    group_by(transcript_id, transcript_position, genome_locus, region) %>% 
    summarise(median_mod_ratio_group1 = median(mod_ratio),
              mean_mod_ratio_group1 = mean(mod_ratio),
              sd_mod_ratio_group1 = sd(mod_ratio),
              mean_n_reads_group1 = mean(n_reads),
              dist_exon_junc_group1 = max(dist_exon_junc)
    ) %>% ungroup() %>% rename("region" = "region_group1")
  
  merge_cols <- c("transcript_id")
  
  # rename for merging
  df_comp <- df_comp %>% rename("group1" = "transcript_id")
  
  m6a_diff_merged <- merge(df_comp, other_data, by=merge_cols)
  
  # make metadata df for other tx in comparison
  other_data <- df %>% 
    group_by(transcript_id, transcript_position, genome_locus, region) %>% 
    summarise(median_mod_ratio_group2 = median(mod_ratio),
              mean_mod_ratio_group2 = mean(mod_ratio),
              sd_mod_ratio_group2 = sd(mod_ratio),
              mean_n_reads_group2 = mean(n_reads),
              dist_exon_junc_group2 = max(dist_exon_junc)
    ) %>% ungroup() %>% rename("region" = "region_group2")
  
  # rename for merging
  m6a_diff_merged <- m6a_diff_merged %>% rename("transcript_id" = "group1") %>% rename("group2" = "transcript_id")
  
  m6a_diff_metadata <- merge(m6a_diff_merged, other_data, by=merge_cols)
  
  # tidy up output
  m6a_diff_metadata <- m6a_diff_metadata %>%
    rename("transcript_id" = "group2",
           "transcript_position.x" = "transcript_position_group1",
           "transcript_position.y" = "transcript_position_group2",
           "genome_locus.x" = "genome_locus_group1",
           "genome_locus.y" = "genome_locus_group2") %>%
    mutate(mod_ratio_difference = median_mod_ratio_group1 - median_mod_ratio_group2)
  
  return(m6a_diff_metadata)
}

# split df by genomic position
list_mod_sites_for_testing <- split(mod_sites_for_testing, mod_sites_for_testing$genome_locus)

# apply function to all genomic positions
mod_sites_results_list <- map(list_mod_sites_for_testing, perform_glmm_analysis)

message("Tested genomic positions for differences")

# combine results output
mod_sites_results <- bind_rows(mod_sites_results_list, .id = "genome_locus")

# adjust p values
mod_sites_results$adjusted_p_values <- p.adjust(mod_sites_results$p_value, method="bonferroni", n=nrow(mod_sites_results))
mod_sites_results <- mod_sites_results[, c(1, 3, 2, 4:ncol(mod_sites_results))]

# export file
write.csv(mod_sites_results, paste0("glmm_results_", brain_region, ".csv"), row.names=F, quote=F)

message("Complete")

