# execute as:
# Rscript name.R high_confidence_modification_sites_from_m6anet.csv gencode.gtf brain_region

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(reshape2)
  library(rstatix)
  library(purrr)
  library(lme4)
  library(multcomp)
  library(GenomicRanges)
  library(GenomicFeatures)
})
options(dplyr.summarise.inform = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# command line variables
modification_sites <- fread(args[1])
gtf_file <- args[2]
brain_region <- args[3]

#setwd("~/Documents/brain-drs/modification_results/genome_transcriptome_coords/wilcox_test_mod_ratios_of_mod_sites/")
#modification_sites <- fread("~/Documents/brain-drs/modification_results/2023_results/all_modified_sites_0.9.csv")
#modification_sites <- modification_sites %>% dplyr::filter(gene_id == "ENSG00000003056" | gene_id == "ENSG00000001460" | gene_id == "ENSG00000005175")

# import gtf as txdb
#txdb <- makeTxDbFromGFF(file = "~/Documents/brain-drs/gencode.v31.annotation.gtf", format = "gtf")
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")

# remove version numbers
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)
exon.names <- data.frame(names=names(exons)) %>% separate(names, into=c("names", NULL), sep="\\.", extra="drop")
names(exons) <- exon.names$names

cds <- cdsBy(txdb, by="tx", use.names=TRUE)
cds.names <- data.frame(names=names(cds)) %>% separate(names, into=c("names", NULL), sep="\\.", extra="drop")
names(cds) <- cds.names$names

# function to convert tx coords to genome coords and annotate positions
convert_to_genome_coordinates <- function(tx_name, tx_pos, exonsdb, cdsdb) {
  # check if transcript_id is in the exonsdb
  # if(!(tx_name %in% names(exonsdb))) {
  #   warning(paste("Transcript ", tx_name, " not found in exonsdb"))
  #   return(list(genome_pos=NA,chromosome=NA,distance_to_nearest_exon_junction=NA,region=NA))
  # }
  
  # subset to one tx
  tx_exons <- exonsdb[[tx_name]]
  chrom <- as.character(seqnames(tx_exons)[1])

  # check tx has exons
  if(length(tx_exons) == 0) {
    return(list(genome_pos=NA,chromosome=NA,distance_to_nearest_exon_junction=NA,region=NA))
  }
  
  # get cumulative length of the exons
  exon_lengths <- width(tx_exons)
  cum_lengths <- cumsum(exon_lengths)
  
  # get exon where the transcript position is located
  idx <- which(cum_lengths >= tx_pos)[1]

  # check transcript position is within an exon
  if(is.na(idx)) {
    warning(paste("Transcript position", tx_pos, "is out of range for transcript", tx_name))
    return(list(genome_pos=NA,chromosome=NA,distance_to_nearest_exon_junction=NA,region=NA))
  }
  
  # get relative position within the exon
  rel_pos <- cum_lengths[idx] - tx_pos
  
  # convert to genomic coordinates
  strand <- as.character(strand(tx_exons[idx])@values)
  if (strand == "+") {
    genome_pos <- end(tx_exons[idx]) - rel_pos + 1
  } else if (strand == "-") {
    genome_pos <- start(tx_exons[idx]) + rel_pos - 1
  }
  
  # calculate distance to nearest exon junction
  if (idx == 1) { # if first exon
    distance_to_nearest_exon_junction <- rel_pos
  } else if (idx == length(cum_lengths)) { # if last exon
    distance_to_nearest_exon_junction <- exon_lengths[idx] - rel_pos
  } else { # if internal exon
    distance_to_nearest_exon_junction <- min(exon_lengths[idx] - rel_pos, rel_pos - cum_lengths[idx-1])
  }
  
  # get region in a tx (CDS, 5' UTR, 3' UTR)
  cdstx <- cdsdb[[tx_name]]
  
  if (length(cdstx) == 0) { # if no CDS, tx is non-coding
    region <- "ncRNA"
  } else {
    cds_cum_lengths <- cumsum(width(cdstx))
    if (tx_pos <= cds_cum_lengths[1]) {
      region <- "UTR5"
    } else if (tx_pos > cds_cum_lengths[length(cds_cum_lengths)]) {
      region <- "UTR3"
    } else {
      region <- "CDS"
    }
  }
  return(list(genome_pos=as.integer(genome_pos),chromosome=chrom,dist_exon_junc=as.integer(distance_to_nearest_exon_junction),region=region))
}

annotated_modification_sites <- modification_sites %>%
  mutate(output = mapply(convert_to_genome_coordinates, transcript_id, transcript_position, MoreArgs = list(exonsdb = exons, cdsdb = cds), SIMPLIFY = FALSE)) %>%
  tidyr::unnest_wider(output, names_sep = "_")

# remove 'output_' from column names
names(annotated_modification_sites) <- gsub("^output_", "", names(annotated_modification_sites))

message("Annotated modification sites")

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

