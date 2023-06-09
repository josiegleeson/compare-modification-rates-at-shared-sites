
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

# filter the databases for only txs with modification sites
cds_filt <- cds[names(cds) %in% modification_sites$transcript_id]
exons_filt <- exons[names(exons) %in% modification_sites$transcript_id]

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
  mutate(output = mapply(convert_to_genome_coordinates, transcript_id, transcript_position, MoreArgs = list(exonsdb = exons_filt, cdsdb = cds_filt), SIMPLIFY = FALSE)) %>%
  tidyr::unnest_wider(output, names_sep = "_")

# remove 'output_' from column names
names(annotated_modification_sites) <- gsub("^output_", "", names(annotated_modification_sites))

write.csv(annotated_modification_sites, paste0("annotated_modification_sites.csv"), row.names=F, quote=F)

message("Annotated modification sites")

