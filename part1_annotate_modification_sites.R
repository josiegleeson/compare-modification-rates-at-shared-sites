# execute as:
# Rscript name.R high_confidence_modification_sites_from_m6anet.csv gencode.gtf

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(reshape2)
  library(purrr)
  library(GenomicRanges)
  library(GenomicFeatures)
})
options(dplyr.summarise.inform = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# command line variables
modification_sites <- fread(args[1])
gtf_file <- args[2]

#setwd("~/Documents/brain-drs/modification_results/genome_transcriptome_coords/stats_testing_mod_ratios_of_mod_sites/")
#modification_sites <- fread("~/Documents/brain-drs/modification_results/genome_transcriptome_coords/all_possible_mod_sites_0.9_prob.csv")
# tx_test <- c("ENST00000366620", "ENST00000366621")
# tx_pos_test <- c(1335,1534)
#modification_sites <- modification_sites %>% dplyr::filter(transcript_id %in% tx_test) %>% dplyr::filter(transcript_position %in% tx_pos_test)

# import gtf as txdb
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")

rtrack_gtf <- rtracklayer::import(gtf_file)
rtrack_df <- data.frame(gene_name = rtrack_gtf$gene_name,
                        transcript_id = rtrack_gtf$transcript_id,
                        transcript_type = rtrack_gtf$transcript_type) %>% 
  separate(transcript_id, into=c("transcript_id", NULL), sep="\\.", extra="drop") %>% 
  group_by(transcript_id) %>% 
  slice_head(n=1) %>% 
  ungroup()

# remove version numbers
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)
exon.names <- data.frame(names=names(exons)) %>% separate(names, into=c("names", NULL), sep="\\.", extra="drop")
names(exons) <- exon.names$names

tx_features <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)

tx_features <- tx_features %>% separate(tx_name, into=c("transcript_id", NULL), sep="\\.", extra="drop")
tx_features <- tx_features %>% dplyr::filter(!endsWith(gene_id, "Y"))

# filter the reference data for only txs with modification sites
tx_features_filt <- tx_features %>% dplyr::filter(transcript_id %in% modification_sites$transcript_id)
exons_filt <- exons[names(exons) %in% modification_sites$transcript_id]

# function to convert tx coords to genome coords and annotate positions
convert_to_genome_coordinates <- function(
  
  tx_name, tx_pos, exonsdb, txfdb
  
) {
  
  # check current transcript_id is in the exonsdb
  if(!(tx_name %in% names(exonsdb))) {
    warning(paste("Transcript ", tx_name, " not found in exonsdb"))
    return(list(genome_pos=NA,chromosome=NA,distance_to_nearest_exon_junction=NA,region=NA))
  }
  
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
    warning(paste("Transcript position", tx_pos, "is out of range for transcript ", tx_name))
    return(list(genome_pos=NA,chromosome=NA,distance_to_nearest_exon_junction=NA,region=NA))
  }
  
  # get relative position to the end site of exon
  rel_pos <- cum_lengths[idx] - tx_pos
  
  # convert to genomic coordinates
  strand <- as.character(strand(tx_exons[idx])@values)
  if (strand == "+") {
    genome_pos <- end(tx_exons[idx]) - rel_pos + 1
  } else if (strand == "-") {
    genome_pos <- start(tx_exons[idx]) + rel_pos - 1
  }
  
  # calculate distance to nearest exon junction
  # includes start or end of transcript as a junction
  
  if (idx == 1) { # if first exon
    distance_to_upstream_exon_junction <- tx_pos
    distance_to_downstream_exon_junction <- rel_pos
    exon_type <- "first"
  } else if (idx == length(cum_lengths)) { # if last exon
    distance_to_upstream_exon_junction <- exon_lengths[idx] - rel_pos
    distance_to_downstream_exon_junction <- rel_pos
    exon_type <- "last"
  } else { # if internal exon
    distance_to_upstream_exon_junction <- exon_lengths[idx] - rel_pos
    distance_to_downstream_exon_junction <- rel_pos
    exon_type <- "internal"
  }
  
  # get region in a tx (CDS, 5' UTR, 3' UTR)
  tx_feat <- tx_features_filt %>% dplyr::filter(transcript_id %in% tx_name)
  if (tx_feat$cds_len == 0) { # if no CDS, tx is non-coding
    region <- "ncRNA"
  } else {
    if (tx_pos <= tx_feat$utr5_len) {
      region <- "UTR5"
    } else if (tx_pos > tx_feat$utr5_len & tx_pos <= (tx_feat$utr5_len + tx_feat$cds_len)) {
      region <- "CDS"
    } else if (tx_pos > (tx_feat$utr5_len + tx_feat$cds_len) & tx_pos <= tx_feat$tx_len) {
      region <- "UTR3"
    } else {
      region <- NA
    }
  }
  
  return(list(genome_pos=as.integer(genome_pos),chromosome=chrom,transcript_length=as.integer(tx_feat$tx_len),dist_up_exon_junc=as.integer(distance_to_upstream_exon_junction),dist_down_exon_junc=as.integer(distance_to_downstream_exon_junction),exon_type=as.factor(exon_type),region=as.factor(region)))

}

# apply function
annotated_modification_sites <- modification_sites %>%
  mutate(output = mapply(convert_to_genome_coordinates, transcript_id, transcript_position, MoreArgs = list(exonsdb = exons_filt, txfdb = tx_features_filt), SIMPLIFY = FALSE)) %>%
  tidyr::unnest_wider(output, names_sep = "_")

# remove 'output_' from column names
names(annotated_modification_sites) <- gsub("^output_", "", names(annotated_modification_sites))

# add in transcript biotypes
annotated_modification_sites <- merge(annotated_modification_sites, rtrack_df, by="transcript_id", all.x=T, all.y=F)

write.csv(annotated_modification_sites, "annotated_modification_sites.csv", row.names=F, quote=F)

message("Annotated modification sites")
