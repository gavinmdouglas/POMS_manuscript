### Parse Almeida read mapping files to produce tables of MAG abundances across samples.
### Will output two sets based on calling MAGs as positive based on 25% and 50% coverage.
### Also prep the KO table at the end of this file.

rm(list=ls(all.names=TRUE))

library(tidyverse)

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset/mapping_results/")

almeida_coverage <- read_csv("bwa_coverage.csv")
almeida_depth <- read_csv("bwa_depth.csv")

sample_col <- colnames(almeida_coverage)
sample_col <- sample_col[2:length(sample_col)]

almeida_depth_25min <- almeida_depth
almeida_depth_50min <- almeida_depth

for(s in sample_col) {
  missing_MAGs_25min <- almeida_coverage[which(almeida_coverage[, s] < 25), "Genome"]$Genome
  almeida_depth_25min[which(almeida_depth_25min$Genome %in% missing_MAGs_25min), s] <- 0
  
  missing_MAGs_50min <- almeida_coverage[which(almeida_coverage[, s] < 50), "Genome"]$Genome
  almeida_depth_50min[which(almeida_depth_50min$Genome %in% missing_MAGs_50min), s] <- 0
}

write_delim(x = almeida_depth_25min,
            path="../mapping_results/modified/bwa_depth_min25coverage.tsv",
            delim = "\t", col_names = TRUE)

write_delim(x = almeida_depth_50min,
            path="../mapping_results/modified/bwa_depth_min50coverage.tsv",
            delim = "\t", col_names = TRUE)


# Also parse KO table to be in format expected by PICRUSt2.
almedia_ko <- read_csv("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset/functional_analyses/kegg_summary.csv")

colnames(almedia_ko)[1] <- "function"

write_delim(x = almedia_ko,
            path="/home/gavin/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset/functional_analyses/modified/kegg_summary.tsv",
            delim = "\t", col_names = TRUE)
