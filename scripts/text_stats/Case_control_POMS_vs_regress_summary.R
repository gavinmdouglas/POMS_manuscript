rm(list = ls(all.names = TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/")

POMS_obesity_ko_sig <- read.table(file = "/home/gavin/github_repos/POMS_manuscript/display_items/Supp_obesity_sig_KOs_RAW.tsv",
                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

POMS_obesity_higher_sig <- read.table(file = "/home/gavin/github_repos/POMS_manuscript/display_items/Supp_obesity_sig_higher_levels_RAW.tsv",
                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")


regress_Almeida_sig <- read.table("/home/gavin/github_repos/POMS_manuscript/key_datafiles/almeida_regress_results.tsv",
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

regress_Almeida_sig_obesity <- regress_Almeida_sig[grep("Obesity", regress_Almeida_sig$dataset), ]
regress_Almeida_sig_CRC <- regress_Almeida_sig[-grep("Obesity", regress_Almeida_sig$dataset), ]

length(which(POMS_obesity_ko_sig$func %in% regress_Almeida_sig_obesity$func))
POMS_obesity_ko_sig[which(POMS_obesity_ko_sig$func %in% regress_Almeida_sig_obesity$func), ]

length(which(POMS_obesity_higher_sig$func %in% regress_Almeida_sig_obesity$func))
POMS_obesity_higher_sig[which(POMS_obesity_higher_sig$func %in% regress_Almeida_sig_obesity$func), ]

