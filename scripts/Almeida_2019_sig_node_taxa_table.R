rm(list=ls(all.names=TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/projects/POMS/MAGs/Almeida2019/")

almeida_taxa <- readRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/Almeida_2019_POMS_output/taxa_table.rds")
almeida_abun <- read.table(file = "mapping_results/modified/bwa_depth_min25coverage.tsv", header=TRUE, sep="\t", check.names=FALSE,
                           row.names=1, quote="", comment.char="")

almeida_sample_info <- read.table("MGS_samples_info_SuppTable1.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")

ERP002061_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP002061"), ]


ERP002061_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP002061_almeida_sample_info$Run)

ERP002061_almeida_sample_info <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Run %in% colnames(ERP002061_almeida_abun)), ]
ERP002061_group1_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP002061_group2_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP002061_out <- readRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/Almeida_2019_POMS_output/ERP002061_POMS_out.rds")

# Write out example table of predominant taxa on each side of significant nodes.
ERP002061_sig_node_taxa <- sapply(ERP002061_out$sig_nodes,
                                  function(x) {
                                    node_taxa(lhs_features = ERP002061_out$balances_info$features[[x]]$lhs,
                                              rhs_features = ERP002061_out$balances_info$features[[x]]$rhs,
                                              taxa=almeida_taxa)
                                  })


ERP002061_sig_node_direction <- pairwise_mean_direction_and_wilcoxon(in_list = ERP002061_out$balances_info$balances,
                                                                     group1 = ERP002061_group1_samples,
                                                                     group2=ERP002061_group2_samples,
                                                                     skip_wilcoxon = TRUE)$mean_direction


ERP002061_sig_node_taxa_df <- data.frame(nodes=names(ERP002061_sig_node_taxa), taxa=ERP002061_sig_node_taxa,
                                         dir=ERP002061_sig_node_direction[names(ERP002061_sig_node_taxa)], stringsAsFactors = FALSE)

ERP002061_sig_node_taxa_df$dir[which(ERP002061_sig_node_taxa_df$dir == "group1")] <- "Obese"
ERP002061_sig_node_taxa_df$dir[which(ERP002061_sig_node_taxa_df$dir == "group2")] <- "Control"

write.table(x = ERP002061_sig_node_taxa_df, file = "/home/gavin/github_repos/POMS_manuscript/data/Almeida_2019_tables/ERP002061_POMS_sig_node_taxa.tsv",
            sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
