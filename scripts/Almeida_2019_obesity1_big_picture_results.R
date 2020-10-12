rm(list=ls(all.names=TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/")

# Note that plot.margin gives parameters from the top clockwise.

almeida_taxa <- readRDS(file = "Almeida_2019_POMS_output/taxa_table.rds")

ERP002061_out <- readRDS(file = "Almeida_2019_POMS_output/ERP002061_POMS_out.rds")

ERP002061_out_tree_prep <- prep_tree_sig_nodes(in_list = ERP002061_out, taxa_table = almeida_taxa)

# Figure out what direction sig. nodes were in.
almeida_sample_info <- read.table("/home/gavin/projects/POMS/MAGs/Almeida2019/MGS_samples_info_SuppTable1.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")
ERP002061_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP002061"), ]
ERP002061_almeida_sample_info <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Run %in% colnames(ERP002061_almeida_abun)), ]
ERP002061_group1_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP002061_group2_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP002061_group1_samples <- ERP002061_group1_samples[which(ERP002061_group1_samples %in% names(ERP002061_out$balances_info$balances[[1]]))]
ERP002061_group2_samples <- ERP002061_group2_samples[which(ERP002061_group2_samples %in% names(ERP002061_out$balances_info$balances[[1]]))]

ERP002061_node_mean_dir <- pairwise_mean_direction_and_wilcoxon(ERP002061_out$balances_info$balances[ERP002061_out$sig_nodes],
                                                      ERP002061_group1_samples, ERP002061_group2_samples,skip_wilcoxon=TRUE)$mean_direction

ERP002061_out_tree_prep$disease_node <- which(ERP002061_out$tree$node.label %in% names(ERP002061_node_mean_dir)[which(ERP002061_node_mean_dir == "group1")]) + length(ERP002061_out$tree$tip.label)
ERP002061_out_tree_prep$control_node <- which(ERP002061_out$tree$node.label %in% names(ERP002061_node_mean_dir)[which(ERP002061_node_mean_dir == "group2")]) + length(ERP002061_out$tree$tip.label)
ERP002061_nonsig_nodes <- names(ERP002061_out$balances_info$balances)[which(! names(ERP002061_out$balances_info$balances) %in% ERP002061_out$sig_nodes)]
ERP002061_out_tree_prep$nonsig <- which(ERP002061_out$tree$node.label %in% ERP002061_nonsig_nodes) + length(ERP002061_out$tree$tip.label)



ERP002061_sig_nodes_tree <- ggtree(ERP002061_out_tree_prep$prepped_tree, layout = "circular") +
                                    geom_nodepoint(aes(subset=(node %in% ERP002061_out_tree_prep$nonsig)), fill="gold", color="black", alpha=0.75, size=4, pch=21) +
                                    geom_nodepoint(aes(subset=(node %in% ERP002061_out_tree_prep$disease_node)), color="black", fill="red", alpha=0.75, size=4, pch=21) +
                                    geom_nodepoint(aes(subset=(node %in% ERP002061_out_tree_prep$control_node)), color="black", fill="blue", alpha=0.75, size=4, pch=21) +
                                    ggtitle("All significant nodes") +
                                    theme(plot.title = element_text(hjust = 0.5, vjust=-27.5),
                                          plot.margin=unit(c(-35,-10,-25,-25), "mm"))


ERP002061_out_num_enriched <- ERP002061_out$df[, c("num_sig_nodes_pos_enrich", "num_sig_nodes_neg_enrich")]

ERP002061_out_num_enriched <- ddply(ERP002061_out_num_enriched, .(num_sig_nodes_pos_enrich, num_sig_nodes_neg_enrich), nrow)

colnames(ERP002061_out_num_enriched)[3] <- "Count"

ERP002061_out_num_enriched_heatmap <- ggplot(ERP002061_out_num_enriched, aes(x=num_sig_nodes_pos_enrich, y=num_sig_nodes_neg_enrich, fill=log10(Count))) +
                                              geom_tile() +
                                              xlab("No. positively enriched nodes") +
                                              ylab("No. negatively enriched nodes") +
                                              labs(fill=expression("log" [10]*"(Count)")) +
                                              scale_x_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 15)) +
                                              scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 12)) +
                                              theme_bw() +
                                              ggtitle("All KO enrichments") +
                                              theme(legend.position = c(0.865, 0.7),
                                                    legend.box.background = element_rect(colour="black"),
                                                    legend.title = element_text(size = 7), 
                                                    legend.text = element_text(size = 7),
                                                    plot.title = element_text(hjust = 0.5)) +
                                              guides(shape = guide_legend(override.aes = list(size = 0.5)),
                                                     color = guide_legend(override.aes = list(size = 0.5)))

ERP002061_sig_nodes_tree_and_heatmap <- plot_grid(ERP002061_sig_nodes_tree, ERP002061_out_num_enriched_heatmap, nrow=1, labels=c('a', 'b'))

ggsave(filename = "../figures/ERP002061_sig_nodes_tree_and_heatmap.pdf", plot = ERP002061_sig_nodes_tree_and_heatmap,
       width = 7.20472, height=3.5)

