rm(list = ls(all.names = TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/results/Almeida_2019_POMS_output/")

almeida_taxa <- readRDS(file = "../../key_inputs/Almeida2019_dataset/taxonomy/taxa_table.rds")

ERP002061_out <- readRDS(file = "ERP002061_POMS_out.rds")

ERP002061_out_tree_prep <- prep_tree_sig_nodes(in_list = ERP002061_out, taxa_table = almeida_taxa)

almeida_abun <- read.table(file = "../../key_inputs/Almeida2019_dataset/mapping_results/modified/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_sample_info <- read.table("../../key_inputs/Almeida2019_dataset/MGS_samples_info_SuppTable1.txt.gz",
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

# Figure out what direction sig. nodes were in.
almeida_sample_info <- read.table("../../key_inputs/Almeida2019_dataset/MGS_samples_info_SuppTable1.txt.gz", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
ERP002061_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study ==  "ERP002061"), ]

ERP002061_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP002061_almeida_sample_info$Run)


ERP002061_almeida_sample_info <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Run %in% colnames(ERP002061_almeida_abun)), ]
ERP002061_group1_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state ==  "Diseased"), "Run"]
ERP002061_group2_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state ==  "Healthy"), "Run"]

ERP002061_group1_samples <- ERP002061_group1_samples[which(ERP002061_group1_samples %in% names(ERP002061_out$balances_info$balances[[1]]))]
ERP002061_group2_samples <- ERP002061_group2_samples[which(ERP002061_group2_samples %in% names(ERP002061_out$balances_info$balances[[1]]))]

ERP002061_node_mean_dir <- pairwise_mean_direction_and_wilcoxon(ERP002061_out$balances_info$balances[ERP002061_out$sig_nodes],
                                                                ERP002061_group1_samples, ERP002061_group2_samples,skip_wilcoxon = TRUE)$mean_direction

ERP002061_out_tree_prep$disease_node <- which(ERP002061_out$tree$node.label %in% names(ERP002061_node_mean_dir)[which(ERP002061_node_mean_dir ==  "group1")]) + length(ERP002061_out$tree$tip.label)
ERP002061_out_tree_prep$control_node <- which(ERP002061_out$tree$node.label %in% names(ERP002061_node_mean_dir)[which(ERP002061_node_mean_dir ==  "group2")]) + length(ERP002061_out$tree$tip.label)
ERP002061_nonsig_nodes <- names(ERP002061_out$balances_info$balances)[which(!names(ERP002061_out$balances_info$balances) %in% ERP002061_out$sig_nodes)]
ERP002061_out_tree_prep$nonsig <- which(ERP002061_out$tree$node.label %in% ERP002061_nonsig_nodes) + length(ERP002061_out$tree$tip.label)



ERP002061_sig_nodes_tree <- ggtree(ERP002061_out_tree_prep$prepped_tree, layout = "circular") +
                                  geom_nodepoint(aes(subset = (node %in% ERP002061_out_tree_prep$nonsig)), fill = "white", color = "black", alpha = 0.75, size = 4, pch = 21) +
                                  geom_nodepoint(aes(subset = (node %in% c(ERP002061_out_tree_prep$disease_node, ERP002061_out_tree_prep$control_node))), color = "black", fill = "gold", alpha = 0.75, size = 4, pch = 21) +
                                  ggtitle("Obesity 1: All significant nodes") +
                                  theme(plot.title = element_text(hjust = 0.5, vjust = -27.5),
                                        plot.margin = unit(c(-35, -10, -25, -25), "mm"))


ERP002061_out_num_enriched <- ERP002061_out$df[, c("num_sig_nodes_group1_enrich", "num_sig_nodes_group2_enrich")]

ERP002061_out_num_enriched <- ddply(ERP002061_out_num_enriched, .(num_sig_nodes_group1_enrich, num_sig_nodes_group2_enrich), nrow)

colnames(ERP002061_out_num_enriched)[3] <- "Count"

ERP002061_out_num_enriched_heatmap <- ggplot(ERP002061_out_num_enriched, aes(x = num_sig_nodes_group1_enrich, y = num_sig_nodes_group2_enrich, fill = log10(Count))) +
                                            geom_tile() +
                                            xlab("No. obesity-enriched nodes") +
                                            ylab("No. control-enriched nodes") +
                                            labs(fill = expression("log"[10]*"(Count)")) +
                                            scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 20)) +
                                            scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 20)) +
                                            theme_bw() +
                                            ggtitle("Obesity 1: All gene family enrichments") +
                                            theme(legend.position = c(0.865, 0.7),
                                                  legend.box.background = element_rect(colour = "black"),
                                                  legend.title = element_text(size = 7), 
                                                  legend.text = element_text(size = 7),
                                                  plot.title = element_text(hjust = 0.5)) +
                                            guides(shape = guide_legend(override.aes = list(size = 0.5)),
                                                   color = guide_legend(override.aes = list(size = 0.5)))

ERP002061_sig_nodes_tree_and_heatmap <- plot_grid(ERP002061_sig_nodes_tree, ERP002061_out_num_enriched_heatmap, nrow = 1, labels = c('a', 'b'))



ERP003612_out <- readRDS(file = "ERP003612_POMS_out.rds")

ERP003612_out_tree_prep <- prep_tree_sig_nodes(in_list = ERP003612_out, taxa_table = almeida_taxa)

# Figure out what direction sig. nodes were in.
almeida_sample_info <- read.table("../../key_inputs/Almeida2019_dataset/MGS_samples_info_SuppTable1.txt.gz", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
ERP003612_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study ==  "ERP003612"), ]

ERP003612_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP003612_almeida_sample_info$Run)

ERP003612_almeida_sample_info <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Run %in% colnames(ERP003612_almeida_abun)), ]
ERP003612_group1_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state ==  "Diseased"), "Run"]
ERP003612_group2_samples <- ERP003612_almeida_sample_info[which(ERP003612_almeida_sample_info$Health.state ==  "Healthy"), "Run"]

ERP003612_group1_samples <- ERP003612_group1_samples[which(ERP003612_group1_samples %in% names(ERP003612_out$balances_info$balances[[1]]))]
ERP003612_group2_samples <- ERP003612_group2_samples[which(ERP003612_group2_samples %in% names(ERP003612_out$balances_info$balances[[1]]))]

ERP003612_node_mean_dir <- pairwise_mean_direction_and_wilcoxon(ERP003612_out$balances_info$balances[ERP003612_out$sig_nodes],
                                                                ERP003612_group1_samples, ERP003612_group2_samples,skip_wilcoxon = TRUE)$mean_direction

ERP003612_out_tree_prep$disease_node <- which(ERP003612_out$tree$node.label %in% names(ERP003612_node_mean_dir)[which(ERP003612_node_mean_dir ==  "group1")]) + length(ERP003612_out$tree$tip.label)
ERP003612_out_tree_prep$control_node <- which(ERP003612_out$tree$node.label %in% names(ERP003612_node_mean_dir)[which(ERP003612_node_mean_dir ==  "group2")]) + length(ERP003612_out$tree$tip.label)
ERP003612_nonsig_nodes <- names(ERP003612_out$balances_info$balances)[which(!names(ERP003612_out$balances_info$balances) %in% ERP003612_out$sig_nodes)]
ERP003612_out_tree_prep$nonsig <- which(ERP003612_out$tree$node.label %in% ERP003612_nonsig_nodes) + length(ERP003612_out$tree$tip.label)



ERP003612_sig_nodes_tree <- ggtree(ERP003612_out_tree_prep$prepped_tree, layout = "circular") +
                                    geom_nodepoint(aes(subset = (node %in% ERP003612_out_tree_prep$nonsig)), fill = "white", color = "black", alpha = 0.75, size = 4, pch = 21) +
                                    geom_nodepoint(aes(subset = (node %in% c(ERP003612_out_tree_prep$disease_node, ERP003612_out_tree_prep$control_node))), color = "black", fill = "gold", alpha = 0.75, size = 4, pch = 21) +
                                    ggtitle("Obesity 2: All significant nodes") +
                                    theme(plot.title = element_text(hjust = 0.5, vjust = -27.5),
                                          plot.margin = unit(c(-35, -10, -25, -25), "mm"))


ERP003612_out_num_enriched <- ERP003612_out$df[, c("num_sig_nodes_group1_enrich", "num_sig_nodes_group2_enrich")]

ERP003612_out_num_enriched <- ddply(ERP003612_out_num_enriched, .(num_sig_nodes_group1_enrich, num_sig_nodes_group2_enrich), nrow)

colnames(ERP003612_out_num_enriched)[3] <- "Count"

ERP003612_out_num_enriched_heatmap <- ggplot(ERP003612_out_num_enriched, aes(x = num_sig_nodes_group1_enrich, y = num_sig_nodes_group2_enrich, fill = log10(Count))) +
                                              geom_tile() +
                                              xlab("No. obesity-enriched nodes") +
                                              ylab("No. control-enriched nodes") +
                                              labs(fill = expression("log"[10]*"(Count)")) +
                                              scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 20)) +
                                              scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 20)) +
                                              theme_bw() +
                                              ggtitle("Obesity 2: All gene family enrichments") +
                                              theme(legend.position = c(0.865, 0.7),
                                                    legend.box.background = element_rect(colour = "black"),
                                                    legend.title = element_text(size = 7), 
                                                    legend.text = element_text(size = 7),
                                                    plot.title = element_text(hjust = 0.5)) +
                                              guides(shape = guide_legend(override.aes = list(size = 0.5)),
                                                     color = guide_legend(override.aes = list(size = 0.5)))

ERP003612_sig_nodes_tree_and_heatmap <- plot_grid(ERP003612_sig_nodes_tree, ERP003612_out_num_enriched_heatmap, nrow = 1, labels = c('c', 'd'))

# Saving the top and bottom of t

ggsave(filename = "../../../figures/Maintext_obesity1_big_picture_RAW.pdf",
       plot = ERP002061_sig_nodes_tree_and_heatmap,
       width = 8, height = 4, dpi = 600, device = "pdf")

ggsave(filename = "../../../figures/Maintext_obesity2_big_picture_RAW.pdf",
       plot = ERP003612_sig_nodes_tree_and_heatmap,
       width = 8, height = 4, dpi = 600, device = "pdf")
